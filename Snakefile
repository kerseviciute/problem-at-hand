import pandas as pd

configfile: 'config.yml'

samples = pd.read_csv(config['sampleSheet'])

rule all:
  input:
    'output/problem-at-hand/S5/run2/features.csv'

rule download:
  output:
    mat = 'raw_data/{project}/{sid}.mat'
  params:
    url = lambda wildcards: samples[ samples['SID'] == wildcards.sid ]['URL'].iloc[0]
  shell:
    '''
      curl {params.url} --location --silent --output {output.mat}
    '''

rule download_montage:
  output:
    mat = 'montage/montage_left_hemisphere.mat'
  params:
    url = 'https://osf.io/download/uw7hr'
  shell:
    '''
      curl {params.url} --location --silent --output {output.mat}
    '''

rule montage_csv:
  input:
    mat = 'montage/montage_left_hemisphere.mat'
  output:
    csv = 'montage/montage_left_hemisphere.csv'
  conda: 'env/mne.yml'
  script: 'python/montage_to_csv.py'

rule prep:
  input:
    montage = 'montage/montage_left_hemisphere.csv',
    raw = 'raw_data/{project}/{sample}_{run}.mat'
  output:
    raw = 'output/{project}/{sample}/{run}/mne_raw.pkl',
    events = 'output/{project}/{sample}/{run}/events.pkl'
  conda: 'env/mne.yml'
  script: 'python/prepare_data.py'

rule filter:
  input:
    raw = 'output/{project}/{sample}/{run}/mne_raw.pkl'
  output:
    filtered = 'output/{project}/{sample}/{run}/mne_filtered.pkl'
  conda: 'env/mne.yml'
  script: 'python/filter.py'

rule drop_bad_channels:
  input:
    filtered = 'output/{project}/{sample}/{run}/mne_filtered.pkl'
  output:
    good_channels = 'output/{project}/{sample}/{run}/mne_good_channels.pkl'
  conda: 'env/mne.yml'
  script: 'python/drop_bad_channels2.py'

# TODO: terribly needs parallelization and optimization
rule extract_features:
  input:
    final = 'output/{project}/{sample}/{run}/mne_good_channels.pkl',
    events = 'output/{project}/{sample}/{run}/events.pkl'
  output:
    features = 'output/{project}/{sample}/{run}/features.csv'
  conda: 'env/mne.yml'
  script: 'python/extract_features.py'
