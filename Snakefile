import pandas as pd

configfile: 'config.yml'

samples = pd.read_csv(config['sampleSheet'])

rule all:
  input:
    expand('{project}/{page}.html', project = config['project'], page = config['report']['pages'])

rule download:
  output:
    mat = 'raw_data/{project}/{sid}.mat'
  params:
    url = lambda wildcards: samples[samples['SID'] == wildcards.sid]['URL'].iloc[0]
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
    good_channels = 'output/{project}/{sample}/{run}/mne_good_channels.pkl',
    badChannels = 'output/{project}/{sample}/{run}/bad_channels.pkl'
  params:
    dropLowQuality = config['filter']['dropLowQuality']
  conda: 'env/mne.yml'
  script: 'python/drop_bad_channels2.py'

rule extract_features:
  input:
    final = 'output/{project}/{sample}/{run}/mne_good_channels.pkl',
    events = 'output/{project}/{sample}/{run}/events.pkl'
  output:
    features = 'output/{project}/{sample}/{run}/features.csv'
  params:
    feature_length = config['features']['length'],
    min_freq = config['features']['min_freq'],
    max_freq = config['features']['max_freq'],
    freq_step = config['features']['freq_step'],
    nperseg = config['features']['nperseg'],
    relative_frequency = config['features']['relative']
  threads: 10
  conda: 'env/mne.yml'
  script: 'python/extract_features.py'


def all_runs_features(wildcards):
  runs = samples[samples['SampleID'] == wildcards.sample]['Run']
  return [f'output/{{project}}/{{sample}}/{run}/features.csv' for run in runs]


def all_runs_raw_mne(wildcards):
  runs = samples[samples['SampleID'] == wildcards.sample]['Run']
  return [f'output/{{project}}/{{sample}}/{run}/mne_raw.pkl' for run in runs]


def all_runs_filtered_mne(wildcards):
  runs = samples[samples['SampleID'] == wildcards.sample]['Run']
  return [f'output/{{project}}/{{sample}}/{run}/mne_filtered.pkl' for run in runs]


def all_runs_bad_channels(wildcards):
  runs = samples[samples['SampleID'] == wildcards.sample]['Run']
  return [f'output/{{project}}/{{sample}}/{run}/bad_channels.pkl' for run in runs]


rule combine_features:
  input:
    features = all_runs_features
  output:
    feature_matrix = 'output/{project}/{sample}/all/feature_matrix.csv',
    feature_key = 'output/{project}/{sample}/all/key.csv'
  conda: 'env/mne.yml'
  script: 'python/combine_features.py'


###################################################################################################
# REPORTS
###################################################################################################

rule report_summary:
  output:
    report = '{project}/index.html'
  params:
    script = 'reports/index.Rmd'
  conda: 'env/r.yml'
  script: 'R/render.R'

rule report_preprocessing:
  input:
    samples = config['sampleSheet'],
    raw = all_runs_raw_mne,
    filtered = all_runs_filtered_mne,
    bad_channels = all_runs_bad_channels
  output:
    report = '{project}/preprocess_{sample}.html'
  params:
    script = 'reports/preprocess.Rmd'
  threads: 4 # avoid running reports in parallel
  conda: 'env/r.yml'
  script: 'R/render.R'
