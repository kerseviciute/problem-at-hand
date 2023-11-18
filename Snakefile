import pandas as pd

configfile: 'config.yml'

samples = pd.read_csv(config['sampleSheet'])

rule all:
  input:
    expand('output/{project}/report/preprocess_S2_run1.html', project = config['project'])
    # expand('raw_data/{project}/{sid}.mat', project = config['project'], sid = samples['SID'])

rule download:
  output:
    mat = 'raw_data/{project}/{sid}.mat'
  params:
    url = lambda wildcards: samples[ samples['SID'] == wildcards.sid ]['URL'].iloc[0]
  shell:
    '''
      curl {params.url} --location --silent --output {output.mat}
    '''

rule preprocess:
  input:
    notebook = 'preprocess.ipynb',
    mat = 'raw_data/{project}/{sid}.mat'
  output:
    report = 'output/{project}/report/preprocess_{sid}.html'
  conda: 'env/mne.yml'
  script: 'python/process_notebook.py'
