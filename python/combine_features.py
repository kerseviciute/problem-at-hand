import pandas as pd
import re
import numpy as np

pd.to_pickle(snakemake, '.combine_features.py.pkl')
# snakemake = pd.read_pickle('.combine_features.py.pkl')

print('Reading the feature files')

# Read features of all the runs from the same person
features = pd.DataFrame()
for file in snakemake.input:
    print(f'Reading file {file}')
    data = pd.read_csv(file)
    features = pd.concat([features, data], ignore_index = True)

# Make sure run ids were extracted correctly and we have 10 runs
assert len(set(features['Run'])) == 10

# Extract feature matrix
x = pd.pivot_table(features, values = 'Feature', index = 'FeatureID', columns = 'EventID')
x = x.apply(pd.to_numeric)

# Extract event details
key = features[['EventID', 'Run', 'Event']].drop_duplicates()
key = key.reset_index(drop = True)

# Make sure everything matches
assert len(key) == x.shape[1]

x.to_csv(snakemake.output['feature_matrix'])
key.to_csv(snakemake.output['feature_key'], index = False)
