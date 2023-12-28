import pickle
import pandas as pd
import re
import numpy as np

with open('.combine_features.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)


# with open('.combine_features.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

def feature_id(row): return f'{row["Channel"]}_{row["Frequency"]}'

def event_id(row): return f'{row["EventID"]}_{row["Run"]}'


print('Reading the feature files')

# Read features of all the runs from the same person
features = pd.DataFrame()
for file in snakemake.input:
    data = pd.read_csv(snakemake.input[0])
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
