import pickle
import mne
import numpy as np
from scipy import signal

with open('.drop_bad_channels.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)


def get_correlation(data, batch_no = 0, batch_size = 256):
    batch_start = batch_no * batch_size
    batch_end = (batch_no + 1) * batch_size
    correlation_matrix = np.ones((batch_size, batch_size))
    for i in range(batch_start, batch_end - 1):
        for j in range(i + 1, batch_end):
            corr = np.corrcoef(data[i, :], data[j, :])[0, 1]
            correlation_matrix[i % batch_size, j % batch_size] = corr
            correlation_matrix[j % batch_size, i % batch_size] = corr
    return correlation_matrix


def get_low_correlation_electrodes(data, min_batch_corr = 0.9, batch_size = 16):
    bad_electrodes = []
    for i in range(0, int(data.shape[0] / batch_size)):
        correlation_matrix = get_correlation(data, batch_no = i, batch_size = batch_size)
        ave_corr = np.mean(correlation_matrix, axis = 0)
        batch_start = i * batch_size
        batch_end = (i + 1) * batch_size
        electrodes = np.array(range(batch_start, batch_end))[ave_corr <= min_batch_corr]
        bad_electrodes = np.append(bad_electrodes, electrodes)
    return bad_electrodes


config = snakemake.config

with open(snakemake.input['filtered'], 'rb') as file:
    raw = pickle.load(file)

print(f'Searching for electrodes with correlation < {config["filter"]["minBatchCorrelation"]}')
print(f'Batch size: {config["filter"]["batchSize"]}')

data = raw.get_data()
bad_electrodes = get_low_correlation_electrodes(
    data,
    min_batch_corr = config['filter']['minBatchCorrelation'],
    batch_size = config['filter']['batchSize']
)

print(f'Detected {len(bad_electrodes)} electrodes with low correlation: {bad_electrodes}')

if config['filter']['removeUncorrelated']:
    print('Removing uncorrelated electrodes')
    raw.drop_channels([f'X{int(x)}' for x in bad_electrodes])
else:
    print('Uncorrelated electrodes were not removed')

with open(snakemake.output['final'], 'wb') as file:
    pickle.dump(raw, file)
