import pickle
import mne
import numpy as np

with open('.drop_bad_channels.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)


# with open('.drop_bad_channels.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)


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


def get_low_correlation_electrodes(data, max_corr_sd = 2, min_batch_corr = 0.9, batch_size = 16):
    bad_electrodes = []
    for i in range(0, int(data.shape[0] / batch_size)):
        correlation_matrix = get_correlation(data, batch_no = i, batch_size = batch_size)
        ave_corr = np.mean(correlation_matrix, axis = 0)
        scaled_ave_corr = np.abs((ave_corr - np.mean(ave_corr)) / np.std(ave_corr))
        cond1 = scaled_ave_corr >= max_corr_sd
        cond2 = ave_corr < min_batch_corr
        batch_start = i * batch_size
        batch_end = (i + 1) * batch_size
        electrodes = np.array(range(batch_start, batch_end))[cond1 & cond2]
        bad_electrodes = np.append(bad_electrodes, electrodes)
    return np.int_(bad_electrodes)


with open(snakemake.input['filtered'], 'rb') as file:
    raw = pickle.load(file)

config = snakemake.config

max_corr_sd = config['filter']['maxCorrSD']
min_batch_corr = config['filter']['minBatchCorr']
batch_size = config['filter']['batchSize']

print(f'Looking for electrodes with max correlation SD = {max_corr_sd} and min correlation = {min_batch_corr}')
print(f'Batch size: {batch_size}')

drop_channels = get_low_correlation_electrodes(
    raw.get_data(),
    max_corr_sd = max_corr_sd,
    min_batch_corr = min_batch_corr,
    batch_size = batch_size
)

print(f'Detected {len(drop_channels)} electrodes with low correlation: {drop_channels}')

if config['filter']['removeBadElectrodes']:
    print('Removing bad electrodes')
    raw.drop_channels([f'X{int(x)}' for x in drop_channels])
else:
    print('Bad electrodes were not removed')

with open(snakemake.output['good_channels'], 'wb') as file:
    pickle.dump(raw, file)
