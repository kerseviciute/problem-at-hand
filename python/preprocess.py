import pickle
import mne
import numpy as np
from scipy import signal

with open('.preprocess.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)


# with open('.preprocess.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

def scale_channels(mne_data):
    data = np.array(mne_data.get_data()).T
    data = data - np.mean(data, axis = 0)
    mne_data._data = data.T


def apply_filter(data, filter):
    filtered_channels = []
    for channel_i in range(data.shape[0]):
        to_filter = data[channel_i, :]
        filtered = signal.sosfiltfilt(filter, to_filter)
        filtered_channels.append(filtered)
    return np.array(filtered_channels)


def remove_powerline(mne_data, lfreq, hfreq):
    data = mne_data.get_data()
    sample_rate = mne_data.info['sfreq']
    filter = signal.butter(
        N = 40,
        Wn = [lfreq, hfreq],
        btype = 'bandstop',
        output = 'sos',
        fs = sample_rate
    )
    data = apply_filter(data, filter)
    mne_data._data = data


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

with open(snakemake.input['raw'], 'rb') as file:
    raw = pickle.load(file)

print('Scaling the channels by subtracting the mean')
scale_channels(raw)

print('Filtering')
print('Remove power line at 60 Hz')
remove_powerline(raw, lfreq = 50, hfreq = 70)

print('Remove power line at 120 Hz')
remove_powerline(raw, lfreq = 115, hfreq = 125)

print('Remove power line at 180 Hz')
remove_powerline(raw, lfreq = 170, hfreq = 190)

print('Remove power line at 240 Hz')
remove_powerline(raw, lfreq = 235, hfreq = 245)

print('Remove frequencies lower than 4 Hz and higher than 290 Hz')
raw.filter(l_freq = 4, h_freq = 290)

print(f'Searching for electrodes with correlation < {config["filter"]["minBatchCorrelation"]} in batches of {config["filter"]["batchSize"]}')

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


with open(snakemake.output['filtered'], 'wb') as file:
    pickle.dump(raw, file)

