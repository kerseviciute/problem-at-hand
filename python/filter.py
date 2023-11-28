import pickle
import mne
import numpy as np
from scipy import signal

with open('.preprocess.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)


# with open('.preprocess.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)


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


config = snakemake.config

with open(snakemake.input['raw'], 'rb') as file:
    raw = pickle.load(file)

if config['filter']['scale']:
    print('Scaling the channel data')
    data = raw.get_data().transpose()
    data = (data - np.mean(data, axis = 0)) / np.std(data, axis = 0)
    raw._data = data.transpose()
else:
    print('Channels were not scaled')

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

with open(snakemake.output['filtered'], 'wb') as file:
    pickle.dump(raw, file)
