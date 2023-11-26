import mne
from scipy.io import loadmat
import pickle
import pandas as pd
import numpy as np

with open('.prepare_data.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)


def extract_event(events, event, sample_rate):
    present_event = events == event
    change_indices = np.where(np.diff(present_event))[0]
    change_indices = np.insert(change_indices, 0, -1)

    event_data = pd.DataFrame({
        'EventStart': (change_indices[:-1] + 1) / sample_rate,
        'EventEnd': (change_indices[1:] + 1) / sample_rate,
        'Event': events[change_indices[:-1] + 1]
    })

    event_data = event_data[event_data['Event'] == event]

    return event_data


montage = mne.channels.read_custom_montage(snakemake.input['montage'])

data = loadmat(snakemake.input['raw'])['y']

timepoints = data[0, :]
events = data[-1, :]

data = data[1:-1, ]

channels = [f'X{x}' for x in range(0, data.shape[0])]

sfreq = snakemake.config['sampleRate']

ch_types = ['eeg' for x in range(0, data.shape[0])]

info = mne.create_info(ch_names = channels, sfreq = sfreq, ch_types = ch_types)

raw = mne.io.RawArray(data, info)
raw.set_montage(montage)

with open(snakemake.output['raw'], 'wb') as file:
    pickle.dump(raw, file)

extracted_events = pd.concat([
    extract_event(events, 0, sample_rate = sfreq),
    extract_event(events, 1, sample_rate = sfreq),
    extract_event(events, 2, sample_rate = sfreq),
    extract_event(events, 3, sample_rate = sfreq),
    extract_event(events, 4, sample_rate = sfreq),
    extract_event(events, 5, sample_rate = sfreq),
    extract_event(events, 6, sample_rate = sfreq)
], ignore_index = True)

extracted_events = extracted_events.sort_values(by = 'EventStart')
extracted_events = extracted_events.reset_index(drop = True)

with open(snakemake.output['events'], 'wb') as file:
    pickle.dump(extracted_events, file)
