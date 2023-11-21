import mne
from scipy.io import loadmat
import pickle
import pandas as pd
import numpy as np


def extract_event(events, event):
    present_event = events == event
    change_indices = np.where(np.diff(present_event))[ 0 ]

    event_data = pd.DataFrame({
        'EventStart': change_indices[ :-1 ] + 1,
        'EventEnd': change_indices[ 1: ],
        'Event': events[ change_indices[ :-1 ] + 1 ]
    })

    event_data = event_data[ event_data[ 'Event' ] == event ]

    return event_data


montage = mne.channels.read_custom_montage(snakemake.input[ 'montage' ])

data = loadmat(snakemake.input[ 'raw' ])[ 'y' ]

timepoints = data[ 0, : ]
events = data[ -1, : ]

data = data[ 1:-1, ]

channels = [ f'X{x}' for x in range(0, data.shape[ 0 ]) ]

sfreq = snakemake.config[ 'sampleRate' ]

ch_types = [ 'eeg' for x in range(0, data.shape[ 0 ]) ]

info = mne.create_info(ch_names = channels, sfreq = sfreq, ch_types = ch_types)

raw = mne.io.RawArray(data, info)
raw.set_montage(montage)

with open(snakemake.output[ 'raw' ], 'wb') as file:
    pickle.dump(raw, file)

extracted_events = pd.concat([
    extract_event(events, 0),
    extract_event(events, 1),
    extract_event(events, 2),
    extract_event(events, 3),
    extract_event(events, 4),
    extract_event(events, 5)
], ignore_index = True)

with open(snakemake.output[ 'events' ], 'wb') as file:
    pickle.dump(extracted_events, file)
