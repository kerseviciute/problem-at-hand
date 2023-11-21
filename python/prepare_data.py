import mne
from scipy.io import loadmat
import pickle

montage = mne.channels.read_custom_montage(snakemake.input['montage'])

data = loadmat(snakemake.input['raw'])['y']
data = data[ 1:-1, ]

channels = [ f'X{x}' for x in range(0, data.shape[0]) ]

sfreq = snakemake.config['sampleRate']

ch_types = [ 'eeg' for x in range(0, data.shape[0]) ]

info = mne.create_info(ch_names = channels, sfreq = sfreq, ch_types = ch_types)

raw = mne.io.RawArray(data, info)
raw.set_montage(montage)

with open(snakemake.output['pickle'], 'wb') as f:
    pickle.dump(raw, f)
