import pickle
import mne
import numpy as np

with open('.drop_bad_channels2.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.drop_bad_channels2.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

with open(snakemake.input['filtered'], 'rb') as file:
    raw = pickle.load(file)

config = snakemake.config

powers = raw.compute_psd().get_data()
corr = np.corrcoef(powers)
ave_corr = np.mean(corr, axis = 0)

max_n_channels = 25
min_corr = 1
bad_channels = []

for i in range(999, 800, -1):
    min_corr = i / 1000
    bad_channels = np.where(ave_corr < min_corr)[0]

    if len(bad_channels) < max_n_channels:
        break

print(f'Detecting channels with average power correlation less than {min_corr}')
print(f'Detected {len(bad_channels)} channels: {bad_channels}')

if config['filter']['dropLowQuality']:
    print('Removing low quality channels')
    raw.drop_channels([f'X{int(x)}' for x in bad_channels])
else:
    print('Low quality channels were not removed')

with open(snakemake.output['good_channels'], 'wb') as file:
    pickle.dump(raw, file)
