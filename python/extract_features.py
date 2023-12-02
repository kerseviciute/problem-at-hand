import pickle
import mne
import numpy as np
import pandas as pd

with open('.extract_average_bandpower.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.extract_average_bandpower.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

##########################################################################################

# Load feature extraction functions
with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

# Read config file
config = snakemake.config
feature_length = config['features']['length']
min_freq = config['features']['min_freq']
max_freq = config['features']['max_freq']
freq_step = config['features']['freq_step']

# Read input data
data = read_file(snakemake.input['final'])

# Read event data and drop the first two events because they are empty
events = read_file(snakemake.input['events'])
events = events.drop(index = [0, 1]).reset_index(drop = True)
events = normalize_events(events, feature_length = feature_length)

# Define the frequency range
freq_range = np.arange(min_freq, max_freq, freq_step)

features = []

# Loop through the events
# TODO: PARALLELIZE!!!!!!!!!!!!!
for index, event in events.iterrows():
    print(f'Processing event {index}')
    event_data = data.copy().crop(tmin = event['Start'], tmax = event['End'])

    # Loop though the channels
    for channel in event_data.ch_names:

        # Loop through the frequency ranges
        # TODO: OPTIMIZE!!!!!
        for freq in freq_range:
            bp = bandpower(
                data = event_data.get_data(picks = channel)[0],
                sf = event_data.info['sfreq'],
                low = freq,
                high = freq + freq_step,
                min_freq = min_freq,
                max_freq = max_freq,
                nperseg = 1200,
                relative = True
            )

            features.append({
                'Channel': channel,
                'Frequency': freq,
                'Type': 'Relative Band Power',
                'Feature': bp,
                'Event': event['Event']
            })

features = pd.DataFrame(features)

features.to_csv(snakemake.output['features'], index = False)
