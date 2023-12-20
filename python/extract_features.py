import pickle
import mne
import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial


def bandpower(data, sf, low, high, min_freq, max_freq, nperseg = None, relative = False):
    from scipy.signal import welch
    from scipy.integrate import simps
    import numpy as np

    if nperseg is None:
        nperseg = sf

    freqs, psd = welch(data, sf, nperseg = nperseg)
    freq_res = freqs[1] - freqs[0]
    idx_band = np.logical_and(freqs >= low, freqs <= high)
    bp = simps(psd[idx_band], dx = freq_res)

    if relative:
        idx_lim = np.logical_and(freqs >= min_freq, freqs <= max_freq)
        bp = bp / simps(psd[idx_lim], dx = freq_res)

    return bp


def process_event(events, data, freq_range, freq_step, min_freq, max_freq, relative_frequency, index):
    print(f'Processing event {index}')

    features = []
    event = events.iloc[index]

    event_data = data.copy().crop(tmin = event['Start'], tmax = event['End'])

    # Loop though the channels
    for channel in event_data.ch_names:
        # Loop through the frequency ranges
        for freq in freq_range:
            bp = bandpower(
                data = event_data.get_data(picks = channel)[0],
                sf = event_data.info['sfreq'],
                low = freq,
                high = freq + freq_step,
                min_freq = min_freq,
                max_freq = max_freq,
                nperseg = 300,
                relative = relative_frequency
            )

            features.append({
                'Channel': channel,
                'Frequency': freq,
                'Type': 'Relative Band Power',
                'Feature': bp,
                'Event': event['Event'],
                'EventID': event['EventID']
            })

    return pd.DataFrame(features)


if __name__ == '__main__':
    with open('.extract_features.py.pkl', 'wb') as file:
        pickle.dump(snakemake, file)

    # with open('.extract_features.py.pkl', 'rb') as file:
    #     snakemake = pickle.load(file)

    ##########################################################################################

    # Load feature extraction functions
    with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
        exec(file.read())

    # Read processing parameters
    feature_length = snakemake.params['feature_length']
    min_freq = snakemake.params['min_freq']
    max_freq = snakemake.params['max_freq']
    freq_step = snakemake.params['freq_step']
    relative_frequency = snakemake.params['relative_frequency']

    # Read input data
    data = read_file(snakemake.input['final'])

    # Read event data and drop the first two events because they are empty
    events = read_file(snakemake.input['events'])
    events = events.drop(index = [0, 1]).reset_index(drop = True)
    events = normalize_events(events, feature_length = feature_length)

    # Define the frequency range
    freq_range = np.arange(min_freq, max_freq, freq_step)

    # Process the events parallely to preserve my very limited patience
    threads = snakemake.threads
    with Pool(threads) as pool:
        features = pool.map(partial(process_event, events, data, freq_range, freq_step, min_freq, max_freq, relative_frequency),
                            range(len(events)))

    features = pd.concat(features)

    features.to_csv(snakemake.output['features'], index = False)
