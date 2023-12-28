import pickle
import mne
import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial
import re


def bandpower(data, sf, low, high, min_freq, max_freq, nperseg = None, relative = False):
    from scipy.signal import welch
    from scipy.integrate import simps
    import numpy as np

    if nperseg is None:
        nperseg = np.min([sf, len(data) - 1])

    freqs, psd = welch(data, sf, nperseg = nperseg)
    freq_res = freqs[1] - freqs[0]
    idx_band = np.logical_and(freqs >= low, freqs <= high)
    bp = simps(psd[idx_band], dx = freq_res)

    if relative:
        idx_lim = np.logical_and(freqs >= min_freq, freqs <= max_freq)
        bp = bp / simps(psd[idx_lim], dx = freq_res)

    return bp


def process_event(events, data, freq_range, freq_step, min_freq, max_freq, relative_frequency, nperseg, run, index):
    print(f'Processing event {index}')

    features = []
    event = events.iloc[index]

    event_data = data.copy().crop(tmin = event['Start'], tmax = event['End'])

    feature_type = 'Relative Band Power' if relative_frequency else 'Absolute Band Power'

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
                nperseg = nperseg,
                relative = relative_frequency
            )

            features.append({
                'Channel': channel,
                'Frequency': freq,
                'Type': feature_type,
                'Feature': bp,
                'Event': event['Event'],
                'Start': event['Start'],
                'End': event['End'],
                'OriginalEventID': event['OriginalEventID'],
                'FeatureID': f'{channel}_{freq}',
                'Run': run,
                'EventID': f'{event["EventID"]}_{run}'
            })

    return pd.DataFrame(features)


def extract_feature_channel(row):
    return re.sub('X([0-9]+)_[0-9]+', '\\1', row['FeatureID'])


def extract_feature_frequency(row):
    return re.sub('X[0-9]+_([0-9]+)', '\\1', row['FeatureID'])


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
    nperseg = snakemake.params['nperseg']

    # Read input data
    data = read_file(snakemake.input['final'])

    # Read event data and drop the first two events because they are empty
    events = read_file(snakemake.input['events'])
    events = events.drop(index = [0, 1]).reset_index(drop = True)
    events = normalize_events(events, feature_length = feature_length)

    # Define the frequency range
    freq_range = np.arange(min_freq, max_freq, freq_step)

    print(f'Processing {len(events)} events')

    # Process the events parallely to preserve my very limited patience
    threads = snakemake.threads
    with Pool(threads) as pool:
        features = pool.map(partial(process_event,
                                    events, data, freq_range, freq_step, min_freq, max_freq,
                                    relative_frequency, nperseg, snakemake.wildcards['run']),
                            range(len(events)))

    features = pd.concat(features)
    features = features.sort_values('Start')

    #######################
    # Validate that everything is okay
    #######################
    print('Validating the events')

    events = read_file(snakemake.input['events'])
    key = features[['Event', 'Start', 'End', 'OriginalEventID']].drop_duplicates()
    key = key.reset_index(drop = True)

    # Validate
    for index, feature_event in key.iterrows():
        eventId = feature_event['OriginalEventID']
        event = events[events['EventID'] == eventId].iloc[0]

        assert feature_event['Start'] >= event['Start'], \
            f'Start position not matching in event {eventId}: {feature_event["Start"]} < {event["Start"]}'
        assert feature_event['End'] <= event['End'], \
            f'End position not matching in event {eventId}: {feature_event["End"]} > {event["End"]}'
        assert feature_event['Event'] == event['Event'], \
            f'Event type not matching in event {eventId}: {feature_event["Event"]} != {event["Event"]}'

    print('Validation successful!')

    # Extract feature matrix
    x = pd.pivot_table(features, values = 'Feature', index = 'FeatureID', columns = 'EventID')
    x = x.apply(pd.to_numeric)

    # Extract event details
    event_info = features[['EventID', 'Run', 'Event', 'OriginalEventID', 'Start', 'End']]
    event_info = event_info.drop_duplicates().reset_index(drop = True)

    # Extract feature details
    feature_info = features[['FeatureID', 'Channel', 'Frequency']]
    feature_info = feature_info.drop_duplicates().reset_index(drop = True)
    feature_info['Frequency'] = feature_info['Frequency'].astype('int', copy = False)

    # Make sure everything matches
    assert len(event_info) == x.shape[1]
    assert len(feature_info) == x.shape[0]

    x.to_csv(snakemake.output['feature_matrix'])
    event_info.to_csv(snakemake.output['event_key'], index = False)
    feature_info.to_csv(snakemake.output['feature_key'], index = False)
