def read_file(filename):
    import pickle

    with open(filename, 'rb') as file:
        raw = pickle.load(file)
    return raw


def center_event_time(event, event_length = 3.0):
    event = event.copy()
    event_sub = (event['End'] - event['Start'] - event_length) / 2
    event['Start'] = event['Start'] + event_sub
    event['End'] = event['End'] - event_sub
    event['Length'] = event['End'] - event['Start']
    assert round(event['Length'], 2) == round(event_length, 2)
    return event


def normalize_events(events, feature_length):
    import pandas as pd

    mini_events = []

    for index, event in events.iterrows():
        n_events = int(event['Length'] / feature_length)
        centered_event = center_event_time(event, n_events * feature_length)

        for i in range(n_events):
            start = centered_event['Start'] + i * feature_length
            end = start + feature_length

            mini_events.append({
                'Start': start,
                'End': end,
                'Event': centered_event['Event'],
                'Length': feature_length
            })

    mini_events = pd.DataFrame(mini_events)
    mini_events = mini_events.sort_values(by = 'Start')
    mini_events = mini_events.reset_index(drop = True)
    mini_events['Event'] = mini_events['Event'].astype('category', copy = False)
    mini_events['EventID'] = [f'Event{x}' for x in range(0, len(mini_events))]

    return mini_events
