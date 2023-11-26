from scipy.io import loadmat
import pandas as pd

with open('.montage_to_csv.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

data = loadmat(snakemake.input['mat'])['pos_256']
data = pd.DataFrame(data)

channels = [f'X{x}' for x in range(0, 256)]
data.insert(0, 'ChannelName', channels)

data.to_csv(snakemake.output['csv'], index = False)
