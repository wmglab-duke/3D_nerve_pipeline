import os

os.chdir('../..')
# import garbagecollector
import gc
import os

import numpy as np
import pandas as pd

alldata = []

for folder in os.listdir('samples'):
    print(folder)
    try:
        for n_sim in [0, 5]:
            folderpath = os.path.join('samples', folder, f'models/0/sims/333/fibersets/{n_sim}')
            for file in [f for f in os.listdir(folderpath) if f.endswith('.dat')]:
                # read the first line of the file
                with open(os.path.join(folderpath, file)) as f:
                    first_line = f.readline()
                # turn into int
                n_sections = int(first_line.split()[0])
                n_nodes = (n_sections - 1) / 11 + 1
                # make sure the number of nodes is an integer
                assert n_nodes.is_integer()
                alldata.append({'folder': folder, 'n_sim': n_sim, 'file': file, 'n_nodes': n_nodes})
            gc.collect()
    except Exception as e:
        print(folder, e)
        continue

# %%
# mean max and median node counts
import pandas as pd

df = pd.DataFrame(alldata)
print(df.groupby('n_sim')['n_nodes'].agg(['mean', 'max', 'median', 'min']))
# if the third digit of the folder is a 3, it is a 3D sample
df['is_3D'] = df['folder'].apply(lambda x: x[2] == '3')
print(df.groupby(['is_3D', 'n_sim'])['n_nodes'].agg(['mean', 'max', 'median', 'min']))
