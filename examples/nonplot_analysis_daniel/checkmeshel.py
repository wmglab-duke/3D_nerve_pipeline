import os

os.chdir(r'D:\threed_final')
import gc
import json
import os

import numpy as np
import pandas as pd

alldata = []

for folder in os.listdir('samples'):
    if folder.startswith('.'):
        continue
    # skip if third digit is 4
    if len(folder) < 3 or folder[2] in ['4', '3']:
        continue
    try:
        folderpath = os.path.join(os.getcwd(), 'samples', folder, 'models', '0', 'model.json')
        # json load
        with open(folderpath) as f:
            data = json.load(f)
        mesh_elements = data.get('mesh').get('stats').get('number_elements')
        alldata.append({'sample': folder, 'mesh_elements': mesh_elements})
        gc.collect()
    except Exception as e:
        print(folder, e)
        alldata.append({'sample': folder, 'mesh_elements': None})
# %%
# mean max and median node counts
import pandas as pd

df = pd.DataFrame(alldata)
df.to_csv('mesh_elements.csv')
print(df['mesh_elements'].agg(['mean', 'max', 'median', 'min']))
print(df['mesh_elements'].mean())
