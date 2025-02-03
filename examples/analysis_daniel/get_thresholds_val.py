#!/usr/bin/env python3
import json
import os
import sys

import pandas as pd

sys.path.append('.')

os.chdir('../..')
from src.core.query import Query

model = 0
source_sim = 3
alldata = []
for sample_ind in [25399, 253999]:
    simint = 333
    q3 = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [sample_ind], 'model': [model], 'sim': [simint]},
        }
    ).run()
    dat3d = q3.data(
        source_sample=25299,
        tortuosity=False,
        peri_site=False,
        zpos=False,
        source_sim=source_sim,
        # label=nerve_label,
        efib_distance=False,
        oneten=False,
        thresh_only=True,
    )
    alldata.append(dat3d)
df = pd.concat(alldata)
import matplotlib.pyplot as plt
import numpy as np

# %% plot
import seaborn as sns

sns.set(style="whitegrid", context='paper', rc={'figure.dpi': 400})
merged_df = pd.merge(
    df.query("sample == 25399"),
    df.query("sample == 253999"),
    on=["master_fiber_index", "nsim"],
    suffixes=("_coarse", "_fine"),
    how="left",
)

merged_df["percdiff"] = (
    np.abs(100 * (merged_df["threshold_coarse"] - merged_df["threshold_fine"])) / merged_df["threshold_coarse"]
)
sns.violinplot(data=merged_df, y='percdiff', x='nsim', palette='RdPu')
sns.stripplot(data=merged_df, y='percdiff', x='nsim', color='black', jitter=False)

plt.ylabel('Percent difference in fiber threshold (%)')
plt.xlabel('Fiber diameter (μm)')
plt.xticks(range(6), [3, 5, 7, 9, 11, 13])
plt.gcf().set_size_inches(4, 3)
print(merged_df['percdiff'].max())
print(merged_df['percdiff'].median())
plt.axhline(1, color='blue')
plt.annotate('threshold search tolerance', xy=(2.4, 0.94), color='blue')

# plot the raw thresholds vs each other
plt.figure()
sns.scatterplot(data=merged_df, x='threshold_coarse', y='threshold_fine', hue='nsim', palette='RdPu')
plt.xlabel('Coarse threshold (μA)')
plt.ylabel('Fine threshold (μA)')
plt.gcf().set_size_inches(4, 3)
# unity line
plt.plot([0, 4], [0, 4], color='black', linestyle='--')
plt.legend(title='Fiber diameter (μm)')
plt.show()


merged_df.iloc[np.argmax(merged_df['percdiff'])]
