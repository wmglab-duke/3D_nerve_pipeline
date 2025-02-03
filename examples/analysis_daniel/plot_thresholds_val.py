# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

os.chdir('../..')

threshload = pd.read_csv(f"thresh_unmatched_sim333_val.csv")

# %% Plotting
sns.set(style="whitegrid", context='paper')
thisdata = threshload.query('sample in [253000, 252000]')
merged_df = pd.merge(
    thisdata.query("sample == 253000"),
    thisdata.query("sample == 252000"),
    on=["master_fiber_index", "nsim"],
    suffixes=("_3d", "_2d"),
    how="left",
)  # TODO remove this how
merged_df["percdiff"] = (
    np.abs(100 * (merged_df["threshold_3d"] - merged_df["threshold_2d"])) / merged_df["threshold_3d"]
)
sns.violinplot(data=merged_df, y='percdiff', x='nsim', palette='RdPu')
sns.stripplot(data=merged_df, y='percdiff', x='nsim', color='black', jitter=False)

plt.ylabel('Percent difference in fiber threshold (%)')
plt.xlabel('Fiber diameter (μm)')
plt.xticks(range(6), [3, 5, 7, 9, 11, 13])
plt.gcf().set_size_inches(4, 3)
print(merged_df['percdiff'].max())
print(merged_df['percdiff'].median())

# %% now plot whether the mean for each fiber diam and nerve is different
# first take mean thresholds
plt.figure()
mean_thresh = merged_df.groupby(['nsim'])[['threshold_3d', 'threshold_2d']].mean().reset_index()
# now percent difference
mean_thresh['percdiff'] = (
    100 * (mean_thresh['threshold_3d'] - mean_thresh['threshold_2d']) / mean_thresh['threshold_3d']
)
# plot
sns.stripplot(data=mean_thresh, y='percdiff', palette='RdPu', jitter=True, x='nsim')
plt.ylabel('Difference in mean fiber threshold (%)')
plt.xlabel('Fiber diameter (μm)')
plt.xticks(range(6), [3, 5, 7, 9, 11, 13])
plt.gcf().set_size_inches(4, 3)
plt.ylim(0, 3)
plt.axhline(1, label='Search tolerance', color='black', linestyle='--')
plt.legend()

# %%
