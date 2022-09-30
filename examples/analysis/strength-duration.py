import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

os.chdir('../../')
threshdat = pd.read_csv('threshall.csv').query('sim==7')
sns.set_style('whitegrid')
# assign nsims as follow: 0 through 2 have fiber diameter of 3
# 3 through 5 have fiber diameter of 13
# 0 and 3 have pulse width of .1
# 1 and 4 have pulse width of .25
# 3 and 6 have pulse width of .4
threshdat['nsim'] = threshdat['nsim'].astype(int)
threshdat['diameter'] = np.nan
threshdat.loc[threshdat['nsim'] < 3, 'diameter'] = 3
threshdat.loc[threshdat['nsim'] > 2, 'diameter'] = 13
threshdat['pulsewidth'] = np.nan
threshdat.loc[threshdat['nsim'] % 3 == 0, 'pulsewidth'] = 0.1
threshdat.loc[threshdat['nsim'] % 3 == 1, 'pulsewidth'] = 0.25
threshdat.loc[threshdat['nsim'] % 3 == 2, 'pulsewidth'] = 0.4
threshdat.rename(columns={'threshold': '2D', 'threshold3d': '3D'}, inplace=True)
threshdat = threshdat.melt(
    id_vars=['sample', 'nsim', 'diameter', 'pulsewidth'],
    value_vars=['2D', '3D'],
    var_name='type',
    value_name='threshold',
)
# seaborn facet scatterplot with x as pulse width, y as threshold, and column as diameter
sns.set(font_scale=1.5)
g = sns.catplot(
    kind='violin',
    data=threshdat,
    col="diameter",
    col_wrap=2,
    sharey=False,
    x='pulsewidth',
    y='threshold',
    palette='colorblind',
    hue='type',
    dodge=True,
    linewidth=1,
)
plt.subplots_adjust(top=0.85)
plt.suptitle('Threshold vs. Pulse Width')
