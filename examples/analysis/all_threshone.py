# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

os.chdir('../../')
threshdat = pd.read_csv('thresh.csv')
sns.set_style('whitegrid')
sns.scatterplot(data=threshdat, x='threshold3d', y='threshold', hue='nsim', palette='colorblind')
# plot one one line out to max of 3d
plt.plot([0, threshdat.threshold.max()], [0, threshdat.threshold.max()], 'r', linewidth=2)
plt.ylabel('2D Threshold')
plt.xlabel('3D Threshold')
plt.xscale('log')
plt.yscale('log')
# make axes squre
plt.gca().set_aspect('equal', adjustable='box')
# calculate correlation between 2d and 3d and percentage 2d greater than 3d, add to title
r, p = pearsonr(threshdat.threshold, threshdat.threshold3d)
perc = sum(threshdat.threshold > threshdat.threshold3d) / len(threshdat.threshold)
plt.title(
    '2D vs 3D Thresholds for all samples and fiber diameters\n'
    f'Correlation between 2D and 3D thresholds: r={r:.2f}, p={p:.2f}\n'
    f'percentage of 2D thresholds greater than 3D couterpart: {perc:.2f}'
)
# calculate percentage of 2d threshold greater than 3d
plt.show()
# generate boxplot of 2d and 3d thresholds
sns.barplot(data=threshdat[['threshold', 'threshold3d']])
plt.ylabel('Threshold')
plt.xlabel('2D vs 3D')
plt.title('2D vs 3D Thresholds for all samples and fiber diameters')
plt.show()
