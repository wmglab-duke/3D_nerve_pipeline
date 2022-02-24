#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import os
import sys
sys.path.append(r'D:\ASCENT\ascent')
os.chdir(r'D:\ASCENT/ascent')

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import numpy as np

import matplotlib.pyplot as plt
from src.core.query import Query
import pandas as pd
import seaborn as sb


linewidth = 0

# plot_type = 'swarm'
plot_type = 'strip'

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)


dats = []

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': [279],
        'model': [0],
        'sim': [33]
    }
}).run()

# builds heatmaps
# q.barcharts_compare_models(logscale=False,
#                            model_labels=['Model 0: Veltink Epineurium, \n              Veltink Perineurium',
#                                          'Model 1: Veltink Epineurium, \n              Goodall Perineurium',
#                                          'Model 2: Goodall Epineurium, \n              Veltink Perineurium',
#                                          'Model 3: Goodall Epineurium, \n              Goodall Perineurium']
#                            )
dats.append(q.threshdat3d(meanify=False))
dats[0] = dats[0].rename(columns = {'mean':'threshold'})
data = pd.concat(dats)
data = data.rename(columns = {'threshold':'Activation Threshold (mA)'})
data.model = data.model.astype(str)
sb.set_theme(style="whitegrid")
g = sb.catplot(x="model", y='Activation Threshold (mA)',hue='nsim',
                data=data, kind=plot_type, height = 5, aspect=.4,linewidth=linewidth,palette = "colorblind",
                sharey=False)
plt.yscale('log')
axs = g.axes

# plt.gcf().suptitle('Activation thresholds for sample {}'.format(pd.unique(data['sample'])[0]))

# retitle = [u'fiber diam: {}\u03bcm'.format(d) for d in [2,5,8,11,13]]

# for ax,title in zip(axs[0],retitle):
#     ax.set_title(title)

plt.subplots_adjust(top = .88,right=.85)
    
plt.suptitle('Activation thresholds for mri279')

# plt.ylim((0.00534827989266171, 11907.037705670795))
plt.savefig('out/analysis/thresholdmri.png',dpi=300)


