#!/usr/bin/env python3.7

"""The copyrights of this software are owned by Duke University.

Please refer to the LICENSE.txt and README.txt files for licensing
instructions. The source code can be found on the following GitHub
repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

import json
import sys

import matplotlib.pyplot as plt

sys.path.append('.')

from src.core.plotter import datamatch
from src.core.query import Query

#%%
model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)


def plot_correlation():
    global ax
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat2d = q.data()
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [simdex]},
        }
    ).run()
    dat3d = q.data(source_sample=samples2d[0])
    dat2d = datamatch(dat2d, dat3d, 'threshold')
    # %%
    sample_labels = ['rostral contact', 'caudal contact']
    import seaborn as sns
    from scipy.stats import pearsonr

    sns.set_theme()
    sns.set(font_scale=1.5)
    dat2d = dat2d.rename(columns={'sample': 'Slice'})
    # %%
    g = sns.lmplot(
        data=dat2d, x="threshold", y="threshold3d", hue="Slice", height=5, col='nsim', sharey=False, sharex=False
    )
    axs = g.axes.ravel()
    axs[0].set_ylabel('3D threshold (mA)')
    plt.suptitle(f'Activation threshold correlation for sample {nerve_label}', fontsize=25)
    plt.subplots_adjust(top=0.85, right=0.93)
    new_labels = ['Anodic\nLeading', 'Cathodic\nLeading']
    for t, l in zip(g._legend.texts, new_labels):
        t.set_text(l)
    for i, ax in g.axes.ravel():
        # ax.set_title(f'fiber diam: {s}μm')
        corr = {}
        for sample in samples2d:
            thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["Slice"] == sample)]
            corr[sample] = round(pearsonr(thisdat['threshold'], thisdat['threshold3d'])[0], 3)
        ax.legend(labels=["r=" + str(corr[sample]) for sample in samples2d])
        ax.set_xlabel('2D threshold (mA)')
    g.savefig(f'out/analysis/{simdex}/threscorr_{nerve_label}', dpi=400)


for simdex in config['sim_data'].keys():
    for sample_data in config['sample_data']:
        samp3d = sample_data['index3d']
        nerve_label = sample_data['name']
        samples2d = [x['index'] for x in sample_data['exsamples']]
        plot_correlation()
