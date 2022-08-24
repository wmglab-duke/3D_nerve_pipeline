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
sim = 3
model = 0
with open('examples/analysis/plotconfig.json') as f:
    config = json.load(f)
for sample_data in config['sample_data']:
    samp3d = sample_data['index3d']
    nerve_label = sample_data['name']
    samples2d = [x['index'] for x in sample_data['exsamples']]
    #%%
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': samples2d, 'model': [model], 'sim': [sim]},
        }
    ).run()

    dat2d = q.data()

    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [samp3d], 'model': [model], 'sim': [sim]},
        }
    ).run()

    dat3d = q.data(source_sample=samples2d[0])

    dat2d = datamatch(dat2d, dat3d, 'threshold')
    #%%

    sample_labels = ['rostral contact', 'caudal contact']

    import seaborn as sns
    from scipy.stats import pearsonr

    sns.set_theme()
    sns.set(font_scale=1.5)

    dat2d = dat2d.rename(columns={'sample': 'Slice'})
    #%%
    g = sns.lmplot(
        data=dat2d, x="threshold3d", y="threshold", hue="Slice", height=5, col='nsim', sharey=False, sharex=False
    )
    axs = g.axes.ravel()
    axs[0].set_ylabel('2D threshold (mA)')
    plt.suptitle(f'Activation threshold correlation for sample {nerve_label}', fontsize=25)
    plt.subplots_adjust(top=0.85, right=0.93)
    new_labels = ['Anodic\nLeading', 'Cathodic\nLeading']
    for t, l in zip(g._legend.texts, new_labels):
        t.set_text(l)
    for i, ax in enumerate(g.axes.ravel()):
        # ax.set_title(f'fiber diam: {s}μm')
        idat = dat2d.query(f'nsim=={i}')
        min_thresh = min([min(idat.threshold), min(idat.threshold3d)])
        max_thresh = max([max(idat.threshold), max(idat.threshold3d)])
        limits = (min_thresh, max_thresh)
        corr = {}
        for sample in samples2d:
            thisdat = dat2d[(dat2d["nsim"] == i) & (dat2d["Slice"] == sample)]
            corr[sample] = round(pearsonr(thisdat['threshold'], thisdat['threshold3d'])[0], 3)
        ax.legend(labels=["r=" + str(corr[sample]) for sample in samples2d])
        ax.set_xlabel('3D threshold (mA)')
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.plot(limits, limits, color='red')
    g.savefig(f'out/analysis/threscorr_{nerve_label}', dpi=400)
