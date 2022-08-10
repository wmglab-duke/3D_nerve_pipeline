#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

import os
import sys

import pandas as pd

sys.path.append(os.path.sep.join([os.getcwd(), '']))
os.chdir('D:/ascent/ascent')

import matplotlib.pyplot as plt
import numpy as np

from src.core.query import Query

sys.path.append(os.getcwd() + '/subrepos/nd_line')
from nd_line.nd_line import nd_line

os.environ['R_HOME'] = r'C:\Users\dpm42\Anaconda3\envs\ascent\lib\R'
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr

pandas2ri.activate()
import rpy2.robjects.lib.ggplot2 as ggplot2

rprint = robjects.globalenv.find("print")
ggpubr = importr("ggpubr")
import seaborn as sb

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14 * 2]) / 2)

samp3d = 653

samples = [652]

models = [0]

sims = [33]

nsim_count = 5

sample_name = '6L'

cuff_contacts = [21, 29]

bigcomp = {0: 'anodic leading', 1: 'cathodic leading'}


q = Query(
    {'partial_matches': False, 'include_downstream': True, 'indices': {'sample': samples, 'model': models, 'sim': sims}}
).run()


def actual_zpos(dat3d, samp3d, model, sim):
    fiberdir = os.path.join('samples', str(samp3d), 'models', str(model), 'sims', str(sim), '3D_fiberset')
    fibers3d = [x for x in os.listdir(fiberdir) if x.endswith('.dat')]
    for file in fibers3d:
        f_ind = int(os.path.splitext(file)[0])
        coord = np.loadtxt(os.path.join(fiberdir, file), skiprows=1)
        fiberline = nd_line(coord)
        datnew = dat3d[(dat3d["model"] == model) & (dat3d["sim"] == sim) & (dat3d["index"] == f_ind)]
        for index, row in datnew.iterrows():
            actual_zpos = fiberline.interp(row['long_ap_pos'])[2]
            dat3d.loc[index, 'activation_zpos'] = actual_zpos
    return dat3d


def datamatch(dest, dat3d, importval):
    dest[importval + '3d'] = np.nan
    for i in range(len(dest)):
        row = dest.iloc[i, :]
        val = dat3d[
            (dat3d["model"] == row['model'])
            & (dat3d["sim"] == row['sim'])
            & (dat3d["nsim"] == row['nsim'])
            & (dat3d["index"] == row['index'])
        ][importval]
        val = list(val)
        if len(val) != 1:
            sys.exit('issue here')
        dest.iloc[i, -1] = val[0]
    if np.any(dest[importval] == np.nan):
        sys.exit('issue here too')
    return dest


def rename_var(df, di):
    for variable, values in di.items():
        for old, new in values.items():
            df = df.replace(to_replace={variable: old}, value=new)
    return df


#%%
dats = []
for sample in samples:
    for n in range(nsim_count):
        ap_path = os.path.join(
            'samples',
            str(sample),
            'models',
            str(models[0]),
            'sims',
            str(sims[0]),
            'n_sims',
            str(n),
            'data',
            'outputs',
            'AP_info.csv',
        )
        dats.append(pd.read_csv(ap_path))
dat2d = pd.concat(dats)
dats3 = []
for n in range(nsim_count):
    ap_path = os.path.join(
        'samples',
        str(samp3d),
        'models',
        str(models[0]),
        'sims',
        str(sims[0]),
        'n_sims',
        str(n),
        'data',
        'outputs',
        'AP_info.csv',
    )
    dats3.append(pd.read_csv(ap_path))
dat3d = pd.concat(dats3)
dat2d = dat2d.reset_index()
dat3d = dat3d.reset_index()

#%%
dat3z = actual_zpos(dat3d, samp3d, models[0], sims[0])

# datfinal = datamatch(dat2d,dat3z,'activation_zpos')

datsb = pd.concat([dat3z, dat2d])

#%%
redict = {
    # "nsim":{
    #     0:'(0) fiber diameter: 2\u03BCm',
    #     1:'(1) fiber diameter: 5\u03BCm',
    #     2:'(2) fiber diameter: 8\u03BCm',
    #     3:'(3) fiber diameter: 11\u03BCm',
    #     4:'(4) fiber diameter: 13\u03BCm'
    #     },
    "sample": {samples[0]: '2D', samp3d: '3D'}
}
datre = rename_var(datsb, redict)
datre.loc[:, 'activation_zpos'] = datre.loc[:, 'activation_zpos'] / 1000
#%%

# plotdata = dat2d[dat2d['sample']==sample]
# sb.set_theme(style="whitegrid")
g = sb.catplot(
    x="sample",
    y='activation_zpos',
    hue='sample',
    col="nsim",
    data=datre,
    kind='strip',
    height=5,
    aspect=0.4,
    linewidth=0,
    order=['2D', '3D'],
    sharey=True,
)
axs = g.axes
axs[0][0].set_ylabel(u'Activation z-position (mm)')
for i, s in enumerate([2, 5, 8, 11, 13]):
    axs[0][i].set_title(f'fiber diam: {s}μm')
plt.subplots_adjust(top=0.88)
plt.suptitle(f"Activation z-position for sample {sample_name} (2D model: slice at {bigcomp[1]} contact)", fontsize=15)
for ax in axs[0]:
    ln = ax.hlines([cuff_contacts], 0.4, 0.6, color='red')
axs[0][-1].legend([ln], ['Cuff Contacts'], loc='center right')
g.savefig(r'out/analysis/{}_zpos.png'.format(samp3d), dpi=400)
