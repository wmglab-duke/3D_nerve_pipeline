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

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.core.query import Query

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


import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import numpy as np

os.chdir('D:/ASCENT/ascent')

import matplotlib.pyplot as plt

from src.core.query import Query

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)

threed = 473

samples = [470, 472]

sample_name = '4R'

models = [0]

sims = [33]

bigcomp = {0: 'anodic leading', 1: 'cathodic leading'}


q = Query(
    {'partial_matches': False, 'include_downstream': True, 'indices': {'sample': samples, 'model': models, 'sim': sims}}
).run()

# builds heatmaps
# q.barcharts_compare_models(logscale=False,
#                            model_labels=['Model 0: Veltink Epineurium, \n              Veltink Perineurium',
#                                          'Model 1: Veltink Epineurium, \n              Goodall Perineurium',
#                                          'Model 2: Goodall Epineurium, \n              Veltink Perineurium',
#                                          'Model 3: Goodall Epineurium, \n              Goodall Perineurium']
#                            )
dat2d = q.threshdat(sl=False, meanify=False)

q = Query(
    {
        'partial_matches': False,
        'include_downstream': True,
        'indices': {'sample': [threed], 'model': models, 'sim': sims},
    }
).run()

dat3d = q.threshdat3d(meanify=False)


def rename_var(df, di):
    for variable, values in di.items():
        for old, new in values.items():
            df = df.replace(to_replace={variable: old}, value=new)
    return df


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


dat2d = datamatch(dat2d, dat3d, 'threshold')
#%% Renaming
redict = {
    "nsim": {
        0: '(0) fiber diameter: 2\u03BCm',
        1: '(1) fiber diameter: 5\u03BCm',
        2: '(2) fiber diameter: 8\u03BCm',
        3: '(3) fiber diameter: 11\u03BCm',
        4: '(4) fiber diameter: 13\u03BCm',
    }
}
dat2d = dat2d.rename(columns={'threshold': '2D', 'threshold3d': '3D'})
datre = rename_var(dat2d, redict)
# datre = dat2d
#%%
for i, sample in enumerate(samples):
    plotdata = datre[datre['sample'] == sample]

    plot = ggpubr.ggpaired(
        plotdata,
        cond1="2D",
        cond2="3D",
        color="condition",
        line_color="gray",
        line_size=0.5,
        point_size=1.2,
        palette="npg",
        facet_by="nsim",
        xlab=False,
        ylab="Threshold (mA)",
        legend="none",
        scales="free_y",
        title=f"Activation thresholds for sample {sample_name} {bigcomp[i]} contact",
    )

    # plot.plot()

    plot.save(r'out/analysis/{}-{}_{}_{}.svg'.format(threed, sample, models[0], sims[0]), dpi=600)


#%%
for nsim in pd.unique(dat2d['nsim']):

    plotdata = dat2d[dat2d['nsim'] == nsim]

    sample_labels = ['rostral contact', 'caudal contact']

    plot = ggpubr.ggpaired(
        plotdata,
        cond1="2D",
        cond2="3D",
        color="condition",
        line_color="gray",
        line_size=0.5,
        point_size=1.2,
        palette="npg",
        facet_by="sample",
        xlab=False,
        ylab="Activation Threshold (mA)",
        legend="none",
        # scales="free_y",
        title=f"Activation thresh for 2D ex models vs. full-3D model ({(nsim + 1) * 2}um)",
    )

    # plot.plot()

    plot.save(r'out/analysis/{}_{}_{}_{}.svg'.format(threed, models[0], sims[0], nsim), dpi=600)
