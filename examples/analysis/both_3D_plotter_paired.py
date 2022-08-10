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

samples3d = [2, 3]

models = [0]

sims = [33]

dats = []


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


for sample in samples3d:
    q = Query(
        {
            'partial_matches': False,
            'include_downstream': True,
            'indices': {'sample': [sample], 'model': models, 'sim': sims},
        }
    ).run()

    dats.append(q.threshdat3d(meanify=False))

data = datamatch(dats[0], dats[1], 'threshold')


def rename_var(df, di):
    for variable, values in di.items():
        for old, new in values.items():
            df = df.replace(to_replace={variable: old}, value=new)
    return df


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
data = data.rename(columns={'threshold': '2D', 'threshold3d': '3D'})
datre = rename_var(data, redict)
# datre = dat2d

plot = ggpubr.ggpaired(
    data,
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
    title="Activation thresholds for 2D extrusion vs 3D curvilinear axons",
)

# plot.plot()

plot.save(r'out/analysis/curvyvsstraight.png', dpi=600)
