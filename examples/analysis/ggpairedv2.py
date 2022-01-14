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
import os
os.environ['R_HOME'] = r'C:\Users\dpm42\Anaconda3\envs\ascent\lib\R'
import rpy2.robjects as robjects
import rpy2
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
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

threed = 653

samples = [650,652]

models = [0]

sims = [3]

dats = []

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': samples,
        'model': models,
        'sim': sims
    }
}).run()

# builds heatmaps
# q.barcharts_compare_models(logscale=False,
#                            model_labels=['Model 0: Veltink Epineurium, \n              Veltink Perineurium',
#                                          'Model 1: Veltink Epineurium, \n              Goodall Perineurium',
#                                          'Model 2: Goodall Epineurium, \n              Veltink Perineurium',
#                                          'Model 3: Goodall Epineurium, \n              Goodall Perineurium']
#                            )
dat2d = q.threshdat(sl=False,meanify=False)

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': [threed],
        'model': models,
        'sim': sims
    }
}).run()

dat3d = q.threshdat3d(meanify = False)


def datamatch(dest,dat3d,importval):
    dest[importval+'3d'] = np.nan
    for i in range(len(dest)):
        row = dest.iloc[i,:]
        val = dat3d[(dat3d["model"]==row['model']) &
                           (dat3d["sim"]==row['sim']) &
                           (dat3d["nsim"]==row['nsim']) &
                           (dat3d["index"]==row['index'])][importval]
        val=list(val)
        if len(val)!=1:sys.exit('issue here')
        dest.iloc[i,-1]=val[0]
    if np.any(dest[importval]==np.nan):
        sys.exit('issue here too')
    return dest

dat2d = datamatch(dat2d,dat3d,'threshold')

dat2d = dat2d[dat2d['nsim']==2]

sample_labels = ['rostral contact','caudal contact']

plot = ggpubr.ggpaired(dat2d, cond1 = "threshold", cond2 = "threshold3d",
                       color = "condition", 
                       line_color = "gray", 
                       line_size = 0.5, 
                       point_size = 1.2,
                       palette = "npg",
                       facet_by = "sample",
                       xlab = False,
                       ylab = "Threshold (mA)",
                        legend = "none",
                        # scales="free_y",
                       title = "Activation thresholds for 2D extrusion models vs. full-3D model (6um)")

plot.plot()

plot.save(r'D:\Box Sync\Home Folder dpm42\lab_meeting/4L_extrusion_comp.png',dpi=600)
