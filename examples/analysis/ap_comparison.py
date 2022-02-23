#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))
os.chdir('D:/ascent/ascent')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src.core.query import Query
sys.path.append(os.getcwd()+'/subrepos/nd_line')
from nd_line.nd_line import nd_line
sys.path.append(r'D:\ASCENT\ascent')
os.chdir(r'D:\ASCENT/ascent')
sys.path.append(os.path.sep.join([os.getcwd(), '']))
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
os.chdir('D:/ASCENT/ascent')
# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)

samp3d = 673

samples =  [670,672]

models = [0]

sims = [3]

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': samples,
        'model': models,
        'sim': sims
    }
}).run()

def actual_zpos(dat,samp3d,model,sim):
    dat3d = dat.reset_index()
    fiberdir = os.path.join('samples',str(samp3d),'models',str(model),'sims',str(sim),'3D_fiberset')
    fibers3d = [x for x in os.listdir(fiberdir) if x.endswith('.dat')]
    for file in fibers3d:
        f_ind = int(os.path.splitext(file)[0])
        coord = np.loadtxt(os.path.join(fiberdir,file),skiprows=1)
        fiberline = nd_line(coord)
        datnew = dat3d[(dat3d["model"]==model) &
                           (dat3d["sim"]==sim) &
                           (dat3d["index"]==f_ind)]
        for index,row in datnew.iterrows():
            actual_zpos = fiberline.interp(row['long_ap_pos'])[2]
            dat3d.loc[index,'activation_zpos']=actual_zpos
    return dat3d

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

dat2d = []
dat3d = []

for sample in samples:
    for model in models:
        for sim in sims:
            simpath = 'samples/{}/models/{}/sims/{}'.format(sample,model,sim)
            for folder in os.listdir(simpath+'/n_sims'):
                inpath = '{}/n_sims/{}/data/outputs/AP_info.csv'.format(simpath,folder)
                dat2d.append(pd.read_csv(inpath))
for model in models:
    for sim in sims:
        simpath = 'samples/{}/models/{}/sims/{}'.format(samp3d,model,sim)
        for folder in os.listdir(simpath+'/n_sims'):
            inpath = '{}/n_sims/{}/data/outputs/AP_info.csv'.format(simpath,folder)
            dat3d.append(pd.read_csv(inpath))
dat2d = pd.concat(dat2d)
dat3d = pd.concat(dat3d)

#%%
dat3final = actual_zpos(dat3d,samp3d,models[0],sims[0])

dat2d = datamatch(dat2d,dat3final,'activation_zpos')

#%%

for sample in samples:
    plotdata = dat2d[dat2d['sample']==sample]

    plot = ggpubr.ggpaired(plotdata, cond1 = "activation_zpos", cond2 = "activation_zpos3d",
                           color = "condition", 
                           line_color = "gray", 
                           line_size = 0.5, 
                           point_size = 1.2,
                           palette = "npg",
                           facet_by = "nsim",
                           xlab = False,
                           ylab = "Position (um)",
                            legend = "none",
                            scales="free_y",
                           title = "Activation position for 2D ex models vs. full-3D model {}vs{}".format(samp3d,sample))
    
    # plot.plot()
    
    plot.save(r'out/analysis/alsotest.png',dpi=600)

