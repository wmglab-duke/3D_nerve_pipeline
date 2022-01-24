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

import numpy as np

import matplotlib.pyplot as plt
from src.core.query import Query
sys.path.append(os.getcwd()+'/subrepos/nd_line')
from nd_line.nd_line import nd_line
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


# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14*2]) / 2)

samp3d = 473

samples =  [470,472]

models = [0]

sims = [33]

nsim_count = 5

bigcomp = {0:'anodic leading',1:'cathodic leading'}


q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': samples,
        'model': models,
        'sim': sims
    }
}).run()



def actual_zpos(dat3d,samp3d,model,sim):
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
        
   
#%%
dats = []
for sample in samples:
    for n in range(nsim_count):
       ap_path = os.path.join('samples',str(sample),'models',str(models[0]),'sims',str(sims[0]),'n_sims',str(n),'data','outputs','AP_info.csv')
       dats.append(pd.read_csv(ap_path))
dat2d = pd.concat(dats)
dats3 = []
for n in range(nsim_count):
   ap_path = os.path.join('samples',str(samp3d),'models',str(models[0]),'sims',str(sims[0]),'n_sims',str(n),'data','outputs','AP_info.csv')
   dats3.append(pd.read_csv(ap_path))
dat3d = pd.concat(dats3)
dat2d = dat2d.reset_index()
dat3d = dat3d.reset_index()

#%%
dat3z = actual_zpos(dat3d,samp3d,models[0],sims[0])

datfinal = datamatch(dat2d,dat3z,'activation_zpos')

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
                           ylab = "Activation zpos (um)",
                            legend = "none",
                            # scales="free_y",
                           title = "Activation zpos for 2D ex models vs. full-3D model {}vs{}".format(samp3d,sample))
    
    # plot.plot()
    
    plot.save(r'out/analysis/{}-{}_{}_{}_ap.png'.format(samp3d,sample,models[0],sims[0]),dpi=600)

#%%


