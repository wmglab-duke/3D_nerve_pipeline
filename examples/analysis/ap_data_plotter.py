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

import matplotlib.pyplot as plt
from src.core.query import Query
sys.path.append(os.getcwd()+'/subrepos/nd_line')
from nd_line.nd_line import nd_line
# os.chdir('D:/ascent/fresh')

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14*2]) / 2)

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

'''
TODO
-add deletion option for vmt
-add function from ggpaired to match 2d and 3d
-add function to get spatial z position from longitudinal z position
'''

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
     
dat2d = q.ap_data(
    delta_V=60,
    absolute_voltage=False,
    delete_vmtime = False)

dat3d = q.ap_data(
    delta_V=60,
    absolute_voltage=False,
    delete_vmtime = False,
    sample_override = samp3d)
#%%
dat3z = actual_zpos(dat3d,samp3d,models[0],sims[0])

dat2d = datamatch(dat2d,dat3z,'activation_zpos')

#%%
