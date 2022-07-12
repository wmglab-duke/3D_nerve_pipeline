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
import seaborn as sb
from scipy.stats import pearsonr
import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import numpy as np
os.chdir('D:/ASCENT/ascent')

import matplotlib.pyplot as plt
from src.core.query import Query

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)

threed = 253

samples = [250,252]

sample_name = '2L'

models = [0]

sims = [33]

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

def rename_var(df,di):
    for variable,values in di.items():
        for old,new in values.items():
            df = df.replace(to_replace={variable:old},value=new)
    return df

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
#%% Renaming
redict = {
  "nsim":{
      0:'fiber diameter: 2\u03BCm',
      1:'fiber diameter: 5\u03BCm',
      2:'fiber diameter: 8\u03BCm',
      3:'fiber diameter: 11\u03BCm',
      4:'fiber diameter: 13\u03BCm'
      }
}
dat2d = dat2d.rename(columns = {'threshold':'2D','threshold3d':'3D'})
datre = rename_var(dat2d,redict)
dat2dnew = datre.drop(columns = '3D').rename(columns={'2D':'threshold'})
dat2dnew['dataset']='2D'
dat3dnew = datre.drop(columns = '2D').rename(columns={'3D':'threshold'})
dat3dnew['dataset']='3D'
datfinal = pd.concat([dat2dnew,dat3dnew])

# datre = dat2d
#%%
sb.set(font_scale = 1.5)
for i,sample in enumerate(samples):
    plotdata = datfinal[datfinal['sample']==sample]
    g = sb.catplot(data = plotdata,kind = 'swarm',col='nsim',hue='inner',y='threshold',x='dataset',sharey=False,palette='colorblind')
    plt.subplots_adjust(top=.85)
    plt.suptitle('Activation thresholds by fascicle (Sample {}, 2D slice: {} contact)'.format(sample_name,bigcomp[i]))
    axs = g.axes.ravel()
    axs[0].set_ylabel('Activation threshold (mA)')
    plt.subplots_adjust(top=.85)
    for j,s in enumerate([2,5,8,11,13]):
        ax = axs[j]
        corr = {}
        thisdat = dat2d[(dat2d["nsim"]==j) & (dat2d["sample"]==sample)]
        corr[sample]=round(pearsonr(thisdat['2D'],thisdat['3D'])[0],3)
        leg = ax.legend(labels = ["r="+str(corr[sample])],handlelength=0, handletextpad=0, fancybox=True, loc = 'lower center')
        for item in leg.legendHandles:
            item.set_visible(False)
        ax.set_xlabel('2D threshold (mA)')

    plt.savefig('out/analysis/colorthresh{}-{}.png'.format(sample_name,bigcomp[i]),dpi=500)
  