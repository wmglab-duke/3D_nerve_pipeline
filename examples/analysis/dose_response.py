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
import seaborn as sb

threed = 473

samples = [470]

samplename = '4R'

models = [0]

sims = [33]

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

data = pd.concat([dat2d,dat3d])
data.reset_index(inplace=True)

# fig,axs = plt.subplots(ncols=5)
def rename_var(df,di):
    for variable,values in di.items():
        for old,new in values.items():
            df = df.replace(to_replace={variable:old},value=new)
    return df

redict = {
  # "nsim":{
  #     0:'(0) fiber diameter: 2\u03BCm',
  #     1:'(1) fiber diameter: 5\u03BCm',
  #     2:'(2) fiber diameter: 8\u03BCm',
  #     3:'(3) fiber diameter: 11\u03BCm',
  #     4:'(4) fiber diameter: 13\u03BCm'
  #     },
  "sample":{
      samples[0]:'2D',
      threed:'3D'
      }
}
data = rename_var(data,redict)


for i in range(5):
    # ax = axs[i]
    sb.color_palette("tab10")
    # plt.figure()
    plotdata = data[data.nsim==i]
    plotdata = plotdata[plotdata['sample']!=670]
    import seaborn as sb
    sb.ecdfplot(data=plotdata,x='threshold',hue = 'sample')
    plt.xscale('log')
    plt.ylabel('Proportion of Fibers Activated')
    plt.xlabel('Activation Threshold (mA, log scale)')
    plt.title('Threshold eCDF for Sample {} (2D slice - highest r)'.format(samplename))
plt.text(.05,-.25,'Note: fiber diameters (\u03bcm) from left to right: [13, 11, 8, 5, 2]',fontstyle='italic')
plt.savefig(r'out/analysis/{}_ecdf.png'.format(threed),dpi=400,bbox_inches = 'tight')
