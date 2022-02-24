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
import pandas as pd
import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import numpy as np
os.chdir('D:/ASCENT/ascent')

import matplotlib.pyplot as plt
from src.core.query import Query

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)

allthreed = [473,653,673]

allsamples = [[470,472],[650,652],[670,672]]

models = [0]

sims = [33]

allsampname = ['4R','6L','6R']

compdat = {}
compdat_pandas = []

for threed,samples,sampname in zip(allthreed,allsamples,allsampname):
    compdat[sampname]={}
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
    #%%
    
    sample_labels = ['rostral contact','caudal contact']
    act_labels = ['anodic-leading contact','cathodic-leading contact']

    for x in sample_labels:
        compdat[sampname][x] = {} 
    import seaborn as sns
    from scipy.stats import pearsonr
    sns.set_theme()
    sns.set(font_scale=1.5)
    
    dat2d = dat2d.rename(columns={'sample':'Slice'})
    #%%
    # Plot sepal width as a function of sepal_length across days
    g = sns.lmplot(
        data=dat2d,
        x="threshold", y="threshold3d", hue="Slice",
        height=5,col = 'nsim',sharey=False,sharex=False
    )
    
    axs = g.axes.ravel()
    axs[0].set_ylabel('3D threshold (mA)')
    plt.suptitle('Activation threshold correlation for sample {}'.format(sampname),fontsize=25)
    plt.subplots_adjust(top=.85,right=.93)
    new_labels = ['Anodic\nLeading', 'Cathodic\nLeading']
    for t, l in zip(g._legend.texts, new_labels):
        t.set_text(l)
    for i,s in enumerate([2,5,8,11,13]):
        ax = axs[i]
        ax.set_title('fiber diam: {}\u03BCm'.format(s))
        corr = {}
        for label_ind,sample in enumerate(samples):
            thisdat = dat2d[(dat2d["nsim"]==i) & (dat2d["Slice"]==sample)]
            correlation = round(pearsonr(thisdat['threshold'],thisdat['threshold3d'])[0],3)
            corr[sample]=correlation
            compdat[sampname][sample_labels[label_ind]][s]= correlation
            compdat_pandas.append({
                'Sample':sampname,
                '2D slice':act_labels[label_ind],
                'fiber diam':s,
                'correlation':correlation })
        ax.legend(labels = ["r="+str(corr[sample]) for sample in samples])
        ax.set_xlabel('2D threshold (mA)')
    g.savefig('out/analysis/threscorr_{}'.format(threed),dpi=400)
#%%
plt.figure()
allcorr = pd.DataFrame(compdat_pandas)
g = sns.scatterplot(data = allcorr,x='fiber diam',y='correlation',style='2D slice',hue = 'Sample',s=100,palette='colorblind')
g = sns.lineplot(data = allcorr,x='fiber diam',y='correlation',style='2D slice',hue = 'Sample',legend = False, palette='colorblind')
plt.legend(bbox_to_anchor=(1.53, .5), loc='center right', borderaxespad=0)
plt.ylabel('Threshold Correlation')
plt.xlabel('Fiber Diameter (\u03bcm)')
plt.title('Activation Threshold Correlation for All Models')

