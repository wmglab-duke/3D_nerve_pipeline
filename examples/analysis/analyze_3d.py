# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 13:12:46 2021

@author: dpm42
"""

import os
import sys

from statistics import variance

sys.path.append(os.path.sep.join([os.getcwd(), '']))

import numpy as np
import pandas as pd
os.chdir('D:/ASCENT/fresh')

import matplotlib.pyplot as plt
from src.core.query import Query
from scipy.stats import ttest_rel, ttest_ind
from scipy.stats import f as eff

def ftest(x,y):
    F=variance(x)/variance(y)
    df1 = len(x)-1
    df2 = len(y)-1
    p_value = eff.cdf(F, df1, df2)
    return p_value

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': [3070,3150,3230],
        'model': [0],
        'sim': [3000]
    }
}).run()

print('NOTE: assumes that nsim 0 is 2d and nsim 1 is 3d')

data = q.ggpaired_3D()

sample_labels = ['rostral contact','center','caudal contact']

datas = {'rostral':data[0][0],
        'center':data[1][0],
        'caudal':data[2][0],
        '3D':data[0][1]}

summary = {}
for k,v in datas.items():
    summary[k]={
        "mean":np.mean(v),
        "std":np.std(v,ddof=1),}

summ={}
for i,dat in enumerate(data):
    t = ttest_rel(data[i][0],data[i][1])
    tf = ttest_ind(data[i][0],data[i][1])
    f = ftest(data[i][0],data[i][1])
    print('For {} paired t-test p value is {}, unpaired t-test value is {}, and f test p value is {}'.format(sample_labels[i],t.pvalue, tf.pvalue,f))
    summ[sample_labels[i]] = {        
            "mean":np.mean(data[i][0]),
            "std":np.std(data[i][0],ddof=1),
            "paired_3D_t":t.pvalue,
            "up_3D_t":tf.pvalue,
            "3D_F":f
            }
summ["3D"]={        
        "mean":np.mean(data[i][1]),
        "std":np.std(data[i][1],ddof=1),
        "paired_3D_t":np.nan,
        "up_3D_t":np.nan,
        "3D_F":np.nan
        }
        
sumdat = pd.DataFrame(summ)
# sumdat = sumdat.round(decimals=5)
print(sumdat)