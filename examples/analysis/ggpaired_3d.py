#!/usr/bin/env python3.7

"""
The copyrights of this software are owned by Duke University.
Please refer to the LICENSE.txt and README.txt files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent
"""

# RUN THIS FROM REPOSITORY ROOT

"""
Created on Thu Nov 11 11:22:44 2021

@author: dpm42
"""
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
os.chdir('D:/ASCENT/fresh')

import matplotlib.pyplot as plt
from src.core.query import Query

# set default fig size
plt.rcParams['figure.figsize'] = list(np.array([16.8, 10.14]) / 2)


# initialize and run Querys
# q = Query({
#     'partial_matches': True,
#     'include_downstream': True,
#     'indices': {
#         'sample': [3],
#         'model': [0, 1, 2, 3],
#         'sim': [0]
#     }
# }).run()
#
# # builds heatmaps
# q.barcharts_compare_models(model_labels=['Original orientation',
#                                          '90-degree rotation',
#                                          '180-degree rotation',
#                                          '270-degree rotation'],
#                            title= 'Upper lobe activation thresholds',
#                            fascicle_filter_indices=[2, 3, 9, 7, 13, 15, 4, 0, 6, 10, 15, 18, 16, 1, 11, 17, 5, 8, 14, 12, 21, 23, 20, 25, 32],
#                            logscale=True)

# q = Query({
#     'partial_matches': True,
#     'include_downstream': True,
#     'indices': {
#         'samples': [1017],
#         'model': [4, 5, 6, 7],
#         'sim': [1040]
#     }
# }).run()

q = Query({
    'partial_matches': False,
    'include_downstream': True,
    'indices': {
        'sample': [3070,3150,3230],
        'model': [0],
        'sim': [3000]
    }
}).run()

# builds heatmaps
# q.barcharts_compare_models(logscale=False,
#                            model_labels=['Model 0: Veltink Epineurium, \n              Veltink Perineurium',
#                                          'Model 1: Veltink Epineurium, \n              Goodall Perineurium',
#                                          'Model 2: Goodall Epineurium, \n              Veltink Perineurium',
#                                          'Model 3: Goodall Epineurium, \n              Goodall Perineurium']
#                            )
data = q.ggpaired_3D()

sample_labels = ['rostral contact','center','caudal contact']

for i,sampdat in enumerate(data):
    d2 = sampdat[0]
    d3 = sampdat[1]
        
    di = {}
    di[sample_labels[i]] = d2
    di["3D"] = d3

    d = pd.DataFrame(di)
    
    plot = ggpubr.ggpaired(d, cond1 = sample_labels[i], cond2 = "3D",
                           color = "condition", 
                           line_color = "gray", 
                           line_size = 0.5, 
                           point_size = 1.3,
                           palette = "npg",
                           xlab = False,
                           ylab = "Threshold (mA)",
                           legend = None,
                           title = "Activation thresholds for {} extrusion model vs. full-3D model".format(sample_labels[i]))
    
    plot.plot()
    
    plot.save('D:/test.png')
