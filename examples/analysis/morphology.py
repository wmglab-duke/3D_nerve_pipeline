import os

import numpy as np
import pandas as pd
import pickle5 as pickle
import seaborn as sns

# load in slides
from matplotlib import pyplot as plt

os.chdir('../../')
slide_dict = {}
for file in os.listdir('input/slides'):
    if file.endswith('.obj'):
        slide_dict[file[:2]] = pickle.load(open(os.path.join('input/slides', file), 'rb'))
alldata = []
for sample_name, slides in slide_dict.items():
    thisdata = []
    simplebranch = -1
    fascicle_n = 0
    for i, slide in enumerate(slides):
        areas = [i.area() / 1000000 for f in slide.fascicles for i in f.inners]
        periareas = [f.outer.area() / 1000000 for f in slide.fascicles]
        nervearea = slide.nerve.area() / 1000000
        if len(areas) != fascicle_n:
            fascicle_n = len(areas)
            simplebranch += 1
            branched = 1
        else:
            branched = 0
        thisdata.append(
            {
                'sample': sample_name,
                'slide_position': i,
                'fascicle_count': len(areas),
                'total_endo_area': sum(areas),
                'total_nerve_area': nervearea,
                'total_epineurium_area': nervearea - sum(periareas),
                'endo-nerve_ratio': sum(areas) / nervearea,
                'mean_inner_area': np.mean(areas),
                'branch_event': branched,
            }
        )
    print(f'Simple branch count for {sample_name}:{simplebranch}')
    alldata.append(pd.DataFrame(thisdata))
alldata = pd.concat(alldata)
alldata.reset_index(inplace=True)
analyses = [
    'fascicle_count',
    'total_endo_area',
    'total_nerve_area',
    'total_epineurium_area',
    'endo-nerve_ratio',
    'mean_inner_area',
]
# lineplot for each
for analysis in analyses:
    plt.figure()
    sns.lineplot(data=alldata, x='slide_position', y=analysis, hue='sample')
# grouping analysis
grouped = alldata.groupby(['sample'])
analysis = grouped.agg(
    {
        'mean_inner_area': [np.mean, np.std],
        'total_endo_area': [np.mean, np.std],
        'total_nerve_area': [np.mean, np.std],
        'total_epineurium_area': [np.mean, np.std],
        'endo-nerve_ratio': [np.mean, np.std],
        'fascicle_count': [np.mean, np.std, np.min, np.max],
    }
)
analysis.columns = column_analysis = ["_".join(col_name).rstrip('_') for col_name in analysis.columns]
# if sample has L it is left, otherwise right
lr = analysis.index.str.contains('L')
analysis['side'] = ['left' if x else 'right' for x in lr]
final = analysis.round(decimals=2)
# TODO make this script check for merge split events
