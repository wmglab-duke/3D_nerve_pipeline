"""Compare thresholds across models using a boxplot.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

For more controls over how the plotting occurs, see the seaborn documentation on barplot:
https://seaborn.pydata.org/generated/seaborn.swarmplot.html
RUN THIS FROM REPOSITORY ROOT
"""

import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))

os.chdir('../..')

import matplotlib.pyplot as plt
import seaborn as sns
from src.core.query import Query

sns.set_theme()

q = Query(
    {
        'partial_matches': False,
        'include_downstream': True,
        'indices': {'sample': [253999, 253000], 'model': [0], 'sim': [333]},
    }
).run()

data = q.common_data_extraction(data_types=['threshold'], threed=True)
g = sns.swarmplot(data=data, x='sample', y='threshold')
plt.title('Threshold swarmplot comparison')
