"""Returns thresholds from the selected sample/model/sim combos as a dataframe.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

Use argument meanify=True to instead get the mean threshold for each nsim with stats
RUN THIS FROM REPOSITORY ROOT
"""

# %% imports
import os
import sys

sys.path.append(os.path.sep.join([os.getcwd(), '']))
from src.core.query import Query

os.chdir('../..')
# %% metadata
samples = [5701]

models = [0]

sims = [3]

dats = []

# %% run query search
q = Query(
    {
        'partial_matches': False,
        'include_downstream': True,
        'indices': {'sample': samples, 'model': models, 'sim': sims},
    }
).run()

# %% obtain thresholds
data = q.data(ignore_missing=False, tonly=True)
data.to_csv(r'C:\Users\dpm42\Box Sync\Home Folder dpm42\Vp_baby\third_go/fiberdat.csv')
