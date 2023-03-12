#!/usr/bin/env python3.7

"""Generate an excel summary of parameters used in runs of ASCENT.

The copyrights of this software are owned by Duke University.
Please refer to the LICENSE and README.md files for licensing instructions.
The source code can be found on the following GitHub repository: https://github.com/wmglab-duke/ascent.

RUN THIS FROM REPOSITORY ROOT
"""

import os
import sys

os.chdir('../..')

from src.core.query import Query

sys.path.append(os.path.sep.join([os.getcwd(), '']))

# initialize and run Querys
q = Query(
    {
        'partial_matches': True,
        'include_downstream': True,
        'indices': {
            'sample': [
                250,
                252,
                370,
                372,
                570,
                572,
                650,
                652,
                670,
                672,
                2509,
                2529,
                3709,
                3729,
                5709,
                5729,
                6709,
                6729,
            ],
            'model': [0],
        },
    }
).run()

q.excel_output(
    'out.xlsx',
    sample_keys=[],
    model_keys=[['min_radius_enclosing_circle'], ['3d_nerve_override']],
    individual_indices=False,
    config_paths=False,
    console_output=False,
    column_width=20,
)
