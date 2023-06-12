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
        'include_downstream': True,  # NOTE REMOVED 651
        'indices': {
            'sample': [
                250,
                251,
                252,
                270,
                271,
                272,
                370,
                371,
                372,
                570,
                571,
                572,
                650,
                652,
                670,
                671,
                672,
                2501,
                2509,
                2515,
                2520,
                2521,
                2524,
                2526,
                2529,
                2530,
                2531,
                2534,
                2535,
                2536,
                3701,
                3709,
                3721,
                3729,
                5701,
                5709,
                5721,
                5729,
                6701,
                6709,
                6721,
                6729,
            ],
            'model': [0],
            'sim': [3],
        },
    }
).run()

q.excel_output(
    'out.xlsx',
    sample_keys=[['scale', 'shrinkage']],
    individual_indices=False,
    config_paths=False,
    console_output=True,
    column_width=20,
)
