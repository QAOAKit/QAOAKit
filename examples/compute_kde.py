'''
If pre-trained KDE models do not load correctly due to conflict in sklearn
version or other problems, recompute them by running the script below.
Warning! Can take a very long time.

Recommended way to run is using GNU parallel:

#!/bin/bash

parallel \
    --jobs 3 \
    """
    python compute_kde.py {1}
    """ ::: "$(seq 1 3)"
'''

import numpy as np
import pickle

import sys
from QAOAKit.parameter_optimization import train_kde

p = int(sys.argv[1])
outpath = f"../data/pretrained_models/kde_n=9_p={p}_large_bandwidth_range.p"
print("Outpath: ", outpath)
# Adjust number of jobs as needed
median, kde = train_kde(p, 9, n_jobs=20, bandwidth_range=np.logspace(-4, 1, 40))
pickle.dump((median, kde), open(outpath, "wb"))
