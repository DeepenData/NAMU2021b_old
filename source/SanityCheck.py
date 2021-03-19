"""
"""

import numpy  as np
import pandas as pd

local_baseline   = pd.read_csv('glico_resultado_local.csv', index_col=0)
cluster_baseline = pd.read_csv('glic2.csv', index_col=0)

BASELINES = np.allclose(local_baseline.to_numpy(), cluster_baseline.to_numpy(), 
    rtol=1e-05, atol=1e-08, equal_nan=False)

print('Baselines son iguales?:', BASELINES)

local_removed   = pd.read_csv('perturbed_naive_FPGS3m_short.csv', index_col=0)
cluster_removed = pd.read_csv('removed_FPGS3m.csv.csv', index_col=0)

REMOVED = np.allclose(local_baseline.to_numpy(), cluster_baseline.to_numpy(), 
    rtol=1e-05, atol=1e-08, equal_nan=False)

print('Removidos son iguales?:', REMOVED)
