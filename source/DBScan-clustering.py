#!/bin/python
"""
"https://xkcd.com/1838/"
"""
# %% --- IMPORTANDO LIBRERIAS Y DATASETS
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

# %% --- IMPORTANDO EL DATASET

import pickle 

# DataFrames Pandas con los datos metabolicos
infile = open('./tmp/excel_dataset.pandas.pkl',     'rb'); df = pickle.load(infile); infile.close()
#infile = open('./tmp/excel_metabolitos.pandas.pkl', 'rb'); df_metabolitos = pickle.load(infile); infile.close()

df = df.set_index('Muestra') # Define la columna muestra como indice

# Esto es lo mismo que 'df_metabolitos' guardado como .pandas.pkl
metabolites_columns = ['Phe', 'Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA', 'C0',
    'C2', 'C3', 'C4', 'C4OH/C3DC', 'C5-OH/C4DC', 'C5', 'C5DC', 'C5:1', 'C6',
    'C6DC', 'C8', 'C8:1', 'C10', 'C10:1', 'C10:2', 'C12', 'C12:1', 'C14',
    'C14:1', 'C14:2', 'C14OH', 'C16', 'C16OH', 'C16:1', 'C16:1OH', 'C18',
    'C18OH', 'C18:1', 'C18:1OH', 'C18:2']

df_metabolitos = df[metabolites_columns]

# %% --- DBSCAN
# DBSCAN - Density-Based Spatial Clustering of Applications with Noise. 
# Finds core samples of high density and expands clusters from them. 
# Good for data which contains clusters of similar density.

EPSILON      = 0.5
MIN_MUESTRAS = 2

dataset = df_metabolitos.dropna() # Elimina filas con NaNs

dataset = dataset.to_numpy() # Pasa a NumPy

from sklearn.preprocessing import StandardScaler, normalize, minmax_scale

dataset = minmax_scale( dataset )

dataset = normalize( dataset )

dbscan = DBSCAN(eps= EPSILON, min_samples = MIN_MUESTRAS).fit( dataset )

#print(dbscan.labels_)
uniq, frec = np.unique(dbscan.labels_, return_counts = True) 
print('Clusters:', uniq)
print('Distribucion:', frec)

# %% --- Compute DBSCAN
from sklearn import metrics

db = DBSCAN(eps= EPSILON, min_samples = MIN_MUESTRAS).fit( dataset )
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)


# %% --- Plot result
import matplotlib.pyplot as plt

X = dataset

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 0] # Circulo vacio

    class_member_mask = (labels == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=14)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
             markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()
# %% --- 
