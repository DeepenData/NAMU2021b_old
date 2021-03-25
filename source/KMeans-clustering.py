#!/bin/python
"""
https://xkcd.com/1838/ 
    Machine Learning, XKCD
"""
# %% --- IMPORTANDO EL DATASET

import pickle 

# DataFrames Pandas con los datos metabolicos
infile = open('./tmp/excel_dataset.pandas.pkl',     'rb'); df = pickle.load(infile); infile.close()
#infile = open('./tmp/excel_metabolitos.pandas.pkl', 'rb'); df_metabolitos = pickle.load(infile); infile.close()

df = df.set_index('Muestra') # Define la columna muestra como indice

print('Entradas iniciales:', df.shape )

# Esto es lo mismo que 'df_metabolitos' guardado como .pandas.pkl
metabolites_columns = ['Phe', 'Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA', 'C0',
    'C2', 'C3', 'C4', 'C4OH/C3DC', 'C5-OH/C4DC', 'C5', 'C5DC', 'C5:1', 'C6',
    'C6DC', 'C8', 'C8:1', 'C10', 'C10:1', 'C10:2', 'C12', 'C12:1', 'C14',
    'C14:1', 'C14:2', 'C14OH', 'C16', 'C16OH', 'C16:1', 'C16:1OH', 'C18',
    'C18OH', 'C18:1', 'C18:1OH', 'C18:2']

df_metabolitos = df[metabolites_columns]

print('Seleccion solo valores metabolicos:', df_metabolitos.shape )

# %% --- REMOVEDOR DE NaNs

df_metabolitos = df_metabolitos.dropna()             # Elimina filas con NaNs

print('Entradas despues de remocion de NaNs:', df_metabolitos.shape )

# %% --- REMOVEDOR DE OUTLIERS
# remueve todos las filas que tengan un outlier en al menos una columna
# calculando el Z (intercuartil) absoluto, 

import numpy as np
from scipy import stats
no_outliers = df_metabolitos[(np.abs(stats.zscore(df_metabolitos)) < 3).all(axis=1)]

print('Entradas despues de remocion de outliers:', no_outliers.shape )

# %% --- ESTANDARIZACION Y ESCALADO 

from sklearn.preprocessing import RobustScaler
X = no_outliers.values
X = RobustScaler().fit_transform(X) # normalizing the features

print('Data escalada para procesamiento:', X.shape )

# %% --- PARAMETROS DE CLUSTERING

CLUSTERS = 5
DIMENSIONALIDAD = 2

# %% --- Multi Dimensional Scalling

from sklearn.manifold import MDS

mds = MDS(n_components=DIMENSIONALIDAD, 
    random_state=0, 
    n_jobs=16).fit_transform(X) # HARDCODED n_jobs=16

mds.shape

# %% --- Isomap

from sklearn.manifold import Isomap

isomap = Isomap(n_components=DIMENSIONALIDAD,
    n_jobs=-1).fit_transform(X)

isomap.shape


# %% --- T-Sne

from sklearn.manifold import TSNE

tsne = TSNE(n_components=DIMENSIONALIDAD, 
    method='barnes_hut', 
    n_jobs=-1).fit_transform(X)

tsne.shape

# %% --- K-Means Clustering

# from sklearn.cluster import KMeans
# kmeans = KMeans(n_clusters= CLUSTERS, random_state=0).fit(X)

from sklearn.cluster import MiniBatchKMeans

kmeans = MiniBatchKMeans(n_clusters= CLUSTERS, 
    batch_size=100, 
    random_state=0).fit(X)

kmeans.labels_

# %% --- Spectral clustering
# Es lentisimo...

# from sklearn.cluster import SpectralClustering
# 
# spectral = SpectralClustering(n_clusters= CLUSTERS,
#     assign_labels="discretize",
#     random_state=0, 
#     n_jobs=-1).fit(X)
# 
# spectral.labels_

# %% --- EMPACANDO RESULTADOS EN UN DATAFRAME

import pandas as pd

data_plot = pd.DataFrame(
    {   
        'MDS_Componente_1': mds[:,0] ,
        'MDS_Componente_2': mds[:,1] , 
        'tSNE_Componente_1': tsne[:,0] , 
        'tSNE_Componente_2': tsne[:,1] , 
        'isomap_Componente_1': isomap[:,0] , 
        'isomap_Componente_2': isomap[:,1] , 
        'k-labels' : kmeans.labels_  
        #'spctral-labels' : spectral.labels_
    }, 
    index= no_outliers.index
)

data_plot['k-labels'].replace([0,1,2,3,4], ['A','B','C','D','E'], inplace=True) # Reemplazo de categorias
# data_plot['spctral-labels'].replace([0,1,2,3,4], ['A','B','C','D','E'], inplace=True) # Reemplazo de categorias

data_plot.astype( {
    'k-labels' : 'category' #, 'spctral-labels' : 'category'
} )

data_plot.to_csv('./results/dataplot.csv')

# %% --- 
import seaborn as sns

sns.scatterplot(data=data_plot, x="tSNE_Componente_1", y="tSNE_Componente_2", hue="k-labels")
# %% --- 
sns.scatterplot(data=data_plot, x="MDS_Componente_1", y="MDS_Componente_2", hue="k-labels")
# %% --- 
sns.scatterplot(data=data_plot, x="isomap_Componente_1", y="isomap_Componente_2", hue="k-labels")
# %% --- 
