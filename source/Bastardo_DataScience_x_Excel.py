"""
Ese comic de SMBC
"""

"""
1. Importar el grafo de Recon
1. centralidad de cada metabolit0
1. Multiplicar toda la columna por la centralidad
2. Aplicar centralidades a los nodos que estamos evaluando
3. Algo de DataScience por aqui... 
"""

# %% --- IMPORTA LA DATA CON PANDAS

import pandas as pd

dataset_centralidad_base = pd.read_csv('data/recon2_metabolite_centralities_metabolome_weights.csv', index_col= 0 )

import pickle 

# DataFrames Pandas con los datos metabolicos
infile = open('./tmp/excel_dataset.pandas.pkl', 'rb'); df = pickle.load(infile); infile.close()

df = df.set_index('Muestra') # Define la columna muestra como indice

print('Entradas iniciales:', df.shape )

# Esto es lo mismo que 'df_metabolitos' guardado como .pandas.pkl
metabolites_columns = ['Phe', 'Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA', 'C0',
    'C2', 'C3', 'C4', 'C4OH/C3DC', 'C5-OH/C4DC', 'C5', 'C5DC', 'C5:1', 'C6',
    'C6DC', 'C8', 'C8:1', 'C10', 'C10:1', 'C10:2', 'C12', 'C12:1', 'C14',
    'C14:1', 'C14:2', 'C14OH', 'C16', 'C16OH', 'C16:1', 'C16:1OH', 'C18',
    'C18OH', 'C18:1', 'C18:1OH', 'C18:2']

dataset_metabolitos = df[metabolites_columns]

# %% --- DEFINE UNA LISTA DE NODOS A LOS QUE APLICAR LAS CENTRALIDADES
# Lista curada de Nodos de Recon para aplicar la cosa
# coincidimos esto con los que salen de 

nodes_to_aply = {
    'Phe' : ['phe_L_e'], 
    'Met' : ['met_L_e'],
    'Val' : ['val_L_e'], 
    'Leu/Ile' : ['leu_L_e'], #, 'ile_L_e'], 
    'tir' : ['tyr_L_e'], 
    'Pro' : ['pro_L_e'], 
    'Arg' : ['arg_L_e'], 
    'Gly' : ['gly_e'], 
    'Ala' : ['ala_L_e'], 
    'Asp' : ['asp_L_e'], 
    'Glt' : ['glu_L_e'], 
    'Cit' : ['citr_L_e'], 
    'Orn' : ['orn_e'], 
    'SA' : [], 
    'C0' : ['crn_e'], 
    'C2' : [], 
    'C3' : [], 
    'C4' : ['c4crn_e'], 
    'C4OH/C3DC' : ['3bcrn_e'], 
    'C5-OH/C4DC' : [], 
    'C5'    : [], 
    'C5DC'  : ['c5dc_e'], 
    'C5:1'  : ['c51crn_e'], 
    'C6'    : ['c6crn_e'], 
    'C6DC'  : [], 
    'C8'    : ['c8crn_e'], 
    'C8:1'  : ['c81crn_e'], 
    'C10'   : ['c10crn_e'], 
    'C10:1' : ['c101crn_e'], 
    'C10:2' : ['decdicrn_e'], 
    'C12'   : ['c12dc_e'], 
    #'C12:1' : ['dodecenoyl'], 
    'C14'   : [], 
    'C14:1' : ['tetdece1crn_e'], 
    'C14:2' : ['tetdec2crn_e'], 
    'C14OH' : ['3tdcrn_e'], 
    'C16'   : ['3hexdcrn_e'], 
    'C16OH' : ['3hexdcrn_e'], 
    'C16:1' : ['3hdececrn_e'], 
    'C16:1OH' : ['3hdececrn_e'], 
    'C18'   : [], 
    'C18OH' : ['3octdeccrn_e'], 
    'C18:1' : [], 
    'C18:1OH' : ['3octdece1crn_e'], 
    'C18:2' : []
}

# met_e = pd.Series( dataset_centralidad.index )

# met_e[met_e.str.contains('gly', case=False)] 

# %% MULTIPLICA EL DATASET DE CENTRALIDADES POR EL DE METABOLITOS

# Extrae los valores del diccionario para buscar los nodos para ultiplicar
nodes_defined = [] # TODO: posiblemente optimizar esto?
for i in nodes_to_aply: 
    nodes_defined = nodes_defined + nodes_to_aply[i] 

dataset_centralidad = dataset_centralidad_base.loc[nodes_defined,'eigenvector'] # Un DataFrame con las centralidades de las cosas

# %% ESTANDARIZA LA MATRIZ DE METABOLITOS

dataset_metabolitos = dataset_metabolitos.dropna() # Elimina filas con NaNs

from sklearn.preprocessing import RobustScaler
X = dataset_metabolitos.values
X = RobustScaler(with_centering=False, with_scaling=False, quantile_range=(25.0, 75.0), copy=True, unit_variance=False).fit_transform(X)

dataset_metabolitos2 = pd.DataFrame( X, columns= dataset_metabolitos.columns, index= dataset_metabolitos.index )

# Solo valores que tengan correspondencia con nodos
metabolitos = [ a for a in nodes_to_aply if nodes_to_aply[a] != [] ]
dataset_centralidad.index = metabolitos

dataset_metabolitos2 = dataset_metabolitos2[metabolitos]

# %% TRUE MULTIPLICA

dataset_metabolitos2.mul( dataset_centralidad , axis=1 )

# %% --- REMOVEDOR DE OUTLIERS
# remueve todos las filas que tengan un outlier en al menos una columna
# calculando el Z (intercuartil) absoluto, 

import numpy as np
from scipy import stats
no_outliers = dataset_metabolitos2[(np.abs(stats.zscore(dataset_metabolitos2)) < 3).all(axis=1)]

# centralidad_base_reduc = dataset_centralidad_base.loc[nodes_defined,]
# no_outliers = centralidad_base_reduc[(np.abs(stats.zscore(centralidad_base_reduc)) < 3).all(axis=1)]

print('Entradas despues de remocion de outliers:', no_outliers.shape )

# %% --- ESTANDARIZACION Y ESCALADO 

# from sklearn.preprocessing import RobustScaler
# X = no_outliers.values
# X = RobustScaler().fit_transform(X) # normalizing the features
# 
# print('Data escalada para procesamiento:', X.shape )

# %% --- MinMax Scaler

from sklearn.preprocessing import MinMaxScaler
X = no_outliers.values
X = MinMaxScaler(feature_range=(0, 10), copy=True, clip=False).fit_transform(X) 

print('Data MinMax para procesamiento:', X.shape )

# %% --- PARAMETROS DE CLUSTERING

CLUSTERS = 4
DIMENSIONALIDAD = 2

# %% --- Multi Dimensional Scalling

# from sklearn.manifold import MDS
# 
# mds = MDS(n_components=DIMENSIONALIDAD, 
#     random_state=0, 
#     n_jobs=16).fit_transform(X) # HARDCODED n_jobs=16
# 
# mds.shape

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
        #'MDS_Componente_1': mds[:,0] ,
        #'MDS_Componente_2': mds[:,1] , 
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

data_plot.to_csv('./results/dataplot_eigenvector.csv')

# %% --- PARAMETROS DE PLOTTING Y DEMAS

import seaborn as sns

# %% --- PLOT TSNE
plot_tsne = sns.scatterplot(data=data_plot, x="tSNE_Componente_1", y="tSNE_Componente_2", hue="k-labels")
plot_tsne = plot_tsne.get_figure()
plot_tsne.savefig("./doc/img/e_plot_tsne.png")
# %% --- PLOT MDS
# plot_mds = sns.scatterplot(data=data_plot, x="MDS_Componente_1", y="MDS_Componente_2", hue="k-labels")
# plot_mds = plot_mds.get_figure()
# plot_mds.savefig("./doc/img/centralidades_plot_mds.png")
# %% --- PLOT ISOMAP
plot_isomap = sns.scatterplot(data=data_plot, x="isomap_Componente_1", y="isomap_Componente_2", hue="k-labels")
plot_isomap = plot_isomap.get_figure()
plot_isomap.savefig("./doc/img/e_plot_isomap.png")
# %% --- KMEANS POST CLUSTERING

Xpc = data_plot[['tSNE_Componente_1','tSNE_Componente_2']]

kmeans_tsne = MiniBatchKMeans(n_clusters= CLUSTERS, 
    batch_size=100, 
    random_state=0).fit(Xpc.values)

kmeans_tsne.labels_

# %% --- COSAS DE GEPHI Y ESO

# './data/Recon2_metabolite_projection.gpickle'



