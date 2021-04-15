"""Esto es codigo test que calcula los delta y FC, y los mete en un modelo de grafo
Esta conectado a las salidas CSV de centralidades base y perturbadas de Recon2. 
Guarda el resultado del fold-change en el grafo de Recon2 y exporta a Gephi y Graphml"""

# %% --- Importando modulos varios

import networkx as nx
import numpy as np
from numpy.lib.shape_base import column_stack
import pandas as pd

# %% --- IMPORTANDO EL GRAFO

# Cuardado desde el archivo ModelToGraph.py
G = nx.read_gpickle('./data/Recon2_rxn_proyected.gpickle')

# %% --- IMPORTANDO LOS DATASETS

baselines  = pd.read_csv( "./results/baseline.csv" , index_col=0 )
baselines.drop(labels='r0399', axis=0, inplace=True)
perturbado = pd.read_csv( "./results/perturbadas_r0399.csv" , index_col=0 )


perturbado.rename(columns={'katz_centrality_numpy': 'katz_centrality'}, inplace=True)

print(len(baselines), len(perturbado))


# Este es solo para un unico nodo removido

# %% --- REINDEXANDO DATASETS
# Esto es porque me di cuenta que se generaron de forma caotica, así que es 
# necesario ordenar las columnas, y reindexar los indices por nodos removidos

baselines         = baselines.reindex( columns = perturbado.columns.sort_values() )
perturbado        = perturbado.reindex( columns = perturbado.columns.sort_values(), index = baselines.index  )
print(
all(baselines.columns == perturbado.columns),
all(baselines.index == perturbado.index))
print(
all(perturbado.katz_centrality >0),
all(baselines.katz_centrality >0))

baselines.drop(labels='katz_centrality', axis=1, inplace=True)
perturbado.drop(labels='katz_centrality', axis=1, inplace=True)


# %% --- CALCULA FOLD CHANGE

fold_change = np.log2(perturbado / baselines)

fold_change.to_csv("./results/fold_change_pku.csv")

# %% --- APLICA ATRIBUTOS AL GRAFO Y LO EXPORTA
for col in fold_change.columns:
    nx.set_node_attributes( G , dict( fold_change[col] ) , str(col) ) 

nx.write_gexf(   G, "./results/graph_files/recon2_FC.gexf")
nx.write_graphml(G, "./results/graph_files/recon2_FC.graphml")

# %% RANKING DE COSAS MODIFICADAS

#def get_top_fc(fold_change, ascending=False):


top_positive = pd.DataFrame({ col : list( fold_change[col].sort_values(ascending=False).head(20).index ) \
    for col in fold_change.columns })


top_positive.to_csv("./results/PKU_recon2rxn_top_positive_fold_change.csv")


top_negative = pd.DataFrame({ col : list( fold_change[col].sort_values(ascending=True).head(20).index ) \
    for col in fold_change.columns })

top_negative.to_csv("./results/PKU_recon2rxn_top_negative_fold_change.csv")