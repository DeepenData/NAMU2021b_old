"""Esto es codigo test que calcula los delta y FC, y los mete en un modelo de grafo
Esta conectado a las salidas CSV de centralidades base y perturbadas de Recon2. 
Guarda el resultado del fold-change en el grafo de Recon2 y exporta a Gephi y Graphml"""

# %% --- Importando modulos varios

import networkx as nx
import numpy as np
import pandas as pd

# %% --- IMPORTANDO EL GRAFO

# Cuardado desde el archivo ModelToGraph.py
G = nx.read_gpickle('./data/Recon2_rxn_proyected.gpickle')

# %% --- IMPORTANDO LOS DATASETS

baselines  = pd.read_csv( "tmp.baseline.csv" , index_col=0 )

perturbado = pd.read_csv( "tmp.perturbadas_r0399.csv" , index_col=0 )
# Este es solo para un unico nodo removido

# %% --- REINDEXANDO DATASETS
# Esto es porque me di cuenta que se generaron de forma caotica, así que es 
# necesario ordenar las columnas, y reindexar los indices por nodos removidos

baselines = baselines.reindex( columns = baselines.columns.sort_values() )

perturbado = perturbado.reindex( columns = baselines.columns.sort_values() , index = baselines.index )
# Lo del indice es para contar por los nodos removidos

# %% --- CALCULA FOLD CHANGE

fold_change = perturbado / baselines

fold_change.to_csv("tmp.fold_change_pku.csv")

# %% --- APLICA ATRIBUTOS AL GRAFO Y LO EXPORTA
for col in fold_change.columns:
    nx.set_node_attributes( G , dict( fold_change[col] ) , str(col) ) 

nx.write_gexf(   G, "recon2_FC.gexf")
nx.write_graphml(G, "recon2_FC.graphml")
