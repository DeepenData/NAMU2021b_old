#!/bin/python
"""Genera un grafo a partir de un modelo de COBRA en formato JSON. Lo exporta en
un formato apto para NetworkX y Graph-tool"""
# %% --- Carga de librerias base

INPUT_MODEL='./data/stimulated_2021.json' # TODO: hacer esto una variable ambiental

# %% --- Sección interesante

import warnings; warnings.filterwarnings('ignore') # Ignora warnings
from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL)

# print("Iniciando optimización del modelo", time.asctime( time.localtime(time.time()) ) )
# solution_fba = model.optimize() # Optimización del modelo para check
# solution_fba.fluxes # Check si los flujos funcionan

import networkx as nx
from cobra.util.array import create_stoichiometric_matrix

S_matrix = create_stoichiometric_matrix(model)
S_matrix = (abs(S_matrix) )
S_matrix = (S_matrix > 0.0).astype(np.int_)

projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
np.fill_diagonal(projected_S_matrix, 0)

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)

G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("Modelo no completamente conexo. Eligiendo componente más grande.")
    
    largest_component = max(nx.connected_components(G), key=len)
    
    print( len(G.nodes) - len(largest_component), 'nodos removidos.')
    print( len(G.nodes), "nodes finales.")

    G = G.subgraph(largest_component)

node_dict = lambda l : dict( zip( list(G.nodes), l ) )

cursed_dict = node_dict( [reaction.id for reaction in model.reactions] )
G_nx = nx.relabel_nodes(G, cursed_dict, copy=True)

# %% --- Exporta el modelo a un formato común

nx.write_graphml(G, "./tmp/graph.graphml", named_key_ids=True)  
