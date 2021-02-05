#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades

Parametros: MODELO, SUBSISTEMAS
Output: CENTRALIDADES.tsv"""


# %% --- Importa el modelo
from cobra.io import load_json_model, save_json_model

INPUT = 'toy_metabolism_AA.json'
model = load_json_model(INPUT)

solution_fba = model.optimize() # Optimización del modelo para check
solution_fba.fluxes # Check si los flujos funcionan

# %% --- Proyección de la matrix estequimetrica unipartita de reacciones 
"""Usamos reacciones para así poder buscar genes. 
1. Crea una matriz estequimetrica que tiene las reacciones del modelo
2. Binarización de matriz de incidencia para eliminar los pesos y dirección de 
   las reacciones. Los coeficientes son un artefacto, cuando solo nos interesa
   la conectividad
3. Proyección de la matriz de incidencia hacia el espacio de las reacciones, lo
   que genera en matriz de adyacencia de reacciones
4. Ceros a la diagonal de la matriz, para eliminar auto-conexiones
4. Binarización, de nuevo. Solo nos interesa saber si una reaccion comparte 
   pathway con otra, no cuantos metabolitos comparten"""

import numpy as np
import networkx as nx
from   cobra.util.array import create_stoichiometric_matrix

S_matrix = create_stoichiometric_matrix(model) # Crea matriz de incidencia

S_matrix = (S_matrix !=0).astype(int) # Binarización de la matriz de incidencia 
  # para eliminar pesos y dirección

projected_S_matrix = np.matmul(S_matrix.T, S_matrix) # Proyección al espacio de 
  # las reacciones de la matriz en matriz de incidencia

"""La diagonal contiene la cantidad de metablitos que estan en cada reacción.
Si una reacción tiene 1 podemos afirmar que es una reacción de Boundary.  
np.matmul(S_matrix, S_matrix.T) # Genera la matriz de metabolitos
La diagonal de esta indicaria en cuantas reacicones participa cada metabolito"""

np.fill_diagonal(projected_S_matrix, 0) # Elimina la conexion consigo mismo

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int) # Re-binariza, 
    # para eliminar el peso al compartir más de un metabolito entre reacciones
reaction_adjacency_matrix.shape

# %% --- Crea el grafo a partir de la matriz de adyacencia
G = nx.from_numpy_matrix(reaction_adjacency_matrix)

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("This graph is not conected. Picking largest component.")
    largest_component = max(nx.connected_components(G), key=len)
    print ('Removed',len(G.nodes) - len(largest_component), 'nodes.')
    G = G.subgraph(largest_component)
    
G.nodes

# %% --- Reasigna nombres a nodos

from networkx import set_node_attributes
set_node_attributes(G, [reaction.id for reaction in model.reactions], 'name' )

# %% --- Definición de subsistemas de interes

subsystems_1 = []
subsystems_2 = []

# %% --- Calculos de centralidades varios
"""Calculos de centralidades baseline,  
Existen 30+ centralidades, las cuales no aplican en todas las condiciones"""

graph = G.copy()

hc  = nx.harmonic_centrality(graph, nbunch=None, distance=None)
ec  = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
dc  = nx.degree_centrality(graph)
bc  = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
"""↓ Here be dragons ↓"""
cc  = nx.closeness_centrality(graph, distance=None, wf_improved=True)
lc  = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
ic  = nx.information_centrality(graph)
soc = nx.second_order_centrality(graph)