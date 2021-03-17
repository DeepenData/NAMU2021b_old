#!/bin/python
"""Genera un grafo a partir de un modelo de COBRA en formato JSON. Lo exporta en
un formato apto para NetworkX y Graph-tool"""

# %% --- IMPORTA EL MODELO

INPUT_MODEL='./data/stimulated_2021.json' # TODO: hacer esto una variable ambiental

import warnings; warnings.filterwarnings('ignore') # Ignora warnings
from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL)

# %% --- OPTIMIZA EL MODELO, PARA FLUJOS Y ESO

# print("Iniciando optimización del modelo", time.asctime( time.localtime(time.time()) ) )
# solution_fba = model.optimize() # Optimización del modelo para check
# solution_fba.fluxes # Check si los flujos funcionan

# %% --- CONVIERTE EL MODELO A UN GRAFO

import networkx as nx
import numpy    as np
from cobra.util.array import create_stoichiometric_matrix

S_matrix = create_stoichiometric_matrix(model)
S_matrix = (abs(S_matrix) ) # Convertir todas las entradas valores positivos
S_matrix = S_matrix.astype(np.int) # Transformar a enteros

projected_S_matrix = np.matmul(S_matrix.T, S_matrix) # Multiplicacion por la derecha para proyectar en el espacio de las reacciones
print('Dimensiones de la matriz:', projected_S_matrix.shape,    # Dimensiones, deberian ser cuadradas
      'Valores no-cero:', np.count_nonzero(projected_S_matrix)) # Valores no-cero, deberian mantenerse

print('Removiendo la diagonal...')
np.fill_diagonal(projected_S_matrix, 0) # Llena las diagonales de ceros. 

projected_S_matrix = (projected_S_matrix !=0).astype(int) # No require importar SKLearn
#from sklearn.preprocessing import Binarizer
#projected_S_matrix = Binarizer().fit_transform(projected_S_matrix) # Binarización
print('Dimensiones de la matriz:', projected_S_matrix.shape,    # Dimensiones, deberian ser cuadradas
      'Valores no-cero:', np.count_nonzero(projected_S_matrix), # Valores no-cero, deberian mantenerse
      'Valores unicos:',  np.unique(projected_S_matrix))        # Deberian ser [0 1]

G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix ) # Convierte la matriz a numpy

# %% --- AÑADE IDS DE LAS REACCIONES COMO NOMBRE DE LOS NODOS

node_dict = lambda l : dict( zip( list(G.nodes), l ) )
cursed_dict = node_dict( [reaction.id for reaction in model.reactions] )

G_nx = nx.relabel_nodes(G, cursed_dict, copy=True)

# %% --- CREA UN GRAFO TOTALMENTE CONEXO ELIMINANDO NODOS HUERFANOS
# Simplifica algunas centralidades, y los huerfanos usualmente son sinks del
# modelo para que este fluya correctamente. 

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("Modelo no completamente conexo. Eligiendo componente más grande.")
    
    largest_component = max(nx.connected_components(G), key=len)
    
    print( len(G.nodes) - len(largest_component), 'nodos removidos.')
    print( len(G.nodes), "nodes finales.")

    G = G.subgraph(largest_component)


# %% --- EXPORTA EL MODELO A UN FORMATO COMÚN

nx.write_graphml(G, "./tmp/graph.graphml", named_key_ids=True)  
