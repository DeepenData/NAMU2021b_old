# %% --- nombrar nodos

from cobra.io import load_json_model, save_json_model
INPUT = 'toy_metabolism_AA.json'
model = load_json_model(INPUT)
import numpy as np
import networkx as nx
from   cobra.util.array import create_stoichiometric_matrix

S_matrix = create_stoichiometric_matrix(model) # Crea matriz de incidencia

S_matrix = (S_matrix !=0).astype(int) # Binarización de la matriz de incidencia 

projected_S_matrix = np.matmul(S_matrix.T, S_matrix) # Proyección al espacio de 

np.fill_diagonal(projected_S_matrix, 0) # Elimina la conexion consigo mismo

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int) # Re-binariza, 
# %%
from scipy.sparse import csr_matrix

G = nx.from_numpy_matrix(reaction_adjacency_matrix)

G.nodes

# %% --- Reasigna nombres a nodos
from networkx import relabel_nodes 

tmp_maping = dict( zip(
    list(G.nodes) , # Labels antiguos
    [reaction.id for reaction in model.reactions] # Labels nuevos
    ))

nx.relabel_nodes(G, tmp_maping ) # TODO: arregla esto que esta roto
G.nodes(data=False) 

# %%
