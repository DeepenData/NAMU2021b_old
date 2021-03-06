# %%
from   heapq import merge
import cobra
from networkx.algorithms.shortest_paths import weighted  
import numpy  as np
import pandas as pd
from   cobra.util.array import create_stoichiometric_matrix
import warnings 
warnings.filterwarnings("ignore")
import networkx as nx

stimulated = cobra.io.load_json_model(
    "/home/alejandro/PostDoc/human-metnet/data/stimulated_2021.json")

S_matrix = create_stoichiometric_matrix(stimulated)
S_matrix = (abs(S_matrix) )
S_matrix = (S_matrix > 0.0).astype(np.int_)

projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
np.fill_diagonal(projected_S_matrix, 0) 

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)
G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )
node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
ids_dict = node_dict( [reaction.id for reaction in stimulated.reactions] )
print(len(ids_dict), len(G.nodes))
G = nx.relabel_nodes(G, ids_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
largest_component = max(nx.connected_components(G), key=len)
G = G.subgraph(largest_component)
print(len(ids_dict), len(G.nodes))
nx.is_connected(G) # Reemplazar por llamadas de función

# %%
def list2attr(grafo, nodos, nombre, atributos):
    """Toma dos listas: nombres de nodos y atributos; y las añade a un grafo

    Parameters
    ----------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    nodos: list
        Una lista de nombres de nodos
    nombre : str
        Nombre del atributo siendo asignado. Ie. nombre de la columna en Gephi
    atributos : list
        Una lista de valores para los nodos en la lista anterior.
        Ambas listas deben ser del mismo largo. 
    
    Returns
    -------
    grafo: bipartite_graph
        Un grafo bipartito con un nuevo atributo para un set de nodos "nodos". 
    """
    assert len(nodos) == len(atributos), "Ambas listas deben ser del mismo largo."
    tmp_list = { nodos[i] : atributos[i] for i in range(0, len(atributos)) }
    
    from networkx import set_node_attributes
    set_node_attributes(grafo, tmp_list, nombre ) # Añade los atributos
    
    return grafo


# %%
#Assign attributes

Results = pd.read_csv('Results.csv')
Ce      = abs(Results.iloc[:,7:]).sum(axis=1)
Ce_IDs  = Results.ID

weighted_centrality_graph = list2attr(grafo = G, nodos = Ce_IDs, nombre= 'Ce',  atributos = Ce)


data_Networkx = pd.read_csv('data_Networkx.csv')
node_type_IDs   = data_Networkx.ID
node_type_types = data_Networkx.Node_type

weighted_centrality_graph = list2attr(
grafo = weighted_centrality_graph, nodos = node_type_IDs, nombre= 'node_type',  atributos = node_type_types)

weighted_centrality_graph.nodes(data = True)


nx.write_gexf(weighted_centrality_graph, "weighted_centrality_graph.gexf")
