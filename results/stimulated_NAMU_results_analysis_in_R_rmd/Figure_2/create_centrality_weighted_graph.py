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

import networkx as nx
import pandas as pd

path                      = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_2/'


nodes_centrality = pd.read_csv(path + "total_centrality_by_node.csv")





# %% 
reaction_projected_stimnulated = nx.read_gpickle(path  + "reaction_projected_stimnulated.gpickle")

reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = nodes_centrality.ID , \
                                            nombre = "Centrality",   atributos = nodes_centrality.Centrality)


reaction_projected_stimnulated.nodes(data = True)


nx.write_gpickle(reaction_projected_stimnulated, path + "/reaction_projected_stimnulated.gpickle")

nx.write_gexf(reaction_projected_stimnulated,  path + "/reaction_projected_stimnulated.gexf")