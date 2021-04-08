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

path_fig_2                      = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_2/'
path_fig_4                      = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_4/'


reaction_projected_stimnulated = nx.read_gpickle(path_fig_2  + "reaction_projected_stimnulated.gpickle")

nodes_modules = pd.read_csv(path_fig_4 + "modules_for_graph.csv")

nodes_modules
# %% 
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = nodes_modules.ID , \
                                            nombre = "Module",   atributos = nodes_modules.Module)


reaction_projected_stimnulated.nodes(data = True)



nx.write_gpickle(reaction_projected_stimnulated, path_fig_4 + "/reaction_projected_modules.gpickle")

nx.write_gexf(reaction_projected_stimnulated,  path_fig_4 + "/reaction_projected_modules.gexf")