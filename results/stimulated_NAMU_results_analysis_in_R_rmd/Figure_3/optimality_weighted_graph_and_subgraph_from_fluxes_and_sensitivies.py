# %%
from re import T
import cobra
import warnings 
import pandas as pd
from pandas.core.base import DataError
warnings.filterwarnings("ignore")

def cobra_to_networkx_rxn_projection(modelo):
    import networkx as nx
    from   cobra.util.array import create_stoichiometric_matrix
    import numpy as np
    from sklearn.preprocessing import Binarizer
    import warnings 
    warnings.filterwarnings("ignore")

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"
    #extraer matriz estequiométrica
    S_matrix = create_stoichiometric_matrix(modelo)
    #convertir todas las entradas valores positivos
    S_matrix = (abs(S_matrix) )
    #transformar a enteros 
    S_matrix = S_matrix.astype(np.int)
    #Multiplicacion por la derecha para proyectar en el espacio de las reacciones
    projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
    #rellenar diagonal con ceros
    np.fill_diagonal(projected_S_matrix, 0) 
    #binarizar
    projected_S_matrix = Binarizer().fit_transform(projected_S_matrix)
    #crear grafo networkx
    G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
    #hacer diccionario con los nombres de las reacciones
    node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
    reaction_dict = node_dict( [reaction.id for reaction in modelo.reactions] )
    #Renombrar los nodos usando el diccionario
    G = nx.relabel_nodes(G, reaction_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
    #G = get_largest_component(G)
    return G


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

def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G

import networkx as nx


path_figure_1          = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_1'
path_figure_3          = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_3'


reaction_projected_stimnulated = nx.read_gpickle( path_figure_1 + "/reaction_projected_stimnulated.gpickle")

optimality_values = pd.read_csv(path_figure_3 + '/optimality_values.csv')






reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = optimality_values.ID
 , \
                                  nombre = "optimality", atributos = optimality_values.value)


nx.write_gexf(reaction_projected_stimnulated,  path_figure_3 + "/reaction_projected_optimalities.gexf")

# %%
#  Extract optimal subgraph

stimulated    = cobra.io.load_json_model(path_figure_1 + '/stimulated_2021.json') 
nodes_attr_df = pd.read_csv(path_figure_1 + '/node_attributes.csv')
fba_solution  = stimulated.optimize()
fluxes        = fba_solution.fluxes.to_frame().query("abs(fluxes) > 1e-12")
reduced_costs = fba_solution.reduced_costs.to_frame().query("abs(reduced_costs) > 1e-12")
optimal_node_names = list(fluxes.join(reduced_costs, how='outer').index)

H1 = reaction_projected_stimnulated.subgraph(optimal_node_names)
H2 = get_largest_component(H1)
print(
len(H1.nodes(data = True)), len(H2.nodes))

# %%
optimal_subgraph = H2.copy()



nx.write_gexf(optimal_subgraph, path + "/optimal_subgraph.gexf")