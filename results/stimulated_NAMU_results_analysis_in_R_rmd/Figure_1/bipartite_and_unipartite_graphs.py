
# %% 
def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G

def cobra_to_bipartite(model):
    import networkx as nx
    from   cobra.util.array import create_stoichiometric_matrix
    import numpy as np
    from sklearn.preprocessing import Binarizer
    import warnings 
    from scipy.sparse import csr_matrix
    from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix


    warnings.filterwarnings("ignore")
    #extraer matriz estequiométrica
    S_matrix = create_stoichiometric_matrix(model)
    #convertir todas las entradas valores positivos
    S_matrix = abs(S_matrix)
    #binarizar
    S_matrix = Binarizer().fit_transform(S_matrix)

    S_matrix = csr_matrix(S_matrix) # Convierte la matriz a una matriz dispersa
    grafo = from_biadjacency_matrix(S_matrix) # Usa la dispersa para un grafo bipartito
    
    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    metabolites_n     = len(model.metabolites) # Numero de metabolitos
    metabolites_names =    [model.metabolites[i].id for i in range(0, metabolites_n) ]

    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
    reactions_n       = len(model.reactions)   # Numero de reacciones
    reactions_names   =    [model.reactions[i].id   for i in range(0, reactions_n)   ]

    names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
    grafo = nx.relabel_nodes(grafo, names_mapped)
   # grafo = get_largest_component(grafo)
    return grafo


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

def solution2attr(solucion, grafo, estandarizar=False, umbral=1e-7):
    """Docstring

    Parameters
    ----------
    solucion : 
    grafo : 
    estandarizar : bool, default False
        Define si se aplican estandarizaciones al modelo final. Estas son entre
    umbral : float
        El umbral en que un flujo o sensibilidad se consideran 0. 
    """

    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones

    print("Metabolitos:", len(metabolites_nodes), " | Reacciones:", len(reactions_nodes)) # Debug

    from numpy import abs # Permite absolutos vectorizados

    shadow_prices = abs(solucion.shadow_prices).tolist() # De solucion, pasa a lista
    fluxes        = abs(solucion.fluxes).tolist()       # posiblemente es más rapido de usar que desde el modelo
    reduced_costs = abs(solucion.reduced_costs).tolist() # considerando el menor parsing pero si requiere pandas

    print("Shadow prices:", len(shadow_prices),"| Fluxes:", len(fluxes),"| Reduced costs:", len(reduced_costs)) # Debug

    # Neutraliza sub-umbral (por defecto 0.0000001)
    shadow_prices = [ 0 if i < umbral else i for i in shadow_prices]
    fluxes        = [ 0 if i < umbral else i for i in fluxes]
    reduced_costs = [ 0 if i < umbral else i for i in reduced_costs]

    if estandarizar == True:
        def estandariza(tmp):
            from numpy import array, log10, inf, min, max
            tmp = log10(array(tmp))
            tmp[tmp == -inf] = tmp[tmp != -inf].min()
            tmp = (tmp - tmp.min() ) / ( tmp.max() - tmp.min() )
            return tmp
        
        shadow_prices = estandariza(shadow_prices)
        fluxes = estandariza(fluxes)
        reduced_costs = estandariza(reduced_costs)

    from networkx import set_node_attributes

    # ASIGNA Shadow_prices, fluxes, reduced_costs
    set_node_attributes(grafo, { metabolites_nodes[i] : shadow_prices[i] for i in range(0, len(shadow_prices)) } , "Shadow Price") 
    set_node_attributes(grafo, { reactions_nodes[i] : fluxes[i] for i in range(0, len(fluxes)) } , "Flux") 
    set_node_attributes(grafo, { reactions_nodes[i] : reduced_costs[i] for i in range(0, len(reduced_costs)) } , "Reduced cost") 

    #grafo.nodes(data=True) # Debug

    return grafo
# %% 
import cobra
import warnings 
import pandas as pd
warnings.filterwarnings("ignore")


path          = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_1'
stimulated    = cobra.io.load_json_model(path + '/stimulated_2021.json') 
nodes_attr_df = pd.read_csv(path + '/node_attributes.csv')
fba_solution        = stimulated.optimize()
# %%
fba_solution        = stimulated.optimize()
stimulated_graph = cobra_to_bipartite(stimulated)
stimulated_graph    = solution2attr(fba_solution, stimulated_graph, estandarizar=True)
rxns   =  [n for n, d in stimulated_graph.nodes(data=True) if d["bipartite"] == 1] 
stimulated_graph  = list2attr(stimulated_graph, nodos = nodes_attr_df.ID , nombre = "type", atributos = nodes_attr_df.Node)
#stimulated_graph.nodes(data=True)


stimulated_graph = get_largest_component(stimulated_graph)
# %%
import networkx as nx

#nx.write_gexf(stimulated_graph, "stimulated_bipartite_graph.gexf")


# %%
reaction_projected_stimnulated = cobra_to_networkx_rxn_projection(stimulated)

#reaction_projected_stimnulated    = solution2attr(fba_solution, reaction_projected_stimnulated, estandarizar=True)
#rxns   =  [n for n, d in reaction_projected_stimnulated.nodes(data=True) if d["bipartite"] == 1] 

fluxes        = fba_solution.fluxes.to_frame()

reduced_costs = fba_solution.reduced_costs.to_frame()

reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = fluxes.index , \
                                  nombre = "fluxes", atributos = fluxes.fluxes)

reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = reduced_costs.index , \
                                  nombre = "reduced_costs", atributos = reduced_costs.reduced_costs)

reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = nodes_attr_df.ID , \
    nombre = "type", atributos = nodes_attr_df.Node)

reaction_projected_stimnulated = get_largest_component(reaction_projected_stimnulated)


nx.write_gpickle(reaction_projected_stimnulated, path + "/reaction_projected_stimnulated.gpickle")

# %%
