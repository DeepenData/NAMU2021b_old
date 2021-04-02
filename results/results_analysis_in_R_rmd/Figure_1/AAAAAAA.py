
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
    grafo = get_largest_component(grafo)
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
    G = get_largest_component(G)
    return G


# %% 
import cobra
import warnings 
import pandas as pd
warnings.filterwarnings("ignore")

stimulated    = cobra.io.load_json_model("stimulated_2021.json") 
nodes_attr_df = pd.read_csv("node_attributes.csv")

# %%
bipartite = cobra_to_bipartite(stimulated)
bipartite.nodes(data=True)