"""
"""

# %% MM --- 2021-01-18 16:40 --- Funciones utilitarias
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
    S_matrix[S_matrix < 0] = -1
    S_matrix[S_matrix == 0] = 0
    S_matrix[S_matrix > 0] = 1
    S_matrix = csr_matrix(S_matrix)
    grafo = from_biadjacency_matrix(S_matrix, create_using=nx.DiGraph)
    
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
# %% MM --- 2021-01-18 16:40 --- Importa el modelo JSON
import cobra
import warnings 
warnings.filterwarnings("ignore")
path          = '/home/alejandro/PostDoc/human-metnet/data/'
model    = cobra.io.load_json_model(path + 'GEM_Recon2_thermocurated_redHUMAN.json') 
di_graph = cobra_to_bipartite(model)
import networkx as nx
nx.write_gpickle(di_graph, path + 'recon2_directed_bipartite_graph.gpickle')
# %%
