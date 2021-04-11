"""
"""

# %% MM --- 2021-01-18 16:40 --- Funciones utilitarias
def get_largest_component(grafo): 
    import networkx as nx
    undirected_graph = grafo.to_undirected()
    largest_component = max(nx.connected_components(undirected_graph), key=len)
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
#di_graph = get_largest_component(di_graph)
import networkx as nx
nx.write_graphml_lxml (di_graph, path + 'recon2_directed_bipartite_graph.graphml')

# %%
hola = nx.read_graphml(path + 'recon2_directed_bipartite_graph.graphml')

graph = hola.copy()
dc    = nx.degree_centrality(graph) # SANITY CHECK: Deberia ser similar a todo lo demas
hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
ec    = nx.eigenvector_centrality(graph, max_iter=2000, tol=1e-05, nstart=None, weight=None)
bc    = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
lc    = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
kz    = nx.katz_centrality(graph, alpha = 0.00001)
pr    = nx.pagerank(graph, alpha = 0.00001)
grc   = nx.global_reaching_centrality(graph)


# %%
from networkx.algorithms.bipartite.centrality  import betweenness_centrality as bi_bet_cen
import networkx as nx
path          = '/home/alejandro/PostDoc/human-metnet/data/'

hola = nx.read_graphml(path + 'recon2_directed_bipartite_graph.graphml')
graph = hola.copy()
mets = [n for n, d in graph.nodes(data=True) if d["bipartite"] == 0] 
rxns = [n for n, d in graph.nodes(data=True) if d["bipartite"] == 1] 
bi_bet = bi_bet_cen(graph, nodes = rxns)
#not for directed
#ic    = nx.information_centrality(graph) # Requiere scipy
#cbc   = nx.communicability_betweenness_centrality(graph) 
# %%
import numpy as np
import networkx as nx
from networkx.algorithms.bipartite.generators import configuration_model as bipartite_config_model




#sum(bi_bet.values())
# %%


hola_deg = hola.degree
dic_mets = { key : hola_deg[key] for key in mets }
dic_rxns = { key : hola_deg[key] for key in rxns }




alien = bipartite_config_model(aseq = list(dic_rxns.values()), bseq = list(dic_mets.values()), create_using=nx.Graph)
alien = alien.to_directed()
alien = nx.DiGraph(alien)

mets_alien = [n for n, d in alien.nodes(data=True) if d["bipartite"] == 0] 
rxns_alien = [n for n, d in alien.nodes(data=True) if d["bipartite"] == 1] 

from networkx.algorithms.bipartite.matrix import biadjacency_matrix

alien_biad = biadjacency_matrix(alien,row_order= mets_alien).toarray()
np.unique(alien_biad)

from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix
from scipy.sparse import csr_matrix

alien_biad[alien_biad < 0] 
alien_biad[alien_biad == 0] 
alien_biad[alien_biad > 0] = np.random.choice([-1,1], len(alien_biad[alien_biad > 0]), p=[.5,.5])




alien_2 = from_biadjacency_matrix(csr_matrix(alien_biad), create_using=nx.DiGraph)


print(np.unique(alien_biad),
nx.is_bipartite(alien_2),
nx.is_multigraphical(alien_2),
nx.is_directed(alien_2))

# %%

print(
nx.is_bipartite(alien),
nx.is_multigraphical(alien),
nx.is_directed(alien))


# %%
g1 = bi_bet_cen(alien, nodes = rxns)
print(
sum(g1.values()),
pd.Series(g1).max(), pd.Series(g1).min()) 

# %%
import pandas as pd
g2 = nx.betweenness_centrality(alien)
print(
sum(g2.values()),
pd.Series(g2).max(), pd.Series(g2).min()) 