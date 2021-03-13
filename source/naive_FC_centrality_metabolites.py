# %%
from   heapq import merge
import time
import cobra
from networkx.drawing.layout import bipartite_layout  
import numpy  as np
import pandas as pd
from   cobra.util.array import create_stoichiometric_matrix
import warnings 
warnings.filterwarnings("ignore")
import networkx as nx


t0        = time.perf_counter()

stimulated = cobra.io.load_json_model(
"/home/alejandro/PostDoc/human-metnet/data/GEM_Recon2_thermocurated_redHUMAN.json")

S_matrix = create_stoichiometric_matrix(stimulated)
S_matrix = (abs(S_matrix) )
S_matrix = (S_matrix > 0.0).astype(np.int_)

projected_S_matrix = np.matmul(S_matrix, S_matrix.T)
np.fill_diagonal(projected_S_matrix, 0) 

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)

G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )
node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
cursed_dict = node_dict( [met.id for met in stimulated.metabolites] )

G = nx.relabel_nodes(G, cursed_dict, copy=True)

largest_component = max(nx.connected_components(G), key=len)
G = G.subgraph(largest_component)
print(len(cursed_dict), len(G.nodes))

print(nx.is_connected(G), list(G.nodes)[0])


t1 = time.perf_counter()
print(f'Finished in {(t1-t0)/60} minutes')

# %%
graph = G.copy()

t0    = time.perf_counter()
dc    = nx.degree_centrality(graph)
t1    = time.perf_counter()
print(f'degree_centrality in {(t1-t0)/60} minutes')

t0    = time.perf_counter()
hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
t1    = time.perf_counter()
print(f'harmonic_centrality in {(t1-t0)/60} minutes')


t0    = time.perf_counter()
ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
t1    = time.perf_counter()
print(f'eigenvector_centrality in {(t1-t0)/60} minutes')

t0    = time.perf_counter()
cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
t1    = time.perf_counter()
print(f'closeness_centrality in {(t1-t0)/60} minutes')


t0    = time.perf_counter()
bc    = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
t1    = time.perf_counter()
print(f'betweenness_centrality in {(t1-t0)/60} minutes')

t0    = time.perf_counter()
lc = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
t1    = time.perf_counter()
print(f'load_centrality in {(t1-t0)/60} minutes')

t0    = time.perf_counter()
ic      = nx.information_centrality(graph)
t1    = time.perf_counter()
print(f'information_centrality in {(t1-t0)/60} minutes')


t0      = time.perf_counter()
cfcc    = nx.current_flow_closeness_centrality(graph)    
t1    = time.perf_counter()
print(f'current_flow_closeness_centrality in {(t1-t0)/60} minutes')

t0    = time.perf_counter()
acfbc   = nx.approximate_current_flow_betweenness_centrality(graph, solver = "cg", kmax=1500, epsilon=0.8)
t1    = time.perf_counter()
print(f'approximate_current_flow_betweenness_centrality in {(t1-t0)/60} minutes')

#t0    = time.perf_counter()
#soc     = nx.second_order_centrality(graph) 
#t1    = time.perf_counter()
#print(f'second_order_centrality in {(t1-t0)/60} minutes')
#Complexity of this implementation, made to run locally on a single machine, is O(n^3), with n the size of G, 
# which makes it viable only for small graphs.

# %% Too heavy, intractable so far.

t0    = time.perf_counter()
cfbc    = nx.current_flow_betweenness_centrality(graph, dtype = float , solver = 'cg')
t1    = time.perf_counter()
print(f'current_flow_betweenness_centrality in {(t1-t0)/60} minutes')


t0    = time.perf_counter()
cbc     = nx.communicability_betweenness_centrality(graph)  
t1    = time.perf_counter()
print(f'communicability_betweenness_centrality in {(t1-t0)/60} minutes')

# %%
def light_centralities(graph):
    hc = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc = nx.degree_centrality(graph)
    bc = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
    cc = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    lc = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
    ic      = nx.information_centrality(graph)
    cfcc    = nx.current_flow_closeness_centrality(graph)    
    acfbc   = nx.approximate_current_flow_betweenness_centrality(graph, solver = "cg", kmax=3500, epsilon=0.5)

    #####-- make DFs --###################################################################

    cfcc_df = pd.DataFrame.from_dict(cfcc, columns = ['cfcc'], orient='index')
    acfbc_df = pd.DataFrame.from_dict(acfbc, columns = ['acfbc'], orient='index')
    hc_df = pd.DataFrame.from_dict(hc, columns = ['hc'], orient='index')
    ec_df = pd.DataFrame.from_dict(ec, columns = ['ec'], orient='index')
    dc_df = pd.DataFrame.from_dict(dc, columns = ['dc'], orient='index')
    bc_df = pd.DataFrame.from_dict(bc, columns = ['bc'], orient='index')
    cc_df = pd.DataFrame.from_dict(cc, columns = ['cc'], orient='index')
    lc_df = pd.DataFrame.from_dict(lc, columns = ['lc'], orient='index')
    ic_df   = pd.DataFrame.from_dict(ic, columns = ['ic'], orient='index')
    #soc_df  = pd.DataFrame.from_dict(soc, columns = ['soc'], orient='index')   

    from functools import reduce
    dfs = [hc_df, ec_df, dc_df, bc_df, cc_df, lc_df, ic_df, cfcc_df, acfbc_df]

    df_all_centralities = reduce(lambda  left, right: left.join(right, how='outer'), dfs)
    # Lambda que reduce dataframes y hace un join
    return df_all_centralities


# %% Get baseline centralities 

G_unperturbed    = G.copy()

t0        = time.perf_counter()

df_light_baseline_centralities = light_centralities(G_unperturbed)

t1 = time.perf_counter()
print(f'Finished in {(t1-t0)/60} minutes')

# %%
df_light_baseline_centralities.to_csv('df_light_baseline_centralities.csv')
# %% Bipartite graph 

def cobra_to_bipartite_graph(modelo, direccionado=True):
    """Toma un modelo de cobra y genera un grafo bipartito de NetworkX

    Parameters
    ----------
    model : cobra.core.model.Model
        Un modelo no-optimizado de reacciones metabolicas en COBRA
    direccionado : bool
        Selecciona si el objeto creado sera un grafo direccionado o no.
        Por defecto crea grafos direccionados. 

    Returns
    -------
    grafo : 
        un grafo bipartito de NetworkX. Incluye atributos de flujo del modelo y
        sensibilidad de metabolitos y reacciones. 
    
    Notes
    -----
    - nx.write_gexf(grafo, "grafo.gexf") Crea una salida para Gephi
    """

    from cobra.util import create_stoichiometric_matrix
    from networkx import relabel_nodes, nodes
    from networkx.algorithms.bipartite.matrix import biadjacency_matrix, from_biadjacency_matrix
    from scipy.sparse import csr_matrix
    import numpy as np

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"

    grafo = abs(create_stoichiometric_matrix(modelo)) # Crea una matriz
    grafo = (grafo > 0.0).astype(np.int_) # Binarización
    grafo = csr_matrix(grafo) # Convierte la matriz a una matriz dispersa
    grafo = from_biadjacency_matrix(grafo) # Usa la dispersa para un grafo bipartito
    
    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    metabolites_n     = len(modelo.metabolites) # Numero de metabolitos
    metabolites_names = [modelo.metabolites[i].id for i in range(0, metabolites_n) ]

    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
    reactions_n       = len(modelo.reactions)   # Numero de reacciones
    reactions_names   = [modelo.reactions[i].id   for i in range(0, reactions_n)   ]

    names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
    grafo = relabel_nodes(grafo, names_mapped)
    return grafo

# %% Bipartite network
bipartite_recon =  cobra_to_bipartite_graph(stimulated, direccionado = False)

len(bipartite_recon.nodes)


G = bipartite_recon.copy()

largest_component = max(nx.connected_components(G), key=len)
G = G.subgraph(largest_component)


print( len(G.nodes), nx.is_connected(G), list(G.nodes)[100:105])
# %% ç

G_bipartite_unperturbed    = G.copy()

t0        = time.perf_counter()

df_light_baseline_centralities = light_centralities(G_bipartite_unperturbed)

t1 = time.perf_counter()
print(f'Finished in {(t1-t0)/60} minutes')

# %%
df_light_baseline_centralities
# %%
