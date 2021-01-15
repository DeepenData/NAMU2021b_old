"""
Creación de un modelo para cobra y futuro analiis de redes  
[Documentación de COBRA](https://cobrapy.readthedocs.io/en/0.17.0/building_model.html)

Toy reactions:

r1: m1 <- 
r2: m1 -> m2
r3: m2 <-
r4: m2 -> 
r5: m2 -> m3
r6: m1 + m5 -> m3 + m4
r7: m3 + m4 -> m1 + m5
r8: m3 + m5 -> m4 
maximize r4 + r8 =    
obj:     m2 + m3 + m5 -> m4

subject to r1 > 1

Agregar el objetivo 

# change the objective to ATPM
model.objective = "obj"

# The upper bound should be 1000, so that we get the actual optimal value
model.reactions.get_by_id("obj").upper_bound = 1000.0
#chequear
linear_reaction_coefficients(model)

#Agregar el consumo de sustrato

model.reactions.get_by_id("r1").bounds = (-1.33, 0.0)

"""
# %% --- Inicialización del modelo
from cobra import Model, Reaction, Metabolite
#del model, r1, r2, r3, r4, r5, r6, r7, r8, solution

model = Model('toy_fba')
#
# --- Creación de metabolitos
m1 = Metabolite("m1")
m2 = Metabolite("m2")
m3 = Metabolite("m3")
m4 = Metabolite("m4")
m5 = Metabolite("m5")

#  --- Creación de reacciones 
r1 = Reaction("r1")
r1.add_metabolites({m1 : 2.0})
r2 = Reaction("r2")
r2.add_metabolites({m1: -1.0, m2: 1.0})
r3 = Reaction("r3")
r3.add_metabolites({m2 : 1.0})
r4 = Reaction("r4")
r4.add_metabolites({m2 : -1.0})
r5 = Reaction("r5")
r5.add_metabolites({m2: -1.0, m3: 1.0})
r6 = Reaction("r6")
r6.add_metabolites({m1: -1.0, m5: -1.0, m3: 2.0, m4: 1.0})
r7 = Reaction("r7")
r7.add_metabolites({m3: -1.0, m4: -1.0, m1: 1.0, m5:1.0})
r8= Reaction("r8")
r8.add_metabolites({m3: -1, m5: -1, m4: 2.0})
obj = Reaction("obj") # Función objetivo
obj.add_metabolites({m2: -1.5, m3: -1.5, m5: -2.5, m4: 1.0}) 
model.add_reactions([r1, r2, r3, r4, r5, r6, r7, r8, obj])
model.objective = "obj"

## Optimización del modelo
solution = model.optimize()
solution = solution.to_frame()
all_reactions = [model.reactions[i].reaction  for i in range(0,len(model.reactions)) ]
solution["all_reactions"] = all_reactions

solution
### === TERMINA EL TOY ===
# %% MM --- 2021-01-15 12:11 --- Exporta un grafo
def cobra2networkx(modelo):
    """Toma un modelo de cobra y genera un grafo bipartito de NetworkX

    Parameters
    ----------
    modelo : cobra_model
        Un modelo de reacciones metabolicas en COBRA
    
    Returns
    -------
    grafo : 
        un grafo bipartito de NetworkX. Incluye atributos de flujo del modelo y
        sensibilidad de metabolitos y reacciones. 
    
    Notes
    -----
        - nx.write_gexf(grafo, "grafo.gexf") Crea una salida para Gephi
    """
    import networkx as nx
    from cobra.util import create_stoichiometric_matrix

    from networkx.algorithms.bipartite.matrix import biadjacency_matrix      # Extrae matriz de adyacencia
    from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix # Crea el grafo
    from sklearn.preprocessing import binarize
    from scipy.sparse import csr_matrix

    tmp = binarize(create_stoichiometric_matrix(model)) # Crea una matriz
    tmp = csr_matrix(tmp) # Convierte la matriz a una matriz dispersa
    tmp = from_biadjacency_matrix(tmp) # Usa la dispersa para un grafo bipartito

    # total_nodes = range(0 , len(tmp.nodes(data=False))) # Rango de nodos
    # particiones = [tmp.nodes(data=True)[i]["bipartite"] for i in total_nodes] # Tipo nodo

    metabolites_nodes = [n for n, d in tmp.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    metabolites_n     = len(model.metabolites) # Numero de metabolitos
    metabolites_names = [model.metabolites[i].id for i in range(0, metabolites_n) ]

    reactions_nodes   = [n for n, d in tmp.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
    reactions_n       = len(model.reactions)   # Numero de reacciones
    reactions_names   = [model.reactions[i].id   for i in range(0, reactions_n)   ]

    names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
    tmp = nx.relabel_nodes(tmp, names_mapped)
    # En este punto el grafo ya cuenta con los nodos y labels segun IDs

    solution = model.optimize()    # Esto posiblemente no deberia estar aqui
    solution = solution.to_frame() # Ni esto tampoco

    # Crea diccionarios con atributos para el grafo
    nodos = metabolites_nodes + reactions_nodes
    # Diccionario para precios (sensibilidad del metabolito) y nx.asignación
    tmp_prices = solution.shadow_prices.tolist() + (['']*reactions_n) # Shadow Prices
    tmp_prices = { nodos[i] : (tmp_prices)[i] for i in range(0, len(tmp_prices)) } 
    nx.set_node_attributes(tmp, tmp_prices, 'prices' ) # Añade una columna de precios
    
    # Diccionario para flujos (flujo de la reaccion) y nx.asignación
    tmp_fluxes = ['']*metabolites_n + solution.fluxes.tolist()  # Flujos de reacciones
    tmp_fluxes = { nodos[i] : (tmp_fluxes)[i] for i in range(0, len(tmp_prices)) } 
    nx.set_node_attributes(tmp, tmp_fluxes, 'fluxes' ) # Añade una columna de flujos

    grafo = tmp # Asigna tmp al grafo
    return grafo

# %% MM --- 2021-01-15 12:10 --- Función toma dos listas y añade atributos
def list2attr(grafo, nodos, atributos, nombre):
    """Toma dos listas de nombres de nodos y atributos y las añade a un grafo

    Parameters
    ----------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    nodos: list
        Una lista de nombres de nodos
    atributos : list
        Una lista de valores para los nodos en la lista anterior.
        Ambas listas deben ser del mismo largo. 
    nombre : str
        Nombre del atributo siendo asignado. Ie. nombre de la columna en Gephi
    
    Returns
    -------
    grafo: bipartite_graph
        Un grafo bipartito con un nuevo atributo para un set de nodos "nodos". 
    """
    tmp_list = { nodos[i] : (atributos)[i] for i in range(0, len(atributos)) }
    nx.set_node_attributes(grafo, tmp_list, nombre ) # Añade los atributos
    return grafo

# %% AA 15_ene_2021 --- Stoichiometric matrix to directed bipartite graph
from sklearn.preprocessing import binarize
from scipy.sparse import csr_matrix
import networkx as nx
from networkx.algorithms.bipartite.matrix import biadjacency_matrix      # Extrae matriz de adyacencia
from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix # Crea el grafo
from cobra.util import create_stoichiometric_matrix
import numpy as np
import pandas as pd

#graph   = nx.read_gpickle(" .gpickle")

stoichiometric_matrix        =  create_stoichiometric_matrix(model)
stoichiometric_matrix        =  np.where(stoichiometric_matrix < 0, -1, stoichiometric_matrix)
stoichiometric_matrix        =  np.where(stoichiometric_matrix > 0,  1, stoichiometric_matrix)


sparse_stoichiometric_matrix = csr_matrix(stoichiometric_matrix)
directed_bipartite_graph     = from_biadjacency_matrix(sparse_stoichiometric_matrix, create_using= nx.DiGraph) 
incidence_matrix = biadjacency_matrix(directed_bipartite_graph, row_order = range(0,stoichiometric_matrix.shape[0]))

print((incidence_matrix == stoichiometric_matrix).all(),directed_bipartite_graph.is_directed())


reactions_nodes   = [n for n, d in directed_bipartite_graph.nodes(data=True) if d["bipartite"] == 1] # 
metabolites_nodes = [n for n, d in directed_bipartite_graph.nodes(data=True) if d["bipartite"] == 0] # 

# --- Poniendole nombres a metabolitos y reacciones
reaction_names    = ["r1","r2","r3","r4","r5","r6","r7","r8",'obj']
metabolites_names = ["m1","m2","m3","m4","m5"]
# 
names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reaction_names  ))
    # el formato + permite uni listas; como vectores en R; y zip crea tuplas

directed_bipartite_graph_labeled_nodes = nx.relabel_nodes(directed_bipartite_graph, names_mapped)
directed_bipartite_graph_labeled_nodes.nodes(data=True) 

node_labels = list(directed_bipartite_graph_labeled_nodes.nodes)

# %% --- Directed centralities


G    = directed_bipartite_graph.copy()
#calculo de centralidades 
dc   = nx.degree_centrality(G)
indc = nx.in_degree_centrality(G)
oudc = nx.out_degree_centrality(G)
pr   = nx.pagerank(G)
kc   = nx.katz_centrality(G)
#ec   = nx.eigenvector_centrality(G)
cc   = nx.closeness_centrality(G)
bc   = nx.betweenness_centrality(G)
lc   = nx.load_centrality(G)
hc   = nx.harmonic_centrality(G)
#cfcc = nx.current_flow_closeness_centrality(G)
#ic   = nx.information_centrality(G)
#cfbc = nx.current_flow_betweenness_centrality(G)
#cbc  = nx.communicability_betweenness_centrality(G)
#soc  = nx.second_order_centrality(G)

my_centralities = "[dc, indc, oudc, pr, kc, cc, bc, lc, hc]"


def mergeDict(dict1, dict2):
   """Une dos diccioanrios y conserva los values de keys coindicentes
   
   Parameters
   ----------
    dict1 : dict
        Diccionario
    dict1 : dict
        Diccionario
   """ 
   dict3 = {**dict1, **dict2}
   for key, value in dict3.items():
       if key in dict1 and key in dict2:
               dict3[key] = [value , dict1[key]]
   return dict3 


import functools
from   iteration_utilities import deepflatten


merged_dicts           = functools.reduce(mergeDict, eval(my_centralities))
merged_centralities    = [[i for i in deepflatten(merged_dicts[i])] for i in range(0,len(merged_dicts))]
nombres_cols           = my_centralities.strip('][').split(',')
nombres_cols.reverse()

merged_directed_centralities_df = pd.DataFrame(merged_centralities, index= node_labels, columns= nombres_cols)

# %% Undirected centralities 
UNdirected_stoichiometric_matrix        =  binarize(abs(create_stoichiometric_matrix(model)))

sparse_undirected_stoichiometric_matrix = csr_matrix(UNdirected_stoichiometric_matrix)
UNdirected_bipartite_graph              = from_biadjacency_matrix(sparse_undirected_stoichiometric_matrix) 
UNdirected_incidence_matrix = biadjacency_matrix(UNdirected_bipartite_graph, row_order = range(0,UNdirected_stoichiometric_matrix.shape[0])) 
print((UNdirected_incidence_matrix == UNdirected_stoichiometric_matrix).all(),UNdirected_bipartite_graph.is_directed())

G2    = UNdirected_bipartite_graph.copy()
#calculo de centralidades 
undirected_dc  = nx.degree_centrality(G2)
undirected_pr   = nx.pagerank(G2)
undirected_kc   = nx.katz_centrality(G2)
undirected_ec   = nx.eigenvector_centrality(G2)
undirected_cc   = nx.closeness_centrality(G2)
undirected_bc   = nx.betweenness_centrality(G2)
undirected_lc   = nx.load_centrality(G2)
undirected_hc   = nx.harmonic_centrality(G2)
undirected_cfcc = nx.current_flow_closeness_centrality(G2)
undirected_ic   = nx.information_centrality(G2)
undirected_cfbc = nx.current_flow_betweenness_centrality(G2)
undirected_cbc  = nx.communicability_betweenness_centrality(G2)
undirected_soc  = nx.second_order_centrality(G2)

my_UNDIRECTED_centralities = \
"[undirected_dc, undirected_pr, undirected_kc, undirected_ec, undirected_cc, undirected_bc, undirected_lc, undirected_hc,\
     undirected_cfcc, undirected_ic, undirected_cfbc, undirected_cbc, undirected_soc]"

merged_dicts           = functools.reduce(mergeDict, eval(my_UNDIRECTED_centralities))
merged_centralities    = [[i for i in deepflatten(merged_dicts[i])] for i in range(0,len(merged_dicts))]
nombres_cols           = my_UNDIRECTED_centralities.strip('][').split(',')
nombres_cols.reverse()
merged_UNdirected_centralities_df = pd.DataFrame(merged_centralities, index= node_labels, columns= nombres_cols)

# %% Merge directed and undirected centralities
(merged_UNdirected_centralities_df["undirected_dc"] == merged_directed_centralities_df["dc"]).all()
all_centralities_DF = merged_directed_centralities_df.join(merged_UNdirected_centralities_df)
#RESULTADO IMPORTANTE 1
all_centralities_DF
# %%
model.optimize()
# %% SAMPLING THE STEADY-STATE
from cobra.sampling import sample
import seaborn as sns
from cobra.sampling.optgp import OptGPSampler

flux_samples = OptGPSampler(model,  processes=16, thinning=500, nproj=100, seed=23)
flux_samples = flux_samples.sample(n = 1000)
#Se pueden hacer más análisis estadísticos!
samples_means    = flux_samples.T.mean(axis=1)
samples_means_DF = pd.DataFrame(samples_means, columns=["samples_mean"])
#RESULTADO IMPORTANTE 2
samples_means_DF

# %%
#sns.pairplot(flux_samples, corner=True)