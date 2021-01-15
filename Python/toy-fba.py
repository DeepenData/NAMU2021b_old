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
model = Model('toy_fba')

# %% --- Creación de metabolitos
m1 = Metabolite("m1")
m2 = Metabolite("m2")
m3 = Metabolite("m3")
m4 = Metabolite("m4")
m5 = Metabolite("m5")

""" Best practise: SBML (Systems Biology Markup Language) compliant IDs
reaction = Reaction('3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.    # This is the default
reaction.upper_bound = 1000. # This is the default
reaction.add_metabolites({ # Reacciones estequiometricas
    malACP_c: -1.0,
    h_c: -1.0
})"""
# %% --- Creación de reacciones 
r1 = Reaction("r1")
r1.add_metabolites({m1 : 1.0})

r2 = Reaction("r2")
r2.add_metabolites({m1: -1.0, m2: 1.0})

r3 = Reaction("r3")
r3.add_metabolites({m2 : 1.0})

r4 = Reaction("r4")
r4.add_metabolites({m2 : -1.0})

r5 = Reaction("r5")
r5.add_metabolites({m2: -1.0, m3: 1.0})

r6 = Reaction("r6")
r6.add_metabolites({m1: -1.0, m5: -1.0, m3: 1.0, m4: 1.0})

r7 = Reaction("r7")
r7.add_metabolites({m3: -1.0, m4: -1.0, m1: 1.0, m5:1.0})

r8= Reaction("r8")
r8.add_metabolites({m3: -1.0, m5: -1.0, m4:1.0})

# %% --- Chequear reacciones 
import numpy as np

print(
r3.reaction,
r4.reaction,
)
# %% 
from cobra.util import create_stoichiometric_matrix

model.add_reactions([r1, r2, r3, r4, r5, r6, r7, r8])
model.reactions.get_by_id('r3')

# %% --- Crea la reacción para objetivo de optimización

obj = Reaction("obj")
obj.add_metabolites({m2: -1.0, m3: -1.0, m5: -1.0, m4: 1.0}) 

model.add_reactions([obj])
model.reactions.get_by_id("obj")

# %% 

#Optimización del modelo
 
model.objective = "obj"
model.reactions.get_by_id("obj").bounds = (0, 1000)

# %% MM  --- 2021-01-15 10:27 --- Exporta un grafo 
def cobra2networkx(modelo):
    """Toma un modelo de cobra y genera un grafo bipartito de NetworkX


    """    
    return grafo

import networkx as nx
from networkx.algorithms.bipartite.matrix import biadjacency_matrix      # Extrae matriz de adyacencia
from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix # Crea el grafo
from sklearn.preprocessing import binarize
from scipy.sparse import csr_matrix

tmp = binarize(create_stoichiometric_matrix(model))
tmp = csr_matrix(tmp)
tmp = from_biadjacency_matrix(tmp) 

total_nodes = range(0,len(tmp.nodes(data=False))) # Rango de nodos
particiones = [tmp.nodes(data=True)[i]["bipartite"] for i in total_nodes] # Tipo nodo

metabolites_nodes = [n for n, d in tmp.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
metabolites_names = ["m1","m2","m3","m4","m5"] # Añade nombres
reactions_nodes   = [n for n, d in tmp.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
reaction_names    = ["r1","r2","r3","r4","r5","r6","r7","r8"] # Añade nombres

names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reaction_names))
tmp = nx.relabel_nodes(tmp, names_mapped)

# %% MM --- 2021-01-15 10:28 --- Funcción que aplica atributos desde dos vectores
def met_atributos(grafo, name, vector, skip = 0):
    """Toma una lista de atributos y se los asigna a un objeto para Gephi

    Parameters
    ----------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    nombre : str
        Nombre del atributo siendo asignado. Ie. nombre de la columna en Gephi
    vector : list
        Una lista con valores para atributos. Se asigna desde el nodo [0] + skip
    skip : int
        Cantidad de nodos a omitir. Por defecto, los nodos de los metabolitos
        van primero en la estructura de datos. Reacciones=len(model.metabolites)
    Returns
    -------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    """
    grafo.nodes(data=False) # Consigue nombre de los nodos
    nodos = metabolites_names + reaction_names # Lista con nombres de nodos

    tmp_prices =   solution.shadow_prices.tolist() + ['']*8 # Shadow Prices
    tmp_fluxes = ['']*5 + solution.fluxes.tolist()          # Flujos

    tmp_prices = { nodos[i] : (tmp_prices)[i] for i in range(0, len(tmp_prices)) } 
    tmp_fluxes = { nodos[i] : (tmp_fluxes)[i] for i in range(0, len(tmp_prices)) } 

    nx.set_node_attributes(tmp, tmp_prices, 'prices' ) # Añade una columna de precios
    nx.set_node_attributes(tmp, tmp_fluxes, 'fluxes' ) # Añade una columna de flujos
    # ----- END
    # tmp.nodes(data=True) # VIS: muestra como deberia verse el grafo
    nx.write_gexf(tmp, "grafo.gexf") # Salida para Gephi

    return grafo

# TODO: Convertir esto en una función de verdad f(modelo, filename)

# %% MM --- 2021-01-15 10:24 --- Función toma dos listas y añade atributos
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

# %% --- Algo más que no se que hace ---

solution = model.optimize()
solution = solution.to_frame()

all_reactions = [model.reactions[i].reaction  for i in range(0,len(model.reactions)) ]
solution["all_reactions"] = all_reactions

# Función que tome los flux y los pegue como atributo NetworkX


# %% --- CENTRALITY
from networkx import degree_centrality
from networkx import eigenvector_centrality
from networkx import closeness_centrality
from networkx import betweenness_centrality
from networkx import load_centrality
from networkx import harmonic_centrality
from networkx import current_flow_closeness_centrality
from networkx import information_centrality
from networkx import current_flow_betweenness_centrality
from networkx import communicability_betweenness_centrality
from networkx import second_order_centrality

import time
import warnings
warnings.filterwarnings('ignore')
# %% 
graph   = nx.read_gpickle(" .gpickle")

dc      = nx.degree_centrality(graph)
ec      = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
cc      = nx.closeness_centrality(graph, distance=None, wf_improved=True)
bc      = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
lc      = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
hc      = nx.harmonic_centrality(graph, nbunch=None, distance=None)
cfcc    = nx.current_flow_closeness_centrality(graph)
ic      = nx.information_centrality(graph)
cfbc    = nx.current_flow_betweenness_centrality(graph)
cbc     = nx.communicability_betweenness_centrality(graph)
soc     = nx.second_order_centrality(graph)

# %% 

model

# %% SAMPLING THE STEADY-STATE

from cobra.sampling import sample

from cobra.sampling.optgp import OptGPSampler

flux_samples = OptGPSampler(model,  processes=16, thinning=500, nproj=100, seed=23)

flux_samples = flux_samples.sample(n = 1000)
flux_samples

# %% Pairwise correlations between flux distributions

import seaborn as sns

sns.pairplot(flux_samples, corner=True)