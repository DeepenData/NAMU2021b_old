#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades con un resumen de las reacciones, y el log2 
entre la reacción sin modificar y modificada. Genera un grafo gephx con esta información codificada en los nodos. 
Parametros: INPUT_MODEL, INPUT_REMOVED_NODES, [N_WORKERS]
Output: OUTPUT_TENSOR_BASELINES, OUTPUT_TENSOR_PERTURBATED"""

import warnings; warnings.filterwarnings('ignore') # Ignora warnings

import ray
import time
import os
# Conectando a cluster de Ray inicializado con sbatch
ray.init(address='auto', _node_ip_address=os.environ["ip_head"].split(":")[0], _redis_password=os.environ["redis_password"])

# --- Sección que importa modulos

import numpy as np
import networkx as nx
import pandas as pd
import sys

# --- Definición de funciones

def eight_centralities(grafo):
    """Calcula centralidades base de una red, usando ocho metodos de NetworkX. 
    Los ultimos dos metodos son omitidos si el grafo no es totalmente conectado,
    y devuelve diccionarios vacios en cambio. 

    ```   
    ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 
    'betweenness_centrality', 'closeness_centrality', 'load_centrality', 
    'information_centrality', 'second_order_centrality']
    ```

    Returns
    -------
    centralities : list
        Una lista de diccionarios {nodo : centralidad} con las centralidades
        segun la lista descrita arriba."""

    import networkx as nx

    largest_component = max(nx.connected_components(grafo), key=len)
    grafo             = grafo.subgraph(largest_component)
    hc  = nx.harmonic_centrality(grafo, nbunch=None, distance=None)
    ec  = nx.eigenvector_centrality(grafo, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc  = nx.degree_centrality(grafo)
    bc  = nx.betweenness_centrality(grafo, normalized=True, weight=None, endpoints=False, seed=None)
    """↓ Here be dragons ↓"""
    cc  = nx.closeness_centrality(grafo, distance=None, wf_improved=True)
    lc  = nx.load_centrality(grafo, cutoff=None, normalized=True, weight=None)
    """ Requieren un grafo full conected """
    #if not (nx.is_connected(grafo)) : 
    #    ic = {}; soc = {}
    #else:
    #ic  = nx.information_centrality(grafo)
    #soc = nx.second_order_centrality(grafo)
    ic  = nx.information_centrality(grafo)
    soc = nx.second_order_centrality(grafo)

    cfcc    = nx.current_flow_closeness_centrality(grafo)    
    cfbc    = nx.current_flow_betweenness_centrality(grafo)
    acfbc   = nx.approximate_current_flow_betweenness_centrality(grafo)
    cbc     = nx.communicability_betweenness_centrality(grafo)   



    centralities = [hc, ec, dc, bc, cc, lc, ic, soc, cfcc, cfbc, acfbc, cbc]

    return centralities

def dict_list_to_numpy(node, node_list):
    """Convierte una lista de diccionarios de centralidades a un array 2D de Numpy

    Parameters
    ----------
    node : list
        Lista de diccionarios de centralidades
    node_list : list
        Lista de los nodos cuyas centralidades se calcularon

    Returns
    -------
    numpy.array
        Un array 2D de forma (m: nodos , c: centralidades)
    """

    import numpy as np
    import pandas as pd

    node = pd.DataFrame.from_dict( node ) # Dataframe centralidades (c) x nodos (m)
    node = node.T # Transpone para nodos (m) x centralidades (c)

    node = node.reindex( node_list ) # Indexa por nodos
    # En esta parte ya tenemos un dataframe de dimensiones (m=n, c) consistente, incluyendo NaNs
    
    node = node.to_numpy(dtype='float32') # Convierte a array de Numpy
    return node

@ray.remote # Usa el decorador de Ray para definir esta como la paralelizable
def delta_centrality(G, removed_nodes):
    """Toma un grafo y una lista de nodos y calcula la centralidad del grafo al
    remover cada uno de esos nodos.

    Prefix y subfix son información sobre la lista de los nodos removidos. El
    atributo final en el grafo es ≪`PREFIX` REMOVED_nodo[i] `SUBFIX` ≫

    Parameters
    ----------
    grafo : 
        Un grafo de NetworkX 
    removed_nodes :
        Una lista de nodos a remover. Puede ser grafo.nodes

    Return
    ------
    grafo :
        Un grafo de NetworkX con atributos de centralidad por cada nodo removido
    centralities : list
        Una lista de diccionarios con calculos de centralidades
    breaks : list
        Lista de nodos que rompen la continuidad del espacio-tiempo (y la red)
    """

    removed_nodes = list(removed_nodes)

    import networkx as nx
    from copy import deepcopy

    breaks = []; centralities = []
    for node in removed_nodes:
        print( time.asctime( time.localtime(time.time()) ), 'Iterando en:', node )
        delta = deepcopy( G )
        delta.remove_node( str(node) ) # Elimina el nodo de la iteración
        centrality = eight_centralities(delta)
        if (centrality[-1] == {}):
            print( str(node), 'breaks continuity!')
            breaks.append( str(node) )
        centrality = dict_list_to_numpy( centrality , G.nodes ) # Convierte a Numpy Array
        centralities.append( centrality )

    return centralities, breaks

# --- Params

INPUT_MODEL = sys.argv[1] # Modelo base

## --- Modelo de Cobra
t1 = time.time()

import warnings; warnings.filterwarnings('ignore') # Ignora warnings
from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL)

print("Iniciando optimización del modelo", time.asctime( time.localtime(time.time()) ) )
#solution_fba = model.optimize() # Optimización del modelo para check
#solution_fba.fluxes # Check si los flujos funcionan

t2 = time.time(); print('\n','TIME Importando y optimizando modelo:', (t2-t1)*1000 ,'ms')

## --- Modelo a grafo de NetworkX
t1 = time.time() # Contador de esta sección

from cobra.util.array import create_stoichiometric_matrix

S_matrix = create_stoichiometric_matrix(model)
S_matrix = (abs(S_matrix) )
S_matrix = (S_matrix > 0.0).astype(np.int_)

projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
np.fill_diagonal(projected_S_matrix, 0)

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)

G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("This graph is not conected. Picking largest component.")
    
    largest_component = max(nx.connected_components(G), key=len)
    print ('Removed',len(G.nodes) - len(largest_component), 'nodes.')
    print ( len(G.nodes), "nodes remaining.")

    G = G.subgraph(largest_component)

node_dict = lambda l : dict( zip( list(G.nodes), l ) )

cursed_dict = node_dict( [reaction.id for reaction in model.reactions] )
G = nx.relabel_nodes(G, cursed_dict, copy=True)

t2 = time.time(); print('\n','TIME Conversión a grafo:', (t2-t1)*1000 ,'ms .', time.asctime( time.localtime(time.time()) ) )

## --- Centralidades delta (Sección paralelizada)
t1 = time.time()

INPUT_REMOVED_NODES = G.nodes # Todos los nodos
N_WORKERS = abs(int(sys.argv[2])) if sys.argv[2] else 1 # Nucleos para procesamiento paralelo
# split_removed_nodes = np.array_split(INPUT_REMOVED_NODES, min(N_WORKERS, INPUT_REMOVED_NODES//2) ) # Usa min(N_WORKERS, INPUT//2) para evitar sobre-paralelizar
split_removed_nodes = np.array_split(INPUT_REMOVED_NODES, min(N_WORKERS, len(INPUT_REMOVED_NODES)//2) ) # Usa min(N_WORKERS, INPUT//2) para evitar sobre-paralelizar

G_ray = ray.put(G) # Pasa esto al object store antes que pasarlo multiples veces
deltas = [delta_centrality.remote(G_ray, NODES) for NODES in split_removed_nodes] # SECCIÓN QUE MANDA CODIGO PARALELO EN RAY

## --- Interludio de centralidades base mientras Ray calcula las perturbadas
baseline_centralities = eight_centralities( G ) # Calcula centralidades base
baseline_centralities = dict_list_to_numpy( baseline_centralities , G.nodes ) 

deltas = ray.get( deltas ) # SECCIÓN QUE TOMA LAS RESPUESTAS DE VUELTA DEL CODIGO PARALELO

t2 = time.time(); print('\n','TIME Centralidades calculadas en paralelo:', time.asctime( time.localtime(time.time()) ) )

## --- Tensor de centralidades perturbadas 

t1 = time.time()
# Como el resultado generado es un par de output parciales, necesito una función
# que me permita pegarlos y ensamblar un output final. Por eso creo dos list
# comprehensions para separar la lista original. Es más rapido que un for loop.
perturbed_centralities = [ delta[0] for delta in deltas] 
breaks =                 [ delta[1] for delta in deltas] 

t2 = time.time(); print('\n','TIME Merge de deltas:', (t2-t1)*1000 ,'ms')

## --- Creando las tablas de salida

t1 = time.time() # TODO: ver que de esto puedo pasar a la ejecución paralelizada

perturbed_centralities = [ i for ii in perturbed_centralities for i in ii ] # Squash N_WORKERS x NO/DO/S -> NODOS
perturbed_centralities = np.asarray( perturbed_centralities ) # Tensor de dimensiones n x m x c (1000 x 999+NaN x 8) listo 

index_nodes = G.nodes # Crea el indice de nodos como lista

baseline_centralities = np.asarray( [ baseline_centralities ] ) # Convierte en tensor de dimensiones 1 x m x c  (1 x 1000 x 8)

import pickle

outfile1 = open('./tmp/baseline_centralities','wb'); pickle.dump(baseline_centralities,outfile1); outfile1.close()
outfile2 = open('./tmp/perturbed_centralities','wb'); pickle.dump(perturbed_centralities,outfile2); outfile2.close()
outfile3 = open('./tmp/index_nodes','wb'); pickle.dump(index_nodes,outfile3); outfile3.close()
outfile4 = open('./tmp/breaks','wb'); pickle.dump(breaks,outfile4); outfile4.close()

"""
for delta in perturbed_centralities:
    tmp = pd.DataFrame.from_dict( delta ) # Selecciona un grupo de 8 centralidades
    tmp = dict( tmp.mean( axis=0 ) ) # Calcula el promedio de centralidades
    perturbed_centralities.append( tmp ) # Al diccionario de centralidades perturbadas

print("Nodes removed for iterative delta centralities calculation:", len( INPUT_REMOVED_NODES ))
print("Unconected graphs generated:", len(breaks) )

perturbed_centralities = pd.DataFrame.from_dict( perturbed_centralities ) # Tabla por nodo

#cols = list(perturbed_centralities.columns); cols = [cols[-1]] + cols[:-1] # Reordena columnas
#perturbed_centralities = perturbed_centralities[cols] # Así la primera es la primera reacción

perturbed = perturbed_centralities.mean( axis=0 ) # Promedio de perturbadas

# Termina de comprimir las centralidades originales porque antes no cargo pandas

baseline = pd.DataFrame.from_dict( baseline_centralities )
baseline = baseline.transpose() # Por algun motivo queda con la otra orientación ?
baseline = baseline.mean( axis=1 )

log2_contribution = np.log2( baseline / perturbed )

tabla_resultado = pd.DataFrame({
    "ids" :         [rx.id           for rx in model.reactions],
    "formula" :     [rx.reaction     for rx in model.reactions],
    "flux" :        [rx.flux         for rx in model.reactions],
    "sensitivity" : [rx.reduced_cost for rx in model.reactions],
    "baseline_centrality":           baseline,
    "preturbed_centrality":          perturbed,
    "log2_centrality_contribution" : log2_contribution
})

OUTPUT_RESUME_TABLE = 'tmp.reactions_delta_centrality.csv'
tabla_resultado.to_csv( OUTPUT_RESUME_TABLE , index=False )
t2 = time.time(); print('\n','TIME Tabla resumen generada:', (t2-t1)*1000 ,'ms') 

# --- Empacando el grafo de gephi
t1 = time.time()
nx.set_node_attributes(G, node_dict( [rx.reaction     for rx in model.reactions] ), "reaction")
nx.set_node_attributes(G, node_dict( [rx.flux         for rx in model.reactions] ), "flux")
nx.set_node_attributes(G, node_dict( [rx.reduced_cost for rx in model.reactions] ), "sensitivity")
nx.set_node_attributes(G, node_dict( baseline ),          "baseline_centrality")
nx.set_node_attributes(G, node_dict( perturbed ),         "preturbed_centrality")
nx.set_node_attributes(G, node_dict( log2_contribution ), "log2_centrality_contribution")

OUTPUT_GRAPH =  "tmp.reactions_delta_centrality.gexf"
nx.write_gexf( G , OUTPUT_GRAPH )
t2 = time.time(); print('\n','TIME Grafo generado:', (t2-t1)*1000 ,'ms')
"""
