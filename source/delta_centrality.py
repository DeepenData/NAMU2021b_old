#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades con un resumen de las reacciones, y el log2 
entre la reacción sin modificar y modificada. Genera un grafo gephx con esta información codificada en los nodos. 
Parametros: INPUT_MODEL, REMOVED_NODES, [N_WORKERS]
Output: OUTPUT_TENSOR_BASELINES, OUTPUT_TENSOR_PERTURBATED"""

import warnings; warnings.filterwarnings('ignore') # Ignora warnings

import ray
import time
import os
# Conectando a cluster de Ray inicializado con sbatch
# ray.init(address='auto', _node_ip_address=os.environ["ip_head"].split(":")[0], _redis_password=os.environ["redis_password"])
ray.init() # Desactivar esto para correr en el cluster

# --- Sección que importa modulos

import numpy as np
import pandas as pd
import sys

import networkx as nx
import graph_tool.all as gt 

# --- Definición de funciones

def eight_centralities(grafo_nx, grafo_gt, index_nodes):
    """Calcula centralidades base de una red, usando ocho metodos de NetworkX. 
    Los ultimos dos metodos son omitidos si el grafo no es totalmente conectado,
    y devuelve diccionarios vacios en cambio. 
    - grafo_nx : toma un objeto grafo en el formato de nx
    - grafo_gt : toma un objeto grafo en el formato de graph-tools
    - 

    ```   
    ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
    'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality',
    'current_flow_closeness_centrality', 'current_flow_betweenness_centrality',
    'approximate_current_flow_betweenness_centrality', 'communicability_betweenness_centrality']
    ```

    Returns
    -------
    centralities : list
        Una lista de diccionarios {nodo : centralidad} con las centralidades
        segun la lista descrita arriba."""

    import numpy    as np
    import pandas   as pd
    import networkx as nx
    import graph_tool.all as gt 

    # CREA UN INDICE TEMPORAL CON EL NODO REMOVIDO
    index_short = list( Gnx.nodes )

    # CALCULA CENTRALIDADES CON GRAPH-TOOLS
    # realmente son pocos, pero algo es algo
    rm , ec = gt.eigenvector(Ggt); ec = pd.Series( ec.a , index= index_short )
    bc , rm = gt.betweenness(Ggt); bc = pd.Series( bc.a , index= index_short )
    cc      = gt.closeness(Ggt);   cc = pd.Series( cc.a , index= index_short )
    
    # CALCULA CENTRALIDADES CON NETWORX
    # lo cual es leeeeeentisimo, pero funciona
    hc    = nx.harmonic_centrality(grafo_nx, nbunch=None, distance=None)            ; hc = pd.Series( hc )
    dc    = nx.degree_centrality(grafo_nx)                                          ; dc = pd.Series( dc )
    lc    = nx.load_centrality(grafo_nx, cutoff=None, normalized=True, weight=None) ; lc = pd.Series( lc )
    
    """ Requieren un grafo full conected """
    largest_component = max(nx.connected_components(grafo_nx), key=len)
    grafo_nx           = grafo_nx.subgraph(largest_component)
    
    ic    = nx.information_centrality(grafo_nx)                                     ; ic = pd.Series( ic )
    soc   = nx.second_order_centrality(grafo_nx)                                    ; soc = pd.Series( soc )
    """ Son nefastos de computar """
    cfcc  = nx.current_flow_closeness_centrality(grafo_nx)                          ; cfcc = pd.Series( cfcc )
    cfbc  = nx.current_flow_betweenness_centrality(grafo_nx)                        ; cfbc = pd.Series( cfbc )
    acfbc = nx.approximate_current_flow_betweenness_centrality(grafo_nx)            ; acfbc = pd.Series( acfbc )
    cbc   = nx.communicability_betweenness_centrality(grafo_nx)                     ; cbc = pd.Series( cbc )
   
    # REINDEXANDO PARA INCLUIR EL NODO OMITIDO
    centralities = [hc, ec, dc, bc, cc, lc, ic, soc, cfcc, cfbc, acfbc, cbc]
    centralities = [ cent.reindex( index_nodes ) for cent in centralities]

    # CONVIERTE A UN ARRAY DE NUMPY
    centralities = [ cent.to_numpy(dtype='float32') for cent in centralities]
    centralities = np.array( centralities ).T # Array 2D Transpuesto. Dins = (999+NaN, 12)

    # OUTPUT FINAL
    return centralities

@ray.remote # Usa el decorador de Ray para definir esta como la paralelizable
def delta_centrality(Gnx, Ggt, removed_nodes, index_nodes):
    """WARPER PARA EJECUCIÓN (no excesivamente) PARALELA EN RAY
    Toma un grafo y una lista de nodos y calcula la centralidad del grafo al
    remover cada uno de esos nodos.

    Prefix y subfix son información sobre la lista de los nodos removidos. El
    atributo final en el grafo es ≪`PREFIX` REMOVED_nodo[i] `SUBFIX` ≫

    Parameters
    ----------
    grafo_nx : 
        Un grafo de NetworkX 
    grafo_gt :
        Un grafo de Graph-tools
    removed_nodes : list
        Una lista de nodos a remover. Puede ser grafo.nodes, pero la idea es que
        es una sección más corta para mejorar el rendimiento. 
    index_nodes : list
        Una lista con todos los nodos. Es usada para ajustar los NaN cuando se
        remueve un nodo. 

    Return
    ------
    centralities : list
        Una lista de arrays con calculos de centralidades
    breaks : list
        Lista de nodos que rompen la continuidad del espacio-tiempo (y la red)
    """

    import networkx as nx
    import graph_tool.all as gt 
    from copy import deepcopy

    centralities = []
    for node in removed_nodes:
        index = index_nodes.index( node ) # Cual es el nodo que fue removido
        print( time.asctime( time.localtime(time.time()) ), 'Iterando en:', node )
        
        delta_nx = deepcopy( Gnx ); delta_gt = deepcopy( Ggt ) # Copias de trabajo
        
        delta_nx.remove_node( str(node) ) # Elimina el nodo de la iteración
        # TODO: hacer la remoción del nodo para Graph-tools
        
        centrality = eight_centralities( delta_nx, delta_gt, index_nodes )
        centrality = np.array(centrality).T # Deberia ser de forma (999+NaN, 12)
        centralities.append( centrality )

    return centralities

# --- Params

INPUT_MODEL = sys.argv[1] # Modelo base
N_WORKERS = abs(int(sys.argv[2])) if sys.argv[2] else 1 # Nucleos para procesamiento paralelo

# %% --- Modelo de Cobra
# REMOVIDO

# %% --- Modelo a grafo de NetworkX
# REMOVIDO

# %% --- Importa los grafos
# Utiliza ModelTOGrapgh.py, el cual genera un archivo de grafo agnostico
Gnx = nx.read_graphml("./tmp/graph.graphml")
Ggt =   gt.load_graph("./tmp/graph.graphml")

index_nodes = list( Gnx.nodes ) # Crea el indice de nodos como lista

# %% --- Pasa objetos indice al Object Store de Ray para preparar el ambiente

index_nodes_ray = ray.put(index_nodes) # Pasa el indice de los nodos

Gnx_ray = ray.put(Gnx) # Pasa el objeto grafo de NetworkX
Ggt_ray = ray.put(Ggt) # Pasa el objeto grafo de Graph-tools

# %% --- Centralidades delta (Sección paralelizada)

REMOVED_NODES = index_nodes # Todos los nodos. Al final los subsistemas se vuelven redundantes, así que esto es innecesario
# split_removed_nodes = np.array_split(REMOVED_NODES, min(N_WORKERS, REMOVED_NODES//2) ) # Usa min(N_WORKERS, INPUT//2) para evitar sobre-paralelizar
split_removed_nodes = np.array_split(REMOVED_NODES, min(N_WORKERS, len(REMOVED_NODES)//2) ) # Usa min(N_WORKERS, INPUT//2) para evitar sobre-paralelizar

deltas = [delta_centrality.remote(Gnx_ray, Ggt_ray, NODES, REMOVED_NODES) for NODES in split_removed_nodes] # SECCIÓN QUE MANDA CODIGO PARALELO EN RAY

# %% --- Interludio de centralidades base mientras Ray calcula las perturbadas
baseline_centralities = eight_centralities( Gnx, Ggt ) # Calcula centralidades base
baseline_centralities = dict_list_to_numpy( baseline_centralities , G.nodes ) 

deltas = ray.get( deltas ) # SECCIÓN QUE TOMA LAS RESPUESTAS DE VUELTA DEL CODIGO PARALELO

# %% --- Tensor de centralidades perturbadas 
# Como el resultado generado es un par de output parciales, necesito una función
# que me permita pegarlos y ensamblar un output final. Por eso creo dos list
# comprehensions para separar la lista original. Es más rapido que un for loop.
perturbed_centralities = [ delta[0] for delta in deltas] 
breaks =                 [ delta[1] for delta in deltas] 

# %% --- Creando las tablas de salida

perturbed_centralities = [ i for ii in perturbed_centralities for i in ii ] # Squash N_WORKERS x NO/DO/S -> NODOS
perturbed_centralities = np.asarray( perturbed_centralities ) # Tensor de dimensiones n x m x c (1000 x 999+NaN x 8) listo 

baseline_centralities = np.asarray( [ baseline_centralities ] ) # Convierte en tensor de dimensiones 1 x m x c  (1 x 1000 x 8)

import pickle

outfile1 = open('./tmp/baseline_centralities','wb'); pickle.dump(baseline_centralities,outfile1); outfile1.close()
outfile2 = open('./tmp/perturbed_centralities','wb'); pickle.dump(perturbed_centralities,outfile2); outfile2.close()
outfile3 = open('./tmp/index_nodes','wb'); pickle.dump(index_nodes,outfile3); outfile3.close()
#outfile4 = open('./tmp/breaks','wb'); pickle.dump(breaks,outfile4); outfile4.close()
