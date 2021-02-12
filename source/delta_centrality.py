#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades con un resumen de las reacciones, y el log2 
entre la reacción sin modificar y modificada. Genera un grafo gephx con esta información codificada en los nodos. 
Parametros: INPUT_MODEL, INPUT_REMOVED_NODES, [N_WORKERS]
Output: OUTPUT_NODE_DELTA, OUTPUT_RESUME_TABLE, OUTPUT_GRAPH"""

import ray
import time
ray.init(address='auto', _redis_password='5241590000000000')

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

    hc  = nx.harmonic_centrality(grafo, nbunch=None, distance=None)
    ec  = nx.eigenvector_centrality(grafo, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc  = nx.degree_centrality(grafo)
    bc  = nx.betweenness_centrality(grafo, normalized=True, weight=None, endpoints=False, seed=None)

    """↓ Here be dragons ↓"""
    cc  = nx.closeness_centrality(grafo, distance=None, wf_improved=True)
    lc  = nx.load_centrality(grafo, cutoff=None, normalized=True, weight=None)
    
    """ Requieren un grafo full conected """
    if not (nx.is_connected(grafo)) : 
        ic = {}; soc = {}
    else:
        ic  = nx.information_centrality(grafo)
        soc = nx.second_order_centrality(grafo)

    centralities = [hc, ec, dc, bc, cc, lc, ic, soc]

    return centralities

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
        print( 'Iterando en:', node )
        delta = deepcopy( G )
        delta.remove_node( str(node) ) # Elimina el nodo de la iteración
        centrality = eight_centralities(delta)
        if (centrality[-1] == {}):
            print( str(node), 'breaks continuity!')
            breaks.append( str(node) )
        centralities.append( centrality )

    return centralities, breaks

# --- Params

INPUT_MODEL = '../data/toy_metabolism_AA.json'

# --- Sección que importa modulos

import numpy as np
import networkx as nx
import pandas as pd

## --- Modelo de Cobra
t1 = time.time()

from cobra.io import load_json_model, save_json_model
model = load_json_model(INPUT_MODEL)

solution_fba = model.optimize() # Optimización del modelo para check
solution_fba.fluxes # Check si los flujos funcionan

t2 = time.time(); print('\n','TIME Importando y optimizando modelo:', (t2-t1)*1000 ,'ms')

## Modelo a grafo de NetworkX
t1 = time.time() # Contador de esta sección

from cobra.util.array import create_stoichiometric_matrix
from sklearn.preprocessing import binarize # TODO: eliminar dependencia

S_matrix = create_stoichiometric_matrix(model)
S_matrix = binarize(abs(S_matrix) ) # TODO: Usar Numpy

projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
np.fill_diagonal(projected_S_matrix, 0)

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)

G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("This graph is not conected. Picking largest component.")
    
    largest_component = max(nx.connected_components(G), key=len)
    print ('Removed',len(G.nodes) - len(largest_component), 'nodes.')
    
    G = G.subgraph(largest_component)

node_dict = lambda l : dict( zip( list(G.nodes), l ) )

cursed_dict = node_dict( [reaction.id for reaction in model.reactions] )
G = nx.relabel_nodes(G, cursed_dict, copy=True)

t2 = time.time(); print('\n','TIME Conversión a grafo:', (t2-t1)*1000 ,'ms')

## --- Centralidades base
t1 = time.time()
baseline_centralities = eight_centralities(G)
t2 = time.time(); print('\n','TIME Calculo de centralidades base:', (t2-t1)*1000 ,'ms')

## --- Centralidades delta (Sección paralelizada)
t1 = time.time()

INPUT_REMOVED_NODES = G.nodes # Todos los nodos
N_WORKERS = 4 # Nucleos para procesamiento paralelo

split_removed_nodes = np.array_split(INPUT_REMOVED_NODES, N_WORKERS)

deltas = ray.get( 
    [delta_centrality.remote(G, NODES) for NODES in split_removed_nodes] 
    )

t2 = time.time(); print('\n','TIME Centralidades calculadas en paralelo:', (t2-t1)*1000 ,'ms')

## --- Centralidades delta (merge)

t1 = time.time()
# Como el resultado generado es un par de output parciales, necesito una función
# que me permita pegarlos y ensamblar un output final. Por eso creo dos listas 
# que seran el final, y luego .append() cada respupuesta parcial del par"""
delta_centralities = []; breaks = []
for delta in deltas:
    delta_centralities.append( delta[0] )
    breaks.append( delta[1] )

t2 = time.time(); print('\n','TIME Merge de deltas:', (t2-t1)*1000 ,'ms')

## --- Creando las tablas de salida

t1 = time.time()

perturbed_centralities = []

for delta in delta_centralities:
    tmp = pd.DataFrame.from_dict( delta ) # Selecciona un grupo de 8 centralidades
    tmp = dict( tmp.mean( axis=0 ) ) # Calcula el promedio de centralidades
    perturbed_centralities.append( tmp ) # Al diccionario de centralidades perturbadas

print("Nodes removed for iterative delta centralities calculation:", len( INPUT_REMOVED_NODES ))
print("Unconected graphs generated:", len(breaks) )

perturbed_centralities = pd.DataFrame.from_dict( perturbed_centralities ) # Tabla por nodo

# TODO: resolver este re-ordenado del dataframe y exportarlo. 
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
