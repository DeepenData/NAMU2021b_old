#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades

Parametros: MODELO, SUBSISTEMAS
Output: CENTRALIDADES.tsv"""


# %% --- Importa el modelo
import time

from cobra.io import load_json_model, save_json_model

INPUT_MODEL = '../data/toy_metabolism_AA.json'
t1 = time.time() 
model = load_json_model(INPUT_MODEL)

solution_fba = model.optimize() # Optimización del modelo para check
solution_fba.fluxes # Check si los flujos funcionan

t2 = time.time(); print('\n TIME Importando y optimizando modelo:', (t2-t1)*1000 ,'ms')

# %% --- Proyección de la matrix estequimetrica unipartita de reacciones 
"""Usamos reacciones para así poder buscar genes. 
1. Crea una matriz estequimetrica que tiene las reacciones del modelo
2. Binarización de matriz de incidencia para eliminar los pesos y dirección de 
   las reacciones. Los coeficientes son un artefacto, cuando solo nos interesa
   la conectividad
5. Proyección de la matriz de incidencia hacia el espacio de las reacciones, lo
   que genera en matriz de adyacencia de reacciones
6. Ceros a la diagonal de la matriz, para eliminar auto-conexiones
7. Binarización, de nuevo. Solo nos interesa saber si una reaccion comparte 
   pathway con otra, no cuantos metabolitos comparten
4. Convierte la matriz en una sparce matix, para optimizar el uso de memoria en
   la computación de redes más grandes"""

import numpy as np
import networkx as nx

from cobra.util.array import create_stoichiometric_matrix
from sklearn.preprocessing import binarize

t1 = time.time()

S_matrix = create_stoichiometric_matrix(model) # Crea matriz de incidencia
S_matrix = binarize(abs(S_matrix) ) # Binarización de la matriz de incidencia 
    # para eliminar pesos y dirección que no son necesarios en este resultado

projected_S_matrix = np.matmul(S_matrix.T, S_matrix) # Proyección al espacio de 
  # las reacciones de la matriz en matriz de incidencia
  # De momento, numpy no soporta multiplicación de matrices sparce. Usa normal. 

"""La diagonal contiene la cantidad de metablitos que estan en cada reacción.
Si una reacción tiene 1 podemos afirmar que es una reacción de Boundary.  
np.matmul(S_matrix, S_matrix.T) # Genera la matriz de metabolitos
La diagonal de esta indicaria en cuantas reacicones participa cada metabolito"""

np.fill_diagonal(projected_S_matrix, 0) # Elimina la conexion consigo mismo

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int) # Re-binariza, 
    # para eliminar el peso al compartir más de un metabolito entre reacciones

t2 = time.time(); print('\n TIME Matriz proyectada:', (t2-t1)*1000 ,'ms')

# %% --- Crea el grafo a partir de la matriz de adyacencia

from networkx.convert_matrix import from_numpy_matrix
t1 = time.time()
G = from_numpy_matrix( reaction_adjacency_matrix )

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("This graph is not conected. Picking largest component.")
    largest_component = max(nx.connected_components(G), key=len)
    print ('Removed',len(G.nodes) - len(largest_component), 'nodes.')
    G = G.subgraph(largest_component)

t2 = time.time(); print('\n TIME Grafo base generado:', (t2-t1)*1000 ,'ms')
# %% --- Reasigna nombres a nodos

from networkx import relabel_nodes

# Crea un diccionario con una lista para cada nodo, en orden
node_dict = lambda l : dict( zip( list(G.nodes), l ) )

cursed_dict = node_dict( [reaction.id for reaction in model.reactions] )
G = relabel_nodes(G, cursed_dict, copy=True)

# %% --- Calculos de centralidades varios
"""Calculos de centralidades baseline,  
Existen 30+ centralidades, las cuales no aplican en todas las condiciones"""

def eight_centralities(grafo):
    """Calcula centralidades base de una red

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

def assign_eight_centralities(grafo, centralities, prefix='',subfix=''):
    """Asigna una lista de centralidades creada por eight_centralities() como 
    atributos

    Parameters
    ----------
    grafo : NetworkX graph
    centralities : list
        Una lista de centralidades creada por eight_centralities()
        Puede ser directamente eight_centralities(grafo)
    prefix : str, default = ''
        Un prefijo para el atributo de centralidad asignada
    subfix : str, default = ''
        Un subfijo para el atributo de centralidad asignada

    Returns
    -------
    grafo : NetworkX graph
    """
    centralities_names = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality',
    'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality']
   
    centralities_names = [ (prefix + name + subfix) for name in centralities_names]
   
    for i in range(0, 7):
        nx.set_node_attributes(grafo, centralities[i], centralities_names[i])
    
    return grafo

t1 = time.time()

baseline_centralities = eight_centralities(G) # Calcula las centralidades
#G = assign_eight_centralities(G, centralities) # Asigna centralidades a nodos

t2 = time.time(); print('\n TIME Centralidades base calculadas:', (t2-t1)*1000 ,'ms')

# %% --- Calcula centralidades para cada nodo del grafo

def delta_centralities(grafo, removed_nodes, prefix='',subfix=''):
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
    prefix : str, default = ''
    subfix : str, default = ''

    Return
    ------
    grafo :
        Un grafo de NetworkX con atributos de centralidad por cada nodo removido
    centralities : list
        Una lista de diccionarios con calculos de centralidades
    breaks : list
        Lista de nodos que rompen la continuidad del espacio-tiempo (y la red)
    """
    import networkx
    from copy import deepcopy
    breaks = []; centralities = []
    for node in removed_nodes:
        print( 'Iterando en:', node )
        delta = deepcopy( grafo )
        #removed = prefix + '_REMOVED_' + str(node) + subfix
        delta.remove_node( str(node) ) # Elimina el nodo de la iteración
        centrality = eight_centralities(delta)
        if (centrality[-1] == {}):
            print( str(node), 'breaks continuity!')
            breaks.append( str(node) )
        centralities.append( centrality )
        # Asigna centralidades calculadas a los nodos
        # grafo = assign_eight_centralities(grafo, centralities, subfix=removed)

    return centralities, breaks

# Calcula centralidades delta para cada nodo del grafo
# SUBSYSTEM
INPUT_REMOVED_NODES = G.nodes
t1 = time.time()
delta_centralities, breaks = delta_centralities(G, INPUT_REMOVED_NODES) 
t2 = time.time(); print('\n TIME Centralidades delta calculadas:', (t2-t1)*1000 ,'ms')

# TODO: hacer algo con 'braks' ? Deberia ser una tabla de nodos y los desconectados

# %% --- Comprimir deltas en objetos manejables
"""Genera un objeto dataframe de [nodos A] columnas y [nodos removidos B] filas. 
Cada celda contiene el promedio de las ocho centralidades calculadas para ese 
nodo A en la remoción del nodo B correspondiente a la fila. Si el valor es NaN, 
es porque el nodo B removido en la fila es el mismo nodo A (cuya centralidad no 
existe porque fue removido)"""

import pandas as pd

t1=time.time()

perturbed_centralities = []

for delta in delta_centralities:
    tmp = pd.DataFrame.from_dict( delta ) # Selecciona un grupo de 8 centralidades
    tmp = dict( tmp.mean( axis=0 ) ) # Calcula el promedio de centralidades
    perturbed_centralities.append( tmp ) # Al diccionario de centralidades perturbadas

print("Nodes removed for iterative delta centralities calculation:", len( INPUT_REMOVED_NODES ))
print("Unconected graphs generated:", len(breaks) )

perturbed_centralities = pd.DataFrame.from_dict( perturbed_centralities ) # Tabla por nodo

cols = list(perturbed_centralities.columns); cols = [cols[-1]] + cols[:-1] # Reordena columnas
perturbed_centralities = perturbed_centralities[cols] # Así la primera es la primera reacción 
# TODO: exportar esta tabla? 

perturbed = perturbed_centralities.mean( axis=0 ) # Promedio de perturbadas

# Termina de comprimir las centralidades originales porque antes no cargo pandas

baseline = pd.DataFrame.from_dict( baseline_centralities )
baseline = baseline.transpose() # Por algun motivo queda con la otra orientación ?
baseline = baseline.mean( axis=1 )

# %% --- Calcula contribución a centralidades de S
"""Como nos interesa saber cual es la contribución del resto de los nodos al sub
-sistema que estamos estudiando, sacamos el fold-change en base al promedio de 
los nodos cuando son removidos. Ie. control / experimental; en lugar de el más 
comun experimental / control. Basicamente porque aqui la existencia del subsistema
es nuestra condición experimental"""

log2_contribution = np.log2( baseline / perturbed )

# %% --- Empacando todo en un Dataframe que servira como resultados

tabla_resultado = pd.DataFrame({
    "ids" :         [rx.id           for rx in model.reactions],
    "formula" :     [rx.reaction     for rx in model.reactions],
    "flux" :        [rx.flux         for rx in model.reactions],
    "sensitivity" : [rx.reduced_cost for rx in model.reactions],
    "baseline_centrality":           baseline,
    "preturbed_centrality":          perturbed,
    "log2_centrality_contribution" : log2_contribution
})

# Exporta esta tabla de datos intermedio
# index=False porque id contiene lo mismo que el index. 
OUTPUT_TABLE = '../results/reactions_delta_centrality.csv'
tabla_resultado.to_csv( OUTPUT_TABLE , index=False )
t2 = time.time(); print('\n TIME Tabla resumen generada:', (t2-t1)*1000 ,'ms') 
# %% --- Empacando en un grafo de gephi
t1 = time.time()
nx.set_node_attributes(G, node_dict( [rx.reaction     for rx in model.reactions] ), "reaction")
nx.set_node_attributes(G, node_dict( [rx.flux         for rx in model.reactions] ), "flux")
nx.set_node_attributes(G, node_dict( [rx.reduced_cost for rx in model.reactions] ), "sensitivity")
nx.set_node_attributes(G, node_dict( baseline ),          "baseline_centrality")
nx.set_node_attributes(G, node_dict( perturbed ),         "preturbed_centrality")
nx.set_node_attributes(G, node_dict( log2_contribution ), "log2_centrality_contribution")

from networkx import write_gexf

OUTPUT_GRAPH =  "../results/reactions_delta_centrality.gexf"
write_gexf( G , OUTPUT_GRAPH )
t2 = time.time(); print('\n TIME Grafo generado:', (t2-t1)*1000 ,'ms')

#"""gexfs: 
#1.- bipartito // listo en otro archivo
#2.- unipartito // Listo aqui! 
# 1.- unipartito flujos.
# 2.- unipartito sensibilidades.
# 3.- unipartito centralidad total(overall)"""
# %%
