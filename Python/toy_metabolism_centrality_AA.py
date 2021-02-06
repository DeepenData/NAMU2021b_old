#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades

Parametros: MODELO, SUBSISTEMAS
Output: CENTRALIDADES.tsv"""


# %% --- Importa el modelo
from cobra.io import load_json_model, save_json_model

INPUT = 'toy_metabolism_AA.json'
model = load_json_model(INPUT)

solution_fba = model.optimize() # Optimización del modelo para check
solution_fba.fluxes # Check si los flujos funcionan

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

S_matrix = create_stoichiometric_matrix(model) # Crea matriz de incidencia
S_matrix = binarize(abs(S_matrix) ) # Binarización de la matriz de incidencia 
    # para eliminar pesos y dirección que no son necesarios en este resultado

projected_S_matrix = np.matmul(S_matrix.T, S_matrix) # Proyección al espacio de 
  # las reacciones de la matriz en matriz de incidencia

"""La diagonal contiene la cantidad de metablitos que estan en cada reacción.
Si una reacción tiene 1 podemos afirmar que es una reacción de Boundary.  
np.matmul(S_matrix, S_matrix.T) # Genera la matriz de metabolitos
La diagonal de esta indicaria en cuantas reacicones participa cada metabolito"""

np.fill_diagonal(projected_S_matrix, 0) # Elimina la conexion consigo mismo

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int) # Re-binariza, 
    # para eliminar el peso al compartir más de un metabolito entre reacciones

# %% --- Crea el grafo a partir de la matriz de adyacencia
from scipy.sparse import csr_matrix

from networkx.convert_matrix import from_scipy_sparse_matrix
reaction_adjacency_matrix = csr_matrix(reaction_adjacency_matrix)
G = from_scipy_sparse_matrix(reaction_adjacency_matrix)

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("This graph is not conected. Picking largest component.")
    largest_component = max(nx.connected_components(G), key=len)
    print ('Removed',len(G.nodes) - len(largest_component), 'nodes.')
    G = G.subgraph(largest_component)

# TODO: funcionalizar esto?
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
        continue
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

centralities = eight_centralities(G) # Calcula las centralidades
G = assign_eight_centralities(G, centralities) # Asigna centralidades a nodos

# %% --- Calcula centralidades para cada nodo del grafo

def delta_centralities(grafo, list_remove, prefix='',subfix=''):
    """Toma un grafo y una lista de nodos y calcula la centralidad del grafo al
    remover cada uno de esos nodos.

    Prefix y subfix son información sobre la lista de los nodos removidos. El
    atributo final en el grafo es ≪`PREFIX` REMOVED_nodo[i] `SUBFIX` ≫

    Parameters
    ----------
    grafo : 
        Un grafo de NetworkX 
    list_remove :
        Una lista de nodos a remover. Puede ser grafo.nodes
    prefix : str, default = ''
    subfix : str, default = ''

    Return
    ------
    grafo :
        Un grafo de NetworkX con atributos de centralidad por cada nodo removido
    breaks : list
        Lista de nodos que rompen la continuidad del espacio-tiempo (y la red)
    """
    import networkx
    from copy import deepcopy
    breaks = []; cents = []
    for node in list_remove:
        print( 'Iterando en:', node )
        delta = deepcopy( grafo )
        removed = prefix + '_REMOVED_' + str(node) + subfix
        delta.remove_node( str(node) ) # Elimina el nodo de la iteración
        centrality = eight_centralities(delta)
        if (centrality[-1] == {}):
            print( str(node), 'breaks continuity!')
            breaks.append( removed )
        cents.append( centrality )
        # Asigna centralidades calculadas a los nodos
        # grafo = assign_eight_centralities(grafo, centralities, subfix=removed)

    return grafo, cents, breaks

# Calcula centralidades delta para cada nodo del grafo
G, cents, breaks = delta_centralities(G, G.nodes) 

# %% --- Celtralidades delta por subsistema
"""
subsystem_s =  [r1, r2, r3]

nodes_without_subsystem_s = set(nodes) - set(subsystem_s)

baseline centralities: (C_r1), (C_r2), (C_r3)
perturbed centralities: C_r1_node_i_removed, C_r2_node_i_removed, C_r3_node_i_removed
unperturbed: mean(C_r1, C_r2, C_r3)
perturbed:   mean(C_r1_node_i_removed, C_r2_node_i_removed, C_r3_node_i_removed )
node_i_contribution_to_S = np.log2(unperturbed/perturbed)

Parallel

with concurrent.futures.ProcessPoolExecutor() as executor:
     for r in executor.map(calc_centr_rmv, my_range):
         


"""