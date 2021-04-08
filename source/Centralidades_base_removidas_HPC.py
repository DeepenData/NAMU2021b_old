#!/bin/python
"""
"L'impatience va éliminer l'irritabilité, la tension nerveuse et sa conséquence, le surmenage."
    - Ejemplo de uso de "surmenage", termino en frrancés para "fatiga mental"

Este es codigo para ejecutar en el cluster virtual Ray dentro del HPC. No intentes ejecutar esto en casa. 
"""
# %% --- IMPORTA EL GRAFO GENERADO

import networkx as nx
import numpy    as np
import pandas   as pd  

LITE=False # Variable para definir si los calculos deben ser complejos o no

G = nx.read_gpickle('./data/Recon2_rxn_proyected.gpickle')

# %% --- DEFINICIONES DE FUNCIONES

def compute_centralities_short(graph):
    """Computa las centralidades rapidas y devuelve una DataFrame (no reindexado) con estas"""
    # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
    hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    # dc    = nx.degree_centrality(graph)
    cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    ic    = nx.information_centrality(graph) # Requiere scipy

    # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
    centralities = {
        'harmonic_centrality' : hc ,
        'eigenvector_centrality' : ec ,
        'closeness_centrality' : cc ,
        'information_centrality' : ic
    }

    # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
    centralities = pd.DataFrame( centralities )

    return centralities

def compute_centralities(graph, lite=False, alpha=0.005):
    """Computa las doce centralidades y devuelve una DataFrame (no reindexado) con estas"""
    if lite == False:
        # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
        dc    = nx.degree_centrality(graph) # SANITY CHECK: Deberia ser similar a todo lo demas
        hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
        ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
        bc    = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
        cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
        lc    = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
        ic    = nx.information_centrality(graph) # Requiere scipy
        cbc   = nx.communicability_betweenness_centrality(graph) 
        kz    = nx.katz_centrality(graph, alpha = alpha)
        pr    = nx.pagerank(graph, alpha = alpha)

        #centralities = [hc, ec, dc, bc, cc, lc, ic, soc, cfcc, cfbc, acfbc, cbc]
        # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
        centralities = {
            'degree_centrality' : dc ,
            'harmonic_centrality' : hc ,
            'eigenvector_centrality' : ec ,
            'betweenness_centrality' : bc ,
            'closeness_centrality' : cc ,
            'load_centrality' : lc ,
            'information_centrality' : ic ,
            'communicability_betweenness_centrality' : cbc ,
            'katz_centrality_numpy' : kz ,
            'pagerank' : pr
        }
        
        # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
        centralities = pd.DataFrame( centralities )
    
    else: 
        # COMPUTA LA VERSION LOW-COST DE LAS CENTRALIDADES  
        centralities = compute_centralities_short(graph)

    return centralities

def get_alpha(G):
    """"Calcula el parametro alpha para PageRank y Katz"""
    tmp = nx.adjacency_matrix(G).toarray()

    tmp = np.linalg.eigvals(tmp)
    tmp = np.max(tmp)
    tmp = np.real(tmp)    

    tmp = .9*(1/tmp)

    return tmp

# %% --- CALCULO DE CENTRALIDADES CON NODOS REMOVIDOS

import time

def removed_nodes_centrality(graph, node, info=False):
    """Helper que remueve un nodo y calcula la centralidad de este"""
    print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()),
        '--- Iterando en nodo', node, # Sanity check de nodo iterado
        '--- Memoria en uso:', ( ray.available_resources()['memory'] - ray.cluster_resources()['memory'] )) # Uso de memoria
    
    G_removed = graph.copy() # CREA UNA COPIA DE TRABAJO PARA EL GRAFO

    G_removed.remove_node( node )  # ELIMINA EL NODO A ITERAR
    G_removed = get_largest_component(G_removed) # ELIMINA NODOS DISCONEXOS

    assert len(graph.nodes) != len(G_removed.nodes), "Ningun nodo removido."
    nodos_removidos = set(graph.nodes) - set(G_removed.nodes)

    if info:
        # INFORMACION EXTRA PARA LOS CALCULOS Y DEMAS
        print( 'Nodos originales:', len(graph.nodes) )
        print( 'Nodos post-remocion:', len(G_removed.nodes) )
        print( 'Nodos removidos:', len(graph.nodes) - len(G_removed.nodes) )
        print( 'Removidos:', nodos_removidos )

    removed_centrality = compute_centralities(G_removed, lite=False, alpha = 0.001) # CENTRALIDADES

    all_nodes = list( graph.nodes ) # REINDEXANDO PARA INCLUIR REMOVIDO
    removed_centrality = removed_centrality.reindex( all_nodes )

    print( 'SANITY_CHECK:', node, 'DIMENSIONES_DATAFRAME:', removed_centrality.shape )

    removed_centrality.name = str( node ) # RENOMBRADO DEL DATAFRAME

    # SANITY CHECK PARA VER QUE EL DATAFRAME TENGA COMO NAN LOS REMOVIDOS
    for removido in nodos_removidos:
        assert np.isnan( removed_centrality.loc[ removido , 'degree_centrality']), 'Nodo removido tiene un valor no NaN. Error de DataFrame.'

    return node, removed_centrality

# %% --- CALCULO DE CENTRALIDADES CON (ITERATIVOS) NODOS REMOVIDOS

import ray

@ray.remote
def hpc_reemoved_centrality( graph, nodes_to_remove ):
    """Calcula las centralidades para nodos removidos de forma distribuida en un cluster"""

    assert type( nodes_to_remove ) == list, 'El input de remocion no es una lista. Usa list()'

    centralities_removed = [ removed_nodes_centrality(graph, node) for node in nodes_to_remove ]

    return centralities_removed

# %% --- COMPUTA LAS BASELINE DEL SISTEMA

if __name__ == '__main__':

    # INICIALIZA LA CONEXION AL SERVIDOR DE RAY
    import os # Importa variables ambientales
    ray.init(address='auto', _node_ip_address=os.environ["ip_head"].split(":")[0], _redis_password=os.environ["redis_password"])
    #ray.init(address='auto', _redis_password='5241590000000000')

    G_ray = ray.put( G ) # Sube el grafo al object store

    # Importa la lista de subsistemas de un archivo externo
    infile = open('./tmp/subsystems_dict.pkl', 'rb'); subsystems_dict = pickle.load(infile); infile.close()

    nodos_remover = subsystems_dict['Phenylalanine metabolism'] + ['r0399']

    # EVITAR EL SOBRE-PARALELISMO DE UN PROCESO POR NODO
    WORKERS = int( ray.cluster_resources()['CPU'] ) # CANTIDAD DE CPUS DEL CLUSTER
    SPLITS = min(WORKERS, len(nodos_remover) ) # Numero de CPUs en el cluster, o nodos en la red
    
    print('Distribuyendo', len(nodos_remover), 'nodos en', SPLITS, 'procesos.') # Sanity check
    
    import random; random.seed(42) ; nodos_remover = random.sample( nodos_remover , WORKERS ) # TODO: muestra de los nodos

    nodos_remover = np.array_split( nodos_remover , SPLITS )  # Lo separa en listas mas pequenas
    nodos_remover = [ list(chunk) for chunk in nodos_remover ] # Convierte a lista de nuevo

    # CALCULO DE CENTRALIDADES PERTURBADAS (MANDA LA TAREA AL CLUSTER VIRTUAL)
    centralidades_distribuidas = [ hpc_reemoved_centrality.remote( G_ray, chunk ) for chunk in nodos_remover ]

    # INTERMEDIO EN QUE LOCALMENTE CALCULA LAS CENTRALIDADES BASE. MENOS DEMANDANTE QUE EL RESTO
    baseline = compute_centralities(G, lite=LITE, alpha=0.001) # TIRA EL CALCULO DE CENTRALIDADES BASE
    print( baseline.info(verbose=True) ) # Sanity check para el formato de las salidas

    # VENGANZA DE LOS PROCESOS ENVIADOS AL CLUSTER (GET.REMOTES)
    centralidades_perturbadas = ray.get( centralidades_distribuidas )
    centralidades_perturbadas = [ii for i in centralidades_perturbadas for ii in i] # Aplasta la lista

    centralidades_perturbadas = { cent[0] : cent[1] for cent in centralidades_perturbadas }

    # GUARDANDO RESULTADOS
    baseline = {'baseline' : baseline } # Para que queden igual
    
    import pickle
    outfile = open('./tmp/baseline.pkl', 'wb'); pickle.dump( baseline ,outfile); outfile.close()
    outfile = open('./tmp/centralidades_perturbadas.pkl', 'wb'); pickle.dump( centralidades_perturbadas ,outfile); outfile.close()
