#!/bin/python
"""
"L'impatience va éliminer l'irritabilité, la tension nerveuse et sa conséquence, le surmenage."
    - Ejemplo de uso de "surmenage", termino en frrancés para "fatiga mental"
"""
# %% --- IMPORTA EL GRAFO GENERADO

import networkx as nx
import numpy    as np
import pandas   as pd  

LITE=True # Variable para definir si los calculos deben ser complejos o no

#G = nx.read_graphml("./tmp/graph.graphml") # Lee el modelo desde un grafo común

# %% --- ELIMINAR --- DEFINICIONES DE FUNCIONES
# TODO: resolver esto y el error de ModelToGrapgh.py

INPUT_MODEL='./data/stimulated_2021.json' # TODO: hacer esto una variable ambiental

import warnings; warnings.filterwarnings('ignore') # Ignora warnings
from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL)

def cobra_to_networkx_rxn_projection(modelo):
    import networkx as nx
    from   cobra.util.array import create_stoichiometric_matrix
    import numpy as np
    from sklearn.preprocessing import Binarizer
    import warnings 
    warnings.filterwarnings("ignore")

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"
    #extraer matriz estequiométrica
    S_matrix = create_stoichiometric_matrix(modelo)
    #convertir todas las entradas valores positivos
    S_matrix = (abs(S_matrix) )
    #transformar a enteros 
    S_matrix = S_matrix.astype(np.int)
    #Multiplicacion por la derecha para proyectar en el espacio de las reacciones
    projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
    #rellenar diagonal con ceros
    np.fill_diagonal(projected_S_matrix, 0) 
    #binarizar
    projected_S_matrix = Binarizer().fit_transform(projected_S_matrix)
    #crear grafo networkx
    G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
    #hacer diccionario con los nombres de las reacciones
    node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
    reaction_dict = node_dict( [reaction.id for reaction in modelo.reactions] )
    #Renombrar los nodos usando el diccionario
    G = nx.relabel_nodes(G, reaction_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
    return G

G = cobra_to_networkx_rxn_projection(model)

def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G

G = get_largest_component(G) # Elimina otras cosas

# %% --- DEFINICIONES DE FUNCIONES

def compute_centralities_short(graph):
    """Computa las centralidades rapidas y devuelve una DataFrame (no reindexado) con estas"""
    # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
    hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc    = nx.degree_centrality(graph)
    cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    ic    = nx.information_centrality(graph) # Requiere scipy

    # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
    centralities = {
        'harmonic_centrality' : hc ,
        'eigenvector_centrality' : ec ,
        'degree_centrality' : dc ,
        'closeness_centrality' : cc ,
        'information_centrality' : ic
    }

    # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
    centralities = pd.DataFrame( centralities )

    return centralities

def compute_centralities(graph, lite=False):
    """Computa las doce centralidades y devuelve una DataFrame (no reindexado) con estas"""
    if lite == False:
        # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
        hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
        ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
        dc    = nx.degree_centrality(graph)
        bc    = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
        cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
        lc    = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
        ic    = nx.information_centrality(graph) # Requiere scipy
        soc   = nx.second_order_centrality(graph) 
        cfcc  = nx.current_flow_closeness_centrality(graph) # Requiere scipy
        cfbc  = nx.current_flow_betweenness_centrality(graph)
        acfbc = nx.approximate_current_flow_betweenness_centrality(graph)
        cbc   = nx.communicability_betweenness_centrality(graph) 

        #centralities = [hc, ec, dc, bc, cc, lc, ic, soc, cfcc, cfbc, acfbc, cbc]
        # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
        centralities = {
            'harmonic_centrality' : hc ,
            'eigenvector_centrality' : ec ,
            'degree_centrality' : dc ,
            'betweenness_centrality' : bc ,
            'closeness_centrality' : cc ,
            'load_centrality' : lc ,
            'information_centrality' : ic ,
            'second_order_centrality' : soc ,
            'current_flow_closeness_centrality' : cfcc ,
            'current_flow_betweenness_centrality' : cfbc ,
            'approximate_current_flow_betweenness_centrality' : acfbc ,
            'communicability_betweenness_centrality' : cbc 
        }
        
        # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
        centralities = pd.DataFrame( centralities )
    
    else: 
        # COMPUTA LA VERSION LOW-COST DE LAS CENTRALIDADES  
        centralities = compute_centralities_short(graph)

    return centralities

# %% --- CALCULO DE CENTRALIDADES CON NODOS REMOVIDOS

import time

def removed_nodes_centrality(graph, node, info=False):
    """Helper que remueve un nodo y calcula la centralidad de este"""
    print( time.localtime(time.time()), '--- Iterando en nodo', node) # Sanity check de nodo iterado
    
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

    removed_centrality = compute_centralities(G_removed, lite=LITE) # CENTRALIDADES

    all_nodes = list( graph.nodes ) # REINDEXANDO PARA INCLUIR REMOVIDO
    removed_centrality = removed_centrality.reindex( all_nodes )

    removed_centrality.name = str( node ) # RENOMBRADO DEL DATAFRAME

    # SANITY CHECK PARA VER QUE EL DATAFRAME TENGA COMO NAN LOS REMOVIDOS
    for removido in nodos_removidos:
        assert np.isnan( removed_centrality.loc[ removido , 'harmonic_centrality']), 'Nodo removido tiene un valor no NaN. Error de DataFrame.'

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

    G_ray = ray.put( G ) # Sube el grafo al object store

    nodos_remover = list( G.nodes )
    #nodos_remover = Glycolysis_astrocyte # TODO: Prueba a escala con solo glico_astrocitos

    # EVITAR EL SOBRE-PARALELISMO DE UN PROCESO POR NODO
    WORKERS = int( ray.cluster_resources()['CPU'] ) # CANTIDAD DE CPUS DEL CLUSTER
    SPLITS = min(WORKERS, len(nodos_remover) ) # Numero de CPUs en el cluster, o nodos en la red
    
    print('Distribuyendo', len(nodos_remover), 'nodos en', SPLITS, 'procesos.') # Sanity check
    
    nodos_remover = np.array_split( nodos_remover , SPLITS )  # Lo separa en listas mas pequenas
    nodos_remover = [ list(chunk) for chunk in nodos_remover ] # Convierte a lista de nuevo

    # CALCULO DE CENTRALIDADES PERTURBADAS (MANDA LA TAREA AL CLUSTER VIRTUAL)
    centralidades_distribuidas = [ hpc_reemoved_centrality.remote( G_ray, chunk ) for chunk in nodos_remover ]

    # INTERMEDIO EN QUE LOCALMENTE CALCULA LAS CENTRALIDADES BASE. MENOS DEMANDANTE QUE EL RESTO
    baseline = compute_centralities(G, lite=LITE) # TIRA EL CALCULO DE CENTRALIDADES BASE
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
