#!/bin/pypy
"""
"L'impatience va éliminer l'irritabilité, la tension nerveuse et sa conséquence, le surmenage."
    - Ejemplo de uso de "surmenage", termino en frrancés para "fatiga mental"
"""
# %% --- IMPORTA EL GRAFO GENERADO

import networkx as nx
import numpy    as np

G = nx.read_graphml("./tmp/graph.graphml") # Lee el modelo desde un grafo común

# %% --- DEFINICIONES DE FUNCIONES

def compute_centralities(graph):
    """Computa las doce centralidades y devuelve una DataFrame (no reindexado) con estas"""
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

    return centralities

baseline = compute_centralities(G)
print( baseline.info(verbose=True) )
