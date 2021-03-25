#!/bin/python

# %% --- 
def cobra_to_networkx_metabolite_projection(modelo):
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
    projected_S_matrix = np.matmul(S_matrix, S_matrix.T)
    #rellenar diagonal con ceros
    np.fill_diagonal(projected_S_matrix, 0) 
    #binarizar
    projected_S_matrix = Binarizer().fit_transform(projected_S_matrix)
    #crear grafo networkx
    G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
    #hacer diccionario con los nombres de las reacciones
    node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
    name_dict = node_dict( [metabolite.id for metabolite in modelo.metabolites] )
    #Renombrar los nodos usando el diccionario
    G = nx.relabel_nodes(G, name_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
    return G
def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G


def compute_centralities_short(graph):
    import pandas as pd
    import networkx as nx
    """Computa las centralidades rapidas y devuelve una DataFrame (no reindexado) con estas"""
    # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
    hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc    = nx.degree_centrality(graph)
    cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    ic    = nx.information_centrality(graph) # Requiere scipy
    kz    = nx.katz_centrality_numpy(graph)
    pr    = nx.pagerank(graph)
#    bc    = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
#    cbc   = nx.communicability_betweenness_centrality(graph) 


    # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
    centralities = {
        'degree' : dc,
        'harmonic' : hc ,
        'eigenvector' : ec ,
        'closeness' : cc ,
        'information' : ic,
        'katz': kz,
        'pagerank':pr#,
    #    'betweenness': bc,
    #    'communicability': cbc
    }

    # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
    centralities = pd.DataFrame( centralities )

    return centralities
# %% --- 

import cobra  
import warnings 
warnings.filterwarnings("ignore")
recon = cobra.io.load_json_model("data/GEM_Recon2_thermocurated_redHUMAN.json")

# %%
G = cobra_to_networkx_metabolite_projection(recon)

G = get_largest_component(G)

# %%
cents = compute_centralities_short(G)

# %%

cents.loc['atp_c']/cents.loc['mal_L_c']


# %%
cents.to_csv('data/recon2_metabolite_centralities_metabolome_weights.csv')

# %%
