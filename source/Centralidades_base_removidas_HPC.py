#!/bin/python
"""
"L'impatience va éliminer l'irritabilité, la tension nerveuse et sa conséquence, le surmenage."
    - Ejemplo de uso de "surmenage", termino en frrancés para "fatiga mental"

Este es codigo para ejecutar en el cluster virtual Ray dentro del HPC. No intentes ejecutar esto en casa. 
"""
# %%% --- IMPORTING BASE LIBRARIES
import networkx as nx
import numpy    as np
import pandas   as pd

import os
import time 

import pickle

# %% --- USING THE RAY LIBRARY FOR DISTRIBUITED COMPUTING
import ray

try:
    ray.init(address='auto', _node_ip_address=os.environ["ip_head"].split(":")[0], _redis_password=os.environ["redis_password"])
except:
    print("No SLURM Cluster detected. Trying a local cluster with default config.")
    ray.init(address='auto',  dashboard_host = '0.0.0.0', _redis_password='5241590000000000')

print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), '--- Connected to Ray Cluster' )

# %% --- WARPER FUNCTIONS FOR RAY

@ray.remote
def remove_node(graph, node):
    """Generates a new graph whitout the selected node

    Args:
        graph (graph): the base graph
        node (str): the node to remove
    
    Returns:
        node (str): the node to remove as sanity check
        graph_node_removed (graph): the largest component of the graph without the node
    """

    #assert node in list( graph.nodes ), "The selected node is not part of the graph"

    G_removed = graph.copy() # CREA UNA COPIA DE TRABAJO PARA EL GRAFO
    G_removed.remove_node( node )  # ELIMINA EL NODO A ITERAR

    # Selects the largest component
    largest_component = max(nx.connected_components(G_removed), key=len)

    graph_node_removed = G_removed.subgraph(largest_component) # Subsets the largest component

    assert nx.is_connected(graph_node_removed), "The graph is not fully conected"

    return node, graph_node_removed

# As these don't take forever (except closeness centrality), is more efficient
# to group them in a single process for the distribuited computing step. 
@ray.remote
def fast_centralities(graph, node, alpha=0.0001):
    """Eigenvector, Pagerank, Katz, y Global reaching
    
    [0] Pass and return removed node as sanity check  
    [1] Dataframe whith results  
    [2] Compute time  
    """

    tic = time.time()

    ec  = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    pr  = nx.pagerank(graph, alpha = alpha)
    kz  = nx.katz_centrality(graph, alpha = alpha)
    grc = nx.global_reaching_centrality(graph) # Takes foreer

    fast_centralities_df = pd.DataFrame({
        "eigenvector_centrality" : ec ,
        "pagerank" : pr ,
        "katz_centrality" : kz ,
        "global_reaching_centrality" : grc 
    })

    toc = time.time()

    return node, fast_centralities_df , (toc - tic)

# And these take forever, so we distribuite these in individual workers

@ray.remote
def closeness_centrality(graph, node):
    tic = time.time()
    out = nx.closeness_centrality(graph, distance=None, wf_improved=True) # Takes forever
    df_closeness_centrality = pd.DataFrame({
        "closeness_centrality" : out
    })
    toc = time.time()
    return node , df_closeness_centrality,  (toc - tic)

@ray.remote
def harmonic_centrality(graph, node):
    tic = time.time()
    out = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    df_harmonic_centrality = pd.DataFrame({
        "harmonic_centrality" : out
    })
    toc = time.time()
    return node , df_harmonic_centrality,  (toc - tic)

@ray.remote
def information_centrality(graph, node):
    tic = time.time()
    out = nx.information_centrality(graph)
    df_information_centrality = pd.DataFrame({
        "information_centrality" : out
    })
    toc = time.time()
    return node , df_information_centrality, (toc - tic)

@ray.remote
def load_centrality(graph, node):
    tic = time.time()
    out = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
    df_load_centrality = pd.DataFrame({
        "load_centrality" : out
    })
    toc = time.time()
    return node , df_load_centrality , (toc - tic)

@ray.remote
def betweenness_centrality(graph, node):
    tic = time.time()
    out = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
    print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()),
            '--- Resolved', node) # This is the slowest, so every other should be done 
    df_betweenness_centrality = pd.DataFrame({
        "betweenness_centrality" : out
    })
    toc = time.time()
    return node , df_betweenness_centrality , (toc - tic)

# %% --- IMPORTING BASE GRAPH

try: 
    G = nx.read_gpickle('./data/stimulated_2021.gpickle')
except:
    G = nx.read_graphml('./data/stimulated_2021.graphml')

print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), '--- Reading gpickle graph' )

G_ray = ray.put(G) # BASE GRAPH TO OBJECT STORE

# %% --- CREATING THE GRAPHS WITH REMOVED NODES

#infile = open('./tmp/subsystems_dict.pkl', 'rb'); subsystems_dict = pickle.load(infile); infile.close()
#pku_set = pd.read_csv('pku_noise.csv', index_col=0) # Fataframe whith a list of nodes
print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), '--- Reading list to remove' )

# List of nodes to remove
NODES_REMOVED = list( G.nodes )

# This launches Ray actors to generate graphs whitout one node from the list
graph_removed = [ remove_node.remote( G_ray , node ) for node in NODES_REMOVED ]
graph_removed = ray.get( graph_removed )

# This sets the returned list of lists in a dictionay { node : graph_node_removed }
graph_removed = { rm[0] : rm[1] for rm in graph_removed }

print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), '--- Done node removal' )

# %% --- LAUNCHING THE CENTRALITY CALCULATIONS
# This launches to the distribuited process in the Ray cluster
fast = [ fast_centralities.remote( graph_removed[node] , node )      for node in graph_removed] 
harc = [ harmonic_centrality.remote( graph_removed[node] , node )    for node in graph_removed] 
info = [ information_centrality.remote( graph_removed[node] , node ) for node in graph_removed] 
load = [ load_centrality.remote( graph_removed[node] , node )        for node in graph_removed] 
btwn = [ betweenness_centrality.remote( graph_removed[node] , node ) for node in graph_removed] 
clos = [ closeness_centrality.remote( graph_removed[node] , node )   for node in NODES_REMOVED] 


# %% --- RETURNING THE CENTRALITIES
# Lists of [[node, dataframe, time], ...] for NODES_REMOVED
fast = ray.get( fast ) 
harc = ray.get( harc ) 
info = ray.get( info ) 
load = ray.get( load ) 
btwn = ray.get( btwn ) 
clos = ray.get( clos )

# %% --- CENTRALITIES TO DATAFRAMES
# Dictionaries to track node and dataframe

dict_df = lambda list_df : dict( zip( [ node[0] for node in list_df ],[ node[1] for node in list_df ]) )

fast_df = dict_df( fast )
harc_df = dict_df( harc )
info_df = dict_df( info )
load_df = dict_df( load )
btwn_df = dict_df( btwn )
clos_df = dict_df( clos )

# Joining dataframes onto single dataframe for each removed node
def join_centrality_dataframes( node ):
    tmp = fast_df[node] # Loads the first dataframe
    tmp = tmp.join( harc_df[node] , how="inner") # joins the harmonic_centrality
    tmp = tmp.join( info_df[node] , how="inner") # joins the information_centrality
    tmp = tmp.join( load_df[node] , how="inner") # joins the load_centrality
    tmp = tmp.join( btwn_df[node] , how="inner") # joins the betweenness_centrality
    tmp = tmp.join( clos_df[node] , how="inner") # joins the closeness_centrality

    return tmp

centralities_final_dict = { node : join_centrality_dataframes( node ) for node in NODES_REMOVED } 

import os; os.makedirs("./tmp", exist_ok=True) # crea .tmp si no existe
outfile = open('./tmp/centralidades_perturbadas.pkl', 'wb'); pickle.dump( centralities_final_dict , outfile); outfile.close()

print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), '--- Done centrality for removed nodes' )

# %% --- TIMING OF THE CALCULATIONS
dict_time = lambda list_df : dict( zip( [ node[0] for node in list_df ],[ node[2] for node in list_df ]) )

fast_time = dict_time( fast )
harc_time = dict_time( harc )
info_time = dict_time( info )
load_time = dict_time( load )
btwn_time = dict_time( btwn )
clos_time = dict_time( clos )

df_times = pd.DataFrame({
    "fast_and_closeness_centralities" : fast_time ,
    "harmonic_centrality" : harc_time ,
    "information_centrality" : info_time ,
    "load_centrality" : load_time ,
    "betweenness_centrality" : btwn_time ,
    "closeness_centrality" : clos_time 
})

print( df_times ) # Info para la terminal

df_times.to_csv("tmp.time_removed_nodes.csv")

print( time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), '--- Done times for centralities' )
