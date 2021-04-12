
"""Esta es una implementación de codigo Naive para paralelismo del calculo de
centralidades. Como resultado, solo las listadas aqui convergen. Ademas, recomiendo
dividir las centralidades segun los workers al final de esto, para evitar demasiados
procesos IDLE por calculos que terminan primero/despues."""

# %% --- IMPORTA LIBRERIAS Y COSAS VARIAS
import networkx as nx
import numpy    as np
import pandas   as pd  

# %% --- DEFINE FUNCIONES DE CALCULO DE TIEMPO

from time import time

tic = time()
toc = time(); print('Done!')

print( toc - tic )

# %% --- RAY START

import ray
ray.init( 
    dashboard_host = '0.0.0.0' , # Accesible para todos
    num_cpus= 9 # Son 9 centralidades a calcular
    )

# %% --- DEFINCIÓN DE FUNCIONES PARA PARALELISMO 

@ray.remote
def harmonic_centrality(graph):
    tic = time()
    out = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    toc = time(); print('Done! - harmonic_centrality')
    return out , (toc - tic)

@ray.remote
def eigenvector_centrality(graph):
    tic = time()
    out = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    toc = time(); print('Done! - eigenvector_centrality')
    return out , (toc - tic)

@ray.remote
def betweenness_centrality(graph):
    tic = time()
    out = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
    toc = time(); print('Done! - betweenness_centrality')
    return out , (toc - tic)

@ray.remote
def closeness_centrality(graph):
    tic = time()
    out = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    toc = time(); print('Done! - closeness_centrality')
    return out , (toc - tic)

@ray.remote
def load_centrality(graph):
    tic = time()
    out = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
    toc = time(); print('Done! - load_centrality')
    return out , (toc - tic)

@ray.remote
def katz_centrality(graph):
    tic = time()
    out = nx.katz_centrality(graph, alpha = 0.00001)
    toc = time(); print('Done! - katz_centrality')
    return out , (toc - tic)

@ray.remote
def pagerank(graph):
    tic = time()
    out = nx.pagerank(graph, alpha = 0.00001)
    toc = time(); print('Done! - pagerank')
    return out , (toc - tic)

@ray.remote
def communicability_betweenness_centrality(graph):
    tic = time()
    out = nx.communicability_betweenness_centrality(graph)
    toc = time(); print('Done! - communicability_betweenness_centrality')
    return out , (toc - tic)

@ray.remote
def information_centrality(graph):
    tic = time()
    out = nx.information_centrality(graph)
    toc = time(); print('Done! - information_centrality')
    return out , (toc - tic)

@ray.remote
def global_reaching_centrality(graph):
    tic = time()
    out = nx.global_reaching_centrality(graph)
    toc = time(); print('Done! - global_reaching_centrality')
    return out , (toc - tic)

# %% --- IMPORTA LIBRERIAS Y COSAS VARIAS
G = nx.read_gpickle('./data/Recon2_rxn_proyected.gpickle')
# G = nx.read_graphml('./data/recon2_directed_bipartite_graph.graphml')
# G = nx.algorithms.bipartite.generators.random_graph(600, 600, 0.0002048, seed=1234, directed=True)

print( nx.is_connected(G) )

G_ray = ray.put(G)

# %% --- INICIO DEL PROCESO PARALELO

hc  = harmonic_centrality.remote( G_ray )
ec  = eigenvector_centrality.remote( G_ray )
bc  = betweenness_centrality.remote( G_ray )

cc  = closeness_centrality.remote( G_ray )
lc  = load_centrality.remote( G_ray )
kz  = katz_centrality.remote( G_ray )

pr  = pagerank.remote( G_ray )
grc = global_reaching_centrality.remote( G_ray )
ic  = information_centrality.remote( G_ray )

# %% --- GETS DE LOS PROCESOS PARALELOS

hc  = ray.get( hc  )
ec  = ray.get( ec  )
bc  = ray.get( bc  )

cc  = ray.get( cc  )
lc  = ray.get( lc  )
kz  = ray.get( kz  )

pr  = ray.get( pr  )
grc = ray.get( grc )
ic  = ray.get( ic  )

# %% --- POSTPROCESADO

process = [
    hc , ec , bc ,
    cc , lc , kz , 
    pr , grc, ic 
]

centrality_names = [
"harmonic_centrality" , "eigenvector_centrality" , "betweenness_centrality" , 
"closeness_centrality" , "load_centrality" , "katz_centrality_numpy" , 
"pagerank" , "information_centrality" , "global_reaching_centrality"  
]

tiempos = pd.Series([ i[1] for i in process], index=centrality_names)


# %%% --- CENTRALIDADES

centralities = {
    'harmonic_centrality' : hc[0] ,
    'eigenvector_centrality' : ec[0] ,
    'betweenness_centrality' : bc[0] ,
    'closeness_centrality' : cc[0] ,
    'load_centrality' : lc[0] ,
    'katz_centrality_numpy' : kz[0] ,
    'pagerank' : pr[0] ,
    'information_centrality' : ic[0] ,
    'global_reaching_centrality' : grc[0]
}

# CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
centralities = pd.DataFrame( centralities )

centralities.to_csv('./tmp/proyección.csv')

# %% --- PLOTLY PARA TIEMPOS

# Distribución sugerida de calculos. Todos los workers demoran ~140m
centrality_workers = [
"Worker_B" , "Worker_A" , "Worker_F" , 
"Worker_A" , "Worker_D" , "Worker_A" , 
"Worker_A" , "Worker_C" , "Worker_A"  
]

df_tiempos = pd.DataFrame({
    'worker' : centrality_workers,
    'centrality_names' : centrality_names,
    'tiempos' : tiempos/60
})

df_tiempos.to_csv('./tmp/tiempos_computo.csv')

import plotly.express as px

fig = px.treemap(df_tiempos, path=['worker', 'centrality_names'], 
    values='tiempos',
    title="Tiempo de computo [min]",
    )

fig.update_traces(textinfo="label+text+value", selector=dict(type='treemap'))

fig.show()

