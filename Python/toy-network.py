"""
Descripción: un archivo de redes de jugete
"""

# %% --- Importando librerias y cosas
import networkx as nx
from   networkx.algorithms import bipartite # Analisis de redes bipartitas
from   networkx.algorithms.bipartite.matrix import from_biadjacency_matrix # Crea el grafo
from   networkx.algorithms.bipartite.matrix import biadjacency_matrix      # Extrae matriz de adyacencia

import numpy as np # NumPy
from scipy.sparse import csr_matrix # NetworkX necesita una matriz sparce
from sklearn.preprocessing import binarize # Una función util de un paquete que no usaremos
# %% --- Creando la matriz
biadjacency_matrix =binarize(abs( \
 np.matrix('1 -1  0  0  0 -1  1  0; \
            0  1  1 -1 -1  0  0  0; \
            0  0  0  0  1  1 -1 -1; \
            0  0  0  0  0  1 -1  1; \
            0  0  0  0  0 -1  1 -1') ) )
# Por mientras la matriz está sin dirección
biadjacency_matrix # Observando la matriz como una array 2d
# %% --- Crear el objeto de matriz sparce con SciPy
# csr_matrix(biadjacency_matrix) # Para ver una salida
sparse_biad_m   = csr_matrix(biadjacency_matrix) # Matriz de adyacencia sparce
bipartite_metab = from_biadjacency_matrix(sparse_biad_m) # Crea la bipartita
bipartite_graph               = bipartite_metab.copy() # Usamos copy para evitar choques de memoria
# %% --- Mostrando diagrama de grafos
bipartite_graph.nodes(data=True) # Requiere una ventana para verlo
# %% --- Extracción de los nodos del grafo bipartito
""" Para extraer un nodo definido
algunnodo = 12
bipartite_graph.nodes(data=True)[algun_nodo]["bipartite"]
"""

total_nodes = range(0,len(bipartite_graph.nodes(data=False))) # Rango de nodos
particiones = [bipartite_graph.nodes(data=True)[i]["bipartite"] for i in total_nodes] # ListComprehensión de

# %% ---- 
print(
particiones.count(0), # Cuenta cuantos ceros hay en la lista
particiones.count(1),
biadjacency_matrix.shape) # Cuenta cuantos unos  hay en la lista

# MOMENTO EN QUE DECIDIMOS PASAR A USAR PuLP y GLPK como solvers de la matriz que CONDApy odia

# %% --- 
reactions_nodes   = [n for n, d in bipartite_graph.nodes(data=True) if d["bipartite"] == 1] # 
metabolites_nodes = [n for n, d in bipartite_graph.nodes(data=True) if d["bipartite"] == 0] # 

# %% --- Poniendole nombres a metabolitos y reacciones
reaction_names    = ["r1","r2","r3","r4","r5","r6","r7","r8"]
metabolites_names = ["m1","m2","m3","m4","m5"]
# %% --- 

names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reaction_names  ))
    # el formato + permite uni listas; como vectores en R; y zip crea tuplas

bipartite_graph = nx.relabel_nodes(bipartite_graph, names_mapped)
# bipartite_graph.nodes(data=True) # VIS: 

# %% --- Creación de subsistemas
"""
Cuando tenemos una red metabolica tenemos subsistemas metabolicos, como las rutas de glicolisis, gluconeogenesis, etc. 
Esto es util para cre
"""
A_sub = dict.fromkeys( ["r1","r2","r3"],{'subsystem': 'A_sub'})
B_sub = dict.fromkeys( ["r4","r5","r6","r7","r8"],{'subsystem': 'B_sub'})

nx.set_node_attributes(bipartite_graph, {**A_sub,**B_sub} )
    # **A_sub,**B_sub hace una concatenación de diccionarios
bipartite_graph.nodes(data=True)
# %% --- Exportando a Gephi 
nx.write_gexf(bipartite_graph, "bipartite_test_1.gexf") # Guarda a espacio local

"""TODO: 
- Analsis de centralidades
- Codigo para hacer graficos bonitos"""

# --- COSAS DE ANALISIS DE BLUJO USANDO EL MODELO ---
