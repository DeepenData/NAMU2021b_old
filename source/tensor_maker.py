#!/bin/python
"""Convierte los objetos generados en el código anterior en un tensor para las computaciones de El Físico


"""
# %% --- Importando cosas iniciales
import pickle # Importa los elementos generados por delta_centrality.py

outfile1 = open('./tmp/baseline_centralities','rb'); baseline_centralities = pickle.load(outfile1); outfile1.close()   # Tensor n x m x c
outfile2 = open('./tmp/perturbed_centralities','rb'); perturbed_centralities = pickle.load(outfile2); outfile2.close() # Tensor 1 x m x c

outfile3 = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(outfile3); outfile3.close() # Indice de largo m == n (todos los nodos)
# Index de las centralidades calculadas
index_centralities = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality']

# index_nodes : es la lista de las reacciones
outfile4 = open('./tmp/breaks','rb'); breaks = pickle.load(outfile4); outfile4.close() # Lista de exceptciones <n

import numpy as np
import pandas as pd

def rcn(tensor): 
    """Convierte un tensor de dimensiones `(n, m, c)` a uno `(n, c, m)`."""
    import numpy as np
    return  np.transpose( tensor, (0,2,1) )

# R = Nodo removido [1:r]
# N = Todos los nodos del modelo [1:n]. len(R) == len(N)
# C = centrality [1:c]. Son 8 en total

# S = subsystem [1:s]. S.nodes = subset( N )

# %% --- Cosas de tensores

# A_n_c: Centralidades (c = 1:8) de todos los nodos (m=1:1000). Dims = 1000 x 8. // Centralidades base
A_n_c = baseline_centralities[0] # (m,c) 

# A_1_n_c: Centralidades base como tensor de alto 1. Dims = 1 x 1000 x 8
A_1_n_c = baseline_centralities # (1,m,c)

# B_r_n_c: Centralidades (c) de los nodos  (n=1, ... NA(en la posición de nodo r) ...,1000) removiendo el nodo r. 
B_r_n_c = perturbed_centralities # Dims: (r, n, c) ~ 1000 x 999+(NaN) x 8

# C_r_n_c se omite porque parece que C es de centralidades y es confuso y caotico 

# D_r_c_n: Es lo mismo que B_r_n_c pero con las centralidades en la segunda dimensión.
D_r_c_n = rcn( B_r_n_c )

get_n( [] )
get_c( [] )
get_n( [] )

# %%
def subsystem_tensor( subsystem_nodes, tensor=perturbed_centralities, complete_index=index_nodes, axis_to_index=1):
    """Genera un tensor con un subset de nodos basados en el subsistema pedido

    Parameters
    ----------
    subsystem_nodes : list
        Lista de nodos por nombre para el subsistema. 
        Eg. `['OxPhox','GluRxn1']`
    tensor : numpy.ndarray, default perturbed_centralities
        Un tensor de dimensiones `(r,n,c)`, by default perturbed_centralities
    complete_index : list, default index_nodes
        La lista con el indice por nombre de los nodos, tiene que ser del mismo
        largo que la dimensión `n` del tensor 
    axis_to_index : int, default 1
        La dimensión en la que se aplica el indexado. Un tensor normal es de 
        forma `(r,n,c)`, pero puede ser transpuesto de forma que n ya no este en
        la posición 1 

    Returns
    -------
    tensor : numpy.ndarray
        Un tensor de dimensiones (r,ns,c), creado a partir de la lista de nodos
        del subsistema. No es el mismo objeto que el tensor original, así que
        puede ser un problema con la memoria. 
    """
    import numpy as np

    # Index por posición a partir de una lista de reacciones
    subsystem_index = [ complete_index.index(node) for node in subsystem_nodes ]

    return tensor[:, subsystem_index ,:]

import random   # DEMOSTRACIÓN DE SUBSISTEMAS RANDOM
random.seed(23) # DE DISTINTO LAGRGO Y SIMILAR

ns1 = [ index_nodes[n] for n in random.sample(range(0, 1055), 7)  ]
ns1 = [ index_nodes[n] for n in random.sample(range(0, 1055), 11) ]

# A_ns_c: Lo mismo que A_n_c pero indexado para ns.
A_ns1_c = subsystem_tensor( ns1, A_1_n_c )[0]
A_ns2_c = subsystem_tensor( ns2, A_1_n_c )[0]

# B_r_ns_c: Lo mismo que B_r_n_c pero indexado para ns. // B_n_c[ns, *, *]
B_r_ns1_c = subsystem_tensor( ns1, B_r_n_c )
B_r_ns2_c = subsystem_tensor( ns2, B_r_n_c )

# %%
# E_c: Promedio por centralidad de todos los nodes en ms. 
E_s1_c = np.mean( A_ns1_c.T, axis=1, keepdims=True ) # .shape = (c,1)
E_s2_c = np.mean( A_ns2_c.T, axis=1, keepdims=True ) # .shape = (c,1)

#F_c_n:Promedio por centralidad de todos los nodes en ms removiendo n. // E_c.mean(axis='n')
F_s1_c_n = np.mean( B_r_ns1_c , axis = 1).T #, keepdims=True )
F_s2_c_n = np.mean( B_r_ns2_c , axis = 1).T #, keepdims=True )
# Las dimensiones correctas (c,n)

F_s1_c_c_n = np.apply_along_axis( np.diag, 0, F_s1_c_n) # Tensor con centralidades diagonales x nodos (c,c,n)
F_s1_inv_c_c_n = np.linalg.inv(F_s1_c_c_n.T).T # Calcula el inverso de las n matrices (c,c).  shape = (c,c,n)

F_s2_c_c_n = np.apply_along_axis( np.diag, 0, F_s2_c_n) # Tensor con centralidades diagonales x nodos (c,c,n)
F_s2_inv_c_c_n = np.linalg.inv(F_s2_c_c_n.T).T # Calcula el inverso de las n matrices (c,c).  shape = (c,c,n)

#F_inv_c_c_n.shape

capa_i_s1 = lambda capa_i : np.matmul( capa_i , E_s1_c ) 
capa_i_s2 = lambda capa_i : np.matmul( capa_i , E_s2_c ) 

Ratios_s1 = np.apply_along_axis( capa_i, 0, F_s1_inv_c_c_n)
Ratios_s2 = np.apply_along_axis( capa_i, 0, F_s2_inv_c_c_n)

FC_s1  = np.log2(Ratios_s1[0].T)
FC_s2  = np.log2(Ratios_s2[0].T)

FC_final = np.asarray( [FC_s1, FC_s2] )      # Shape (s,n,c)
FC_final = np.transpose( FC_final, (2,1,0) ) # Shape (n,c,s)

subsystems = {
    'name1' : [reaction1, reaction2, reaction3], 
    'name2' : [reaction3, reaction5, reaction9]
}

get_s(): 
 Returs 
    
    dataframe 
        nodenames, centralitnames
