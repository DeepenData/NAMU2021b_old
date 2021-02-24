#!/bin/python
"""Convierte los objetos generados en el código anterior en un tensor para las computaciones de El Físico


"""
# %% --- Importando cosas iniciales
import pickle # Importa los elementos generados por delta_centrality.py

infile1 = open('./tmp/baseline_centralities','rb'); baseline_centralities = pickle.load(infile1); infile1.close()   # Tensor n x m x c
infile2 = open('./tmp/perturbed_centralities','rb'); perturbed_centralities = pickle.load(infile2); infile2.close() # Tensor 1 x m x c

# index_nodes : es la lista de las reacciones
infile3 = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(infile3); infile3.close() # Indice de largo m == n (todos los nodos)
# index_centralities : nombres de las centralidades calculadas
index_centralities = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality']

import numpy as np
import pandas as pd

# S = subsystem [1:s]. S.nodes = subset( N )

# R = Nodo removido [1:r]
# N = Todos los nodos del modelo [1:n]. len(R) == len(N)
# C = centrality [1:c]. Son 8 en total

# %% --- Cosas de tensores

def subsystem_tensor( subsystem_nodes, tensor, complete_index=index_nodes, axis_to_index=1):
    """Genera un tensor con un subset de nodos basados en el subsistema pedido

    Parameters
    ----------
    subsystem_nodes : list
        Lista de nodos por nombre para el subsistema. 
        Eg. `['OxPhox','GluRxn1']`
    tensor : numpy.ndarray
        Un tensor de dimensiones `(r,n,c)`, con r > 0
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

# %% --- Tensor final en si
def contrib_to_centrality( ns ): 

    ns_baselines = subsystem_tensor( ns, baseline_centralities )[0] # Set aplanado del baseline. Dims = (ns,c)
    ns_perturbed = subsystem_tensor( ns, perturbed_centralities )   # Tensor de perturbated. Dims = (r,ns,c)

    ns_mean_baselines = np.mean( ns_baselines.T, axis=1) # Dims = (c,1)
    ns_mean_perturbed = np.mean( ns_perturbed.T, axis=1) # Dims = (c,r)

    ns_perturbed_square = np.apply_along_axis( np.diag, 0, ns_mean_perturbed) # Tensor con centralidades diagonales. Dims = (c,c,r)
    ns_perturbed_square_inv = np.linalg.inv( ns_perturbed_square.T ).T # Calcula el inverso de las r matrices (c,c). Dims = (c,c,r)

    baselintes_matmul = lambda square_inv : np.matmul( square_inv, ns_mean_baselines ) # Un helper para iterativos
    ns_ratios = np.apply_along_axis( baselintes_matmul, 0, ns_perturbed_square_inv) # Ratios en primera fase
    ns_ratios = np.log2( ns_ratios ) # Fold change de baseline/perturbed. Ie, controbución a centralidad. Dims = (c,r)

    return ns_ratios

# %% --- Ejemplos de subsistemas

import random   # DEMOSTRACIÓN DE SUBSISTEMAS RANDOM
random.seed(23) # DE DISTINTO LAGRGO Y SIMILAR

# EJEMPLO DUMMY DE SUBSISTEMAS
subsystems = {
    'ns1' : [ index_nodes[n] for n in random.sample(range(0, 1055), 7)  ] ,
    'ns2' : [ index_nodes[n] for n in random.sample(range(0, 1055), 11) ]
}

# %% --- Crea el tensor final de contribución a centralidad. 

FC_final = [ contrib_to_centrality( ns ) for ns in subsystems.values() ]

FC_final = np.asarray( FC_final ) # A tensor. Dims = (s,c,n)
FC_final = FC_final.T        # Transposición. Dims = (n,c,s)

outfile1 = open('./tmp/subsystems_centrality','wb'); pickle.dump(FC_final,outfile1); outfile1.close()

# %% --- Función que genera dataframes a pedido
def fx(subsystem_name): 

    subsystem_keys  = list( subsystems.keys() )
    subsystem_index = subsystem_keys.index( subsystem_name )
 
    import pandas as pd

    subsystem_frame = FC_final[:,:, subsystem_index ]
    subsystem_frame = pd.DataFrame( subsystem_frame )

    subsystem_frame.columns = index_centralities
    subsystem_frame.index   = index_nodes

    return subsystem_frame