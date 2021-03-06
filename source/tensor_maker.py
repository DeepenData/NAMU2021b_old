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
    # TODO excepción en caso de que algo no exista en la lista, por la remoción y eso

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

# %% --- Importa subsistemas

# Los subsistemas son creados como un diccionario 'subsistema' : [nodo1, nodo2, ...] en otro
# archivo, a partir de las definiciones del modelo que estamos usando. Un problema de los modelos
# es que hay que realizar curaciones manuales, así que por eso es preferible que la lista de 
# reacciones por ID este en otra parte

infile = open('./tmp/subsystems_dict','rb'); subsystems = pickle.load(infile); infile.close()

# %% --- Crea el tensor final de contribución a centralidad. 

FC_final = [ contrib_to_centrality( ns ) for ns in subsystems.values() ]

FC_final = np.asarray( FC_final ) # A tensor. Dims = (s,c,n)
FC_final = FC_final.T        # Transposición. Dims = (n,c,s)

outfile1 = open('./tmp/node_controbution_to_subsystems','wb'); pickle.dump(FC_final,outfile1); outfile1.close()

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

# %% --- 
# tensor (1000,8,2) -> tensor_plano (1000,8*2)
# %T% TMAGÍA: calculo de correlaciones (1000,8*2) -> (1000,1000)
# Suma de absolutos por columna (1000,1) // Suma centralidades

# TODO: Calculo de las correlaciones entre nodos y (centralidades * subsistema)
"""
# Import required libraries 
from scipy.stats import kendalltau 

# Taking values from the above example in Lists 
X = [1, 2, 3, 4, 5, 6, 7] 
Y = [1, 3, 6, 2, 7, 4, 5] 

# Calculating Kendall Rank correlation 
corr, _ = kendalltau(X, Y) 
print('Kendall Rank correlation: %.5f' % corr) 

# This code is contributed by Amiya Rout 
"""
# Dims: (n,(s*c))
#               s1_harmonic_centrality s2_harmonic_centrality 
# DM_10fthf5glu               0.000878               0.000878

cumulative_centrality = np.sum( np.sum( abs(FC_final) , axis=1), axis=1) # Dims =  (1054) # elimina las otras dos dimensiones

# %% --- TABLA DE SALIDA CON RESUMEN DE CENTRALIDAD AGREGADA
# [reaction_name, reaction_formula ,flux, reduced_cost]) = final_results

INPUT_MODEL = './data/stimulated_2021.json'
import warnings; warnings.filterwarnings('ignore') # Ignora warnings
from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL) # Genera infinitos warning
model.optimize()

tabla_resultado = pd.DataFrame({
    "ids" :         [rx.id           for rx in model.reactions],
    "formula" :     [rx.reaction     for rx in model.reactions],
    "flux" :        [rx.flux         for rx in model.reactions],
    "sensitivity" : [rx.reduced_cost for rx in model.reactions]
    })

tabla_resultado = tabla_resultado[tabla_resultado['ids'].isin( index_nodes )] # Elimina nodos no-conectados
tabla_resultado['cumulative_centtrality'] = cumulative_centrality # Añade centralidad agregada

OUTPUT_RESUME_TABLE = './tmp/cumulative_centrality_FBA_table.csv'
tabla_resultado.to_csv( OUTPUT_RESUME_TABLE , index=False )
