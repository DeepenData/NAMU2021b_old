#!/bin/python
"""
"An Ixian machine? You defy the Jihad!" 
"There's a lesson in that, too. What do such machines really do? They increase the number of things we can do without thinking. Things we do without thinking — there's the real danger."
    - God-Emperor Leto II 
"""

# %% --- Importar resultados
import pickle 

# Esto asume que está en 'human-metnet/'
infile = open('./tmp/centralidades_perturbadas.pkl','rb'); perturbed_centralities = pickle.load(infile); infile.close()
infile = open('./tmp/baseline.pkl' ,'rb'); baseline_centralities  = pickle.load(infile); infile.close()
infile = open('./tmp/subsystems_dict','rb'); subsystems = pickle.load(infile); infile.close()

# %% --- Importando librerias utiles
import pandas   as pd
import numpy    as np
import networkx as nx

# %% --- SANITY CHECK Reordena el diccionario de centralidades perturbadas

orden_reacciones = baseline_centralities['baseline'].index # Orden como la salida de los nodos
perturbed_centralities = { index : perturbed_centralities[index] for index in orden_reacciones }

# %% --- FULL TENSOR

def full_tensor(diccionarios):
    """Genera un tensor desde el diccionario"""
    subsys = [ diccionarios[ KEY ].to_numpy() for KEY in list(diccionarios.keys()) ]
    subsys = np.asarray( subsys )
    return subsys

baseline_centralities_tensor  = full_tensor(baseline_centralities)
perturbed_centralities_tensor = full_tensor(perturbed_centralities)

import pickle
outfile = open('./tmp/baseline_tensor.pkl', 'wb'); pickle.dump( baseline_centralities_tensor ,outfile); outfile.close()
outfile = open('./tmp/centralidades_perturbadas_tensor.pkl', 'wb'); pickle.dump( perturbed_centralities_tensor ,outfile); outfile.close()

# %% --- Cosas de los subsistemas ???

def subsystem_tensor( subsystem_nodes, diccionarios):
    """Genera un tensor con un subset de nodos basados en el subsistema pedido"""
    subsys = [ diccionarios[ KEY ].loc[ subsystem_nodes ,].to_numpy() for KEY in list(diccionarios.keys()) ]
    subsys = np.asarray( subsys )
    return subsys

# TODO: eliminar esto? 
# Une dos subsistemas
subsystem = subsystems['glicolisis_astros']

perturbed_subsystem = subsystem_tensor( subsystem, perturbed_centralities )
baseline_subsystem  = subsystem_tensor( subsystem, baseline_centralities   )

# %% --- Colapsando la segunda dimensión... (nodos totales)
# Esto calcula 4 promedios de las centralidades para los nodos

from numba import jit

# Un helper porque estas funciones no pueden lidiar con NaNs porque son superiores a los NaN o algo así
noNaN_perturbed_subsystem = np.nan_to_num( perturbed_subsystem , copy=True, nan=0.0, posinf=None, neginf=None)

# Es el unico que no require lambdas raros :D
aritmetic_perturbed = np.nanmean( noNaN_perturbed_subsystem, axis= 1)
aritmetic_baseline  = np.nanmean(        baseline_subsystem, axis= 1)

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
@jit(nopython=True)
def geometric_mean(array):
    array = array[~np.isnan(array)] # Filtra los NaN
    fcorr = len(array[array > 0])/len(array) # Factor de corrección
    array = array[array > 0] # Selecciona solo mayores a cero
    try: gmean = ( np.prod( array )**(1 / len(array)) ) * fcorr # Geometrica por factor de corrección
    except: return np.nan
    else : return gmean

geometric_perturbed = np.apply_along_axis( geometric_mean, 1, noNaN_perturbed_subsystem)
geometric_baseline  = np.apply_along_axis( geometric_mean, 1,        baseline_subsystem)

# Define una función helper para cuadraticos (Root Square Mean). No en scipy base
@jit(nopython=True)
def quadratic_mean(array):
    mean =  np.sqrt( np.nanmean( array*array ) )
    return mean

quadratic_perturbed = np.apply_along_axis( quadratic_mean, 1, noNaN_perturbed_subsystem)
quadratic_baseline  = np.apply_along_axis( quadratic_mean, 1,        baseline_subsystem)

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
@jit(nopython=True)
def harmonic_mean(array):
    array = array[~np.isnan(array)] # Filtra los NaN
    fcorr = len(array[array > 0])/len(array) # Factor de corrección
    array = array[array > 0] # Selecciona solo mayores a cero
    try : hmean = len(array)/np.sum(1 / array ) * fcorr # Geometrica por factor de corrección
    except: return np.nan
    else : return hmean

harmonic_perturbed  = np.apply_along_axis( harmonic_mean,  1, noNaN_perturbed_subsystem)
harmonic_baseline   = np.apply_along_axis( harmonic_mean,  1,        baseline_subsystem)

# %% --- Divisiones! Calculando como un nodo afecta conectividad de otro
# Calcula las metricas para definir los efectos en los nodos

fold_change_aritmetic = np.log2( aritmetic_baseline / aritmetic_perturbed )
fold_change_geometric = np.log2( geometric_baseline / geometric_perturbed )
fold_change_quadratic = np.log2( quadratic_baseline / quadratic_perturbed )
fold_change_harmonic  = np.log2( harmonic_baseline  / harmonic_perturbed  )

delta_aritmetic = ( aritmetic_baseline - aritmetic_perturbed ) / aritmetic_baseline
delta_geometric = ( geometric_baseline - geometric_perturbed ) / geometric_baseline
delta_quadratic = ( quadratic_baseline - quadratic_perturbed ) / quadratic_baseline
delta_harmonic  = ( harmonic_baseline  - harmonic_perturbed  ) / harmonic_baseline 

# %% --- Convierte a Dataframes

index_centralities = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality',
      'current_flow_closeness_centrality', 'current_flow_betweenness_centrality',
       'approximate_current_flow_betweenness_centrality', 'communicability_betweenness_centrality']

def crear_dataframe( matriz ):
    df = pd.DataFrame(
        matriz,
        index = baseline_centralities['baseline'].index,
        columns = index_centralities
    )
    return df

# %% --- Exporta a R

# TODO: eliminar este código estúpido que no deberia existir
# me da asco verlo y solo existe por la capacidad multi-cursor de VS Code
#   - Manu, 2021-03-13 23:04

exportar = [ 
    fold_change_aritmetic , fold_change_geometric , 
    fold_change_quadratic , fold_change_harmonic , 
    delta_aritmetic , delta_geometric , 
    delta_quadratic , delta_harmonic ]

exportar_como = [ 
    'fold_change_aritmetic' , 'fold_change_geometric' , 
    'fold_change_quadratic' , 'fold_change_harmonic' , 
    'delta_aritmetic' , 'delta_geometric' , 
    'delta_quadratic' , 'delta_harmonic' ]

for i in range(len(exportar)):
    filename = './results/' + exportar_como[i] + '.csv'
    frame = crear_dataframe( exportar[i] )
    frame.to_csv( filename )

# %% --- 
