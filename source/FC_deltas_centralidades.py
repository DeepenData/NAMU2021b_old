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

# %% --- Importando librerias utiles
import pandas   as pd
import numpy    as np
import networkx as nx

# %% --- SANITY CHECK Reordena el diccionario de centralidades perturbadas

INDEX_NODES = baseline_centralities['baseline'].index # Orden como la salida de los nodos
perturbed_centralities = { index : perturbed_centralities[index] for index in INDEX_NODES }

# %% --- FULL TENSOR

def full_tensor(diccionarios):
    """Genera un tensor desde el diccionario"""
    subsys = [ diccionarios[ KEY ].to_numpy() for KEY in list(diccionarios.keys()) ]
    subsys = np.asarray( subsys )
    return subsys

baseline_centralities_tensor  = full_tensor(baseline_centralities)
perturbed_centralities_tensor = full_tensor(perturbed_centralities)

import pickle
import os; os.makedirs("./tmp", exist_ok=True) # crea .tmp si no existe
outfile = open('./tmp/baseline_tensor.pkl', 'wb'); pickle.dump( baseline_centralities_tensor ,outfile); outfile.close()
outfile = open('./tmp/centralidades_perturbadas_tensor.pkl', 'wb'); pickle.dump( perturbed_centralities_tensor ,outfile); outfile.close()

# %% --- Cosas de los subsistemas

def subsystem_tensor( subsystem_nodes, diccionarios):
    """Genera un tensor con un subset de nodos basados en el subsistema pedido"""
    subsys = [ diccionarios[ KEY ].loc[ subsystem_nodes ,].to_numpy() for KEY in list(diccionarios.keys()) ]
    subsys = np.asarray( subsys )
    return subsys

# %% --- DEFINE COSAS DE PROMEDIOS

from numba import jit

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
@jit(nopython=True)
def geometric_mean(array):
    array = array[~np.isnan(array)] # Filtra los NaN
    fcorr = len(array[array > 0])/len(array) # Factor de corrección
    array = array[array > 0] # Selecciona solo mayores a cero
    try: gmean = ( np.prod( array )**(1 / len(array)) ) * fcorr # Geometrica por factor de corrección
    except: return np.nan
    else : return gmean

# Define una función helper para cuadraticos (Root Square Mean). No en scipy base
@jit(nopython=True)
def quadratic_mean(array):
    mean =  np.sqrt( np.nanmean( array*array ) )
    return mean

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
@jit(nopython=True)
def harmonic_mean(array):
    array = array[~np.isnan(array)] # Filtra los NaN
    fcorr = len(array[array > 0])/len(array) # Factor de corrección
    array = array[array > 0] # Selecciona solo mayores a cero
    try : hmean = len(array)/np.sum(1 / array ) * fcorr # Geometrica por factor de corrección
    except: return np.nan
    else : return hmean

# %% --- Convierte a Dataframes

index_centralities = [
    'degree_centrality', 'harmonic_centrality', 'eigenvector_centrality', 'betweenness_centrality', 
    'closeness_centrality', 'load_centrality', 'information_centrality', 'communicability_betweenness_centrality', 
    'katz_centrality', 'pagerank'
]

def crear_dataframe( matriz ):
    df = pd.DataFrame(
        matriz,
        index = INDEX_NODES,
        columns = index_centralities
    )
    return df

# %% --- INICIA EL LOOP ITERATIVO DE SUBSISTEMAS

# Importando el diccionario de {subsistema : [nodos] }
# infile = open('./tmp/subsystems_dict','rb'); subsystems = pickle.load(infile); infile.close()

subsystems = {
    'Glycolysis_astrocyte' : ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m'] ,
    'Glycolysis_neuron' : ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron', 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron'] ,
    'ETC_neuron'    : ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron'] ,
    'ETC_astrocyte' :  ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']
}

subsystems

# %% ---  CALCULOS TEDIOSOS

for sub in subsystems: 
    subsystem = subsystems[sub]

    perturbed_subsystem = subsystem_tensor( subsystem, perturbed_centralities )
    baseline_subsystem  = subsystem_tensor( subsystem, baseline_centralities  )

    # --- Colapsando la segunda dimensión... (nodos totales)
    #     Esto calcula 4 promedios de las centralidades para los nodos

    # Un helper porque estas funciones no pueden lidiar con NaNs porque son superiores a los NaN o algo así
    noNaN_perturbed_subsystem = np.nan_to_num( perturbed_subsystem , copy=True, nan=0.0, posinf=None, neginf=None)

    # Es el unico que no require lambdas raros :D
    aritmetic_perturbed = np.nanmean( noNaN_perturbed_subsystem, axis= 1)
    aritmetic_baseline  = np.nanmean(        baseline_subsystem, axis= 1)

    geometric_perturbed = np.apply_along_axis( geometric_mean, 1, noNaN_perturbed_subsystem)
    geometric_baseline  = np.apply_along_axis( geometric_mean, 1,        baseline_subsystem)

    quadratic_perturbed = np.apply_along_axis( quadratic_mean, 1, noNaN_perturbed_subsystem)
    quadratic_baseline  = np.apply_along_axis( quadratic_mean, 1,        baseline_subsystem)

    harmonic_perturbed  = np.apply_along_axis( harmonic_mean,  1, noNaN_perturbed_subsystem)
    harmonic_baseline   = np.apply_along_axis( harmonic_mean,  1,        baseline_subsystem)

    # --- Divisiones! Calculando como un nodo afecta conectividad de otro
    #     Calcula las metricas para definir los efectos en los nodos

    fold_change_aritmetic = np.log2( aritmetic_baseline / aritmetic_perturbed )
    fold_change_geometric = np.log2( geometric_baseline / geometric_perturbed )
    fold_change_quadratic = np.log2( quadratic_baseline / quadratic_perturbed )
    fold_change_harmonic  = np.log2( harmonic_baseline  / harmonic_perturbed  )

    delta_aritmetic = ( aritmetic_baseline - aritmetic_perturbed ) / aritmetic_baseline
    delta_geometric = ( geometric_baseline - geometric_perturbed ) / geometric_baseline
    delta_quadratic = ( quadratic_baseline - quadratic_perturbed ) / quadratic_baseline
    delta_harmonic  = ( harmonic_baseline  - harmonic_perturbed  ) / harmonic_baseline 

    # --- Exporta a R (via Excel)

    # todo: eliminar este código estúpido que no deberia existir
    # me da asco verlo y solo existe por la capacidad multi-cursor de VS Code
    #   - Manu, 2021-03-13 23:04

    # Es feo pero funciona, y realmente cuando esta en un loop se ve casi decente
    #   - Manu, 2021-03-20 21:07

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

    #for i in range(len(exportar)):
    #    filename = './results/' + exportar_como[i] + '.csv'
    #    frame = crear_dataframe( exportar[i] )
    #    frame.to_csv( filename )

    # Define un creador del Excel
    SALIDA = './results/' + str(sub) + '_centralies_2021-03-31.xlsx' # Nombre del archivo, subsistema_centralities.xsls
    salida_Excel = pd.ExcelWriter( SALIDA , engine="xlsxwriter")   # Crea un onjeto para guardar el Excel

    exportar_dframes = [ crear_dataframe( cent ) for cent in exportar ] # Convierte todo a DataFrames (de nuevo)

    # Loop dentro de la lista de hojas a exportar
    for i, cent in enumerate (exportar_dframes):
        cent.to_excel(salida_Excel, sheet_name=exportar_como[i] , index=True)
    
    # Guarda el Excel
    salida_Excel.save()

# %% --- 
