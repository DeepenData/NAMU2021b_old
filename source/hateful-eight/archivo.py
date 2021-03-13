#!/bin/python

# %% --- Importar resultados
import pickle 

# Esto asume que está en 'human-metnet/'
infile = open('./tmp/perturbed_centralities','rb'); perturbed_centralities = pickle.load(infile); infile.close()
infile = open('./tmp/baseline_centralities' ,'rb'); baseline_centralities  = pickle.load(infile); infile.close()
infile = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(infile); infile.close()
infile = open('./tmp/subsystems_dict','rb'); subsystems = pickle.load(infile); infile.close()

# %% --- Importando librerias utiles
import pandas   as pd
import numpy    as np
import networkx as nx

# %% --- Algo informativo sobre el numero de NaNs

# count_nan = lambda promedios : np.count_nonzero(np.isnan( promedios ))
# def nan_info(promedios):
#     """Devuelve un print con los datos totales, numero de NaNs, y porcentaje"""
#     print ( pd.Series(
#         [
#             promedios.size, 
#             count_nan( promedios ), 
#             round(100*count_nan( promedios )/promedios.size, 0)
#         ], 
#         index=['Datos totales', 'NaNs', 'Porcentaje NaN'], 
#         dtype='uint'
#         ) )
# 
# print( 'baseline_centralities' )
# nan_info( baseline_centralities )
# print( 'perturbed_centralities' )
# nan_info( perturbed_centralities )




# %% --- Colapsando la tercera dimensión... (centralidades)
# Esto calcula 4 promedios de las centralidades para los nodos

from scipy import stats

# Un helper porque estas funciones no pueden lidiar con NaNs porque son superiores a los NaN o algo así
noNaN_perturbed_centralities = np.nan_to_num( perturbed_centralities , copy=True, nan=0.0, posinf=None, neginf=None)

# Es el unico que no require lambdas raros :D
aritmetic_perturbed = np.nanmean( noNaN_perturbed_centralities, axis= 1)
aritmetic_baseline  = np.nanmean(        baseline_centralities, axis= 1)

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
geometric_mean = lambda array :  stats.gmean( array[array != 0] ) * len(array[array != 0])/len(array)
geometric_perturbed = np.apply_along_axis( geometric_mean, 1, noNaN_perturbed_centralities)
geometric_baseline  = np.apply_along_axis( geometric_mean, 1,        baseline_centralities)

# Define una función helper para cuadraticos (Root Square Mean). No en scipy base
quadratic_mean = lambda array : np.sqrt( np.nanmean( array*array ) )
quadratic_perturbed = np.apply_along_axis( quadratic_mean, 1, noNaN_perturbed_centralities)
quadratic_baseline  = np.apply_along_axis( quadratic_mean, 1,        baseline_centralities)

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
# solo funciona para valores superiores a cero, por lo que decidimos ignorar los negativos
harmonic_mean = lambda array : stats.hmean( array[array > 0] ) * len(array[array > 0])/len(array)
harmonic_perturbed  = np.apply_along_axis( harmonic_mean,  1, noNaN_perturbed_centralities)
harmonic_baseline   = np.apply_along_axis( harmonic_mean,  1,        baseline_centralities)

# TODO: probar un JIT para la velocidad de estas lacras

# %% --- Divisiones! Calculando como un nodo afecta conectividad de otro
# 

log2_aritmetic = np.log2( aritmetic_baseline / aritmetic_perturbed )
log2_geometric = np.log2( geometric_baseline / geometric_perturbed )
log2_quadratic = np.log2( quadratic_baseline / quadratic_perturbed )
log2_harmonic  = np.log2( harmonic_baseline  / harmonic_perturbed  )

ratio_aritmetic = ( aritmetic_baseline - aritmetic_perturbed ) / aritmetic_baseline
ratio_geometric = ( geometric_baseline - geometric_perturbed ) / geometric_baseline
ratio_quadratic = ( quadratic_baseline - quadratic_perturbed ) / quadratic_baseline
ratio_harmonic  = ( harmonic_baseline  - harmonic_perturbed  ) / harmonic_baseline 

# %% --- 
# 1056 x 12
# 
# delta_centralidad_c = ratio(baseline/pertubed)
# 
# nodos_removido x delta_centralidad_c # (1000,12) ~ (1,12)/(1000,12)
# nodos_removido x delta_centralidad_c
# 
# nodos_removido x delta_centralidad_c
# nodos_removido x delta_centralidad_c
# 
# 
# nodos_removido x delta_centralidad_c_del_sudsystema_agregada
