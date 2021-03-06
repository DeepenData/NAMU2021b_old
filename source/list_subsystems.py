"""
import random   # DEMOSTRACIÓN DE SUBSISTEMAS RANDOM
random.seed(23) # DE DISTINTO LAGRGO Y SIMILAR

# EJEMPLO DUMMY DE SUBSISTEMAS
subsystems = {
    'ns1' : [ index_nodes[n] for n in random.sample(range(0, 1055), 7)  ] ,
    'ns2' : [ index_nodes[n] for n in random.sample(range(0, 1055), 11) ]
}"""

import cobra 

INPUT_MODEL = './data/stimulated_2021.json'

import warnings; warnings.filterwarnings('ignore') # Ignora warnings
from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL) # Genera infinitos warning

# %% --- CREA UN DATAFRAME PARA SLECCIONAR REACCIONES POR SUBSISTEMA
import pandas as pd
tabla_primaria = pd.DataFrame({
    "subsystem" :   [rx.subsystem for rx in model.reactions],
    "ids" :         [rx.id           for rx in model.reactions],
    "name" :        [rx.name         for rx in model.reactions],
    "formula" :     [rx.reaction     for rx in model.reactions]
    })

# %% --- REACCIONES EXPLORATORIAS PARA GLICOLISIS

glicolisis = tabla_primaria[tabla_primaria["subsystem"].str.contains('Glycolysis', case=False)==True] # Dataframe exploratorio
# Esta algo contaminado con reacciones que no corresponden a glicolisis, por lo que hare una curación manual para seleccionar
# solo las relevantes a glicolisis y gluconeogenesis

glicolisis_ids = ['ENO', 'FBA', 'GAPD', 'HEX1', 'PFK', 'PGI', 'PGK', 'PGM', 'PYK', 'TPI'] # Curado manual
# PEPCK no existe en el grafo porque no es fully-conected
glicolisis_ids = glicolisis_ids + [ ids + '_Neuron' for ids in glicolisis_ids ] # Stimulated tiene "IDs duplicados" por ID_Neuron

glicolisis = tabla_primaria[tabla_primaria['ids'].isin( glicolisis_ids )] # Dataframe con reacciones de glicolisis curadas
glicolisis.to_csv( './tmp/glicolisis.csv' , index=False ) # Salida a CSV para legilibilidad y eso

glicolisis_ids = list( glicolisis['ids'] ) # Es necesario que queden absolutamente iguales

# %% --- REACCIONES EXPLORATORIAS PARA GLICOLISIS

oxphox = tabla_primaria[tabla_primaria['subsystem'].str.contains('Oxidative Phosphorylation', case=False)==True] # Dataframe exploratorio

#oxphox_ids = ['ATPS4m','CYOOm2','CYOR-u10m','P45011A1m','SUCCt2m'] # Curación manual no necesaria
oxphox.to_csv( './tmp/oxphox.csv' , index=False ) # Salida a CSV para legilibilidad y eso

oxphox_ids = list( oxphox['ids'] ) # Es necesario que queden absolutamente iguales

# %% --- DICCIONARIO DE SUBSISTEMAS MANUALMENTE CURADOS
# Este es para pasarlo a otros analisis, como un objeto Python

subsystems_dict = {
    'glicolisis' : glicolisis_ids, 
    'oxphox'     : oxphox_ids
}

import pickle 
outfile = open('./tmp/subsystems_dict', 'wb'); pickle.dump( subsystems_dict , outfile ); outfile.close()
