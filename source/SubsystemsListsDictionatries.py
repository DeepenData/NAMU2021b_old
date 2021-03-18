#!/bin/python
"""
"Of the Horcrux, wickedest of magical inventions, we shall not speak nor give direction"
    — Magick Moste Evile
"""
# %% --- IMPORTA EL MODELO Y COSAS
import cobra 

# TODO: que use una variable ambiental? Aunque igual es manualmente curado
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

glicolisis_ids_astros = ['ENO', 'FBA', 'GAPD', 'HEX1', 'PFK', 'PGI', 'PGK', 'PGM', 'PYK', 'TPI'] # Curado manual
# PEPCK no existe en el grafo porque no es fully-conected
glicolisis_ids_neuron = [ ids + '_Neuron' for ids in glicolisis_ids_astros ] # Stimulated tiene "IDs duplicados" por ID_Neuron

glicolisis_astros = tabla_primaria[tabla_primaria['ids'].isin( glicolisis_ids_astros )] # Dataframe con reacciones de glicolisis curadas
glicolisis_astros.to_csv( './tmp/glicolisis_astros.csv' , index=False ) # Salida a CSV para legilibilidad y eso

glicolisis_neuron = tabla_primaria[tabla_primaria['ids'].isin( glicolisis_ids_neuron )] # Dataframe con reacciones de glicolisis curadas
glicolisis_neuron.to_csv( './tmp/glicolisis_neuron.csv' , index=False ) # Salida a CSV para legilibilidad y eso

# SANITIY CHECK DE IDS CORRESPONDIENDO AL MODELO
glicolisis_ids_astros  = list( glicolisis_astros['ids'] ) # Es necesario que queden absolutamente iguales
glicolisis_ids_neuron = list( glicolisis_neuron['ids'] ) # same 

# %% --- REACCIONES EXPLORATORIAS PARA FOSFORILACIÓN OXIDATIVA

oxphox = tabla_primaria[tabla_primaria['subsystem'].str.contains('Oxidative Phosphorylation', case=False)==True] # Dataframe exploratorio

#oxphox_ids = ['ATPS4m','CYOOm2','CYOR-u10m','P45011A1m','SUCCt2m'] # Curación manual no necesaria
oxphox.to_csv( './tmp/oxphox.csv' , index=False ) # Salida a CSV para legilibilidad y eso

oxphox_ids = list( oxphox['ids'] ) # Es necesario que queden absolutamente iguales

oxphox_ids_astros = [ id for id in oxphox_ids if not id.endswith('_Neuron')]
oxphox_ids_neuron = [ id for id in oxphox_ids if id.endswith('_Neuron')]

# %% --- REACCIONES PARA TRANSPORTE DE ELECTRONES
# Manualmente curados, posiblemente requiere algo más automatico

# TODO: revisar esto
etc_astros = ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']
etc_neuron = ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron']

# %% --- DICCIONARIO DE SUBSISTEMAS MANUALMENTE CURADOS
# Este es para pasarlo a otros analisis, como un objeto Python

subsystems_dict = {
    'glicolisis_astros' : glicolisis_ids_astros, 
    'etc_astros'     : etc_astros,
    'glicolisis_neuron' : glicolisis_ids_neuron, 
    'etc_neuron'     : etc_neuron
}

import pickle 
outfile = open('./tmp/subsystems_dict', 'wb'); pickle.dump( subsystems_dict , outfile ); outfile.close()
