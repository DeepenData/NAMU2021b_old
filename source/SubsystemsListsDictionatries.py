#!/bin/python
"""
"Of the Horcrux, wickedest of magical inventions, we shall not speak nor give direction"
    — Magick Moste Evile
"""
# %% --- IMPORTA EL MODELO Y COSAS
import cobra 

# TODO: que use una variable ambiental? Aunque igual es manualmente curado
INPUT_MODEL = './data/GEM_Recon2_thermocurated_redHUMAN.json'

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

# Devuelve un filtrado por subsistema
df_subsystem = lambda subsys : tabla_primaria[ tabla_primaria['subsystem'] == subsys ] 

# Lista de subsistemas unicos en el modelo
unique_subsystems = pd.unique( tabla_primaria['subsystem'] )

# List of dataframes by subsystems
# dfs_by_subsys = [ df_subsystem( sub ) for sub in unique_subsystems ]

# %% --- DICCIONARIO DE SUBSISTEMAS MANUALMENTE CURADOS
# Este es para pasarlo a otros analisis, como un objeto Python

subsystems_dict = { sub : list(df_subsystem( sub )['ids']) for sub in unique_subsystems }

import pickle 
import os; os.makedirs("./tmp", exist_ok=True) # crea .tmp si no existe
outfile = open('./tmp/subsystems_dict.pkl', 'wb'); pickle.dump( subsystems_dict , outfile ); outfile.close()


