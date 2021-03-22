#!/bin/python
"""
# TODO: "Una cita interesante"
"""
# %% --- IMPORTANDO LIBRERIAS Y DATASETS
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

# %% --- OPCIONES DE CONFIGURACION  

import pickle 

# DataFrames Pandas con los datos metabolicos
infile = open('./tmp/excel_dataset.pandas.pkl',     'rb'); df_metabolitos = pickle.load(infile); infile.close()
#infile = open('./tmp/excel_metabolitos.pandas.pkl', 'rb'); df_metabolitos = pickle.load(infile); infile.close()
# %% --- 
