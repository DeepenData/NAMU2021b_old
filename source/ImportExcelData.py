#!/bin/python

import sys
import pandas as pd

#EXCEL_FILE = sys.argv[1] # Ubicación del archivo Excel.
EXCEL_FILE = './data/Compilado_Anual_2017 (referencia).xlsx' # Path de muestra

# %% --- Importa la tabla, remueve cosas innecesarias

df = pd.read_excel(EXCEL_FILE, 
    sheet_name=0 ,
    skiprows=3 ,
    skipfooter= 5,
    header=0 ,
    engine='openpyxl')

print('Importando excel. Filas de datos:', df.shape[0] )

# %% --- Rework dtypes, convierte a formatos adecuados

import numpy as np

# Sección que elimina indeseables y los reemplaza por NaN
df.replace('ND ', np.NaN, inplace=True)
df.replace('ND', np.NaN, inplace=True)
df.replace('RN', np.NaN, inplace=True)

# Este codigo disminuye en cerca de un 45% el almacenamiento necesario
# y hace el procesamiento más limpio de correr
df = df.astype({
    'Fecha de nacimiento' : 'datetime64[ns]', 'Fecha entrega informe' : 'datetime64[ns]',
    'Edad (días calculados)' : 'float32',
    'RN/No_RN' : 'category', 'Muestra' : 'string', 
    'Nombre' : 'string', 'Apellido' : 'string',
    'Sexo' : 'category', 'RUT' : 'string', 
    'Médico' : 'string', 'Centro derivador' : 'category', 
    'Procedencia' : 'category',
    'Fecha obtención de muestra' : 'datetime64[ns]', 'Hora obtención de muestra' : 'object',
    'Fecha recepción de laboratorio' : 'datetime64[ns]', 'Hora recepción en laboratorio' : 'object',
    'Tipo de muestra' : 'category', 
    'Método' : 'category',
    'Phe'   : 'float32', 'Met'   : 'float32', 'Val'     : 'float32', 'Leu/Ile'   : 'float32', 'tir'        : 'float32',
    'Pro'   : 'float32', 'Arg'   : 'float32', 'Gly'     : 'float32', 'Ala'       : 'float32', 'Asp'        : 'float32',
    'Glt'   : 'float32', 'Cit'   : 'float32', 'Orn'     : 'float32', 'SA'        : 'float32', 'C0'         : 'float32',
    'C2'    : 'float32', 'C3'    : 'float32', 'C4'      : 'float32', 'C4OH/C3DC' : 'float32', 'C5-OH/C4DC' : 'float32', 
    'C5'    : 'float32', 'C5DC'  : 'float32', 'C5:1'    : 'float32', 'C6'        : 'float32', 'C6DC'       : 'float32', 
    'C8'    : 'float32', 'C8:1'  : 'float32', 'C10'     : 'float32', 'C10:1'     : 'float32', 'C10:2'      : 'float32', 
    'C12'   : 'float32', 'C12:1' : 'float32', 'C14'     : 'float32', 'C14:1'     : 'float32', 'C14:2'      : 'float32', 
    'C14OH' : 'float32', 'C16'   : 'float32', 'C16OH'   : 'float32', 'C16:1'     : 'float32', 'C16:1OH'    : 'float32', 
    'C18'   : 'float32', 'C18OH' : 'float32', 'C18:1'   : 'float32', 'C18:1OH'   : 'float32', 'C18:2'      : 'float32'
})
print( df.info() )

# %% --- Sección que crea una subtabla solo con los metabolitos e IDS de muestra

metabolites_columns = ['Phe', 'Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA', 'C0',
    'C2', 'C3', 'C4', 'C4OH/C3DC', 'C5-OH/C4DC', 'C5', 'C5DC', 'C5:1', 'C6',
    'C6DC', 'C8', 'C8:1', 'C10', 'C10:1', 'C10:2', 'C12', 'C12:1', 'C14',
    'C14:1', 'C14:2', 'C14OH', 'C16', 'C16OH', 'C16:1', 'C16:1OH', 'C18',
    'C18OH', 'C18:1', 'C18:1OH', 'C18:2']

df_metabolites = df[['Muestra'] + metabolites_columns]
df_metabolites.set_index('Muestra', inplace=True)

# %% --- Exportando dataframes 

import pickle 

outfile = open('./tmp/excel_dataset', 'wb'); pickle.dump( df , outfile ); outfile.close()
outfile = open('./tmp/excel_metabolitos', 'wb'); pickle.dump( df_metabolites , outfile ); outfile.close()

print("Excel a pickles terminado")
