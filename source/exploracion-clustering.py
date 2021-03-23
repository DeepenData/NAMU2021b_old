#!/bin/python
"""
# TODO: "Una cita interesante"
"""
# %% --- IMPORTANDO EL DATASET

import pickle 

# DataFrames Pandas con los datos metabolicos
infile = open('./tmp/excel_dataset.pandas.pkl',     'rb'); df = pickle.load(infile); infile.close()
#infile = open('./tmp/excel_metabolitos.pandas.pkl', 'rb'); df_metabolitos = pickle.load(infile); infile.close()

df = df.set_index('Muestra') # Define la columna muestra como indice
df = df.dropna()             # Elimina filas con NaNs


# Esto es lo mismo que 'df_metabolitos' guardado como .pandas.pkl
metabolites_columns = ['Phe', 'Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA', 'C0',
    'C2', 'C3', 'C4', 'C4OH/C3DC', 'C5-OH/C4DC', 'C5', 'C5DC', 'C5:1', 'C6',
    'C6DC', 'C8', 'C8:1', 'C10', 'C10:1', 'C10:2', 'C12', 'C12:1', 'C14',
    'C14:1', 'C14:2', 'C14OH', 'C16', 'C16OH', 'C16:1', 'C16:1OH', 'C18',
    'C18OH', 'C18:1', 'C18:1OH', 'C18:2']

df_metabolitos = df[metabolites_columns]

# %% --- 
import seaborn as sns
from sklearn.preprocessing import StandardScaler, normalize, minmax_scale
# %% --- 

metabolites_columns2 = ['Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA']

for col in metabolites_columns2:
    print( df_metabolitos[col] )
    #set = df_metabolitos[col].values
    #set = minmax_scale( set )
    #sns.histplot(data= set , bins=30)


# %% --- PAIR PLOT --- Demora una cantidad de siglos... 

import matplotlib.pyplot as plt 

metabolitos_pairplot = sns.pairplot(df_metabolitos, diag_kind="kde")
metabolitos_pairplot.map_lower(sns.kdeplot, levels=4, color=".2")

metabolitos_pairplot.savefig("./doc/img/metabolitos_pairplot.png")

plt.clf() # Clean parirplot figure from sns 

# %% --- 
