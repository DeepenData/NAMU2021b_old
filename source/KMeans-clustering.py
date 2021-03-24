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

# %% --- ESTANDARIZACION Y ESCALADO 

from sklearn.preprocessing import StandardScaler
x = df_metabolitos.values
x = StandardScaler().fit_transform(x) # normalizing the features


# %% --- PCA

from sklearn.decomposition import PCA
pca_metabolitos = PCA(n_components=2)
principalComponents = pca_metabolitos.fit_transform(x)

# %% --- EMPAQUE EN UN DATAFRAME
# tengo codigo R que hace esto, pero es una oportunidad de aprender Python

PCA_dataframe = pd.DataFrame(data = principalComponents, columns = ['PC_1', 'PC_2' ]) #, 'PC_3', 'PC_4', 'PC_5'])

PCA_dataframe

# %% --- 

print('Explained variation per principal component: {}'.format(pca_metabolitos.explained_variance_ratio_))

# %% --- 
