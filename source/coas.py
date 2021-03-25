# To add a new cell, type '# %% --- '
# To add a new markdown cell, type '# %% ---  [markdown]'
# %% ---  Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-toolsai.jupyter added
#import os
#try:
#	os.chdir(os.path.join(os.getcwd(), '../../../../tmp/4c3f94c9-da94-496f-ad37-e20ec2a79e7e'))
#	print(os.getcwd())
#except:
#	pass
# %% --- 
import pickle

# %% --- 
infile = open('./data/Recon2_metabolite_projection.gpickle', 'rb'); g = pickle.load(infile); infile.close()

# %% --- 
type(g)

# %% --- 
import networkx as nx

# %% --- 
import pandas as pd

# %% --- 
df = pd.read_csv('results/dataplot.csv', index_col=0); df
# %% --- 
df_A = df[ df['k-labels'] == 'A' ]
df_B = df[ df['k-labels'] == 'B' ]
df_C = df[ df['k-labels'] == 'C' ]
df_D = df[ df['k-labels'] == 'D' ]
df_E = df[ df['k-labels'] == 'E' ]

# %% --- 

g.nodes

# %% --- 


dataset_centralidad_base = pd.read_csv('data/recon2_metabolite_centralities_metabolome_weights.csv', index_col= 0 )

import pickle 

# DataFrames Pandas con los datos metabolicos
infile = open('./tmp/excel_dataset.pandas.pkl', 'rb'); df = pickle.load(infile); infile.close()

df = df.set_index('Muestra') # Define la columna muestra como indice

print('Entradas iniciales:', df.shape )

# Esto es lo mismo que 'df_metabolitos' guardado como .pandas.pkl
metabolites_columns = ['Phe', 'Met', 'Val', 'Leu/Ile', 'tir',
    'Pro', 'Arg', 'Gly', 'Ala', 'Asp', 'Glt', 'Cit', 'Orn', 'SA', 'C0',
    'C2', 'C3', 'C4', 'C4OH/C3DC', 'C5-OH/C4DC', 'C5', 'C5DC', 'C5:1', 'C6',
    'C6DC', 'C8', 'C8:1', 'C10', 'C10:1', 'C10:2', 'C12', 'C12:1', 'C14',
    'C14:1', 'C14:2', 'C14OH', 'C16', 'C16OH', 'C16:1', 'C16:1OH', 'C18',
    'C18OH', 'C18:1', 'C18:1OH', 'C18:2']

dataset_metabolitos = df[metabolites_columns]

dataset_metabolitos = dataset_metabolitos.dropna()


