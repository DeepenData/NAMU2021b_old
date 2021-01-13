""" # Cosas para instalar con conda. Son magic de Jupyter, posiblemente fallen
!conda update conda
!pip install update pip
!pip install update cobra
!pip install update pip
!pip install update  joblib 
!pip install update  multiprocess
!pip install update networkx
!pip install update  pandas 
!pip install update numpy 
!pip install update scipy    
"""
# %% --- 0.1 Importación de librerias
import cobra  
import numpy  as np
import pandas as pd

import warnings 
warnings.filterwarnings("ignore")

# %% --- 0.2 Carga el modelo
normal = cobra.io.load_json_model("GEM_Recon2_thermocurated_redHUMAN.json")

# %% --- 1. Optimizar el modelo para biomasa
solution_normal = normal.optimize() # Optimiza la cosa de FBA
normal.summary()         # El modelo se optimiza para biomasa

# %% - Función metabolitos a reacciones
def metabolites2reactions(modelo, metabolitos): # Toma metabolitos
  """Crea un agregado de las reacciones segun un listado de metabolitos. 

  Parameters
  ----------
  modelo: model
      El modelo inicial del cual se toman las reacciones
  metabolitos : list
      IDs de un metabolito. La función acepta multiples para crear un joint. 

  Returns
  -------
  df
      crea un dataframe pandas agregada de ids, nombre de reacciones, y reacciones en si
  """
  frozen = modelo.metabolites.get_by_id(metabolitos[0]).reactions # Crea el inicio
  for met in metabolitos: # Añade metabolistos al listado de reacciones
    frozen = frozen.union(modelo.metabolites.get_by_id(met).reactions)
  
  ids   = [list(frozen)[i].id       for i in range(0,len(frozen))] # List Comprehesion
  names = [list(frozen)[i].name     for i in range(0,len(frozen))]
  rxns  = [list(frozen)[i].reaction for i in range(0,len(frozen))]

  list_of_tuples = list(zip(ids, names, rxns))     # DataFrame en Python; crea triplets
  columns        = ["ids","names","rxns"]   # Titulos de columnas de la tripleta
  df             =  pd.DataFrame(list_of_tuples, columns = columns) # Titulos asignados a la lista de tuplas
  
  return df

# %% - Lista de metabolitos
# Nos interesa tirosina y fenilalanina porque son targets de mutaciones en enfermedades metabolicas
# así que extraemos reacciones asociadas a Tyr y Phe

glucose_oxygen = metabolites2reactions(normal, ["o2_e", "glc_D_e"])
glucose_oxygen[glucose_oxygen.names.str.contains('exchange')]

# %% --- Extacción de IDs y busqueda de cosas relacionadas a glucosa

rxn_reactions = [normal.reactions[i].reaction for i in range(0, len(normal.reactions)) ]
rxn_names =     [normal.reactions[i].name     for i in range(0, len(normal.reactions)) ]
rxn_ids =       [normal.reactions[i].id       for i in range(0, len(normal.reactions)) ]

model_df = pd.DataFrame( list(zip(rxn_ids, rxn_names, rxn_reactions)), columns = ["rxn_ids", "rxn_names", "rxn_reactions"] )

model_df[model_df.rxn_names.str.contains('Glucose', case = False)]

# %% --- Formas de glucosa en el modelo
rxn_names      = [normal.reactions[i].name for i in range(0, len(normal.reactions))]
rxn_ids        = [normal.reactions[i].id for i in range(0, len(normal.reactions))]
rxn_reaction   = [normal.reactions[i].reaction for i in range(0, len(normal.reactions))]

model_df_rxns       = pd.DataFrame(list(zip(rxn_ids,rxn_names, rxn_reaction)), columns = ['rxn_ids','rxn_names',"rxn_reaction"])


met_names      = [normal.metabolites[i].name for i in range(0, len(normal.metabolites))]
met_ids        = [normal.metabolites[i].id for i in range(0, len(normal.metabolites))]
met_formula    = [normal.metabolites[i].formula for i in range(0, len(normal.metabolites))]

model_df_mets       = pd.DataFrame(list(zip(met_ids,met_names, met_formula)), columns = ['met_ids','met_names',"met_formula"])

model_df_mets[model_df_mets.met_ids.str.contains('glc_D_')]
#model_df_mets[model_df_mets.met_ids.str.contains('^o2_')]

# %% --- Identificando reacciones de intercambio de glucosa y oxigeno
glucose_oxygen = metabolites2reactions(normal, ["o2_e", "glc_D_e"])
glucose_oxygen[glucose_oxygen.names.str.contains('exchange')]

print(
"Intercambio de Oxigeno:", normal.reactions.get_by_id("EX_o2_e").bounds  ,
"Intercambio de Glucosa:", normal.reactions.get_by_id("EX_glc_e").bounds )

upkates = ["EX_o2_e","EX_glc_e"]
solution_normal[upkates]
# %% --- Elimina flujos

normal2 = normal # Crea una copia del modelo

kill = ["EX_asn_L_e", "EX_dag_hs_e", "EX_dgsn_e", "EX_din_e","EX_gln_L_e",
         "EX_ins_e","EX_ser_L_e","EX_thymd_e","EX_tyr_L_e","EX_citr_L_e",
         "EX_alaala_e","EX_glypro_e","EX_leuleu_e"]

for i in kill: 
  normal.reactions.get_by_id(i).bounds = (0,0) # Desactiva metabolismo de glucosa

normal.reactions.get_by_id("EX_o2_e" ).bounds = (-100,100)
normal.reactions.get_by_id("EX_glc_e").bounds = (-100,100)

# %% --- Rutas metabolicas alteradas
solution_normal2 = normal2.optimize()

shadows = solution_normal2.shadow_prices
shadows = shadows.to_frame()

shadows = shadows.query('abs(shadow_prices) > 5e-17') # Elimina cambios menores
shadows = shadows.query('shadow_prices > -3') # ELIMINAR uno que se va de la regla
shadows = shadows.sort_values('shadow_prices', ascending = False)

shadows = shadows.rename_axis('metabolites').reset_index() 
shadows[shadows.metabolites.str.contains('_e$', case = False)] # Metabolitos _extracelulares

# %% --- Plot de rutas metabolicas alteradas
from matplotlib import pyplot as plt

y = shadows.shadow_prices # Sensibilidad
x = shadows.metabolites   # Metabolito

plt.plot(x, y, '-')
plt.fill_between(x, y, 200, where = (y > 195), facecolor='g', alpha=0.6)

plt.title("Sensibilidad a metabolitos")
plt.show()
# %% --- Cambio numerico del uptake
upkates = ["EX_o2_e","EX_glc_e"]
solution_normal2[upkates]
