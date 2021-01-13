#-------------------------------------------------------------------------------
# 
#  Modelo por Optimización termodinamica
# 
#-------------------------------------------------------------------------------
#
#!conda update conda
#!pip install update pip
#!pip install update cobra
#!pip install update pip
#!pip install update  joblib 
#!pip install update  multiprocess
#!pip install update networkx
#!pip install update  pandas 
#!pip install update numpy 
#!pip install update scipy    
# 
#-------------------------------------------------------------------------------

# %% --- Importación de librerias
import cobra  as cb
import pandas as pd

import warnings 
warnings.filterwarnings("ignore")

# %% --- Carga los modelos y datos preexistentes
normal = cobra.io.load_json_model("GEM_Recon2_thermocurated_redHUMAN.json")
normal

# %% - 
print(
normal.reactions.get_by_id("PHETHPTOX2").bounds,
normal.reactions.get_by_id("r0399").bounds)

# %% - Lista de metabolitos

frozen       = normal.metabolites.get_by_id("tyr_L_e").reactions.union(normal.metabolites.get_by_id("tyr_L_c").reactions).union(normal.metabolites.get_by_id("tyr_L_m").reactions)
tyr_L_ids    = [list(frozen)[i].id for i in range(0,len(frozen))]
tyr_L_names  = [list(frozen)[i].name for i in range(0,len(frozen))]
tyr_L_rxns   = [list(frozen)[i].reaction for i in range(0,len(frozen))]
list_of_tuples   = list(zip(tyr_L_ids, tyr_L_names, tyr_L_rxns))  
columns          = ["tyr_L_e_ids","tyr_L_e_names","tyr_L_e_rxns"]
tyr_L          =  pd.DataFrame(list_of_tuples, columns = columns)  

frozen         = normal.metabolites.get_by_id("phe_L_e").reactions.union(normal.metabolites.get_by_id("phe_L_c").reactions).union(normal.metabolites.get_by_id("phe_L_m").reactions)
phe_L_ids    = [list(frozen)[i].id for i in range(0,len(frozen))]
phe_L_names  = [list(frozen)[i].name for i in range(0,len(frozen))]
phe_L_rxns   = [list(frozen)[i].reaction for i in range(0,len(frozen))]

list_of_tuples   = list(zip(phe_L_ids,phe_L_names, phe_L_rxns))  
columns          = ["phe_L_e_ids","phe_L_e_names","phe_L_e_rxns"]
phe_L          =  pd.DataFrame(list_of_tuples, columns = columns)

# %% - FBA

solution_normal = normal.optimize()
phe_L_fluxes   = solution_normal.fluxes[phe_L_ids].to_frame() #
phe_L_fluxes.query('abs(fluxes) > 1e-10', inplace=True)
phe_L_fluxes.rename_axis('phe_L_e_ids').reset_index()
phe_L_fluxes = phe_L_fluxes.rename_axis('phe_L_e_ids').reset_index()
pd.merge(phe_L_fluxes, phe_L)

# %% - 

tyr_L_fluxes   = sol_normal.fluxes[tyr_L_ids].to_frame()  #.query('abs(fluxes) > 1e-7', inplace=True)
tyr_L_fluxes.query('abs(fluxes) > 1e-10', inplace=True)
tyr_L_fluxes = tyr_L_fluxes.rename_axis('tyr_L_e_ids').reset_index()
pd.merge(tyr_L_fluxes, tyr_L)

# %% - Fenotipo - Finding glucose and oxugen uptakes

#normal.
#[list(frozen)[i].reaction for i in range(0,len(frozen))]

rxn_names      = [normal.reactions[i].name for i in range(0, len(normal.reactions))]
rxn_ids        = [normal.reactions[i].id for i in range(0, len(normal.reactions))]
rxn_reaction   = [normal.reactions[i].reaction for i in range(0, len(normal.reactions))]

model_df_rxns       = pd.DataFrame(list(zip(rxn_ids,rxn_names, rxn_reaction)), columns = ['rxn_ids','rxn_names',"rxn_reaction"])


met_names      = [normal.metabolites[i].name for i in range(0, len(normal.metabolites))]
met_ids        = [normal.metabolites[i].id for i in range(0, len(normal.metabolites))]
met_formula   = [normal.metabolites[i].formula for i in range(0, len(normal.metabolites))]

model_df_mets       = pd.DataFrame(list(zip(met_ids,met_names, met_formula)), columns = ['met_ids','met_names',"met_formula"])

model_df_mets[model_df_mets.met_ids.str.contains('glc_D_')]
#model_df_mets[model_df_mets.met_ids.str.contains('^o2_')]

# %% - Función metabolitos a reacciones
# ToDo: cambiar esto por import module.py

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

# %% - Observa glucosa y oxigeno

glucose_oxygen = metabolites2reactions(normal, ["o2_e", "glc_D_e"])
glucose_oxygen[glucose_oxygen.names.str.contains('exchange')]

# %% - Uptakes y bounds para glucosa y oxigeno

normal.reactions.get_by_id("EX_o2_e").reaction
#normal.reactions.get_by_id("EX_glc_e").reaction

# %% 

upkates = ["EX_o2_e","EX_glc_e"]
solution_normal[upkates]

# %%

pd.set_option('display.max_rows', None)

solution_normal_df          = pd.merge(solution_normal.to_frame().rename_axis('rxn_ids').reset_index(),model_df_rxns)
solution_normal_df_filtered = solution_normal_df.query('abs(fluxes) > 1e-10') 
df_EXs                      = solution_normal_df_filtered[solution_normal_df_filtered.rxn_ids.str.contains('EX')]
df_EXs.query('(fluxes) < 0') 

# %% --- Make the model glucose-oxygent dependent

solution_normal2.shadow_prices.to_frame()

# %% --- Parsimonious FBA

pfba_solution = cobra.flux_analysis.pfba(normal)
sol_normal = pfba_solution

phe_L_fluxes   = sol_normal.fluxes[phe_L_ids].to_frame() #.query('abs(fluxes) > 1e-7', inplace=True)
phe_L_fluxes.query('abs(fluxes) > 1e-7', inplace=True)
phe_L_fluxes = phe_L_fluxes.rename_axis('phe_L_e_ids').reset_index()
pd.merge(phe_L_fluxes, phe_L)

tyr_L_fluxes   = sol_normal.fluxes[tyr_L_ids].to_frame()  #.query('abs(fluxes) > 1e-7', inplace=True)
tyr_L_fluxes.query('abs(fluxes) > 1e-7', inplace=True)
tyr_L_fluxes = tyr_L_fluxes.rename_axis('tyr_L_e_ids').reset_index()
pd.merge(tyr_L_fluxes, tyr_L)

# %% --- geometric FBA

geometric_fba_sol = cobra.flux_analysis.geometric_fba(normal, epsilon=1e-01, max_tries=50)
sol_normal = geometric_fba_sol
phe_L_fluxes   = sol_normal.fluxes[phe_L_e_ids].to_frame() #.query('abs(fluxes) > 1e-7', inplace=True)
phe_L_fluxes.query('abs(fluxes) > 1e-7', inplace=True)
phe_L_fluxes = phe_L_fluxes.rename_axis('phe_L_e_ids').reset_index()
pd.merge(phe_L_fluxes, phe_L)