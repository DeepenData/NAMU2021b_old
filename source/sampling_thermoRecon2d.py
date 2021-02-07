
# %% --- 

import cobra
import pandas as pd
import numpy as np
import warnings 
warnings.filterwarnings("ignore")
#GEM_Recon2_thermocurated_redHUMAN = cobra.io.load_matlab_model('GEM_Recon2_thermocurated_redHUMAN.mat')
#cobra.io.save_json_model(GEM_Recon2_thermocurated_redHUMAN, "GEM_Recon2_thermocurated_redHUMAN.json")
normal = cobra.io.load_json_model("GEM_Recon2_thermocurated_redHUMAN.json")
print(
normal.reactions.get_by_id("PHETHPTOX2").bounds,
normal.reactions.get_by_id("r0399").bounds)
frozen         = normal.metabolites.get_by_id("tyr_L_e").reactions.union(normal.metabolites.get_by_id("tyr_L_c").reactions).union(normal.metabolites.get_by_id("tyr_L_m").reactions)
tyr_L_ids    = [list(frozen)[i].id for i in range(0,len(frozen))]
tyr_L_names  = [list(frozen)[i].name for i in range(0,len(frozen))]
tyr_L_rxns   = [list(frozen)[i].reaction for i in range(0,len(frozen))]

list_of_tuples   = list(zip(tyr_L_ids,tyr_L_names, tyr_L_rxns))  
columns          = ["tyr_L_ids","tyr_L_names","tyr_L_rxns"]
tyr_L          =  pd.DataFrame(list_of_tuples, columns = columns)  


frozen         = normal.metabolites.get_by_id("phe_L_e").reactions.union(normal.metabolites.get_by_id("phe_L_c").reactions).union(normal.metabolites.get_by_id("phe_L_m").reactions)
phe_L_ids    = [list(frozen)[i].id for i in range(0,len(frozen))]
phe_L_names  = [list(frozen)[i].name for i in range(0,len(frozen))]
phe_L_rxns   = [list(frozen)[i].reaction for i in range(0,len(frozen))]

list_of_tuples   = list(zip(phe_L_ids,phe_L_names, phe_L_rxns))  
columns          = ["phe_L_ids","phe_L_names","phe_L_rxns"]
phe_L          =  pd.DataFrame(list_of_tuples, columns = columns)

# %% ---

from cobra.sampling import OptGPSampler, ACHRSampler, sample
s = sample(normal, n = 100, method='optgp',   seed=23, processes=16, thinning=100)

# %% ---

s[phe_L.phe_L_ids].T

# %% ---

s[phe_L.phe_L_ids].T
