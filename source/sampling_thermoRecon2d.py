
# %% --- 

import cobra
import pandas as pd
import numpy as np
from cobra.sampling import OptGPSampler, ACHRSampler, sample
import warnings 
warnings.filterwarnings("ignore")
#GEM_Recon2_thermocurated_redHUMAN = cobra.io.load_matlab_model('GEM_Recon2_thermocurated_redHUMAN.mat')
#cobra.io.save_json_model(GEM_Recon2_thermocurated_redHUMAN, "GEM_Recon2_thermocurated_redHUMAN.json")
path          = '/home/alejandro/PostDoc/human-metnet/data/'
model    = cobra.io.load_json_model(path + 'GEM_Recon2_thermocurated_redHUMAN.json') 



#s = sample(model, n = 100, method='optgp',   seed=23, processes=16, thinning=100)



# %% 

#perturbed_model  = model.copy()
print(
model.reactions.get_by_id("PHETHPTOX2").bounds,
model.reactions.get_by_id("r0399").bounds)



def block_a_flux(model, rxn):
    hola = model.copy()
    hola.reactions.get_by_id(rxn).bounds = 0, 0
    return hola

perturbed_model = block_a_flux(model,"PHETHPTOX2" )


print(
perturbed_model.reactions.get_by_id("PHETHPTOX2").bounds,
perturbed_model.reactions.get_by_id("r0399").bounds)


perturbed_model = block_a_flux(model,"r0399" )


print(
perturbed_model.reactions.get_by_id("PHETHPTOX2").bounds,
perturbed_model.reactions.get_by_id("r0399").bounds)


# %%

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



# %% ---

s[phe_L.phe_L_ids].T

# %% ---

s[phe_L.phe_L_ids].T
