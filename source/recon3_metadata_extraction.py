# %% --- 
import cobra
import pandas as pd
import numpy as np
from cobra.sampling import OptGPSampler, ACHRSampler, sample
import warnings 
warnings.filterwarnings("ignore")
#GEM_Recon2_thermocurated_redHUMAN = cobra.io.load_matlab_model('GEM_Recon2_thermocurated_redHUMAN.mat')
#cobra.io.save_json_model(GEM_Recon2_thermocurated_redHUMAN, "GEM_Recon2_thermocurated_redHUMAN.json")

path          = '/home/alejandro/AA_PostDoc/human-metnet/data/'
model    = cobra.io.load_json_model(path + 'GEM_Recon3_thermocurated_redHUMAN.json') 
import pandas as pd 

import pandas as pd 

rxn_ids        = [model.reactions[i].id for i in range(len(model.reactions))]
rxn_names      = [model.reactions[i].name  for i in range(len(model.reactions))]
reactions      = [model.reactions[i].reaction for i in range(len(model.reactions))]
rxn_bounds     = [model.reactions[i].bounds for i in range(len(model.reactions))]
rxn_subsystems = [model.reactions[i].subsystem for i in range(len(model.reactions))]

met_ids         = [i.id          for i in model.metabolites]
met_names       = [i.name        for i in model.metabolites]
met_formulas    = [i.formula     for i in model.metabolites]
met_compartment = [i.compartment for i in model.metabolites]
met_charge      = [i.charge      for i in model.metabolites]

gene_ids     = [model.genes[i].id  for i in range(len(model.genes))]


mets_df = pd.DataFrame({
    'ID' : met_ids,
    'Name' : met_names,
    'Formula' : met_formulas,
    'Compartment' : met_compartment,
    'Charge' : met_charge
})


genes_rxns   = []
for i in range(len(model.genes)):
    rxn_set    = list(model.genes[i].reactions)
    gene_rxns = [rxn_set[j].id for j in range(len(rxn_set))]
    genes_rxns.append(gene_rxns)


genes_and_rxns = pd.DataFrame({
    'Gene_ID' : gene_ids,
    'Reactions' : genes_rxns
})


df_rxns = pd.DataFrame({
    'rxn_ids' : rxn_ids,
    'rxn_names' : rxn_names,
    'reactions'  : reactions,
    'rxn_bounds':rxn_bounds,
    'rxn_subsystems': rxn_subsystems
})

mets_df.to_csv("./data/recon3_metabolite_metadata.csv")
genes_and_rxns.to_csv("./data/recon3_genes_to_reactionIDs.csv")
df_rxns.to_csv("./data/recon3_reactions_metadata.csv")

# %%
