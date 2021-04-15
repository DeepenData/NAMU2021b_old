# %%
import cobra
import warnings 
warnings.filterwarnings("ignore")
path          = '/home/alejandro/PostDoc/human-metnet/data/'
model    = cobra.io.load_json_model(path + 'GEM_Recon2_thermocurated_redHUMAN.json') 
# %%
#biomassi
import pandas as pd 

rxn_ids      = [model.reactions[i].id for i in range(len(model.reactions))]
rxn_names    = [model.reactions[i].name  for i in range(len(model.reactions))]
reactions    = [model.reactions[i].reaction for i in range(len(model.reactions))]
rxn_bounds    = [model.reactions[i].bounds for i in range(len(model.reactions))]


met_ids      = [model.metabolites[i].id       for i in range(len(model.metabolites))]
met_names    = [model.metabolites[i].name     for i in range(len(model.metabolites))]
met_formulas = [model.metabolites[i].formula  for i in range(len(model.metabolites))]

gene_ids = [model.genes[i].id  for i in range(len(model.genes))]

genes_rxns = []
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
    'rxn_bounds':rxn_bounds
})

genes_and_rxns.to_csv("./data/recon2_genes_to_rxns.csv")
df_rxns.to_csv("./data/recon2_rxn_info.csv")

# %%
import re
import numpy as np

def filter_list_by_regex(mylist, pattern):
    r = re.compile(pattern)
    newlist = list(filter(r.match, mylist)) # Read Note
    return newlist 

filter_list_by_regex(met_ids, 'phe_L_c')

#df_rxns[df_rxns['reactions'].str.contains('.* $', case=False, regex=True)]


l1 = list(model.metabolites.get_by_id('phe_L_e').reactions)

[l1[i].reaction for i in range(len(l1))]

# %%    

df_tyr_exch = df_rxns[df_rxns['reactions'].str.contains('tyr_L_e', case=False, regex=True)]

tyr_uptakes = df_tyr_exch[df_tyr_exch['reactions'].str.contains('<.*tyr_L_e|<=>|tyr_L_e.*>', case=False, regex=True)]
tyr_uptakes_ids = list(tyr_uptakes.rxn_ids)
[model.reactions.get_by_id(tyr_uptakes_ids[i]).bounds for i in  range(len(tyr_uptakes_ids))]

# %%  

#model.reactions.get_by_id('10FTHF5GLUtm')
# %% 
#SOLO DEJARLOS SALIR 


tyr_uptakes_ids = list(tyr_uptakes.rxn_ids)

for i in range(len(tyr_uptakes_ids)):
    model.reactions.get_by_id(tyr_uptakes_ids[i]).bounds = 0,0

[model.reactions.get_by_id(tyr_uptakes_ids[i]).bounds for i in  range(len(tyr_uptakes_ids))]

# %% 
fba = model.optimize().to_frame()

fba.loc[['PHETHPTOX2','r0399'],]
# %% 

objective = df_rxns[df_rxns['rxn_names'].str.contains('.*mass.*', case=False, regex=True)]
pd.options.display.max_colwidth = 1000
objective

biomass_mets = list(model.reactions.get_by_id('biomass').metabolites)

[biomass_mets[i].id for i in range(len(biomass_mets))]