#!/bin/python
"""
"What do such machines really do? They increase the number of things we can do without thinking. Things we do without thinking — there's the real danger."
    - God-Emperor Leto II 
"""

# %% --- 
from   heapq import merge
import cobra  
import numpy  as np
import pandas as pd
from   cobra.util.array import create_stoichiometric_matrix
import warnings 
warnings.filterwarnings("ignore")
import networkx as nx
 
model_path = '/home/alejandro/NAMU_in_progress/Figures_creation/Code/Results_Acevedo_et_al_2021/01_phpps_bipartite_fluxes_sensis_optimality/stimulated_2021.json'

stimulated = cobra.io.load_json_model(model_path)



    

def get_cobra_model_annotations(model):

    model_rxns = model.reactions
    model_mets = model.metabolites

    rxn_reactions  = [model_rxns[i].reaction for i in range(len(model_rxns)) ]
    rxn_names      = [model_rxns[i].name     for i in range(len(model_rxns)) ]
    rxn_ids        = [model_rxns[i].id       for i in range(len(model_rxns)) ]
    met_names      = [model_mets[i].name for i in range(len(model_mets))]
    met_ids        = [model_mets[i].id for i in range(len(model_mets))]
    met_formula    = [model_mets[i].formula for i in range(len(model_mets))]

    return rxn_reactions, rxn_names, rxn_ids, met_names, met_ids, met_formula

rxn_reactions, rxn_names, rxn_ids, met_names, met_ids, met_formula = get_cobra_model_annotations(stimulated)

fba =  stimulated.optimize()


reaction_info_FBA = pd.DataFrame(list(zip( rxn_ids, rxn_names, rxn_reactions, fba.fluxes, fba.reduced_costs)),
 columns = ['rxn_ids','rxn_names',"rxn_reactions", "fluxes", "reduced_costs"])
#reaction_info.index = reaction_info.rxn_ids

reaction_info_FBA.to_csv('reaction_info_FBA.csv')





