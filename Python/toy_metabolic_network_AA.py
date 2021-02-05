"""
A_Exchange:               A[e]<-                         
A_Uptake:                 A[e]-> A[c]
r1:             A[c] + atp[c] -> B[c]
r2:                     B[c]  -> 2 atp[c] + 3 nadh[c] + C[c]
r3:                      C[c] -> 2 nadh[c] + 0.2 C[c]
r5:          C[c] + 2 nadh[c] -> 3 E[c]
OxPhox:       nadh[c] + o2[c] -> 2 atp[c]
ATP_demand:            atp[c] ->
C_sink:                  C[c] ->
O2_EX:                   o2[e]<- 
O2_Uptake:               o2[e]-> o2[c]
E_Exchange              E[e] <=> 
E_Uptake                E[e] <=> E[c]

"""

# %% --- 
from cobra import Model, Reaction, Metabolite

model = Model('toy_metabolism')
# --- Metabolites creation
A_e    = Metabolite("A[e]")
A_c    = Metabolite("A[c]")
B_c    = Metabolite("B[c]")
atp_c  = Metabolite("atp[c]")
nadh_c = Metabolite("nadh[c]")
C_c    = Metabolite("C[c]")
E_c    = Metabolite("E[c]")
o2_e   = Metabolite("o2[e]")
o2_c   = Metabolite("o2[c]")
E_e    = Metabolite("E[e]")
# --- Reactions creation
A_Exchange  = Reaction("A_Exchange")
A_Uptake    = Reaction("A_Uptake")
r1          = Reaction("r1")
r2          = Reaction("r2")
r3          = Reaction("r3")
r5          = Reaction("r5")
OxPhox      = Reaction("OxPhox")
ATP_demand  = Reaction("ATP_demand")
C_sink      = Reaction("C_sink")
O2_Exchange = Reaction("O2_Exchange")
O2_Uptake   = Reaction("O2_Uptake")
E_Exchange  = Reaction("E_Exchange")
E_Uptake    = Reaction("E_Uptake")

A_Exchange.add_metabolites({A_e: -1})
A_Uptake.add_metabolites({ A_e: -1, A_c: 1})
r1.add_metabolites({A_c: -1, atp_c: -1,B_c: 1})
r2.add_metabolites({B_c: -1, atp_c: 2, nadh_c: .5, C_c: 1})
r3.add_metabolites({C_c: -1, nadh_c: 4})
r5.add_metabolites({C_c: -1, nadh_c: -2, E_c: 3})
OxPhox.add_metabolites({nadh_c: -6, o2_c: -1, atp_c: 2})
ATP_demand.add_metabolites({atp_c: -1})
C_sink.add_metabolites({C_c: -1, atp_c: -1})
O2_Exchange.add_metabolites({o2_e: -1})
O2_Uptake.add_metabolites({ o2_e: -1, o2_c: 1})
E_Exchange.add_metabolites({ E_e: -1})
E_Uptake.add_metabolites({E_e: -1, E_c: 1})

model.add_reactions([A_Exchange, A_Uptake, r1, r2, r3, r5, OxPhox, \
                     ATP_demand, C_sink, O2_Exchange, O2_Uptake, \
                     E_Exchange, E_Uptake])




model.reactions.get_by_id("A_Exchange").bounds = (-1000, 0)       
model.reactions.get_by_id("A_Uptake").bounds = (0, 3.3)    
model.reactions.get_by_id("r1").bounds = (0, 1000)     
model.reactions.get_by_id("r2").bounds = (0, 1000)     
model.reactions.get_by_id("r3").bounds = (0, 1000)    
model.reactions.get_by_id("r5").bounds =(-1000, 1000)     
model.reactions.get_by_id("OxPhox").bounds = (0, 1000)     
model.reactions.get_by_id("ATP_demand").bounds = (3, 1000)   
model.reactions.get_by_id("C_sink").bounds = (0, 1000)   
model.reactions.get_by_id("O2_Exchange").bounds = (-1000, 0)    
model.reactions.get_by_id("O2_Uptake").bounds = (0, 10)
model.reactions.get_by_id("E_Exchange").bounds = (-10000, 1000)     
model.reactions.get_by_id("E_Uptake").bounds = (-1000, 2.2)     
import pandas as pd

#pd.DataFrame(list(zip([model.reactions[i].id for i in range(0,len(model.reactions))],
#[model.reactions[i].reaction for i in range(0,len(model.reactions))])))   

# %%
model.objective = "C_sink"
sol = model.optimize()
sol.fluxes

# %%
from cobra.io import save_json_model

save_json_model(model, "toy_metabolism_AA.json")