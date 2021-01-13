"""
Creación de un modelo para cobra y futuro analiis de redes  
[Documentación de COBRA](https://cobrapy.readthedocs.io/en/0.17.0/building_model.html)

Toy reactions:

r1: m1 <- 
r2: m1 -> m2
r3: m2 <-
r4: m2 -> 
r5: m2 -> m3
r6: m1 + m5 -> m3 + m4
r7: m3 + m4 -> m1 + m5
r8: m3 + m5 -> m4 
maximize r4 + r8 =    
obj:     m2 + m3 + m5 -> m4

subject to r1 > 1

Agregar el objetivo 

# change the objective to ATPM
model.objective = "obj"

# The upper bound should be 1000, so that we get the actual optimal value
model.reactions.get_by_id("obj").upper_bound = 1000.0
#chequear
linear_reaction_coefficients(model)

#Agregar el consumo de sustrato

model.reactions.get_by_id("r1").bounds = (-1.33, 0.0)




"""
# %% --- Inicialización del modelo
from cobra import Model, Reaction, Metabolite
model = Model('toy_fba')

# %% --- Creación de metabolitos
m1 = Metabolite("m1")
m2 = Metabolite("m2")
m3 = Metabolite("m3")
m4 = Metabolite("m4")
m5 = Metabolite("m5")

""" Best practise: SBML (Systems Biology Markup Language) compliant IDs
reaction = Reaction('3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0.    # This is the default
reaction.upper_bound = 1000. # This is the default
reaction.add_metabolites({ # Reacciones estequiometricas
    malACP_c: -1.0,
    h_c: -1.0
})"""
# %% --- Creación de reacciones 
r1 = Reaction("r1")
r1.add_metabolites({m1 : 1.0})

r2 = Reaction("r2")
r2.add_metabolites({m1: -1.0, m2: 1.0})

r3 = Reaction("r3")
r3.add_metabolites({m2 : 1.0})

r4 = Reaction("r4")
r4.add_metabolites({m2 : -1.0})

r5 = Reaction("r5")
r5.add_metabolites({m2: -1.0, m3: 1.0})

r6 = Reaction("r6")
r6.add_metabolites({m1: -1.0, m5: -1.0, m3: 1.0, m4: 1.0})

r7 = Reaction("r7")
r7.add_metabolites({m3: -1.0, m4: -1.0, m1: 1.0, m5:1.0})

r8= Reaction("r8")
r8.add_metabolites({m3: -1.0, m5: -1.0, m4:1.0})

# %% --- Chequear reacciones 
import numpy as np

print(
r3.reaction,
r4.reaction,
)
# %% 
from cobra.util import create_stoichiometric_matrix

model.add_reactions([r1, r2, r3, r4, r5, r6, r7, r8])
model.reactions.get_by_id('r3')

create_stoichiometric_matrix(model)


# %% --- Crea la reacción para objetivo de optimización

obj = Reaction("obj")
obj.add_metabolites({m2: -1.0, m3: -1.0, m5: -1.0, m4: 1.0}) 


#ptimización del modelo
 
model.obj
ective = "obj"
model.reactions.get_by_id("obj").upper_bound = 1000.0
# %%

m
odel.reactions
