#!/usr/bin/env python
"""Un programa que importa modelos dados en ./data/ y crea un agragado
Por ahora, se limita a duplicar modelos
Parametros: modelo_base, *modelos_suma
"""
# %% --- Importar modelo
from cobra.io import load_json_model, save_json_model
model = load_json_model('toy_model.json')

# Lista de modelos a agregar
# TODO: convertir esto en input del programa

models = [model]

print ('%i modelos importados' % len(models)) # Sanity check

# %% --- Revisar conflictos de ID
"""Al tener varios modelos en un mismo espacio, los nombres de metabolitos y
reacciones chocarán. Por eso, el nombre de estos idealmente vendria sobrepuesto
por el ID del modelo. Pero si dos modelos tienen el mismo ID, como el caso de mi
intento de duplicación, necesito una forma de resolver eso idealmente añadiendo
un contador al ID del modelo. Eg. [mod_a, mod_a] -> [mod_a_1, mod_a_2]"""

i = 1 # Parte en 1 porque 0 es el modelo base
for model in models:
    model.id = str(i) + '_' + model.id
    i = i + 1 # Contador de ID en el array de modelos

# %% --- Renombar reacciones y metabolitos
"""Esta sección se encarga de obtener las reacciones en una lista de modelos,
cambiar el nombre de las reacciones y metabolitos para incluir el ID del modelo.
Consideraciones:
- Reacciones boundry pueden comportarse de forma extraña
- Algunas reacciones boundry pueden importar directamente a la célula; que pasa
si chocan con otras celulas que importan directamente?
- NO ES NECESARIO RENOMBRAR COMPARTIMIENTOS, porque aunque esten en el mismo
  lugar los metabolitos no interactuarian entre si, por las nuevas reacciones"""

from cobra import Model, Reaction, Metabolite

# Omitimos los metabolitos externos, dado que eventualmente caen a un pool comun
EXTERNAL = ['e'] # Valores para compartimientos no-internos
is_external = lambda metabolite : ((metabolite.id.endswith('[e]')) or (metabolite.compartment in EXTERNAL))

def subpend_id(model):
    """Renames IDs for reactions and metabolites, model.id"""
    for metabolite in model.metabolites:
        if (is_external(metabolite)): # Metabolito es exterior
            continue
        metabolite.id = model.id + '_' + metabolite.id
    model.repair() # Se que esta duplicado pero mejor asegurarme
    
    for reaction in model.reactions:
        # Posiblemente una linea para eliminar los boundary ?
        reaction.id   = model.id + '_' + reaction.id
    model.repair()

for model in models:
    subpend_id(model)
    # Posiblemente remover boundary

# %% --- Inicialización del agragado
"""Carga el modelo base que se utiliza para añadir el resto de los modelos. 
Una ventaja de usar un modelo base en lugar de uno vacio es que así contamos con
el medio, una reacción a optmizar, y en general los atributos que nos sirven 
para continuar con el modelado"""

BASE = 'base.json'

aggregated = load_json_model(BASE)
aggregated.id = '0_' + aggregated.id # Nuevo ID inicializanco en cero

subpend_id(aggregated) # Cambia el nombre de las reacciones y metabolitos

# %% --- Agregador super-agregado
for model in models: 
    aggregated.add_metabolites( model.metabolites )
    aggregated.add_reactions( model.reactions )

print('%i reacciones  finales' % len(model.reactions))
print('%i metabolitos finales' % len(model.metabolites))

# %% --- Sopa
"""Esta sección se encarga de exportar un JSON finla agregado, que es una sopa"""

OUTPUT = 'soup.json'
save_json_model(aggregated, OUTPUT)