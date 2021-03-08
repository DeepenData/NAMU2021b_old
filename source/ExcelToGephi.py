#!/bin/python

# %% --- Importa dataset de metabolitos

import pickle

# Importa el dataset pandas de metabolitos
infile = open('../tmp/excel_metabolitos','rb'); df_metabolitos = pickle.load(infile); infile.close()

# %% --- Importando modulos varios

import numpy    as np
import pandas   as pd
import networkx as nx

# %% --- Calcula promedios para asignar a los nodos de interes

from scipy import stats

# Es el unico que no require lambdas raros :D
aritmetic = df_metabolitos.mean()

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
geometric_mean = lambda array :  stats.gmean( array[array > 0] ) * len(array[array != 0])/len(array)
geometric = df_metabolitos.apply( geometric_mean, axis=0, raw=True )

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
harmonic_mean = lambda array : stats.hmean( array[array > 0] ) * len(array[array != 0])/len(array)
harmonic  = df_metabolitos.apply( harmonic_mean, axis=0, raw=True )

# Define una función helper para cuadraticos (Root Square Mean). No en scipy base
quadratic_mean = lambda array : np.sqrt( np.nanmean( array*array ) )
quadratic = df_metabolitos.apply( quadratic_mean, axis=0, raw=True )

# Ensambla un DataFrame

df_means = pd.DataFrame({
    'aritmetic' : aritmetic,
    'geometric' : geometric,
    'harmonic'  : harmonic,
    'quadratic' : quadratic
}, dtype = 'float32')

# %% --- ELIMINAR --- Definde funciones de Grafos
# TODO: Cambiar nombres de archivos para permitir importación de modulos

def cobra_to_bipartite_graph(modelo, direccionado=True):
    """Toma un modelo de cobra y genera un grafo bipartito de NetworkX

    Parameters
    ----------
    model : cobra.core.model.Model
        Un modelo no-optimizado de reacciones metabolicas en COBRA
    direccionado : bool
        Selecciona si el objeto creado sera un grafo direccionado o no.
        Por defecto crea grafos direccionados. 

    Returns
    -------
    grafo : 
        un grafo bipartito de NetworkX. Incluye atributos de flujo del modelo y
        sensibilidad de metabolitos y reacciones. 
    
    Notes
    -----
    - nx.write_gexf(grafo, "grafo.gexf") Crea una salida para Gephi
    """

    from cobra.util import create_stoichiometric_matrix
    from networkx import relabel_nodes, nodes
    from networkx.algorithms.bipartite.matrix import biadjacency_matrix, from_biadjacency_matrix
    from scipy.sparse import csr_matrix
    import numpy as np

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"

    grafo = abs(create_stoichiometric_matrix(modelo)) # Crea una matriz
    grafo = (grafo > 0.0).astype(np.int_) # Binarización
    grafo = csr_matrix(grafo) # Convierte la matriz a una matriz dispersa
    grafo = from_biadjacency_matrix(grafo) # Usa la dispersa para un grafo bipartito
    
    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    metabolites_n     = len(modelo.metabolites) # Numero de metabolitos
    metabolites_names = [modelo.metabolites[i].id for i in range(0, metabolites_n) ]

    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
    reactions_n       = len(modelo.reactions)   # Numero de reacciones
    reactions_names   = [modelo.reactions[i].id   for i in range(0, reactions_n)   ]

    names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
    grafo = relabel_nodes(grafo, names_mapped)
    return grafo

# %% --- Importa modelo Cobra y hace el grafo bipartito

from cobra.io import load_json_model

import warnings; warnings.filterwarnings('ignore') # Ignora warnings
model = load_json_model('../data/stimulated_2021.json') # MODELO INPUT

grafo = cobra_to_bipartite_graph(model) # Crea el bipartito

warnings.filterwarnings('ignore') # Reestablece warnings

# %% --- Exploratorio del modelo para asignaciones

# met = pd.Series(
#     data  = [ met.name for met in model.metabolites],
#     index = [ met.id   for met in model.metabolites],
# )
# 
# met.sort_values(inplace = True) # Alfabetico para legibilidad
# 
# met[met.str.contains('carnitine', case=False)] # Busquedas

# %% --- Diccionario de cosas para el grafo de cobra 

# Esto es una lista curada. 
# El modelo no cuenta con la mayor parte de las carnitinas consideradas en las
# pruebas bioquimicas, por lo que conviene tener más información de las pruebas 
# para la asignación al grupo más cercano

nodes_to_aply = {
    'Phe' : ['phe-L[mA]', 'phe-L[mN]'],
    'Met' : [],
    'Val' : ['val-L[I]', 'val-L[cA]', 'val-L[cN]', 'val-L[e]', 'val-L[mA]', 'val-L[mN]'],
    'Leu/Ile' : [ 'leu-L[I]', 'leu-L[cA]', 'leu-L[cN]', 'leu-L[e]', 'leu-L[mA]', 'leu-L[mN]', 
    'ile-L[I]', 'ile-L[cA]', 'ile-L[cN]', 'ile-L[e]', 'ile-L[mA]', 'ile-L[mN]'],
    'tir' : ['tyr-L[mA]', 'tyr-L[mN]'],
    'Pro' : ['pro-L[I]', 'pro-L[cA]', 'pro-L[cN]', 'pro-L[e]', 'pro-L[mA]', 'pro-L[mN]'],
    'Arg' : [ 'arg-L[I]', 'arg-L[cA]', 'arg-L[cN]', 'arg-L[e]', 'arg-L[mA]', 'arg-L[mN]'],
    'Gly' : ['gly[I]', 'gly[cA]', 'gly[cN]', 'gly[e]', 'gly[mA]', 'gly[mN]'], 
    'Ala' : ['ala-B[I]', 'ala-B[cA]', 'ala-B[cN]', 'ala-B[e]', 'ala-B[mA]', 'ala-B[mN]', 
    'ala-L[I]', 'ala-L[cA]', 'ala-L[cN]', 'ala-L[e]', 'ala-L[mA]', 'ala-L[mN]'],
    'Asp' : ['asp-L[cA]', 'asp-L[cN]', 'asp-L[mA]', 'asp-L[mN]'],
    'Glt' : ['glu-L[I]', 'glu-L[cA]', 'glu-L[cN]', 'glu-L[e]', 'glu-L[mA]', 'glu-L[mN]'],
    'Cit' : ['cit[cA]', 'cit[cN]', 'cit[mA]', 'cit[mN]'], 
    'Orn' : ['orn[I]', 'orn[cA]', 'orn[cN]', 'orn[e]', 'orn[mA]', 'orn[mN]'],
    'SA' : [],
    'C0' : [ 'crn[I]', 'crn[cA]', 'crn[cN]', 'crn[e]', 'crn[mA]', 'crn[mN]'],
    'C2' : ['acrn[cA]', 'acrn[cN]', 'acrn[mA]', 'acrn[mN]'],
    'C3' : ['pcrn[mA]', 'pcrn[mN]'], 
    'C4' : [], 
    'C4OH/C3DC' : [], 
    'C5-OH/C4DC' : [], 
    'C5'    : [], 
    'C5DC'  : [], 
    'C5:1'  : [], 
    'C6'    : [], 
    'C6DC'  : [], 
    'C8'    : [], 
    'C8:1'  : [], 
    'C10'   : [], 
    'C10:1' : [], 
    'C10:2' : [], 
    'C12'   : [], 
    'C12:1' : [], 
    'C14'   : [], 
    'C14:1' : [], 
    'C14:2' : [], 
    'C14OH' : [], 
    'C16'   : [], 
    'C16OH' : [], 
    'C16:1' : [], 
    'C16:1OH' : [], 
    'C18'   : [], 
    'C18OH' : [], 
    'C18:1' : [], 
    'C18:1OH' : [], 
    'C18:2' : []
}

# %% --- Aplica atributos en un For Loop

for i in df_means.index:
    # Crea un diccionario con la lista de nodos, y una lista de los valores repetirdos de los promedios
    d1 = { ii : df_means.loc[ i ,'aritmetic'] for ii in nodes_to_aply[ i ] }
    d2 = { ii : df_means.loc[ i ,'geometric'] for ii in nodes_to_aply[ i ] }
    d3 = { ii : df_means.loc[ i ,'harmonic' ] for ii in nodes_to_aply[ i ] }
    d4 = { ii : df_means.loc[ i ,'quadratic'] for ii in nodes_to_aply[ i ] }

    # Empieza a aplicar los atributos
    nx.set_node_attributes( grafo , d1 , 'aritmetic') 
    nx.set_node_attributes( grafo , d2 , 'geometric')
    nx.set_node_attributes( grafo , d3 , 'harmonic' )
    nx.set_node_attributes( grafo , d4 , 'quadratic')

# %% --- Exporta a Gephhi

from networkx import write_gexf
write_gexf(grafo, "stimulatedasd.gexf")
print('done!')

# %% --- Genera histogramas de distribución de promedios
# TODO: Generar histogramas de distribución de promedios
