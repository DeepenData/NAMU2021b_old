"""Prueba de las funciones previamente definidas en un grafo usando un JSON con
las reacciones del metabolismo humano
"""
# %% MM --- 2021-01-18 16:40 --- Importa el modelo JSON
import cobra 
import networkx as nx
import pandas as pd

# %% MM --- 2021-01-18 16:40 --- Funciones utilitarias
def cobra2networkx(modelo, direccionado=True):
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
    from sklearn.preprocessing import binarize
    from scipy.sparse import csr_matrix

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"

    grafo = binarize(abs(create_stoichiometric_matrix(modelo))) # Crea una matriz
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

def list2attr(grafo, nodos, nombre, atributos):
    """Toma dos listas: nombres de nodos y atributos; y las añade a un grafo

    Parameters
    ----------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    nodos: list
        Una lista de nombres de nodos
    nombre : str
        Nombre del atributo siendo asignado. Ie. nombre de la columna en Gephi
    atributos : list
        Una lista de valores para los nodos en la lista anterior.
        Ambas listas deben ser del mismo largo. 
    
    Returns
    -------
    grafo: bipartite_graph
        Un grafo bipartito con un nuevo atributo para un set de nodos "nodos". 
    """
    assert len(nodos) == len(atributos), "Ambas listas deben ser del mismo largo."
    tmp_list = { nodos[i] : atributos[i] for i in range(0, len(atributos)) }
    
    from networkx import set_node_attributes
    set_node_attributes(grafo, tmp_list, nombre ) # Añade los atributos
    
    return grafo

def attr2partition(grafo, lista, nombre, asignar=0):
    """Añade atributos a todo un grafo, solo metabolitos, o solo reacciones.
    Parameters
    ----------
    grafo:
        Un objeto de NetworkX.
    lista: list
        Una lista de atributos a asignar a un subset de nodos. 
        Debe ser del largo del subset de metabolitos, reacciones, o todo.
    nombre: str
        El nombre bajo el que se asignaran los atributos
    asignar: int, default 0
        A que subset se esta asignando. Por defecto todos los nodos. 
        0 = todos, 1 = metabolitos, 2 = reacciones
    Returns
    -------
    grafo:
        Un objeto de NetworkX con nuevos atributos. 
    """
    assert asignar in [0,1,2], "La asignación debe ser 0 = todos, 1 = metabolitos, o 2 = reacciones"

    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones

    print("Metabolitos:", len(metabolites_nodes), " | Reacciones:", len(reactions_nodes)," | Z:", len(metabolites_nodes) + len(reactions_nodes)) # Debug

    if   asignar == 0: 
        nodos = metabolites_nodes + reactions_nodes  # 0: todos
        #assert len(lista) == len(nodos),             "El largo debe ser len(grafo.nodes)"
    elif asignar == 1: 
        nodos = metabolites_nodes                    # 1: metabolitos
        #assert len(lista) == len(metabolites_nodes), "El largo no coincide con los metabolitos"
    elif asignar == 2: 
        nodos = reactions_nodes                      # 2: reacciones
        #assert len(lista) == len(reactions_nodes),   "El largo no coincide con las reacciones"

    from networkx import set_node_attributes
    set_node_attributes(grafo, {nodos[i] : lista[i] for i in range(0, len(lista))}, nombre)
    
    return grafo

def solution2attr(solucion, grafo):
    """Docstring
    """

    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones

    print("Metabolitos:", len(metabolites_nodes), " | Reacciones:", len(reactions_nodes)) # Debug

    shadow_prices = solucion.shadow_prices.tolist() # De solucion, pasa a lista
    fluxes        = solucion.fluxes.tolist()        # posiblemente es más rapido de usar que desde el modelo
    reduced_costs = solucion.reduced_costs.tolist() # considerando el menor parsing pero si requiere pandas

    print("Shadow prices:", len(shadow_prices),"| Fluxes:", len(fluxes),"| Reduced costs:", len(reduced_costs)) # Debug

    from networkx import set_node_attributes

    # ASIGNA Shadow_prices, fluxes, reduced_costs
    set_node_attributes(grafo, { metabolites_nodes[i] : shadow_prices[i] for i in range(0, len(shadow_prices)) } , "Shadow Price") 
    set_node_attributes(grafo, { reactions_nodes[i] : fluxes[i] for i in range(0, len(fluxes)) } , "Flux") 
    set_node_attributes(grafo, { reactions_nodes[i] : reduced_costs[i] for i in range(0, len(reduced_costs)) } , "Reduced cost") 

    #grafo.nodes(data=True) # Debug

    return grafo

# %% MM --- 2021-01-18 16:44 --- Empieza cosas interesantes (human)

human = cobra.io.load_json_model("GEM_Recon2_thermocurated_redHUMAN.json") # INPUT

grafo = cobra2networkx(human) # Función interesante 1

solucion = human.optimize()

grafo = solution2attr(solucion, grafo) # Función interesante 2

from networkx import write_gexf
write_gexf(grafo, "human_thermo2.gexf")

# %% AA --- 2021-01-19 16:47 --- Cosas interesantes con modelo de cerebro

brain = 

grafo_brain = cobra2networkx(brain)
solucion_brain = brain.optimize()
grafo_brain = solution2attr(solucion_brain, grafo_brain)

# %% MM --- 2021-01-19 17:20 --- Calculos de centralidades (funcionalizados!)