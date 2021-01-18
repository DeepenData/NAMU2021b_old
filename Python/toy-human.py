"""Prueba de las funciones previamente definidas en un grafo usando un JSON con
las reacciones del metabolismo humano
"""
# %% MM --- 2021-01-18 16:40 --- Importa el modelo JSON
import cobra 
import networkx as nx
import pandas as pd

human = cobra.io.load_json_model("GEM_Recon2_thermocurated_redHUMAN.json")

# %% MM --- 2021-01-18 16:40 --- Funciones utilitarias
def cobra2networkx(model, direccionado=True):
    """Toma un modelo de cobra y genera un grafo bipartito de NetworkX

    Parameters
    ----------
    model : cobra_model
        Un modelo de reacciones metabolicas en COBRA
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
    import networkx
    from cobra.util import create_stoichiometric_matrix

    from networkx.algorithms.bipartite.matrix import biadjacency_matrix      # Extrae matriz de adyacencia
    from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix # Crea el grafo
    from sklearn.preprocessing import binarize
    from scipy.sparse import csr_matrix

    tmp = binarize(abs(create_stoichiometric_matrix(model))) # Crea una matriz
    tmp = csr_matrix(tmp) # Convierte la matriz a una matriz dispersa
    tmp = from_biadjacency_matrix(tmp) # Usa la dispersa para un grafo bipartito
    # Eventualmente hacer gtrafos direccionoas y no-direccionados
    # total_nodes = range(0 , len(tmp.nodes(data=False))) # Rango de nodos
    # particiones = [tmp.nodes(data=True)[i]["bipartite"] for i in total_nodes] # Tipo nodo

    metabolites_nodes = [n for n, d in tmp.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    metabolites_n     = len(model.metabolites) # Numero de metabolitos
    metabolites_names = [model.metabolites[i].id for i in range(0, metabolites_n) ]

    reactions_nodes   = [n for n, d in tmp.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
    reactions_n       = len(model.reactions)   # Numero de reacciones
    reactions_names   = [model.reactions[i].id   for i in range(0, reactions_n)   ]

    names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
    tmp = networkx.relabel_nodes(tmp, names_mapped)
    
    grafo = tmp # Asigna tmp al grafo
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
    import networkx as nx
    assert len(nodos) == len(atributos), "Ambas listas deben ser del mismo largo."
    tmp_list = { nodos[i] : atributos[i] for i in range(0, len(atributos)) }
    nx.set_node_attributes(grafo, tmp_list, nombre ) # Añade los atributos
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
    asignar: int
        A que subset se esta asignando. Por defecto todos los nodos. 
        0 = todos, 1 = metabolitos, 2 = reacciones
    Returns
    -------
    grafo:
        Un objeto de NetworkX con nuevos atributos. 
    """
    assert asignar in [0,1,2], "La asignación debe ser 0 = todos, 1 = metabolitos, o 2 = reacciones"
    
    #total_nodes = range(0 , len(grafo.nodes(data=False))) # Rango de nodos
    #particiones = [grafo.nodes(data=True)[i]["bipartite"] for i in total_nodes] # Tipo nodo (0:met, 1:rxn)

    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones

    if   asignar == 0: 
        nodos = metabolites_nodes + reactions_nodes  # 0: todos
        assert len(lista) == len(nodos),             "El largo debe ser len(grafo.nodes)"
    elif asignar == 1: 
        nodos = metabolites_nodes                    # 1: metabolitos
        assert len(lista) == len(metabolites_nodes), "El largo no coincide con los metabolitos"
    elif asignar == 2: 
        nodos = reactions_nodes                      # 2: reacciones
        assert len(lista) == len(reactions_nodes),   "El largo no coincide con las reacciones"

    grafo = list2attr(grafo, nodos, lista, nombre)
    return grafo

def fba_solutions(modelo):
    """Extrae flujos y sensibilidades de un modelo optimizado

    Parameters
    ----------
    modelo : Solution
        Un modelo COBRA optimizado. La función fallara sin los parametros 
        creados por la optimización. Utiliza modelo.optimize().
    
    Returns
    -------
    shadow_prices : list
        Sensibilidad de los metabolitos. 
    fluxes : list
        Flujo de metabolitos en cada reacción.
    reduced_cost : list
        Sensibilidad de las reacciones. 
    """
    #assert ,"El modelo no esta optimizado. Utiliza modelo.optimize()"

    shadow_prices = [modelo.metabolites[i].shadow_price for i in range(0, len(modelo.metabolites)) ]
    fluxes        = [modelo.reactions[i].flux           for i in range(0, len(modelo.reactions  )) ]
    reduced_costs = [modelo.reactions[i].reduced_cost   for i in range(0, len(modelo.reactions  )) ]

    return shadow_prices, fluxes, reduced_costs

# %% MM --- 2021-01-18 16:44 --- Empieza cosas interesantes

human = human.optimize() # Resolución del FBA

grafo = cobra2networkx(human) # Crea el bipartito

# %% MM --- 2021-01-18 16:51 --- Añade atributos al grafo

shadow_prices, fluxes, reduced_costs = fba_solutions(human)

grafo  = attr2partition(grafo, shadow_prices, "Shadow Price", 1) # Asigna sensibilidad (met)
grafo  = attr2partition(grafo, fluxes, "Flux", 2) # Asigna flujos 
grafo  = attr2partition(grafo, reduced_costs, "Reduced cost", 2) # Asigna sensibilidad (rxn)

# %% MM --- 2021-01-18 16:51 --- Exporta a Gephi

nx.write_gexf(grafo, "human_thermo2.gexf") # Crea una salida para Gephi