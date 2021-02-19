#!/usr/bin/env python
"""Un programa que importa un modelo y genera una tabla de centralidades con un resumen de las reacciones, y el log2 
entre la reacción sin modificar y modificada. Genera un grafo gephx con esta información codificada en los nodos. 
Parametros: INPUT_MODEL, INPUT_REMOVED_NODES, [N_WORKERS]
Output: OUTPUT_NODE_DELTA, OUTPUT_RESUME_TABLE, OUTPUT_GRAPH


Matrix_rows_cols_width

A_m_c: Centralidades (c = 1:8) de todos los nodos (m=1:1000). Dims = 1000 x 8. // Centralidades base
tA_1_m_c: Centralidades base como tensor de alto 1. 

B_n_m_c: Centralidades (c) de los nodos  (m=1, ... NA(en la posición de nodo n) ...,1000) removiendo el nodo n. Dims: 1000 x 999+(NaN) x 8

D_m_c_n: Es lo mismo que B_n_m_c pero con las centralidades en la segunda dimensión.

ms: Nodos que pertenecen al subsistema s. Dims = length({ nodos E s}) x 1.

A_ms_c: Lo mismo que A_m_c pero indexado para ms.
B_n_ms_c: Lo mismo que B_n_m_c pero indexado para ms. // B_m_c[ms, *, *]
D_ms_c_n: Lo mismo que D_m_c_n pero indexado para ms.

E_c: Promedio por centralidad de todos los nodes en ms. // mean(axis='c')
F_c_n:Promedio por centralidad de todos los nodes en ms removiendo n. // E_c.mean(axis='n')

## Operaciones:

Para los nodos del conjunto ms1, tal que  ms1 pertenece el subsistema s1.

l_ms1 = length(ms1) : el número de nodos en el subsistema s1. // len(ms1)

- A_ms1_c (l_ms1 x 8) => calcular el promedio de cada columna ->  (1 x 8) -> transponer => E_c (8 x 1) // by_columns.means()
- E_c = columns_means(A_ms1_c)
   
- D_ms1_c_n (l_ms1 x 8 x n, donde n = son los nodos removidos = 1000)  => calcular el promedio de cada columna (c = [1:8]) => 
  G( 1 x 8 x n) => colapsar la primera dimensión  => matriz( 8 x n) = (c x n) => F_c_n // cosa = cosa[0] 
  pseudocode : means_by_cols(D_ms1_c_n) = G( 1 x 8 x n) => esto es un arreglo de largo n de vectores fila de tamaño 1 x 8 
              o un tensor aplanado.
              G => cada vector fila (1 x 8) apendiarlo en una nueva matriz M tal que M tenga dims 8 x n = F_c_n.

- F_c_n (8 x n) => extraer cada columna (n) por separado, formar una matriz diagonal (8 x 8) desde cada columna e invertir (elevar a -1).
            Guardar todas las matrices diagonales en un tensor (o ndarray) => F_inv_c_c_n (8 x 8 x 1000).
            pseudocode : F_inv_c_c_n = 
                F_c_n.T.aslist() // lista len() = 1000 de arrays len() = 8 // 
                for n in F_c_n : diagonalize(n)
                F_c_n.asarray() // ahora vuelve a ser un array de 1000 x 8 x 8...
                transpose(?) ... array 8 x 8 x 1000                 

- F_inv_c_c_n => multiplicar cada capa (matrices diagonales 8 x 8) del tensor por el vector E_c (8 x 1) =>
                cada unos de los vectores resultantes (8 x 1, estas son las razones) apendiarlos en una matrix. => (8 x n) => 
                transponer => Ratios (n x 8).
                pseudocode : F_inv_c_c_n(...,...,i) x E_c = r_i
                    Ratios = transpose(matrix([r_1, r_2, ... r_i]))

- Ratios => aplicar log2 a cada entrada de la matriz => FC (1000 x 8). pseudocode: FC_s1 = log2(Ratios).
  ***FC_s1 (1000 x 8) : row = nodo removido, col= fold change de una medida de centralidad (para el subsistema s1).
  *** interpretación de FC_s1: el efecto (contribución) de cada nodo sobre las distintas centralidades (c=[1:8]) del subsistema s1.

- FC_s1 => repetir todo lo anterior pero con los nodos (ms2) del siguiente subsistema (s2) => FC_s2 => 
  iterar s (número de subsistemas) veces =>
  FC_s1, FC_s2, ..., FC_ss (j = 1,..., s)

- Construir el tensor final (tFC) apendiando las capas: FC_s1, FC_s2, ..., FC_ss (j = 1,..., s) => tFC (1000 x 8 x s).
  pseudocode: tFC = tensor(FC_s1, FC_s2, ..., FC_ss).

git add source/delta_centrality.py
git commit -m "Sync lacra"
git pam

"""

import ray
import time
import os
# Conectando a cluster de Ray inicializado con sbatch
ray.init(address='auto', _node_ip_address=os.environ["ip_head"].split(":")[0], _redis_password=os.environ["redis_password"])

# --- Definición de funciones

def eight_centralities(grafo):
    """Calcula centralidades base de una red, usando ocho metodos de NetworkX. 
    Los ultimos dos metodos son omitidos si el grafo no es totalmente conectado,
    y devuelve diccionarios vacios en cambio. 

    ```   
    ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 
    'betweenness_centrality', 'closeness_centrality', 'load_centrality', 
    'information_centrality', 'second_order_centrality']
    ```

    Returns
    -------
    centralities : list
        Una lista de diccionarios {nodo : centralidad} con las centralidades
        segun la lista descrita arriba."""

    import networkx as nx

    hc  = nx.harmonic_centrality(grafo, nbunch=None, distance=None)
    ec  = nx.eigenvector_centrality(grafo, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc  = nx.degree_centrality(grafo)
    bc  = nx.betweenness_centrality(grafo, normalized=True, weight=None, endpoints=False, seed=None)

    """↓ Here be dragons ↓"""
    cc  = nx.closeness_centrality(grafo, distance=None, wf_improved=True)
    lc  = nx.load_centrality(grafo, cutoff=None, normalized=True, weight=None)

    """ Requieren un grafo full conected """
    if not (nx.is_connected(grafo)) : 
        ic = {}; soc = {}
    else:
        ic  = nx.information_centrality(grafo)
        soc = nx.second_order_centrality(grafo)

    centralities = [hc, ec, dc, bc, cc, lc, ic, soc]

    return centralities

@ray.remote # Usa el decorador de Ray para definir esta como la paralelizable
def delta_centrality(G, removed_nodes):
    """Toma un grafo y una lista de nodos y calcula la centralidad del grafo al
    remover cada uno de esos nodos.

    Prefix y subfix son información sobre la lista de los nodos removidos. El
    atributo final en el grafo es ≪`PREFIX` REMOVED_nodo[i] `SUBFIX` ≫

    Parameters
    ----------
    grafo : 
        Un grafo de NetworkX 
    removed_nodes :
        Una lista de nodos a remover. Puede ser grafo.nodes

    Return
    ------
    grafo :
        Un grafo de NetworkX con atributos de centralidad por cada nodo removido
    centralities : list
        Una lista de diccionarios con calculos de centralidades
    breaks : list
        Lista de nodos que rompen la continuidad del espacio-tiempo (y la red)
    """

    removed_nodes = list(removed_nodes)

    import networkx as nx
    from copy import deepcopy

    breaks = []; centralities = []
    for node in removed_nodes:
        print( 'Iterando en:', node )
        delta = deepcopy( G )
        delta.remove_node( str(node) ) # Elimina el nodo de la iteración
        centrality = eight_centralities(delta)
        if (centrality[-1] == {}):
            print( str(node), 'breaks continuity!')
            breaks.append( str(node) )
        centralities.append( centrality )

    return centralities, breaks

# --- Params

import sys

INPUT_MODEL = sys.argv[1] # Modelo base

# --- Sección que importa modulos

import numpy as np
import networkx as nx
import pandas as pd

## --- Modelo de Cobra
t1 = time.time()

from cobra.io import load_json_model
model = load_json_model(INPUT_MODEL)

print("Iniciando optimización del modelo")
#solution_fba = model.optimize() # Optimización del modelo para check
#solution_fba.fluxes # Check si los flujos funcionan

t2 = time.time(); print('\n','TIME Importando y optimizando modelo:', (t2-t1)*1000 ,'ms')

## Modelo a grafo de NetworkX
t1 = time.time() # Contador de esta sección

from cobra.util.array import create_stoichiometric_matrix

S_matrix = create_stoichiometric_matrix(model)
S_matrix = (abs(S_matrix) )
S_matrix = (S_matrix > 0.0).astype(np.int_)

projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
np.fill_diagonal(projected_S_matrix, 0)

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)

G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )

if not (nx.is_connected(G)):
    # Check si este grafo tiene nodos huerfanos
    print("This graph is not conected. Picking largest component.")
    
    largest_component = max(nx.connected_components(G), key=len)
    print ('Removed',len(G.nodes) - len(largest_component), 'nodes.')
    print ( len(G.nodes), "nodes remaining.")

    G = G.subgraph(largest_component)

node_dict = lambda l : dict( zip( list(G.nodes), l ) )

cursed_dict = node_dict( [reaction.id for reaction in model.reimporactions] )
G = nx.relabel_nodes(G, cursed_dict, copy=True)

t2 = time.time(); print('\n','TIME Conversión a grafo:', (t2-t1)*1000 ,'ms')

## --- Centralidades base
t1 = time.time()
baseline_centralities = eight_centralities(G)
t2 = time.time(); print('\n','TIME Calculo de centralidades base:', (t2-t1)*1000 ,'ms')

## --- Centralidades delta (Sección paralelizada)
t1 = time.time()

INPUT_REMOVED_NODES = G.nodes # Todos los nodos
N_WORKERS = abs(int(sys.argv[2])) if sys.argv[2] else 1 # Nucleos para procesamiento paralelo

G_ray = ray.put(G) # Pasa esto al object store antes que pasarlo multiples veces

# Usa min(N_WORKERS, INPUT//2) para evitar sobre-paralelizar
split_removed_nodes = np.array_split(INPUT_REMOVED_NODES, min(N_WORKERS, INPUT_REMOVED_NODES//2) )

deltas = ray.get( 
    # SECCIÓN QUE ESTA CORRIENDO CODIGO PARALELO EN RAY
    [delta_centrality.remote(G_ray, NODES) for NODES in split_removed_nodes] 
    )

t2 = time.time(); print('\n','TIME Centralidades calculadas en paralelo:', (t2-t1)*1000 ,'ms')

## --- Centralidades perturbadas (merge)

t1 = time.time()
# Como el resultado generado es un par de output parciales, necesito una función
# que me permita pegarlos y ensamblar un output final. Por eso creo dos list
# comprehensions para separar la lista original. Es más rapido que un for loop.
perturbed_centralities = [ delta[0] for delta in deltas] 
breaks =                 [ delta[1] for delta in deltas] 

t2 = time.time(); print('\n','TIME Merge de deltas:', (t2-t1)*1000 ,'ms')

## --- Creando las tablas de salida

t1 = time.time()
# TODO: ver que de esto puedo pasar a la ejecución paralelizada
import pandas as pd

perturbed_centralities = [ i for ii in perturbed_centralities for i in ii ] # Squash N_WORKERS x NO/DO/S -> NODOS
perturbed_centralities = [ pd.DataFrame.from_dict( node ) for node in perturbed_centralities ] # Dict a Dataframes

centralities = [
    'harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
    'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality'
]

for node in perturbed_centralities: node.index = centralities # Indexa las centralidades

perturbed_centralities = [ node.T for node in perturbed_centralities ] # Transposición de dataframes, cols=centralities
perturbed_centralities = [ node.reindex( G.nodes ) for node in perturbed_centralities ] # Indexa por nodos

# En esta parte ya tenemos un set de dataframes de dimensiones (m=n, c) consistentes, incluyendo NaN


"""
n = nodos removidos por interación; m = resto de los nodos del grafo + 1 NaN; c = centralidades (8)
nota que len(n) == len(m)

A_m_c : Array 2D << Dataframe << [list {centralidades} ] # Un array de tamaño (m,c) 
B_n_m_c : Array 3D << [list  A_m_c ] # un tensor de tamaño (n, m, c) 
"""

perturbed_centralities = [ node.to_numpy() for node in perturbed_centralities ] # Convierte a Numpy
# TODO: vale la pena ahorrar memoria con .to_numpy( dtype='float32' ) ?
perturbed_centralities = perturbed_centralities.as_array()
# TODO: añadir los index al tensor u objeto creado
# TODO: VS Code desde el Browser (Clarktech) 

import pickle

outfile1 = open('./tmp/baseline_centralities','wb'); pickle.dump(baseline_centralities,outfile1); outfile1.close()
# TODO: tensorizar esto
outfile2 = open('./tmp/perturbed_centralities','wb'); pickle.dump(perturbed_centralities,outfile2); outfile2.close()
outfile3 = open('./tmp/breaks','wb'); pickle.dump(breaks,outfile3); outfile3.close()

"""
perturbed_centralities = []


for delta in perturbed_centralities:
    tmp = pd.DataFrame.from_dict( delta ) # Selecciona un grupo de 8 centralidades
    tmp = dict( tmp.mean( axis=0 ) ) # Calcula el promedio de centralidades
    perturbed_centralities.append( tmp ) # Al diccionario de centralidades perturbadas

print("Nodes removed for iterative delta centralities calculation:", len( INPUT_REMOVED_NODES ))
print("Unconected graphs generated:", len(breaks) )

perturbed_centralities = pd.DataFrame.from_dict( perturbed_centralities ) # Tabla por nodo

# TODO: resolver este re-ordenado del dataframe y exportarlo.
#cols = list(perturbed_centralities.columns); cols = [cols[-1]] + cols[:-1] # Reordena columnas
#perturbed_centralities = perturbed_centralities[cols] # Así la primera es la primera reacción

perturbed = perturbed_centralities.mean( axis=0 ) # Promedio de perturbadas

# Termina de comprimir las centralidades originales porque antes no cargo pandas

baseline = pd.DataFrame.from_dict( baseline_centralities )
baseline = baseline.transpose() # Por algun motivo queda con la otra orientación ?
baseline = baseline.mean( axis=1 )

log2_contribution = np.log2( baseline / perturbed )

tabla_resultado = pd.DataFrame({
    "ids" :         [rx.id           for rx in model.reactions],
    "formula" :     [rx.reaction     for rx in model.reactions],
    "flux" :        [rx.flux         for rx in model.reactions],
    "sensitivity" : [rx.reduced_cost for rx in model.reactions],
    "baseline_centrality":           baseline,
    "preturbed_centrality":          perturbed,
    "log2_centrality_contribution" : log2_contribution
})

OUTPUT_RESUME_TABLE = 'tmp.reactions_delta_centrality.csv'
tabla_resultado.to_csv( OUTPUT_RESUME_TABLE , index=False )
t2 = time.time(); print('\n','TIME Tabla resumen generada:', (t2-t1)*1000 ,'ms') 

# --- Empacando el grafo de gephi
t1 = time.time()
nx.set_node_attributes(G, node_dict( [rx.reaction     for rx in model.reactions] ), "reaction")
nx.set_node_attributes(G, node_dict( [rx.flux         for rx in model.reactions] ), "flux")
nx.set_node_attributes(G, node_dict( [rx.reduced_cost for rx in model.reactions] ), "sensitivity")
nx.set_node_attributes(G, node_dict( baseline ),          "baseline_centrality")
nx.set_node_attributes(G, node_dict( perturbed ),         "preturbed_centrality")
nx.set_node_attributes(G, node_dict( log2_contribution ), "log2_centrality_contribution")

OUTPUT_GRAPH =  "tmp.reactions_delta_centrality.gexf"
nx.write_gexf( G , OUTPUT_GRAPH )
t2 = time.time(); print('\n','TIME Grafo generado:', (t2-t1)*1000 ,'ms')
"""
