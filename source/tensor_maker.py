#!/bin/python
"""Convierte los objetos generados en el código anterior en un tensor para las computaciones de El Físico


"""
# %% --- Importando cosas iniciales
import pickle # Importa los elementos generados por delta_centrality.py

outfile1 = open('./tmp/baseline_centralities','rb'); baseline_centralities = pickle.load(outfile1); outfile1.close()   # Tensor n x m x c
outfile2 = open('./tmp/perturbed_centralities','rb'); perturbed_centralities = pickle.load(outfile2); outfile2.close() # Tensor 1 x m x c
outfile3 = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(outfile3); outfile3.close() # Indice de largo m == n (todos los nodos)
outfile4 = open('./tmp/breaks','rb'); breaks = pickle.load(outfile4); outfile4.close() # Lista de exceptciones <n

# Index de las centralidades calculadas
ct = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality']

import numpy as np

def ncm(tensor): 
    """Convierte un tensor de dimensiones `(n, m, c)` a uno `(n, c, m)`."""
    import numpy as np
    return  np.transpose( tensor, (0,2,1) )


# N = Nodo removido [1:n]
# M = Todos los nodos del modelo [1:m]. len(n) == len(m)
# C = centrality [1:c]. Son 8 en total

# S = subsystem [1:s]. S.nodes = subset( N )

# %% --- Cosas de tensores

# tA_1_m_c: Centralidades base como tensor de alto 1. Dims = 1 x 1000 x 8
tA_m_c = baseline_centralities

# A_m_c: Centralidades (c = 1:8) de todos los nodos (m=1:1000). Dims = 1000 x 8. // Centralidades base
A_m_c = baseline_centralities[0]

# B_n_m_c: Centralidades (c) de los nodos  (m=1, ... NA(en la posición de nodo n) ...,1000) removiendo el nodo n. Dims: 1000 x 999+(NaN) x 8
B_n_m_c = perturbed_centralities

# D_n_c_m: Es lo mismo que B_n_m_c pero con las centralidades en la segunda dimensión.
D_n_c_m = ncm( B_n_m_c )

# TODO ms: Nodos que pertenecen al subsistema s. Dims = length({ nodos E s}) x 1.

# A_ms_c: Lo mismo que A_m_c pero indexado para ms.
A_ms_c = A_m_c[1:25,:]

# B_n_ms_c: Lo mismo que B_n_m_c pero indexado para ms. // B_m_c[ms, *, *]
B_n_ms_c = B_n_m_c[:,1:25,:]

# D_n_c_ms: Lo mismo que D_n_c_ms pero indexado para ms.
D_n_c_ms = ncm( B_n_ms_c )
# TODO: abajo se contradice y dice que deberia ser shape=(ms,c,n)
D_ms_c_n = np.transpose( B_n_ms_c, axes=(1,2,0) )

# E_c: Promedio por centralidad de todos los nodes en ms. 
E_c = np.mean( A_ms_c.T, axis=1, keepdims=True ) # .shape = (c,1)

#F_c_n:Promedio por centralidad de todos los nodes en ms removiendo n. // E_c.mean(axis='n')
F_c_n = np.mean( B_n_ms_c , axis = 1) #, keepdims=True )
F_c_n = F_c_n.T # Las dimensiones correctas (c,n)

F_c_c_n = np.apply_along_axis( np.diag, 0, F_c_n.T) # Tensor con centralidades diagonales x nodos (c,c,n)
F_inv_c_c_n = np.linalg.inv(F_c_c_n.T).T # Calcula el inverso de las n matrices (c,c).    shape = (c,c,n)

Ratios = F_inv_c_c_n * E_c # .shape = (c,c,n) * (1,c)

# %% --- Tensores extraños
"""## Operaciones:

Para los nodos del conjunto ms1, tal que  ms1 pertenece el subsistema s1.

l_ms1 = length(ms1) : el número de nodos en el subsistema s1. // len(ms1)

- A_ms1_c (l_ms1 x 8) => calcular el promedio de cada columna ->  (1 x 8) -> transponer => E_c (8 x 1) // by_columns.means()

# TODO ?? D_ms1_c_n (l_ms1 x 8 x n, donde n = son los nodos removidos = 1000)  => calcular el promedio de cada columna (c = [1:8]) => 

  G( 1 x 8 x n) => colapsar la primera dimensión  => matriz( 8 x n) = (c x n) => F_c_n // cosa = cosa[0] 
  # TODO: G es una matriz intermedia innecesaria

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

"""
# %%

Ratios = baseline_centralities / perturbed_centralities
Ratios = np.log2( Ratios )