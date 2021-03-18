# %%
from   heapq import merge
import cobra  
import numpy  as np
import pandas as pd
from   cobra.util.array import create_stoichiometric_matrix
import warnings 
warnings.filterwarnings("ignore")
import networkx as nx

# creacion el grafo en networkx
#
#importar modelo cobra que està en formato json jason
stimulated = cobra.io.load_json_model(
    "/home/alejandro/PostDoc/human-metnet/data/stimulated_2021.json")
#extraer matriz estequiométrica
S_matrix = create_stoichiometric_matrix(stimulated)
#convertir todas las entradas valores positivos
S_matrix = (abs(S_matrix) )
#transformar a enteros 
S_matrix = S_matrix.astype(np.int)
#Multiplicacion por la derecha para proyectar en el espacio de las reacciones
projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
print(projected_S_matrix.shape, np.count_nonzero(projected_S_matrix))
np.fill_diagonal(projected_S_matrix, 0) 
#binarizar
from sklearn.preprocessing import Binarizer
projected_S_matrix = Binarizer().fit_transform(projected_S_matrix)
#chequear
print(projected_S_matrix.shape, np.count_nonzero(projected_S_matrix), np.unique(projected_S_matrix))
#crear grafo networkx
G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
#hacer diccionario con los nombres de las reacciones
node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
cursed_dict = node_dict( [reaction.id for reaction in stimulated.reactions] )
#chequear connectividad, los tamaños y que los nombres aun NO existen
print(nx.is_connected(G) , len(cursed_dict), len(G.nodes), 'DPGM_Neuron' in list(G.nodes) )

# %%
#Renombrar los nodos usando el diccionario
G = nx.relabel_nodes(G, cursed_dict, copy=True) # Revisar que esto este antes de la remoción del grafo

largest_component = max(nx.connected_components(G), key=len)
G = G.subgraph(largest_component)

print(nx.is_connected(G) , len(cursed_dict), len(G.nodes), 'DPGM_Neuron' in list(G.nodes) )

Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']
Glycolysis_neuron = ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron',
 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron']

ETC_neuron    = ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron']
ETC_astrocyte =  ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']
# %% Perturbed centralities for TYRTAm

# Calculated Centralities Removed
G_with_a_removal = G.copy() # Sus
G_with_a_removal.remove_node('TYRTAm')
# Generate connected components and select the largest:
largest_component       = max(nx.connected_components(G_with_a_removal), key=len)
# Create a subgraph of G consisting only of this component:
graph_largest_component   = G_with_a_removal.subgraph(largest_component)

###### compute centralities ##########
G_with_a_removal  = graph_largest_component.copy()

grafo = G_with_a_removal.copy()

print(len(G.nodes), len(grafo.nodes))
# %% 

hc  = nx.harmonic_centrality(grafo, nbunch=None, distance=None)
ec  = nx.eigenvector_centrality(grafo, max_iter=1000, tol=1e-05, nstart=None, weight=None)
dc  = nx.degree_centrality(grafo)
bc  = nx.betweenness_centrality(grafo, normalized=True, weight=None, endpoints=False, seed=None)
"""↓ Here be dragons ↓"""
cc  = nx.closeness_centrality(grafo, distance=None, wf_improved=True)
lc  = nx.load_centrality(grafo, cutoff=None, normalized=True, weight=None)
""" Requieren un grafo full conected """
#if not (nx.is_connected(grafo)) : 
#    ic = {}; soc = {}
#else:
#ic  = nx.information_centrality(grafo)
#soc = nx.second_order_centrality(grafo)
ic  = nx.information_centrality(grafo)
soc = nx.second_order_centrality(grafo)

cfcc    = nx.current_flow_closeness_centrality(grafo)    
cfbc    = nx.current_flow_betweenness_centrality(grafo)
acfbc   = nx.approximate_current_flow_betweenness_centrality(grafo)
cbc     = nx.communicability_betweenness_centrality(grafo)   

#####-- make DFs --###################################################################

cfcc_df = pd.DataFrame.from_dict(cfcc, columns = ['cfcc'], orient='index')
cfbc_df = pd.DataFrame.from_dict(cfbc, columns = ['cfbc'], orient='index')
acfbc_df = pd.DataFrame.from_dict(acfbc, columns = ['acfbc'], orient='index')
cbc_df = pd.DataFrame.from_dict(cbc, columns = ['cbc'], orient='index')



hc_df = pd.DataFrame.from_dict(hc, columns = ['hc'], orient='index')
ec_df = pd.DataFrame.from_dict(ec, columns = ['ec'], orient='index')
dc_df = pd.DataFrame.from_dict(dc, columns = ['dc'], orient='index')
bc_df = pd.DataFrame.from_dict(bc, columns = ['bc'], orient='index')
cc_df = pd.DataFrame.from_dict(cc, columns = ['cc'], orient='index')
lc_df = pd.DataFrame.from_dict(lc, columns = ['lc'], orient='index')
ic_df   = pd.DataFrame.from_dict(ic, columns = ['ic'], orient='index')
soc_df  = pd.DataFrame.from_dict(soc, columns = ['soc'], orient='index')   

from functools import reduce
dfs = [hc_df, ec_df, dc_df, bc_df, cc_df, lc_df, ic_df, soc_df, cfcc_df, cfbc_df, acfbc_df,cbc_df]
df_all_perturbed_centralities = reduce(lambda  left, right: left.join(right, how='outer'), dfs)

# %%

numeric_all_perturbed_centralities = np.array(df_all_perturbed_centralities)
print(
np.count_nonzero(np.isnan(numeric_all_perturbed_centralities)))
df_all_perturbed_centralities_from_HPC \
     = pd.read_csv( "/home/alejandro/PostDoc/human-metnet/source/hateful-eight/TYRTAm_df.csv",index_col=0)

numeric_all_perturbed_centralities_from_HPC = np.array(df_all_perturbed_centralities_from_HPC)
print(
np.count_nonzero(np.isnan(numeric_all_perturbed_centralities_from_HPC)))
# %%
#Promedio aritmetico MANU
# Un helper porque estas funciones no pueden lidiar con NaNs porque son superiores a los NaN o algo así
noNaN_perturbed_subsystem = np.nan_to_num( numeric_all_perturbed_centralities , copy=True, nan=0.0, posinf=None, neginf=None)
# Es el unico que no require lambdas raros :D
aritmetic_perturbed = np.nanmean( noNaN_perturbed_subsystem, axis= 1)
print(
aritmetic_perturbed)

noNaN_perturbed_subsystem = np.nan_to_num( numeric_all_perturbed_centralities_from_HPC , copy=True, nan=0.0, posinf=None, neginf=None)
# Es el unico que no require lambdas raros :D
aritmetic_perturbed = np.nanmean( noNaN_perturbed_subsystem, axis= 0)


print(aritmetic_perturbed.shape,
aritmetic_perturbed)

pd.DataFrame(aritmetic_perturbed)

# %%
#Promedio aritmetico NAIVE

rxn_to_remove = 'TYRTAm'

pd.DataFrame(df_all_perturbed_centralities.mean(axis=0), columns={rxn_to_remove})




# %%
df_all_perturbed_centralities.loc[Glycolysis_astrocyte,:].mean(axis=0)
# %%
Glycolysis_astrocyte_idxs = [417, 139, 414, 415, 445, 324, 242, 479, 412, 252, 281, 241, 258, 278]
df_all_perturbed_centralities_from_HPC.loc[Glycolysis_astrocyte_idxs,:].mean(axis=0)
# %%
