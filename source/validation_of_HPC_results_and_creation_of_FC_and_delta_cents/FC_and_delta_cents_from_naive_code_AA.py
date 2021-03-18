# naive_FC_centrality
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
#Renombrar los nodos usando el diccionario
G = nx.relabel_nodes(G, cursed_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
largest_component = max(nx.connected_components(G), key=len)
G = G.subgraph(largest_component)
print(nx.is_connected(G) , len(cursed_dict), len(G.nodes), 'DPGM_Neuron' in list(G.nodes) )


# %% --- Extrae subsistemas con una lista

Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']
Glycolysis_neuron = ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron',
 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron']
ETC_neuron    = ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron']
ETC_astrocyte =  ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']

# %%
def compute_centralities(graph):
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
    df_all_centralities = reduce(lambda  left, right: left.join(right, how='outer'), dfs)
    # Lambda que reduce dataframes y hace un join
    return df_all_centralities


def calc_centr_rmv(rxn_to_remove):
    # Calculated Centralities Removed
    G_with_a_removal = G.copy() # Sus
    G_with_a_removal.remove_node(rxn_to_remove)

    # Generate connected components and select the largest:
    largest_component       = max(nx.connected_components(G_with_a_removal), key=len)
    # Create a subgraph of G consisting only of this component:
    graph_largest_component   = G_with_a_removal.subgraph(largest_component)
    # Hacer una extracción del componente más grande; calcular centralidades con los otros dos que antes se ignoraban
    # Tirar lista que indica cuantos componentes tiene el grafo; idealmente [1079, 1]
    #nodes_in_small_components = list(set(G_with_a_removal.nodes) - set(graph_largest_component.nodes))
    

    ###### compute centralities ##########
    graph  = graph_largest_component.copy()
    df_all_perturbed_centralities = compute_centralities(graph)

    # Función que extrae solo los subsistemas para hacer un mean más corto
    def get_centrality_mean_by_subsys(a_subsystem):
        # Linea importante que saca los promedios. 
        # TODO: Usar otras formas de promedio. 
        df = pd.DataFrame(df_all_perturbed_centralities.loc[eval(a_subsystem)].mean(axis=0), columns={rxn_to_remove})
        new_row_names = [s + str('_'+a_subsystem) for s in list(df.index.values)]
        df.index = new_row_names
        return df

    df_1 = get_centrality_mean_by_subsys('ETC_astrocyte')
    df_2 = get_centrality_mean_by_subsys('Glycolysis_astrocyte')
    df_3 = get_centrality_mean_by_subsys('ETC_neuron')
    df_4 = get_centrality_mean_by_subsys('Glycolysis_neuron')

    from functools import reduce
    dfs = [df_1, df_2, df_3, df_4]
    df_all_subsystems = reduce(lambda  left, right: left.append(right, ignore_index=False), dfs)
    return df_all_subsystems

# %% --- Resultados naive con dos reacicones
L_LACt2r = calc_centr_rmv('L-LACt2r')

# Dataframes distintas que se juntan y hacen cosas raras
# %%

# Genera un dataframe de una columna/vector con nombre de columna de la reacción
# Los rownames corresponden a la centralidad más el subsistema más la célula
# más abajo estos dataframes/vectores se ensamblan en uno con una columna por reacción

# %% 
from functools import reduce
df_perturbed_centralities  = reduce(lambda  left, right: left.join(right, how='outer'), [UPP3S_Neuron, TYRTAm])
#df_perturbed_centralities
# %% Get baseline centralities 

# Desde aqui obtenemos las centralidades de linea base
# esto hace ¿un dataframe... que hace?

G_unperturbed    = G.copy()
df_all_baseline_centralities = compute_centralities(G_unperturbed)

def get_centrality_mean_by_subsys(a_subsystem):
    df = pd.DataFrame(df_all_baseline_centralities.loc[eval(a_subsystem)].mean(axis=0), columns={'baseline'})
    new_row_names = [s + str('_'+a_subsystem) for s in list(df.index.values)]
    df.index = new_row_names
    return df # 

df_1 = get_centrality_mean_by_subsys('ETC_astrocyte')
df_2 = get_centrality_mean_by_subsys('Glycolysis_astrocyte')
df_3 = get_centrality_mean_by_subsys('ETC_neuron')
df_4 = get_centrality_mean_by_subsys('Glycolysis_neuron')

from functools import reduce
dfs = [df_1, df_2, df_3, df_4]
df_all_subsystems_baseline = reduce(lambda  left, right: left.append(right, ignore_index=False), dfs)

df_all_subsystems_baseline
# Esto genera un dataframe con 
# %%
df_all_subsystems_baseline
# %% Get fold changes

#en el caso de que ninguno cambie todos son 0, por tanto el promedio aritmético también

FC_UPP3S_Neuron =  np.log2(df_all_subsystems_baseline.baseline/df_perturbed_centralities.UPP3S_Neuron)
FC_TYRTAm       = np.log2(df_all_subsystems_baseline.baseline/df_perturbed_centralities.TYRTAm)

data  = {'FC_TYRTAm': FC_TYRTAm, 'FC_UPP3S_Neuron': FC_UPP3S_Neuron}
final_FC = pd.DataFrame(data)

final_FC # Tiene nombres de las columnas FC_{{rxn}}_{{cell}}
# TODO: transponer esto en una matriz larga en lugar de matriz ancha
# TODO: deberia coincidir con los valores que ya hemos procesado

# %% --- Histograma(s)

# TODO: hacer esto con distitnas tecnicas de promedio
# TODO: extraer valores de los dataframe para hacer un histograma
#np.array(final.values).flatten() 
# %>% histograma() 
# %>% plot()
# %% Calculating delta centralities
#en el caso de que ninguno cambie todos son 0, por tanto el promedio aritmético también

def get_delta_centrality(unperturbed, perturbed):
    delta = (unperturbed-perturbed)/unperturbed
    return delta

delta_UPP3S_Neuron = get_delta_centrality(df_all_subsystems_baseline.baseline, df_perturbed_centralities.UPP3S_Neuron)
delta_TYRTAm       = get_delta_centrality(df_all_subsystems_baseline.baseline, df_perturbed_centralities.TYRTAm)


data  = {'delta_TYRTAm': delta_TYRTAm, 'delta_UPP3S_Neuron': delta_UPP3S_Neuron}
final_delta = pd.DataFrame(data)
# %%

final_FC


# %%
