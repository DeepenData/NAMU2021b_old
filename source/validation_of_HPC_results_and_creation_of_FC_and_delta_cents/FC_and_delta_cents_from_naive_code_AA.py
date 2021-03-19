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
G0 = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
#hacer diccionario con los nombres de las reacciones
node_dict   = lambda l : dict( zip( list(G0.nodes), l ) )
cursed_dict = node_dict( [reaction.id for reaction in stimulated.reactions] )
#chequear connectividad, los tamaños y que los nombres aun NO existen
print(nx.is_connected(G0) , len(cursed_dict), len(G0.nodes), 'DPGM_Neuron' in list(G0.nodes) )
#Renombrar los nodos usando el diccionario
G0 = nx.relabel_nodes(G0, cursed_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
largest_component = max(nx.connected_components(G0), key=len)


G = G0.subgraph(largest_component)
print(nx.is_connected(G) , len(cursed_dict), len(G.nodes), 'DPGM_Neuron' in list(G.nodes) )


# %% --- Extrae subsistemas con una lista

Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']
Glycolysis_neuron = ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron',
 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron']
ETC_neuron    = ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron']
ETC_astrocyte =  ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']

# %%
def compute_centralities(grafo):
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

def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G

import pickle 

infile      = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(infile); infile.close()
index_nodes = list(index_nodes)
Glycolysis_astrocyte_idxs =  [ index_nodes.index(node) for node in Glycolysis_astrocyte ]



def compute_centralities_short(graph):
    """Computa las centralidades rapidas y devuelve una DataFrame (no reindexado) con estas"""
    # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
    hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc    = nx.degree_centrality(graph)
    cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    ic    = nx.information_centrality(graph) # Requiere scipy

    # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
    centralities = {
        'harmonic_centrality' : hc ,
        'eigenvector_centrality' : ec ,
        'degree_centrality' : dc ,
        'closeness_centrality' : cc ,
        'information_centrality' : ic
    }

    # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
    centralities = pd.DataFrame( centralities )

    return centralities



# %% baseline

G_unperturbed    = G.copy()
df_all_baseline_centralities = compute_centralities(G_unperturbed)

baseline_Glycolysis_astrocyte_df = df_all_baseline_centralities.loc[Glycolysis_astrocyte]


# %% --- Resultados naive con dos reacicones
baseline_Glycolysis_astrocyte_df.to_csv(
    "/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/Naive_baseline_Glycolysis_astrocyte.csv")

# %%
G_with_a_removal = G.copy() 
G_with_a_removal.remove_node('FPGS3m')
G_with_a_removal = get_largest_component(G_with_a_removal)
perturbed_naive_L_LACt2r = compute_centralities_short(G_with_a_removal)
# %%


perturbed_naive_FPGS3m_short = perturbed_naive_L_LACt2r


perturbed_naive_FPGS3m_short.to_csv('perturbed_naive_FPGS3m_short.csv')
# %%

len(set(G.nodes) - set(G_with_a_removal.nodes))

perturbed_naive_FPGS3m_short.shape
# %%
#L_LACt2r     =  [ index_nodes.index(node) for node in ['L-LACt2r'] ]
#perturbed_naive_L_LACt2r_df = \
perturbed_naive__Glycolysis_astrocyte_L_LACt2r_df = \
perturbed_naive_L_LACt2r.loc[Glycolysis_astrocyte]
perturbed_naive__Glycolysis_astrocyte_L_LACt2r_df
# %% FC y delta cents
def aritmetic_mean(df):
    return np.nanmean( df, axis= 0)

aritmetic_mean_perturbed_naive__Glycolysis_astrocyte_L_LACt2r = \
    aritmetic_mean(perturbed_naive__Glycolysis_astrocyte_L_LACt2r_df)


# %% 

path = \
'/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_perturbed_L_LACt2r_df.csv'

left_nodes =  list(set(G.nodes) - set(G_with_a_removal.nodes))
left_nodes_idxs =  [ index_nodes.index(node) for node in left_nodes ]
hpc_perturbed_L_LACt2r_df = pd.read_csv(path, index_col=0)
hpc_perturbed_L_LACt2r_df.drop(left_nodes_idxs, axis=0, inplace = True)

100*(np.linalg.norm(perturbed_naive_L_LACt2r)/ \
np.linalg.norm(hpc_perturbed_L_LACt2r_df))

# %% 


path = \
'/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_perturbed_Glycolysis_astrocyte_L_LACt2r_df.csv'

hpc_perturbed_Glycolysis_astrocyte_L_LACt2r_df = pd.read_csv(path, index_col=0)
#hpc_perturbed_Glycolysis_astrocyte_L_LACt2r_df.drop(left_nodes_idxs, axis=0, inplace = True)

100*(np.linalg.norm(perturbed_naive__Glycolysis_astrocyte_L_LACt2r_df)/ \
np.linalg.norm(hpc_perturbed_Glycolysis_astrocyte_L_LACt2r_df))

[index_nodes[index] for index in hpc_perturbed_Glycolysis_astrocyte_L_LACt2r_df.index]

# %% 

path = \
'/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_baseline_Glycolysis_astrocyte_df.csv'

hpc_baseline_Glycolysis_astrocyte_df = pd.read_csv(path, index_col=0)

#Glycolysis_df_all_baseline_centralities
hpc_baseline_Glycolysis_astrocyte_df
# %% 
from functools import reduce
#df_perturbed_centralities  = reduce(lambda  left, right: left.join(right, how='outer'), [UPP3S_Neuron, TYRTAm])
#df_perturbed_centralities Get baseline centralities 

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
'''
#en el caso de que ninguno cambie todos son 0, por tanto el promedio aritmético también

#FC_UPP3S_Neuron =  np.log2(df_all_subsystems_baseline.baseline/df_perturbed_centralities.UPP3S_Neuron)
#FC_TYRTAm       = np.log2(df_all_subsystems_baseline.baseline/df_perturbed_centralities.TYRTAm)

#data  = {'FC_TYRTAm': FC_TYRTAm, 'FC_UPP3S_Neuron': FC_UPP3S_Neuron}
final_FC = pd.DataFrame(data)

final_FC # Tiene nombres de las columnas FC_{{rxn}}_{{cell}}
# TODO: transponer esto en una matriz larga en lugar de matriz ancha
# TODO: deberia coincidir con los valores que ya hemos procesado'''

# %% --- Histograma(s)

# TODO: hacer esto con distitnas tecnicas de promedio
# TODO: extraer valores de los dataframe para hacer un histograma
#np.array(final.values).flatten() 
# %>% histograma() 
# %>% plot()
# %% Calculating delta centralities
#en el caso de que ninguno cambie todos son 0, por tanto el promedio aritmético también
'''
def get_delta_centrality(unperturbed, perturbed):
    delta = (unperturbed-perturbed)/unperturbed
    return delta

#delta_UPP3S_Neuron = get_delta_centrality(df_all_subsystems_baseline.baseline, df_perturbed_centralities.UPP3S_Neuron)
#delta_TYRTAm       = get_delta_centrality(df_all_subsystems_baseline.baseline, df_perturbed_centralities.TYRTAm)


data  = {'delta_TYRTAm': delta_TYRTAm, 'delta_UPP3S_Neuron': delta_UPP3S_Neuron}
final_delta = pd.DataFrame(data)'''
# %%

'''final_FC

'''
# %%
