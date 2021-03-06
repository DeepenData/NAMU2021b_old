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

stimulated = cobra.io.load_json_model(
    "/home/alejandro/PostDoc/Metabolic_diseases_lab/proyecto-redes-metabolicas/data/stimulated_2021.json")


S_matrix = create_stoichiometric_matrix(stimulated)
S_matrix = (abs(S_matrix) )
S_matrix = (S_matrix > 0.0).astype(np.int_)

projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
np.fill_diagonal(projected_S_matrix, 0) 

reaction_adjacency_matrix = (projected_S_matrix !=0).astype(int)

G = nx.convert_matrix.from_numpy_matrix( reaction_adjacency_matrix )
node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
cursed_dict = node_dict( [reaction.id for reaction in stimulated.reactions] )
print(len(cursed_dict), len(G.nodes))

# %%
G = nx.relabel_nodes(G, cursed_dict, copy=True) # Revisar que esto este antes de la remoción del grafo

largest_component = max(nx.connected_components(G), key=len)
G = G.subgraph(largest_component)
print(len(cursed_dict), len(G.nodes))

nx.is_connected(G) # Reemplazar por llamadas de función
# %% --- Extrae subsistemas
import re

graph_subsystems = [stimulated.reactions.get_by_id(a_node_name).subsystem for a_node_name in list(G.nodes)]

data = {'node_name': list(G.nodes), 'subsystem': graph_subsystems}
df   = pd.DataFrame(data)

Glycolysis = df[df['subsystem'].str.contains("Glycolysis")].node_name
ETC        = df[df['subsystem'].str.contains("Oxidative")].node_name


# %% --- Extrae subsistemas con una lista
Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']

Glycolysis_neuron = ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron',
 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron', 'PEPCK_Neuron']

subs="Neuron"
#Glycolysis_neuron = [x for x in Glycolysis if re.search(subs, x)] 
ETC_neuron           = [x for x in ETC if re.search(subs, x)] 
#Glycolysis_astrocyte = list(set(Glycolysis) - set(Glycolysis_neuron))
ETC_astrocyte        = list(set(ETC) - set(ETC_neuron))

# Separar subsistemas por astro/neurona




# %%
def compute_centralities(graph):
    hc = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    dc = nx.degree_centrality(graph)
    bc = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
    cc = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    lc = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
    # 
#   cfcc    = current_flow_closeness_centrality(graph)
    ic      = nx.information_centrality(graph)
#   cfbc    = current_flow_betweenness_centrality(graph)
#   cbc     = communicability_betweenness_centrality(graph)
    soc     = nx.second_order_centrality(graph) 
    #####-- make DFs --###################################################################3
    hc_df = pd.DataFrame.from_dict(hc, columns = ['hc'], orient='index')
    ec_df = pd.DataFrame.from_dict(ec, columns = ['ec'], orient='index')
    dc_df = pd.DataFrame.from_dict(dc, columns = ['dc'], orient='index')
    bc_df = pd.DataFrame.from_dict(bc, columns = ['bc'], orient='index')
    cc_df = pd.DataFrame.from_dict(cc, columns = ['cc'], orient='index')
    lc_df = pd.DataFrame.from_dict(lc, columns = ['lc'], orient='index')
    # 5 new
#   cfcc_df = pd.DataFrame.from_dict(cfcc, columns = ['cfcc'], orient='index')
    ic_df   = pd.DataFrame.from_dict(ic, columns = ['ic'], orient='index')
#   cfbc_df = pd.DataFrame.from_dict(cfbc, columns = ['cfbc'], orient='index')
#   cbc_df  = pd.DataFrame.from_dict(cbc, columns = ['cbc'], orient='index')
    soc_df  = pd.DataFrame.from_dict(soc, columns = ['soc'], orient='index')   
    from functools import reduce
    dfs = [hc_df, ec_df, dc_df, bc_df, cc_df, lc_df, ic_df, soc_df]
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
    nodes_in_small_components = list(set(G_with_a_removal.nodes) - set(graph_largest_component.nodes))
    

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
UPP3S_Neuron = calc_centr_rmv('UPP3S_Neuron')
TYRTAm       = calc_centr_rmv('TYRTAm')
# Dataframes distintas que se juntan y hacen cosas raras
# %%
UPP3S_Neuron
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
# Esto genera un dataframe con 
# %%
df_all_subsystems_baseline
# %% Get fold changes
FC_UPP3S_Neuron =  np.log2(df_all_subsystems_baseline.baseline/df_perturbed_centralities.UPP3S_Neuron)
FC_TYRTAm       = np.log2(df_all_subsystems_baseline.baseline/df_perturbed_centralities.TYRTAm)

data  = {'FC_TYRTAm': FC_TYRTAm, 'FC_UPP3S_Neuron': FC_UPP3S_Neuron}
final = pd.DataFrame(data)

final # Tiene nombres de las columnas FC_{{rxn}}_{{cell}}
# TODO: transponer esto en una matriz larga en lugar de matriz ancha
# TODO: deberia coincidir con los valores que ya hemos procesado

# %% --- Histograma(s)

# TODO: hacer esto con distitnas tecnicas de promedio
# TODO: extraer valores de los dataframe para hacer un histograma
np.array(final.values).flatten() 
# %>% histograma() 
# %>% plot()
