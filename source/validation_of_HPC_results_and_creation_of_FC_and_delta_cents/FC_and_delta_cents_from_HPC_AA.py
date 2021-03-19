#!/bin/python
"""
"An Ixian machine? You defy the Jihad!" 
"There's a lesson in that, too. What do such machines really do? They increase the number of things we can do without thinking. Things we do without thinking — there's the real danger."
    - God-Emperor Leto II 
"""
# %% --- Importar resultados
from networkx.algorithms.isomorphism.ismags import intersect
import pandas   as pd
import numpy    as np
import networkx as nx
import pickle 
# Esto asume que está en 'human-metnet/'
#%cd ..
#%cd ..
infile = open('./tmp/perturbed_centralities','rb'); perturbed_centralities = pickle.load(infile); infile.close()
infile = open('./tmp/baseline_centralities' ,'rb'); baseline_centralities  = pickle.load(infile); infile.close()
infile = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(infile); infile.close()
infile = open('./tmp/subsystems_dict','rb'); subsystems = pickle.load(infile); infile.close() 

index_nodes = list(index_nodes)

Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']
Glycolysis_neuron = ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron',
 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron']
ETC_neuron    = ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron']
ETC_astrocyte =  ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']

index_centralities = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality',
      'current_flow_closeness_centrality', 'current_flow_betweenness_centrality',
       'approximate_current_flow_betweenness_centrality', 'communicability_betweenness_centrality']
# %% 
def cobra_to_networkx_rxn_projection(modelo):
    import networkx as nx
    from   cobra.util.array import create_stoichiometric_matrix
    import numpy as np
    from sklearn.preprocessing import Binarizer
    import warnings 
    warnings.filterwarnings("ignore")

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"
    #extraer matriz estequiométrica
    S_matrix = create_stoichiometric_matrix(modelo)
    #convertir todas las entradas valores positivos
    S_matrix = (abs(S_matrix) )
    #transformar a enteros 
    S_matrix = S_matrix.astype(np.int)
    #Multiplicacion por la derecha para proyectar en el espacio de las reacciones
    projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
    #rellenar diagonal con ceros
    np.fill_diagonal(projected_S_matrix, 0) 
    #binarizar
    projected_S_matrix = Binarizer().fit_transform(projected_S_matrix)
    #crear grafo networkx
    G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
    #hacer diccionario con los nombres de las reacciones
    node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
    reaction_dict = node_dict( [reaction.id for reaction in modelo.reactions] )
    #Renombrar los nodos usando el diccionario
    G = nx.relabel_nodes(G, reaction_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
    return G


def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G
# %%  generar grafo
import cobra  
import warnings 
warnings.filterwarnings("ignore")

stimulated = cobra.io.load_json_model(
    "/home/alejandro/PostDoc/human-metnet/data/stimulated_2021.json")
    
G0 = cobra_to_networkx_rxn_projection(stimulated)
G  = get_largest_component(G0)

# %% Encontrar todos los nodos que son 'bridges' del componente más grande del grafo
#nx.has_bridges(G) 

bridges = set(np.array(list(nx.bridges(G))).flatten())
degrees =  dict(nx.degree(G))
degrees_one = {key:value for (key, value) in degrees.items() if value == 1}
degrees_one = set(list(degrees_one.keys()))
disjoining_nodes =  list(set.difference(bridges, degrees_one))
#check 
G_with_a_removal = G.copy()
G_with_a_removal.remove_node(disjoining_nodes[0])
print(disjoining_nodes[0],nx.is_connected(G_with_a_removal))




# %% Baseline

Glycolysis_astrocyte_idxs          =  [ index_nodes.index(node) for node in Glycolysis_astrocyte ]


reshaped_tensor = baseline_centralities[ :,Glycolysis_astrocyte_idxs,:].reshape(14,12)

baseline_Glycolysis_astrocyte_df = pd.DataFrame(reshaped_tensor, index= Glycolysis_astrocyte, columns= index_centralities)
baseline_Glycolysis_astrocyte_df.to_csv(
    "/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_baseline_Glycolysis_astrocyte_df.csv")

baseline_Glycolysis_astrocyte_df
# %%
#perturbed_hpc_Glycolysis_astrocyte_L_LACt2r


L_LACt2r     =  [ index_nodes.index(node) for node in ['L-LACt2r'] ]

perturbed_hpc_L_LACt2r_df  =  pd.DataFrame(perturbed_centralities[L_LACt2r,: ,:].reshape(1056,12))

perturbed_hpc_L_LACt2r_df.to_csv(
    "/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_perturbed_L_LACt2r_df.csv")


perturbed_hpc_Glycolysis_astrocyte_L_LACt2r_df = perturbed_hpc_L_LACt2r_df.loc[Glycolysis_astrocyte_idxs,:]


perturbed_hpc_Glycolysis_astrocyte_L_LACt2r_df.to_csv(
    "/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_perturbed_Glycolysis_astrocyte_L_LACt2r_df.csv")
# %%
# %%





# %% 

def aritmetic_mean(df):
    return np.nanmean( df, axis= 0)






hpc_aritmetic_mean_perturbed_Glycolysis_astrocyte_L_LACt2r = pd.DataFrame(aritmetic_mean(perturbed_hpc_Glycolysis_astrocyte_L_LACt2r_df))
hpc_aritmetic_mean_perturbed_Glycolysis_astrocyte_L_LACt2r.to_csv(
    "/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_aritmetic_mean_perturbed_Glycolysis_astrocyte_L_LACt2r.csv")


hpc_aritmetic_mean_baseline_Glycolysis_astrocyte = pd.DataFrame(aritmetic_mean(baseline_Glycolysis_astrocyte_df))
hpc_aritmetic_mean_baseline_Glycolysis_astrocyte.to_csv(
    "/home/alejandro/PostDoc/human-metnet/source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/hpc_aritmetic_mean_baseline_Glycolysis_astrocyte.csv")


# %% FC y delta cents



ratio_from_arit_mean_L_LACt2r  = aritmetic_mean(baseline_Glycolysis_astrocyte_df)/aritmetic_mean(perturbed_hpc_Glycolysis_astrocyte_L_LACt2r_df)
FC_from__arit_mean_L_LACt2r = np.log2(ratio_from_arit_mean_L_LACt2r)
delta_from_arit_mean_L_LACt2r  = (aritmetic_mean(baseline_Glycolysis_astrocyte_df) - aritmetic_mean(perturbed_hpc_Glycolysis_astrocyte_L_LACt2r_df))/aritmetic_mean(baseline_Glycolysis_astrocyte_df)


# %%


# %% --- Cosas de los subsistemas ???

# TODO: reemplazar esto por una importación desde Tensor Maker
def subsystem_tensor( subsystem_nodes, tensor):
    """Genera un tensor con un subset de nodos basados en el subsistema pedido"""
    subsystem_index = [ index_nodes.index(node) for node in subsystem_nodes ]
    return tensor[:, subsystem_index ,:]


perturbed_Glycolysis_astrocyte = subsystem_tensor( Glycolysis_astrocyte, perturbed_centralities )
baseline_Glycolysis_astrocyte  = subsystem_tensor( Glycolysis_astrocyte, baseline_centralities   )

perturbed_Glycolysis_neuron = subsystem_tensor( Glycolysis_neuron, perturbed_centralities )
baseline_Glycolysis_neuron  = subsystem_tensor( Glycolysis_neuron, baseline_centralities   )

perturbed_ETC_neuron = subsystem_tensor( ETC_neuron, perturbed_centralities )
baseline_ETC_neuron  = subsystem_tensor( ETC_neuron, baseline_centralities   )

perturbed_ETC_astrocyte= subsystem_tensor( ETC_astrocyte, perturbed_centralities )
baseline_ETC_astrocyte  = subsystem_tensor( ETC_astrocyte, baseline_centralities   )



#np.concatenate([perturbed_ETC_astrocyte, perturbed_Glycolysis_astrocyte], axis = 0).shape

# %% --- Colapsando la segunda dimensión... (nodos totales)
# Esto calcula 4 promedios de las centralidades para los nodos

def get_8_aggregations(baseline_subsystem, perturbed_subsystem):

    from numba import jit

    # Un helper porque estas funciones no pueden lidiar con NaNs porque son superiores a los NaN o algo así
    noNaN_perturbed_subsystem = np.nan_to_num( perturbed_subsystem , copy=True, nan=0.0, posinf=None, neginf=None)

    # Es el unico que no require lambdas raros :D
    aritmetic_perturbed = np.nanmean( noNaN_perturbed_subsystem, axis= 1)
    aritmetic_baseline  = np.nanmean(        baseline_subsystem, axis= 1)

    # Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
    @jit(nopython=True)
    def geometric_mean(array):
        array = array[~np.isnan(array)] # Filtra los NaN
        fcorr = len(array[array > 0])/len(array) # Factor de corrección
        array = array[array > 0] # Selecciona solo mayores a cero
        try: gmean = ( np.prod( array )**(1 / len(array)) ) * fcorr # Geometrica por factor de corrección
        except: return np.nan
        else : return gmean

    geometric_perturbed = np.apply_along_axis( geometric_mean, 1, noNaN_perturbed_subsystem)
    geometric_baseline  = np.apply_along_axis( geometric_mean, 1,        baseline_subsystem)

    # Define una función helper para cuadraticos (Root Square Mean). No en scipy base
    @jit(nopython=True)
    def quadratic_mean(array):
        mean =  np.sqrt( np.nanmean( array*array ) )
        return mean

    quadratic_perturbed = np.apply_along_axis( quadratic_mean, 1, noNaN_perturbed_subsystem)
    quadratic_baseline  = np.apply_along_axis( quadratic_mean, 1,        baseline_subsystem)

    # Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
    @jit(nopython=True)
    def harmonic_mean(array):
        array = array[~np.isnan(array)] # Filtra los NaN
        fcorr = len(array[array > 0])/len(array) # Factor de corrección
        array = array[array > 0] # Selecciona solo mayores a cero
        try : hmean = len(array)/np.sum(1 / array ) * fcorr # Geometrica por factor de corrección
        except: return np.nan
        else : return hmean

    harmonic_perturbed  = np.apply_along_axis( harmonic_mean,  1, noNaN_perturbed_subsystem)
    harmonic_baseline   = np.apply_along_axis( harmonic_mean,  1,        baseline_subsystem)
    ####################################
    fold_change_aritmetic = np.log2( aritmetic_baseline / aritmetic_perturbed )
    fold_change_geometric = np.log2( geometric_baseline / geometric_perturbed )
    fold_change_quadratic = np.log2( quadratic_baseline / quadratic_perturbed )
    fold_change_harmonic  = np.log2( harmonic_baseline  / harmonic_perturbed  )

    delta_aritmetic = ( aritmetic_baseline - aritmetic_perturbed ) / aritmetic_baseline
    delta_geometric = ( geometric_baseline - geometric_perturbed ) / geometric_baseline
    delta_quadratic = ( quadratic_baseline - quadratic_perturbed ) / quadratic_baseline
    delta_harmonic  = ( harmonic_baseline  - harmonic_perturbed  ) / harmonic_baseline 

    return fold_change_aritmetic, fold_change_geometric, fold_change_quadratic, fold_change_harmonic, delta_aritmetic, delta_geometric, delta_quadratic, delta_harmonic

# %% --- Divisiones! Calculando como un nodo afecta conectividad de otro
# Calcula las metricas para definir los efectos en los nodos





aggregation = [ 
    'fold_change_aritmetic' , 'fold_change_geometric' , 
    'fold_change_quadratic' , 'fold_change_harmonic' , 
    'delta_aritmetic' , 'delta_geometric' , 
    'delta_quadratic' , 'delta_harmonic' ]



def get_df_8_aggregations_x_subsys(baseline, perturbed, name_string):

    my_8_aggregation = get_8_aggregations(baseline, perturbed)
    growing_df       = pd.DataFrame({'' : []})

    for i in range(len(my_8_aggregation)):
        col_names         =  [s + str('_'+aggregation[i] + '_' + name_string) for s in index_centralities ]

        df      = pd.DataFrame(
                my_8_aggregation[i],
                index = index_nodes,
                columns = col_names)
        new_df     = df
        growing_df = pd.concat([growing_df, new_df], axis=1)

    return growing_df  


 # %% 
perturbed_ETC_astrocyte.shape

  # %% 

ETC_astrocyte_full        = get_df_8_aggregations_x_subsys(baseline_ETC_astrocyte, perturbed_ETC_astrocyte, 'ETC_astrocyte')
ETC_neuron_full           = get_df_8_aggregations_x_subsys(baseline_ETC_neuron, perturbed_ETC_neuron, 'ETC_neuron')
Glycolysis_astrocyte_full = get_df_8_aggregations_x_subsys(baseline_Glycolysis_astrocyte, perturbed_Glycolysis_astrocyte, 'Glycolysis_astrocyte')
Glycolysis_neuron_full    = get_df_8_aggregations_x_subsys(baseline_Glycolysis_neuron, perturbed_Glycolysis_neuron, 'Glycolysis_neuron')

 # %% 

def get_by_aggregation(srint_aggregation): 
    def filter_cols_by_regex(df, patten):
        return df.loc[:, df.columns.str.contains(patten, regex=True)]

    srint_aggregation 
    agg = aggregation[aggregation.index(srint_aggregation)]

    df =  pd.concat([
    filter_cols_by_regex(ETC_astrocyte_full,agg).T,
    filter_cols_by_regex(ETC_neuron_full,agg).T,
    filter_cols_by_regex(Glycolysis_astrocyte_full,agg).T,
    filter_cols_by_regex(Glycolysis_neuron_full,agg).T,], axis=0).T
    return df


fold_change_aritmetic = get_by_aggregation('fold_change_aritmetic')
fold_change_geometric = get_by_aggregation('fold_change_geometric')
fold_change_quadratic = get_by_aggregation('fold_change_quadratic')
fold_change_harmonic = get_by_aggregation('fold_change_harmonic')
delta_aritmetic = get_by_aggregation('delta_aritmetic')
delta_geometric = get_by_aggregation('delta_geometric')
delta_quadratic = get_by_aggregation('delta_quadratic')
delta_harmonic = get_by_aggregation('delta_harmonic')
# %% Proof
print(
np.array(
fold_change_aritmetic.filter(regex = r'aritmetic_Glycolysis_astrocyte', axis=1).loc['L-LACt2r',:]) \
     == FC_from__arit_mean_L_LACt2r , \

np.array(
delta_aritmetic.filter(regex = r'aritmetic_Glycolysis_astrocyte', axis=1).loc['L-LACt2r',:]) \
     == delta_from_arit_mean_L_LACt2r )
######   FIN  ######