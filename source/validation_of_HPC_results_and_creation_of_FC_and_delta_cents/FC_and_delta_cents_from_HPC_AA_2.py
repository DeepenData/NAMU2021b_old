# %% 

from networkx.algorithms.isomorphism.ismags import intersect
import pandas   as pd
import numpy    as np
import networkx as nx
import pickle 
# Esto asume que está en 'human-metnet/'
#%cd ..
#%cd ..

infile = open('./tmp/index_nodes','rb');     index_nodes = pickle.load(infile); infile.close()
#infile = open('./tmp/subsystems_dict','rb'); subsystems = pickle.load(infile); infile.close() 

infile = open('./tmp/centralidades_perturbadas.pkl','rb'); centralidades_perturbadas = pickle.load(infile); infile.close()
infile = open('./tmp/baseline.pkl' ,'rb'); baseline  = pickle.load(infile); infile.close()


index_nodes = list(index_nodes)

Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']
index_centralities = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality',
      'current_flow_closeness_centrality', 'current_flow_betweenness_centrality',
       'approximate_current_flow_betweenness_centrality', 'communicability_betweenness_centrality']

# %% Baseline

baseline['baseline'].loc[Glycolysis_astrocyte]

# %% Perturbed 


centralidades_perturbadas['FPGS3m'].loc[Glycolysis_astrocyte]


# %% FC y delta cents
def aritmetic_mean(df):
    return np.nanmean( df, axis= 0)



baseline_Glycolysis_astrocyte_df             = baseline['baseline'].loc[Glycolysis_astrocyte]
perturbed_hpc_Glycolysis_astrocyte_FPGS3m_df = centralidades_perturbadas['FPGS3m'].loc[Glycolysis_astrocyte]

ratio_from_arit_mean_FPGS3m  = aritmetic_mean(baseline_Glycolysis_astrocyte_df)/aritmetic_mean(perturbed_hpc_Glycolysis_astrocyte_FPGS3m_df)


FC_from__arit_mean_FPGS3m    = np.log2(ratio_from_arit_mean_FPGS3m)
delta_from_arit_mean_FPGS3m  = (aritmetic_mean(baseline_Glycolysis_astrocyte_df) - aritmetic_mean(perturbed_hpc_Glycolysis_astrocyte_FPGS3m_df))/aritmetic_mean(baseline_Glycolysis_astrocyte_df)
print(
FC_from__arit_mean_FPGS3m,
delta_from_arit_mean_FPGS3m)

# %% --- Importar los tensores con los resultados HPC
import pandas   as pd
import numpy    as np
import networkx as nx
import pickle 


infile = open('./tmp/centralidades_perturbadas_tensor.pkl','rb'); perturbed_centralities = pickle.load(infile); infile.close()
infile = open('./tmp/baseline_tensor.pkl' ,'rb'); baseline_centralities  = pickle.load(infile); infile.close()


#infile = open('./tmp/index_nodes','rb'); index_nodes = pickle.load(infile); infile.close()



#index_nodes = list(index_nodes)

index_nodes = list( baseline['baseline'].index )

# %% 

# TODO: reemplazar esto por una importación desde Tensor Maker
def subsystem_tensor( subsystem_nodes, tensor):
    """Genera un tensor con un subset de nodos basados en el subsistema pedido"""
    subsystem_index = [ index_nodes.index(node) for node in subsystem_nodes ]
    return tensor[:, subsystem_index ,:]

Glycolysis_astrocyte = ['PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m']
Glycolysis_neuron = ['ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron',
 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron']
ETC_neuron    = ['ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron']
ETC_astrocyte =  ['PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA']

perturbed_Glycolysis_astrocyte = subsystem_tensor( Glycolysis_astrocyte, perturbed_centralities )
baseline_Glycolysis_astrocyte  = subsystem_tensor( Glycolysis_astrocyte, baseline_centralities   )

perturbed_Glycolysis_neuron = subsystem_tensor( Glycolysis_neuron, perturbed_centralities )
baseline_Glycolysis_neuron  = subsystem_tensor( Glycolysis_neuron, baseline_centralities   )

perturbed_ETC_neuron = subsystem_tensor( ETC_neuron, perturbed_centralities )
baseline_ETC_neuron  = subsystem_tensor( ETC_neuron, baseline_centralities   )

perturbed_ETC_astrocyte= subsystem_tensor( ETC_astrocyte, perturbed_centralities )
baseline_ETC_astrocyte  = subsystem_tensor( ETC_astrocyte, baseline_centralities   )

index_centralities = ['harmonic_centrality', 'eigenvector_centrality', 'degree_centrality', 'betweenness_centrality', 
      'closeness_centrality', 'load_centrality', 'information_centrality', 'second_order_centrality',
      'current_flow_closeness_centrality', 'current_flow_betweenness_centrality',
       'approximate_current_flow_betweenness_centrality', 'communicability_betweenness_centrality']


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
# %% --- EXPORTA ESTO EN UN EXCEL POR SHEETS

dflist = [
    fold_change_aritmetic, fold_change_geometric, fold_change_quadratic, fold_change_harmonic,
    delta_aritmetic, delta_geometric, delta_quadratic, delta_harmonic 
    ]

dflist_names = [
    'fold_change_aritmetic', 'fold_change_geometric', 'fold_change_quadratic', 'fold_change_harmonic',
    'delta_aritmetic', 'delta_geometric', 'delta_quadratic', 'delta_harmonic'
]

# Define un creador del Excel
SALIDA = './source/validation_of_HPC_results_and_creation_of_FC_and_delta_cents/delta_fc_centralidades.xlsx'
Excelwriter = pd.ExcelWriter( SALIDA ,engine="xlsxwriter")

#We now loop process the list of dataframes
for i, df in enumerate (dflist):
    df.to_excel(Excelwriter, sheet_name=dflist_names[i] , index=True)
#And finally save the file
Excelwriter.save()


# %% Proof


FC_tensor = np.array(fold_change_aritmetic.filter(regex = r'aritmetic_Glycolysis_astrocyte', axis=1).loc['FPGS3m',:])
FC_naive  =  FC_from__arit_mean_FPGS3m



delta_tensor = np.array(delta_aritmetic.filter(regex = r'aritmetic_Glycolysis_astrocyte', axis=1).loc['FPGS3m',:])
delta_naive  =  delta_from_arit_mean_FPGS3m


print(
np.around( FC_tensor, 5) == np.around(FC_naive , 5),
np.around( delta_tensor, 5) == np.around(delta_naive , 5))





######   FIN  ######


# %%
