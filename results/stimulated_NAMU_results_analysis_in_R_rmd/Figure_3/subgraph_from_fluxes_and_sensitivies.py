# %%
from re import T
import cobra
import warnings 
import pandas as pd
from pandas.core.base import DataError
warnings.filterwarnings("ignore")
def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G


path          = '/home/alejandro/PostDoc/human-metnet/results/stimulated_NAMU_results_analysis_in_R_rmd/Figure_1'
import networkx as nx
reaction_projected_stimnulated = nx.read_gpickle( path + "/reaction_projected_stimnulated.gpickle")

stimulated    = cobra.io.load_json_model(path + '/stimulated_2021.json') 
nodes_attr_df = pd.read_csv(path + '/node_attributes.csv')
fba_solution  = stimulated.optimize()
fluxes        = fba_solution.fluxes.to_frame().query("abs(fluxes) > 1e-12")
reduced_costs = fba_solution.reduced_costs.to_frame().query("abs(reduced_costs) > 1e-12")
optimal_node_names = list(fluxes.join(reduced_costs, how='outer').index)


#  Extract optimal subgraph

H1 = reaction_projected_stimnulated.subgraph(optimal_node_names)
H2 = get_largest_component(H1)
print(
len(H1.nodes(data = True)), len(H2.nodes))

# %%
optimal_subgraph = H2.copy()



nx.write_gexf(optimal_subgraph, path + "/optimal_subgraph.gexf")