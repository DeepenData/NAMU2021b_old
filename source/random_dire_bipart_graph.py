# %%
def str_match(node_list, pattern):
    import re
    hh  = [re.findall(pattern, a_node) for a_node in node_list]
    nodes = [item for sublist in hh for item in sublist]
    return nodes

def generate_rbd_graph(A, B, p):
    from networkx.algorithms.bipartite.generators import random_graph
    import networkx as nx
    import random as rd
    import numpy as np
    rbd_graph = random_graph(A,B,p,directed=False)
    rbd_graph = nx.DiGraph(rbd_graph).to_directed()
    ebunch    = rd.sample(list(rbd_graph.edges),B)
    rbd_graph.remove_edges_from(ebunch)

    while np.invert(nx.is_connected(nx.Graph(rbd_graph))):
        rbd_graph = random_graph(A,B,p,directed=False)
        rbd_graph = nx.DiGraph(rbd_graph).to_directed()
        ebunch    = rd.sample(list(rbd_graph.edges),B)
        rbd_graph.remove_edges_from(ebunch)
    
    print("Bipartite:",
    nx.is_bipartite(rbd_graph),"\nDirected:",
    nx.is_directed(rbd_graph),"\nConnected:",
    nx.is_connected(nx.Graph(rbd_graph)))

    return rbd_graph


def plot_bipartite(G, node_color = "grey"):
    from networkx import bipartite
    import networkx as nx
    import matplotlib.pyplot as plt

    X, Y = bipartite.sets(G)
    pos = dict()
    pos.update( (n, (1, i)) for i, n in enumerate(X) ) # put nodes from X at x=1
    pos.update( (n, (2, i)) for i, n in enumerate(Y) ) # put nodes from Y at x=2
    nx.draw(G,  node_color= node_color,  pos=pos, cmap=plt.get_cmap('viridis'), with_labels=True, font_color='red')
    return plt.show()

def set_node_names(G):
    import networkx as nx
    from networkx import bipartite

    A, B = bipartite.sets(G)

    rxns      = ["rxn" + str(i) for i in range(0,len(A))]
    mets      = ["met" + str(i) for i in range(0,len(B))]
    mapping_A = dict(zip(A, rxns ))
    mapping_B = dict(zip(B,  mets ))

    G = nx.relabel_nodes(G, mapping_A)
    G = nx.relabel_nodes(G, mapping_B)
    return G

def make_union(G1, H1):
    import networkx as nx
    from networkx.algorithms.operators.binary import union
    G1_H1 = union(G1,H1, rename=('G-','H-') )
    G1_H1 = nx.DiGraph(G1_H1).to_directed()
    return G1_H1

def make_bridges(G1_H1):
    node_list =  list(G1_H1.nodes)
    G_rxn = str_match(node_list, "G.rxn.*")
    G_met = str_match(node_list, "G.met.*")
    H_rxn = str_match(node_list, "H.rxn.*")
    H_met = str_match(node_list, "H.met.*")

    import random as rd

    G_rxn_samples = rd.sample(G_rxn,3)
    G_met_samples = rd.sample(G_met,3)

    H_rxn_samples = rd.sample(H_rxn,3)
    H_met_samples = rd.sample(H_met,3)

    bridges = list()
    for G_r, G_m, H_r, H_m in zip(G_rxn_samples, G_met_samples, H_rxn_samples, H_met_samples):
        
        e1 = (G_r, H_m,{})
        e2 = (H_m, G_r,{})
        e3 = (H_r, G_m,{})
        e4 = (G_m, H_r,{})
        bridges.append([e1, e2, e3, e4])
        
    bridges = [item for sublist in bridges for item in sublist]
    G1_H1.add_edges_from(bridges)  

    import networkx as nx

    print("Bipartite:",
    nx.is_bipartite(G1_H1),"\nDirected:",
    nx.is_directed(G1_H1),"\nConnected:",
    nx.is_connected(nx.Graph(G1_H1)))
    return G1_H1

def plot_bipartite_subgraphs(G1_H1):
    from networkx import bipartite
    import networkx as nx
    import matplotlib.pyplot as plt

    X, Y = bipartite.sets(G1_H1)

    GH=str_match(list(G1_H1.nodes), "G|H") 

    import re

    GH = [re.sub('G','green', gh) for gh in GH ]
    GH = [re.sub('H','blue', gh) for gh in GH ]

    pos = dict()
    pos.update( (n, (1, i)) for i, n in enumerate(X) ) # put nodes from X at x=1
    pos.update( (n, (2, i)) for i, n in enumerate(Y) ) # put nodes from Y at x=2
    nx.draw_networkx(G1_H1, arrowsize = 5, width = .6, node_size =30, alpha = 0.6, \
        node_color= GH,  pos=pos, edge_color = "grey",  \
        with_labels=False, font_color='red')
    return plt.show()

def generate_pa_graph(A, p):
    from networkx.algorithms.bipartite.generators import preferential_attachment_graph
    import networkx as nx
    import random as rd
    import numpy as np
    from networkx import bipartite
    pa_graph = preferential_attachment_graph(aseq = A_degree_sequence,  p = p)
    pa_graph = nx.DiGraph(pa_graph).to_directed()
    ebunch    = rd.sample(list(pa_graph.edges),len(A))
    pa_graph.remove_edges_from(ebunch)
    

    while  np.invert(nx.is_connected(nx.Graph(pa_graph))) :
        pa_graph = preferential_attachment_graph(aseq = A_degree_sequence,  p = p)
        pa_graph = nx.DiGraph(pa_graph).to_directed()
        ebunch    = rd.sample(list(pa_graph.edges),len(A)-2)
        pa_graph.remove_edges_from(ebunch)
        #A2, B2 = bipartite.sets(pa_graph)
    A2, B2 = bipartite.sets(pa_graph)
    print("Bipartite:",
    nx.is_bipartite(pa_graph),"\nDirected:",
    nx.is_directed(pa_graph),"\nConnected:",
    nx.is_connected(nx.Graph(pa_graph)), "\nlen(A2) > len(B2)" ,(len(A2) > len(B2)))

    return pa_graph

def plot_degree_distribution(G):
    import collections
    import matplotlib.pyplot as plt
    import networkx as nx

    G = nx.Graph(G)

    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color="b")

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)

    # draw graph in inset
    plt.axes([0.4, 0.4, 0.5, 0.5])
    Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    pos = nx.spring_layout(G)
    plt.axis("off")
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    return plt.show()

# %%
#Random
#Creación de dos grafos random (generate_rbd_graph): G1 y H1. 
G1    =  generate_rbd_graph(A = 16, B = 12, p = 0.11) #A y B, son la cantidad de nodos en cada partición. p: probabilidad de edge
#Se asignan nombres de pseudo-reacciones y metabolitos (set_node_names).
G1    =  set_node_names(G1)
H1    =  generate_rbd_graph(A = 16, B = 12, p = 0.11)
H1    =  set_node_names(H1)
#La union (make_union) de de los grafos y conexión (make_bridges) de los dos grafos (G1 y H1), 
# cada nombre de nodo tendrá un prefijo de su grafo de origen (G- o H-). 
# Los grafos de origen son las subredes.
G1_H1 =  make_union(G1, H1)
G1_H1 =  make_bridges(G1_H1)

# %%
#Free-scale
#Generar una secuencia de grados (A_degree_sequence) en base al grafo anterior G1 para crear una grafo con "Prefential Attachment" 
#que da origen a la topología de 'escala-libre'.
#crear A_degree_sequence:
from networkx import bipartite
degree_sequence = [d for n, d in G1.degree()]
A, B = bipartite.sets(G1)
A_degree_sequence = degree_sequence[0:len(A)]

G2    =  generate_pa_graph(A = A_degree_sequence, p = 0.2)
G2    =  set_node_names(G2)
H2    =  generate_pa_graph(A = A_degree_sequence, p = 0.2)
H2    =  set_node_names(H2)
G2_H2 =  make_union(G2, H2)
G2_H2 =  make_bridges(G2_H2)
# %%
plot_bipartite_subgraphs(G1_H1)
plot_degree_distribution(G1_H1)
plot_bipartite_subgraphs(G2_H2)
plot_degree_distribution(G2_H2)
# %%
