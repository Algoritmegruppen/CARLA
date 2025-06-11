from itertools import combinations
from collections import defaultdict, namedtuple
import numpy as np
import pandas as pd
import networkx as nx
from bentley_ottman_crossings import bentley_ottmann


Instance = namedtuple(
    "Instance", ["G", "graph_name", "node_to_community", "communities", "alpha"]
)


# DATA LOADING
def get_weighted_community_graph(G, community_alg="louvain_communities"):
    """Returns a graph where the nodes are the communities and the edges are
    the connections between the communities and a dictionary that maps the
    nodes to the communities.  The weight of the edges is the number of
    connections between the communities.

    """
    com_graph = nx.Graph()
    com_edge_weights = defaultdict(int)
    node_to_community = {}
    if community_alg == "louvain_communities":
        communities = nx.algorithms.community.louvain_communities(G)
    elif community_alg == "by_country":
        communities = defaultdict(list)
        for n, d in G.nodes(data=True):
            country = d.get("country")
            if country is not None:
                communities[country].append(n)
        communities = list(communities.values())
    elif community_alg == "label":
        communities = defaultdict(list)
        for n, d in G.nodes(data=True):
            class_ = d.get("label")
            if class_ is not None:
                communities[class_].append(n)
        communities = list(communities.values())
    else:
        ValueError("Invalid community_alg")

    for i, c in enumerate(communities):
        for n in c:
            node_to_community[n] = i

    for u, v, d in G.edges(data=True):
        new_u, new_v = node_to_community[u], node_to_community[v]
        if new_u != new_v:
            com_edge_weights[tuple(sorted((new_u, new_v)))] += 1

    max_weight = max(com_edge_weights.values())
    for (u, v), w in com_edge_weights.items():
        com_edge_weights[(u, v)] = round(w / max_weight, 3)

    for (u, v), w in com_edge_weights.items():
        com_graph.add_edge(u, v, weight=w)

    return com_graph, node_to_community, communities


def load_graph_with_com_based_weights(
    file_name, alpha=0.3, community_alg="louvain_communities"
):
    """Loads a graph from a file and computes community-based edge weights.
    Args:
    - file_name: path of the file containing the graph.
    - alpha: A scaling factor for edge weights.
    - type: Type of the file (e.g., 'gml', 'gexf', 'csv').
    - community_alg: Algorithm to use for community detection.
    Returns:
    - A dictionary containing the graph, community graph, node-to-community
    mapping, communities, graph name, and alpha.

    """
    filetype = file_name.split(".")[-1].lower().strip()
    G = nx.read_gexf(file_name)
    if filetype == "gml":
        G = nx.read_gml(file_name)
    elif filetype == "gexf":
        G = nx.read_gexf(file_name)
    elif filetype == "csv":
        df = pd.read_csv(file_name)
        G = nx.from_pandas_edgelist(
            df, source="source", target="target", edge_attr=True
        )
    else:
        print("No file", file_name)
        exit()

    # need a connected graph, so we take the biggest component
    biggest_component = max(nx.connected_components(G.to_undirected()), key=len)
    G = G.subgraph(biggest_component).copy()

    _, node_to_community, communities = get_weighted_community_graph(
        G, community_alg=community_alg
    )

    for u, v, d in G.edges(data=True):
        if node_to_community[u] == node_to_community[v]:
            d["edge_weight"] = 1 / alpha
        else:
            d["edge_weight"] = alpha
            # Might improve the layout by adding the weight of the community edge
    return Instance(
        G=G,
        graph_name=file_name.split(".")[0],
        node_to_community=node_to_community,
        communities=communities,
        alpha=alpha,
    )


# RECURSIVE INITIALIZATION FUNCTION
def recursive_community_layout(
    G, communities, center=(0, 0), scale=1, community_size=10, layer=0
):
    """Recursive function to detect communities in a graph, construct a
    community quotient graph, and process larger communities recursively.

    Args:
    - G: The input graph (networkx Graph)
    - center: Tuple (x, y) specifying the center of the layout for the current graph
    - scale: Scale factor for the layout
    - community_size: The maximum size of a community before recursive subdivision

    Returns:
    - A layout dictionary for the graph

    """

    total_nodes = len(G.nodes())

    if len(communities) == 1:
        community_subgraph = G.subgraph(communities[0])
        community_pos = nx.spring_layout(community_subgraph, center=center, scale=scale)
        layout = {node: community_pos[node] for node in communities[0]}
        return layout

    community_graph = nx.Graph()
    community_mapping = {
        node: idx for idx, comm in enumerate(communities) for node in comm
    }
    community_layout = {}

    for u, v in G.edges():
        cu, cv = community_mapping[u], community_mapping[v]
        if cu != cv:
            if community_graph.has_edge(cu, cv):
                community_graph[cu][cv]["weight"] += 1
            else:
                community_graph.add_edge(cu, cv, weight=1)
                community_layout = nx.spring_layout(
                    community_graph, center=center, scale=scale
                )
                layout = {}
    for idx, comm in enumerate(communities):
        comm_center = community_layout[idx]
        comm_scale = scale * (len(comm) / total_nodes) / 1

        if len(comm) > community_size:
            subgraph = G.subgraph(comm)
            sub_communities = nx.algorithms.community.louvain_communities(subgraph)
            sub_layout = recursive_community_layout(
                subgraph,
                sub_communities,
                center=comm_center,
                scale=comm_scale,
                community_size=community_size,
                layer=layer + 1,
            )
            layout.update(sub_layout)
        else:
            community_subgraph = G.subgraph(comm)
            community_pos = nx.spring_layout(
                community_subgraph, center=comm_center, scale=comm_scale
            )
            for node in comm:
                layout[node] = community_pos[node]

    return layout


# MEASURING FUNCTIONS
def atedge_length(G, pos):
    """
    Returns the total edge length of the graph G.
    """
    return sum(
        np.linalg.norm(np.array(pos[u]) - np.array(pos[v])) for u, v in G.edges()
    )


def scaled_atedge_length(G, pos):
    """
    Returns the scaled edge length of the graph G.
    sum of the edge lengths / total distance between all pairs of nodes
    """
    atedge_len = atedge_length(G, pos)
    vertex_len = sum(
        np.linalg.norm(np.array(pos[u]) - np.array(pos[v]))
        for u, v in combinations(G.nodes(), 2)
    )
    return atedge_len / vertex_len


def norm_atedge_lenght(G, pos):
    return scaled_atedge_length(G, pos) / ((len(G.edges) / (len(G.nodes) ** 2)))


def crossings(G, pos):
    """
    Returns the number of crossings in the graph G.
    """
    return bentley_ottmann(G, pos)


def plot(G, graph_name, layout):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6, 6))
    nx.draw_networkx_labels(
        G, layout, font_size=6, font_color="black", font_family="sans-serif"
    )
    nx.draw_networkx_edges(G, layout, alpha=0.1, width=0.6)
    nx.draw_networkx_nodes(G, layout, node_size=200, node_color="white", edgecolors="blue")
    title = f"Network: {graph_name}. "
    title += f"Crossings: {crossings(G, layout)}. "
    nel = 100*norm_atedge_lenght(G, layout)
    title += f"NEL: {nel:.1f}%."
    plt.title(title)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    import sys

    ITERATIONS = 1000
    ALPHA = 0.3
    if len(sys.argv) != 2:
        exit("Usage: python carla.py fname")
    fname = sys.argv[1]
    graph_data = load_graph_with_com_based_weights(fname, alpha=ALPHA)

    init_layout = recursive_community_layout(
        graph_data.G,
        graph_data.communities,
        center=(0, 0),
        scale=1,
        community_size=15,
    )
    layout = nx.spring_layout(
        graph_data.G,
        iterations=ITERATIONS,
        weight="edge_weight",
        pos=init_layout,
    )
    plot(graph_data.G, graph_data.graph_name, layout)
