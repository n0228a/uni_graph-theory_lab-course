"""Computes a reaction signature using a Weisfeiler-Lehman-based graph kernel."""
from synkit.IO import rsmi_to_its, rsmi_to_graph
from synkit.Vis import GraphVisualizer
import matplotlib.pyplot as plt
import networkx as nx
import hashlib
import pandas as pd

def get_hash(data: str):
    """Returns a stable 64-bit integer from a string."""
    return int(hashlib.blake2b(data.encode(), digest_size=8).hexdigest(), 16)

# Example reaction SMILES string
rsmi = '[CH3:1][CH:2]=[O:3].[CH:4]([H:7])([H:8])[CH:5]=[O:6]>>[CH3:1][CH:2]=[CH:4][CH:5]=[O:6].[O:3]([H:7])([H:8])'


# Parse SMILES into educt and product graphs
educt_graph, product_graph = rsmi_to_graph(rsmi)
its_graph = rsmi_to_its(rsmi)


def getWL(graph, h_max):
    """
    Computes node, edge, and shortest-path features for a graph using a WL-like algorithm.

    Args:
        graph (nx.Graph): Input molecular graph.
        h_max (int): Number of WL iterations.

    Returns:
        tuple[set, set, set]: Node, edge, and shortest-path feature sets.
    """
    feature_setN = set()
    feature_setE = set()
    feature_setSP = set()
    shortest_pathsL = dict(nx.shortest_path_length(graph))
    shortest_paths = dict(nx.shortest_path(graph))
    # Initialize labels with element types
    labels = {n: str(graph.nodes[n].get('element')) for n in graph.nodes()}
    print(labels)
    for label in labels.values():
        feature_setN.add(get_hash(label))
    for h in range(h_max+1):
        # Generate edge features
        for u,v,d in graph.edges(data=True):
            l_a = labels[u]
            l_b = labels[v]
            l_ab = str(d.get('order'))

            node_pair = sorted([l_a, l_b])
            triplet = f"{node_pair[0]}{l_ab}{node_pair[1]}"
            #print(triplet)
            feature_setE.add(get_hash(triplet))

        new_labels = {}
        # Generate shortest-path features
        for n in graph.nodes():
            for m in shortest_pathsL[n]:
                if n != m:
                    # Keep the original Element-Distance-Element feature
                    l_n = labels[n]
                    l_m = labels[m]
                    distance = shortest_pathsL[n][m]
                    node_pair = sorted([l_m, l_n])
                    sp_feature_dist = f"{node_pair[0]}-{distance}-{node_pair[1]}"
                    #print(sp_feature_dist)
                    feature_setSP.add(get_hash(sp_feature_dist))

                    # Concatenated labels along the shortest path
                    path_nodes = shortest_paths[n][m]
                    path_labels = [labels[node] for node in path_nodes]
                    sp_feature_forward = "".join(path_labels)
                    sp_feature_backward = "".join(reversed(path_labels))
                    # directional invariance of the paths
                    sp_feature_path = min(sp_feature_forward, sp_feature_backward)
                    #print(sp_feature_path)
                    feature_setSP.add(get_hash(sp_feature_path))

            # Update node labels for the next iteration (WL aggregation)
            if h<h_max:
                current = labels[n]
                #print(current)
                neighbor_labels = []
                for neighbor in graph.neighbors(n):
                    neighbor_labels.append(labels[neighbor])
                    #print(neighbor_labels)
                neighbor_labels.sort()
                combined = current + "".join(neighbor_labels)
                #print(combined)
                new_labels[n] = combined

        # Update labels and add new features
        if h<h_max:
            labels = new_labels
            #print(labels)
            for label in labels.values():
                feature_setN.add(get_hash(label))
        #print(labels)
        print(len(feature_setSP))

    print(f"Final feature set size: {len(feature_setN)}")
    return feature_setN, feature_setSP, feature_setE

# --- Main Execution ---
# Generate features for educt and product graphs
a1, a2, a3 = getWL(educt_graph, 4)
b1, b2, b3 = getWL(product_graph, 4)
c1, c2, c3 = getWL(its_graph, 4)

# The reaction signature is the symmetric difference of the node features.
signature1 = a1.symmetric_difference(b1)
signature2 = a2.symmetric_difference(b2)
signature3 = a3.symmetric_difference(b3)

# Create dataframe consisting of 6 columns:
df = pd.DataFrame(columns=[
    "DRF Edges",
    "DRF Shortest Paths",
    "ITS Nodes",
    "ITS Edges",
    "ITS Shortest Paths"
])

if __name__ == "__main__":
    # For each SMILES execute the WL algorithm and add to dataframe
    df_schneider = pd.read_csv("schneider50k_clean.tsv", sep="\t")
    for index, row in df_schneider.iterrows():
        rsmi = row["clean_rxn"]
        try:
            educt_graph, product_graph = rsmi_to_graph(rsmi)
            its_graph = rsmi_to_its(rsmi)
            nodes_e, sps_e, edges_e = getWL(educt_graph, 4)
            nodes_p, sps_p, edges_p = getWL(product_graph, 4)
            nodes_its, sps_its, edges_its = getWL(its_graph, 4)

            drf_nodes = nodes_e.symmetric_difference(nodes_p)
            drf_sps = sps_e.symmetric_difference(sps_p)
            drf_edges = edges_e.symmetric_difference(edges_p)
        except Exception as e:
            print(f"Error processing reaction {index}: {e}")
            drf_nodes = drf_edges = drf_sps = nodes_its = edges_its = sps_its = set()

        df = pd.concat([df, pd.DataFrame([{
            "DRF Nodes": drf_nodes,
            "DRF Edges": drf_edges,
            "DRF Shortest Paths": drf_sps,
            "ITS Nodes": nodes_its,
            "ITS Edges": edges_its,
            "ITS Shortest Paths": sps_its
        }])], ignore_index=True)

    # Write to Excel
    df.to_excel("pre-computed-feature_sets.xlsx", index=False)