from phi_transformation import *
from networkx import Graph
from synkit.IO import rsmi_to_graph
from multiset import Multiset
import logging
import sys

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)s - %(name)s - %(message)s',
    stream=sys.stdout,
    force=True
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def wl_vertex(graph, h_max):
    """
    Weisfeiler-Lehman algorithm using vertex-based phi transformation.
    Iteratively refines node labels based on neighborhood information.
    
    Args:
        graph: NetworkX graph object
        h_max: Maximum number of iterations
        
    Returns:
        Multiset of all feature labels collected across all iterations
    """
    # Initialize feature dictionary with vertex labels (no hashing)
    feature_set = phi_vertex_dict_graph_no_hash(graph)
    
    # Store all features across iterations
    all_features = Multiset()
    
    # Add initial features
    for node, labels in feature_set.items():
        for label in labels:
            all_features.add(label)
    
    logger.debug(f"wl_vertex - Initial feature_set: {feature_set}")
    logger.debug(f"wl_vertex - Initial all_features: {all_features}")
    
    # Iterate h_max times to refine labels
    for h in range(h_max):
        new_feature_set = {}
        
        for node, labels in feature_set.items():
            # Get immediate neighbors (1-hop)
            neighbors = list(graph.neighbors(node))
            
            # Collect neighbor labels
            neighbor_labels = Multiset()
            for neighbor in neighbors:
                neighbor_labels.update(feature_set[neighbor])
            
            # Combine current label with sorted neighbor labels
            combined = Multiset(labels)
            combined.update(neighbor_labels)
            
            # Create new label by concatenating sorted labels
            new_label = ''.join(sorted(combined))
            new_feature_set[node] = [new_label]
            
            # Add to global feature collection
            all_features.add(new_label)
            
            logger.debug(f"Iteration {h}, Node {node}: labels={labels}, neighbors={neighbors}, neighbor_labels={neighbor_labels}, new_label={new_label}")
        
        # Update feature set for next iteration
        feature_set = new_feature_set
    
    logger.debug(f"wl_vertex - Final all_features: {all_features}")
    return all_features


def wl_vertex_hashed(graph, h_max):
    """
    Weisfeiler-Lehman algorithm using vertex-based phi transformation with hashing.
    
    Args:
        graph: NetworkX graph object
        h_max: Maximum number of iterations
        
    Returns:
        Set of all hashed feature labels collected across all iterations
    """
    # Initialize feature dictionary with hashed vertex labels
    feature_set = phi_vertex_dict_graph(graph)
    
    # Store all features across iterations
    all_features = set()
    
    # Add initial features
    for node, labels in feature_set.items():
        for label in labels:
            all_features.add(label)
    
    logger.debug(f"wl_vertex_hashed - Initial feature_set: {feature_set}")
    
    # Iterate h_max times to refine labels
    for h in range(h_max):
        new_feature_set = {}
        
        for node, labels in feature_set.items():
            # Get immediate neighbors (1-hop)
            neighbors = list(graph.neighbors(node))
            
            # Collect neighbor labels
            neighbor_labels = []
            for neighbor in neighbors:
                neighbor_labels.extend(feature_set[neighbor])
            
            # Combine current label with sorted neighbor labels
            combined = labels + neighbor_labels
            combined_sorted = sorted(combined)
            
            # Hash the combined label
            combined_str = ''.join(combined_sorted)
            new_label_hash = get_hash(combined_str)
            new_feature_set[node] = [new_label_hash]
            
            # Add to global feature collection
            all_features.add(new_label_hash)
            
            logger.debug(f"Iteration {h}, Node {node}: combined={combined_sorted[:3]}..., new_hash={new_label_hash[:16]}...")
        
        # Update feature set for next iteration
        feature_set = new_feature_set
    
    logger.debug(f"wl_vertex_hashed - Total features collected: {len(all_features)}")
    return all_features


def wl_edge(graph, h_max):
    """
    Weisfeiler-Lehman algorithm using edge-based phi transformation.
    Refines edge labels based on incident edges.
    
    Args:
        graph: NetworkX graph object
        h_max: Maximum number of iterations
        
    Returns:
        Set of all edge feature labels collected across all iterations
    """
    # Get initial edge features (no hashing)
    edge_features = phi_edge_graph_no_hash(graph)
    
    # Store all features across iterations
    all_features = set(edge_features)
    
    logger.debug(f"wl_edge - Initial edge_features: {edge_features}")
    
    # Create edge dictionary for easier lookup
    edge_dict = {}
    for u, v, d in graph.edges(data=True):
        u, v = sorted([u, v])
        edge_key = (u, v)
        edge_dict[edge_key] = f"{u}{d['order']}{v}"
    
    # Iterate h_max times
    for h in range(h_max):
        new_edge_dict = {}
        
        for (u, v), edge_label in edge_dict.items():
            # Find all edges incident to this edge (connected via u or v)
            incident_edges = []
            
            # Edges connected to u
            for neighbor in graph.neighbors(u):
                edge_key = tuple(sorted([u, neighbor]))
                if edge_key in edge_dict and edge_key != (u, v):
                    incident_edges.append(edge_dict[edge_key])
            
            # Edges connected to v
            for neighbor in graph.neighbors(v):
                edge_key = tuple(sorted([v, neighbor]))
                if edge_key in edge_dict and edge_key != (u, v):
                    incident_edges.append(edge_dict[edge_key])
            
            # Combine current edge label with sorted incident edge labels
            combined = edge_label + ''.join(sorted(incident_edges))
            new_edge_dict[(u, v)] = combined
            
            # Add to global feature collection
            all_features.add(combined)
            
            logger.debug(f"Iteration {h}, Edge ({u},{v}): label={edge_label}, incidents={len(incident_edges)}, new_label={combined[:50]}...")
        
        # Update edge dictionary for next iteration
        edge_dict = new_edge_dict
    
    logger.debug(f"wl_edge - Total features collected: {len(all_features)}")
    return all_features


def wl_edge_hashed(graph, h_max):
    """
    Weisfeiler-Lehman algorithm using edge-based phi transformation with hashing.
    
    Args:
        graph: NetworkX graph object
        h_max: Maximum number of iterations
        
    Returns:
        Set of all hashed edge feature labels collected across all iterations
    """
    # Get initial edge features (hashed)
    edge_features = phi_edge_graph(graph)
    
    # Store all features across iterations
    all_features = set(edge_features)
    
    logger.debug(f"wl_edge_hashed - Initial edge_features count: {len(edge_features)}")
    
    # Create edge dictionary with hashed labels
    edge_dict = {}
    for u, v, d in graph.edges(data=True):
        u, v = sorted([u, v])
        edge_key = (u, v)
        edge_label = f"{u}{d['order']}{v}"
        edge_dict[edge_key] = get_hash(edge_label)
    
    # Iterate h_max times
    for h in range(h_max):
        new_edge_dict = {}
        
        for (u, v), edge_hash in edge_dict.items():
            # Find all edges incident to this edge
            incident_hashes = []
            
            # Edges connected to u
            for neighbor in graph.neighbors(u):
                edge_key = tuple(sorted([u, neighbor]))
                if edge_key in edge_dict and edge_key != (u, v):
                    incident_hashes.append(edge_dict[edge_key])
            
            # Edges connected to v
            for neighbor in graph.neighbors(v):
                edge_key = tuple(sorted([v, neighbor]))
                if edge_key in edge_dict and edge_key != (u, v):
                    incident_hashes.append(edge_dict[edge_key])
            
            # Combine current edge hash with sorted incident hashes
            combined_str = edge_hash + ''.join(sorted(incident_hashes))
            new_hash = get_hash(combined_str)
            new_edge_dict[(u, v)] = new_hash
            
            # Add to global feature collection
            all_features.add(new_hash)
            
            logger.debug(f"Iteration {h}, Edge ({u},{v}): incidents={len(incident_hashes)}, new_hash={new_hash[:16]}...")
        
        # Update edge dictionary for next iteration
        edge_dict = new_edge_dict
    
    logger.debug(f"wl_edge_hashed - Total features collected: {len(all_features)}")
    return all_features


def wl_shortest_path(graph, h_max=1):
    """
    Weisfeiler-Lehman algorithm using shortest path-based phi transformation.
    Uses all shortest paths as features.
    
    Note: h_max is typically 1 for shortest path variant since paths already
    capture neighborhood information at various distances.
    
    Args:
        graph: NetworkX graph object
        h_max: Maximum number of iterations (default 1)
        
    Returns:
        Set of all shortest path feature labels
    """
    # Get shortest path features (no hashing)
    path_features = phi_shortest_path_graph_no_hash(graph)
    
    logger.debug(f"wl_shortest_path - Path features count: {len(path_features)}")
    logger.debug(f"wl_shortest_path - Sample paths: {list(path_features)[:5]}")
    
    # For shortest path variant, we typically don't iterate
    # as the paths already encode structural information
    all_features = set(path_features)
    
    # Optional: iterate to refine features
    if h_max > 1:
        current_features = path_features
        for h in range(1, h_max):
            new_features = set()
            # Combine existing paths to create new features
            for path in current_features:
                # Create variations by combining with other paths
                for other_path in list(current_features)[:10]:  # Limit to avoid explosion
                    combined = ''.join(sorted([path, other_path]))
                    new_features.add(combined)
            
            all_features.update(new_features)
            current_features = new_features
            logger.debug(f"Iteration {h}, new features: {len(new_features)}")
    
    logger.debug(f"wl_shortest_path - Total features collected: {len(all_features)}")
    return all_features


def wl_shortest_path_hashed(graph, h_max=1):
    """
    Weisfeiler-Lehman algorithm using shortest path-based phi transformation with hashing.
    
    Args:
        graph: NetworkX graph object
        h_max: Maximum number of iterations (default 1)
        
    Returns:
        Set of all hashed shortest path feature labels
    """
    # Get shortest path features (hashed)
    path_features = phi_shortest_path_graph(graph)
    
    logger.debug(f"wl_shortest_path_hashed - Path features count: {len(path_features)}")
    
    all_features = set(path_features)
    
    # Optional: iterate to refine features
    if h_max > 1:
        current_features = path_features
        for h in range(1, h_max):
            new_features = set()
            # Combine existing path hashes
            for path_hash in current_features:
                for other_hash in list(current_features)[:10]:  # Limit combinations
                    combined = ''.join(sorted([path_hash, other_hash]))
                    new_hash = get_hash(combined)
                    new_features.add(new_hash)
            
            all_features.update(new_features)
            current_features = new_features
            logger.debug(f"Iteration {h}, new features: {len(new_features)}")
    
    logger.debug(f"wl_shortest_path_hashed - Total features collected: {len(all_features)}")
    return all_features


def wl_kernel(graph1, graph2, h_max, method='vertex', use_hash=False):
    """
    Compute Weisfeiler-Lehman graph kernel between two graphs.
    
    Args:
        graph1: First NetworkX graph
        graph2: Second NetworkX graph
        h_max: Maximum number of iterations
        method: 'vertex', 'edge', or 'shortest_path'
        use_hash: Whether to use hashed versions
        
    Returns:
        Kernel value (number of common features)
    """
    # Select appropriate WL function
    if method == 'vertex':
        wl_func = wl_vertex_hashed if use_hash else wl_vertex
    elif method == 'edge':
        wl_func = wl_edge_hashed if use_hash else wl_edge
    elif method == 'shortest_path':
        wl_func = wl_shortest_path_hashed if use_hash else wl_shortest_path
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Get features for both graphs
    features1 = wl_func(graph1, h_max)
    features2 = wl_func(graph2, h_max)
    
    # Compute kernel as intersection size
    if isinstance(features1, Multiset) and isinstance(features2, Multiset):
        # For multisets, compute the sum of minimum counts
        common = features1 & features2
        kernel_value = sum(common.values())
    else:
        # For sets, compute intersection size
        kernel_value = len(features1 & features2)
    
    logger.debug(f"wl_kernel - Method: {method}, Hash: {use_hash}")
    logger.debug(f"wl_kernel - Graph1 features: {len(features1)}, Graph2 features: {len(features2)}, Common: {kernel_value}")
    
    return kernel_value


# Test helper functions
def get_test_graph():
    """Get a test graph from a sample reaction SMILES"""
    g, _ = rsmi_to_graph("[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][N:8](Cc2ccccc2)[CH2:13][CH2:12]1>>[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1")
    return g


def test_all_methods():
    """Test all WL variants on a sample graph"""
    graph = get_test_graph()
    h_max = 3
    
    print("=" * 60)
    print("Testing Weisfeiler-Lehman variants")
    print("=" * 60)
    
    print("\n1. WL Vertex (no hash):")
    features = wl_vertex(graph, h_max)
    print(f"   Features collected: {len(features)}")
    print(f"   Sample: {list(features)[:3]}")
    
    print("\n2. WL Vertex (hashed):")
    features = wl_vertex_hashed(graph, h_max)
    print(f"   Features collected: {len(features)}")
    print(f"   Sample: {list(features)[:3]}")
    
    print("\n3. WL Edge (no hash):")
    features = wl_edge(graph, h_max)
    print(f"   Features collected: {len(features)}")
    print(f"   Sample: {list(features)[:3]}")
    
    print("\n4. WL Edge (hashed):")
    features = wl_edge_hashed(graph, h_max)
    print(f"   Features collected: {len(features)}")
    print(f"   Sample: {list(features)[:3]}")
    
    print("\n5. WL Shortest Path (no hash):")
    features = wl_shortest_path(graph, h_max=1)
    print(f"   Features collected: {len(features)}")
    print(f"   Sample: {list(features)[:3]}")
    
    print("\n6. WL Shortest Path (hashed):")
    features = wl_shortest_path_hashed(graph, h_max=1)
    print(f"   Features collected: {len(features)}")
    print(f"   Sample: {list(features)[:3]}")
    
    print("\n" + "=" * 60)
    print("Testing WL Kernel")
    print("=" * 60)
    
    # Test kernel with same graph
    kernel_val = wl_kernel(graph, graph, h_max=2, method='vertex', use_hash=False)
    print(f"\nKernel (same graph, vertex): {kernel_val}")
    
    kernel_val = wl_kernel(graph, graph, h_max=2, method='edge', use_hash=True)
    print(f"Kernel (same graph, edge, hashed): {kernel_val}")


if __name__ == "__main__":
    test_all_methods()
