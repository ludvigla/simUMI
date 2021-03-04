import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

def edit_dist(a, b):
    """

    Edit distance

    Returns the edit distance/hamming distances between 
    two strings of equal length.

    Parameters
    ----------
    a, b : str
    	UMI string. Lengths of a and b have to be equal

    Returns
    ----------
    int
        hamming distance

    """

    try:
        assert len(a) == len(b)
    except AssertionError:
        print("strings are of unequal lengths")
    dist = sum([not a == b for a, b in zip(a, b)])
    return dist


def dedup_naive(molecular_barcodes, mismatches=1):
    """

    Naive duplicate removal

    Sorts UMIs and creates a new group if the edit distance 
    between a UMI and the preceeding UMI is larger than 
    the number of allowed mismatches.

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values
    mismatched : int
        number of allowed mismatches

    Returns
    ----------
    int
        Number of UMI groups/clusters

    """

    clusters_dict = {}
    nclusters = 0
    for i, molecular_barcode in enumerate(sorted(molecular_barcodes.keys())):
        if i == 0:
            clusters_dict[nclusters] = [molecular_barcode]
        else:
            # compare distant of previous molecular barcodes and new one
            # if distance is between threshold we add it to the cluster 
            # otherwise we create a new cluster
            if edit_dist(clusters_dict[nclusters][-1], molecular_barcode) <= mismatches:
                clusters_dict[nclusters].append(molecular_barcode)
            else:
                nclusters += 1
                clusters_dict[nclusters] = [molecular_barcode]
    return len(clusters_dict)


def dedup_hierarchical(molecular_barcodes, mismatches=1, method="single"):
    """

    Hierarchical duplicate removal

    Runs a hierarchical clustering on the edit distance matrix
    computed for all paris of UMIs.

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values
    mismatches : int
        number of allowed mismatches
    method : str
        method to be used ["single", "complete", "ward"]

    Returns
    ----------
    int
        Number of UMI groups/clusters

    """

    molecular_barcodes = list(molecular_barcodes.keys())

    def d(coord):
        i, j = coord
        return edit_dist(molecular_barcodes[i], molecular_barcodes[j])

    # Create hierarchical clustering and obtain flat clusters at the distance given
    indices = np.triu_indices(len(molecular_barcodes), 1)
    distance_matrix = np.apply_along_axis(d, 0, indices)
    linkage_cluster = linkage(distance_matrix, method=method)
    flat_clusters = fcluster(linkage_cluster, mismatches, criterion='distance')
    return len(set(flat_clusters))


def dedup_unique(molecular_barcodes):
    """
    
    Unique duplicate removal

    Count all unique UMIs.

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values

    Returns
    ----------
    int
        Number of UMI groups

    """
    return len(molecular_barcodes.keys())


def dedup_percentile(molecular_barcodes):
    """
    
    Percentile duplicate removal

    Count all UMIs with a total count higher than the 
    average count across all UMIs divided by 100.

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values

    Returns
    ----------
    int
        Number of UMI groups

    """
    threshold = np.mean(list(molecular_barcodes.values())) / 100
    return len([umi for umi in list(molecular_barcodes.keys()) if molecular_barcodes[umi] > threshold])


def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)

    return found


def dedup_graph(molecular_barcodes, mismatches=1):
    """
    
    Graph duplicate removal

    Count connected components of graph.

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values

    Returns
    ----------
    int
        Number of UMI groups

    """
    def get_adj_list_cluster(umis):
        return {umi: [umi2 for umi2 in umis if edit_dist(umi, umi2) <= mismatches] for umi in umis}

    def get_connected_components_cluster(graph, molecular_barcodes):
        found = list()
        components = list()
        for node in sorted(graph, key=lambda x: molecular_barcodes[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)
        return components

    adj_list = get_adj_list_cluster(list(molecular_barcodes.keys()))
    clusters = get_connected_components_cluster(adj_list, molecular_barcodes)
    return len(clusters)


def dedup_adj(molecular_barcodes, mismatches=1):
    """
    
    Adjecncy duplicate removal

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values

    Returns
    ----------
    int
        Number of UMI groups

    """
    def get_adj_list_adjacency(umis):
        return {umi: [umi2 for umi2 in umis if edit_dist(umi, umi2) <= mismatches] for umi in umis}

    def get_connected_components_adjacency(graph, molecular_barcodes):
        found = list()
        components = list()
        for node in sorted(graph, key=lambda x: molecular_barcodes[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)
        return components

    def remove_umis(adj_list, cluster, nodes):
        '''removes the specified nodes from the cluster and returns
        the remaining nodes '''
        # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
        nodes_to_remove = set([node
                               for x in nodes
                               for node in adj_list[x]] + nodes)
        return cluster - nodes_to_remove

    def get_best_adjacency(cluster, adj_list, counts):
        if len(cluster) == 1:
            return list(cluster)
        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)
        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i + 1])) == 0:
                return sorted_nodes[:i + 1]

    def reduce_clusters_adjacency(adj_list, clusters, counts):
        # TS - the "adjacency" variant of this function requires an adjacency
        # list to identify the best umi, whereas the other variants don't
        # As temporary solution, pass adj_list to all variants
        n = 0
        for cluster in clusters:
            parent_umis = get_best_adjacency(cluster, adj_list, counts)
            n += len(parent_umis)
        return n

    adj_list = get_adj_list_adjacency(list(molecular_barcodes.keys()))
    clusters = get_connected_components_adjacency(adj_list, molecular_barcodes)
    count = reduce_clusters_adjacency(adj_list, clusters, molecular_barcodes)
    return count


def dedup_dir_adj(molecular_barcodes, mismatches=1):
    """
    
    Directed adjecncy duplicate removal

    Parameters
    ----------
    molecular_barcodes : dict
        dictionary with UMIs as keys and UMI counts as values

    Returns
    ----------
    int
        Number of UMI groups

    """
    def get_adj_list_directional_adjacency(umis, counts):
        return {umi: [umi2 for umi2 in umis if edit_dist(umi, umi2) <= mismatches and
                      counts[umi] >= (counts[umi2] * 2) - 1] for umi in umis}

    def get_connected_components_adjacency(graph, molecular_barcodes):
        found = list()
        components = list()
        for node in sorted(graph, key=lambda x: molecular_barcodes[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)
        return components

    def remove_umis(adj_list, cluster, nodes):
        '''removes the specified nodes from the cluster and returns
        the remaining nodes '''
        # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
        nodes_to_remove = set([node
                               for x in nodes
                               for node in adj_list[x]] + nodes)
        return cluster - nodes_to_remove

    def reduce_clusters_directional_adjacency(adj_list, clusters, counts):
        n = 0
        for cluster in clusters:
            n += 1
        return n

    adj_list = get_adj_list_directional_adjacency(list(molecular_barcodes.keys()), molecular_barcodes)
    clusters = get_connected_components_adjacency(adj_list, molecular_barcodes)
    count = reduce_clusters_directional_adjacency(adj_list, clusters, molecular_barcodes)
    return count