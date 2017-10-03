def bfs(G, s):
    """Perform Breadth-first search.

    Parameters
    ----------
    G : dict
        a graph with nodes as keys with values as a set of destination nodes
    s : tuple
        node to start BFS

    Returns
    -------
    seen : set
        set of seen nodes while using BFS
    """
    seen = set([s])
    current = G[s]
    while current:
        tmp_current = set()
        for node in current:
            if node not in seen:
                seen.add(node)
                tmp_current |= G[node] - seen
        current = tmp_current
    return seen


def connected_components(G):
    """Find the connected components of a graph.

    Parameters
    ----------
    G : dict
        graph of residue neighbors

    Returns
    -------
    component_list : list
        list of connnected components in G
    """
    all_nodes = G.keys()
    seen = set()
    component_list = []
    for node in all_nodes:
        if node in seen:
            continue
        else:
            component = bfs(G, node)
            component_list.append(component)
            seen |= component
    return component_list
