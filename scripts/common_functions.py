# Standard Library imports
from __future__ import annotations
from collections import defaultdict
from typing import TypeVar, Hashable, cast

# External Imports
import numpy as np
import pandas as pd
from scipy import stats

# Local Imports


# Define the type for the nodes
T = TypeVar("T", bound=Hashable)


def find_connected_components[T](
    node_list: list[T], edge_list: list[tuple[T, T]]
) -> list[set[T]]:
    """
    Find the connected components of a graph

    Parameters
    ----------
    node_list : list of Hashable
        List of the nodes which are in the graph (nodes must be hashable)
    edge_list : list of tuple of Hashable
        List of edges in the graph

    Returns
    -------
    connected_components : list of sets of Hashable
        List of the connected components of the graph
    """
    # Create a parents dict
    parents = {node: node for node in node_list}
    # For each edge, ensure that the nodes have a common parent
    for n1, n2 in edge_list:
        try:
            p1 = _get_parent(n1, parents)
            p2 = _get_parent(n2, parents)
        except KeyError as e:
            raise ValueError(
                "Edge list contains nodes not found in node list!"
            ) from e
        if p1 != p2:
            common_parent = min(p1, p2)  # Using min to select a single parent
            parents[n1] = common_parent
            parents[n2] = common_parent
    # The parents for every node can now be evaluated
    components: defaultdict[T, set[T]] = defaultdict(set)
    for node in node_list:
        p = _get_parent(node, parents)
        components[p].add(node)
    # Return the list of the component sets
    return list(components.values())


def _get_parent(node: T, parents: dict[T, T]) -> T:
    """
    Get the (root) parent of a node

    Parameters
    ----------
    node : str
        The node to determine the parent of
    parents : dict of str to str
        Dictionary describing the parent of each node
    """
    if node == parents[node]:
        return node
    else:
        parents[node] = _get_parent(parents[node], parents)
        return parents[node]


def find_degree(
    node_list: list[T], edge_list: list[tuple[T, T]]
) -> dict[T, int]:
    """
    Find the degree of nodes within a graph

    Parameters
    ----------
    node_list : list of Hashable
        List of the nodes which are in the graph (nodes must be hashable)
    edge_list : list of tuple of Hashable
        List of edges in the graph

    Returns
    -------
    degree_dict : dict of Hashable to int
        Dictionary keyed by node, with value corresponding to degree of
        the node
    """
    degree: defaultdict[T, int] = defaultdict(int)
    for n1, n2 in edge_list:
        degree[n1] += 1
        degree[n2] += 1
    for n in node_list:
        degree[n] += 0
    return degree


# Create a function for getting the neighbors of nodes in the graph
def find_neighbors(
    node_list: list[T], edge_list: list[tuple[T, T]]
) -> dict[T, set[T]]:
    """
    Find the neighbors of nodes within a graph

    Parameters
    ----------
    node_list : list of Hashable
        List of the nodes which are in the graph (nodes must be hashable)
    edge_list : list of tuple of Hashable
        List of edges in the graph

    Returns
    -------
    neighors : dict of Hashable to set of Hashable
        A dictionary keyed by node, with values of sets of neighbors to
        the key node
    """
    neighbors: dict[T, set[T]] = defaultdict(set)
    for n1, n2 in edge_list:
        neighbors[n1].add(n2)
        neighbors[n2].add(n1)
    for n in node_list:
        # This just forces an empty set to be created if
        # it wasn't already created during edge iteration
        _ = neighbors[n]
    return neighbors


# Create a function for finding the representative nodes in a graph
# with connected components
def find_representative_nodes(
    node_list: list[T], edge_list: list[tuple[T, T]]
) -> dict[T, set[T]]:
    """
    Find representative nodes in a graph by iteratively selecting nodes of
    high degree from seperate components
    """
    # Create a dict to represent the representative nodes
    representative_nodes: dict[T, set[T]] = defaultdict(set)
    # Create a set of nodes to process,
    # the process will continue until all of the
    # nodes are assigned to a representative node
    nodes_to_process: set[T] = set(node_list)
    # Get a copy of the edges list to update so that
    # the original list isn't modified
    # This can be shallow, as the items aren't changed, only the
    # list itself
    edge_list: list[tuple[T, T]] = edge_list.copy()
    while nodes_to_process:
        node_list = list(nodes_to_process)
        edge_list = cast(list[tuple[T, T]], edge_list)
        components = find_connected_components(node_list, edge_list)
        neighbors = find_neighbors(node_list, edge_list)
        # Use len to get degree so that we aren't iterating twice
        # through the edges list
        degree = {
            node: len(neighboring_nodes)
            for node, neighboring_nodes in neighbors.items()
        }
        # Create a set of nodes which will be removed from the nodes to process set
        nodes_to_remove: set[T] = set()
        # For each component, choose a representative node
        # representing the node with the highest degree in the
        # component
        for cur_component in components:
            rep_node = max(cur_component, key=lambda it: degree[it])
            representative_nodes[rep_node] = neighbors[rep_node]
            nodes_to_remove.add(rep_node)
            nodes_to_remove.update(neighbors[rep_node])
        # Remove the nodes to remove from the nodes_to_process set
        nodes_to_process -= nodes_to_remove
        # Remove edges which include a node to remove from the edges list
        edge_list = list(
            filter(
                lambda e: (
                    not ((e[0] in nodes_to_remove) or (e[1] in nodes_to_remove))
                ),
                edge_list,
            )
        )
    return representative_nodes


def rpkm_to_tpm(rpkm: pd.DataFrame):
    """Convert RPKM normalized counts to TPM normalized counts

    Parameters
    ----------
    rpkm : pd.DataFrame
        RPKM normalized count data, with genes as columns and samples as
        rows

    Returns
    -------
    pd.DataFrame
        TPM normalized counts
    """
    return rpkm.divide(rpkm.sum(axis=1), axis=0) * 1.0e6


def fdr_with_nan(ps: pd.Series, **kwargs) -> pd.Series:
    """
    Perform false discovery correction on a series of p-values which may contain
    missing values

    Parameters
    ----------
    ps : pd.Series
        Pandas series containing p-values to correct for false discovery
    kwargs
        Key word arguments passed to scipy.stats.false_discovery_control

    Returns
    -------
    corrected_ps : pd.Series
        Series with corrected ps, with the same index as the input series
    """
    corrected_ps = pd.Series(np.nan, index=ps.index)
    no_nan_series = ps.dropna(inplace=False)
    corrected_ps[no_nan_series.index] = stats.false_discovery_control(
        no_nan_series, **kwargs
    )
    return corrected_ps
