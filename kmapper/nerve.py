import itertools
from collections import defaultdict
import numpy as np

__all__ = ["GraphNerve", "SimplicialNerve"]


class Nerve:
    """Base class for implementations of a nerve finder to build a Mapper complex."""

    def __init__(self):
        pass

    def compute(self, nodes, links):
        raise NotImplementedError()


class GraphNerve(Nerve):
    """Creates the 1-skeleton of the Mapper complex.

    Parameters
    -----------

    min_intersection: int, default is 1
        Minimum intersection considered when computing the nerve. An edge will be created only when the intersection between two nodes is greater than or equal to `min_intersection`
    """

    def __init__(self, min_intersection=1):
        self.min_intersection = min_intersection

    def __repr__(self):
        return "GraphNerve(min_intersection={})".format(self.min_intersection)

    def compute(self, nodes):
        """Helper function to find edges of the overlapping clusters.

        Parameters
        ----------
        nodes:
            A dictionary with entires `{node id}:{list of ids in node}`

        Returns
        -------
        edges:
            A 1-skeleton of the nerve (intersecting  nodes)

        simplicies:
            Complete list of simplices

        """

        result = defaultdict(list)

        # Create links when clusters from different hypercubes have members with the same sample id.
        candidates = itertools.combinations(nodes.keys(), 2)
        for candidate in candidates:
            # if there are non-unique members in the union
            if (
                len(set(nodes[candidate[0]]).intersection(nodes[candidate[1]]))
                >= self.min_intersection
            ):
                result[candidate[0]].append(candidate[1])

        edges = [[x, end] for x in result for end in result[x]]
        simplices = [[n] for n in nodes] + edges
        return result, simplices


class SimplicialNerve(Nerve):
    """Creates the entire Cech complex of the covering defined by the nodes.

    Parameters
    -----------

    min_intersection: int, default is 1
        Minimum intersection considered when computing the nerve. An edge will be created only when the intersection between two nodes is greater than or equal to `min_intersection`   

    max_dimension: int, default is 1
        Maximum dimension of simplicies generated. A simplex will only be created if it has dimension less than or equal to `max_dimension`

    cubes_shape: tuple, default is none
        Optional parameter showing the arrangement of cubes in the cover, enabling more efficient simplex construction.

    """
    def __init__(self, min_intersection=1, max_dimension=1, cubes_shape: tuple=None):
        self.min_intersection = min_intersection
        self.max_dimension = max_dimension
        self.cubes_shape = cubes_shape

    def __repr__(self):
        return "SimplicialNerve(min_intersection={}, max_dimension={})".format(self.min_intersection, self.max_dimension)
    
    def _compute_cluster_intersection(self, nodes: dict, clusters: list):
        intersection = set(nodes[clusters[0]])
        for i in range(1,len(clusters)):
            intersection = intersection.intersection(nodes[clusters[i]])
        return intersection
    
    def _compute_neighbor_cubes(self, cube: int):
        neighbors = self._compute_neighbors_recurse(cube, self.cubes_shape)
        return neighbors
    
    def _compute_neighbors_recurse(self, cube:int, cubes_shape:dict):
        result = [cube]
        for idx,size in enumerate(cubes_shape):
            result += self._compute_neighbors_recurse(cube+size**idx, cubes_shape[:idx])
            result += self._compute_neighbors_recurse(cube-size**idx, cubes_shape[:idx])
        return result

    def _get_nodes_in_surrounding_of_cube(self, cube_id:int, nodes_in_cubes: list):
        num_cubes = self._num_cubes()
        neighbors = self._compute_neighbor_cubes(cube_id)
        comb_nodes = []
        for neighbor in neighbors:
            if(neighbor < num_cubes and neighbor >= 0):
                comb_nodes += nodes_in_cubes[neighbor]
        return comb_nodes

    def _add_applicable_candidates(self, simplices, candidates, nodes):
        for candidate in candidates:
            intersection = self._compute_cluster_intersection(nodes, candidate)
            if(len(intersection) >= self.min_intersection):
                simplices.append(candidate)

    def _compute_n_dimensional_simplices_sectioned(self, dimension: int,nodes:dict, nodes_in_cubes: list):
        simplices = []
        num_simplex_nodes = dimension + 1
        for cube_id, cube in enumerate(nodes_in_cubes):
            comb_nodes = self._get_nodes_in_surrounding_of_cube(cube_id, nodes_in_cubes)
            candidates = itertools.combinations(comb_nodes, num_simplex_nodes)
            self._add_applicable_candidates(simplices, candidates, nodes)
        return simplices
    
    def _compute_n_dimensional_simplices_trivial(self, dimension:int, nodes: dict):
        simplices = []
        num_simplex_nodes = dimension + 1
        candidates = itertools.combinations(nodes.keys(), num_simplex_nodes)
        self._add_applicable_candidates(simplices, candidates, nodes)
                
        return simplices
       

    def _compute_n_dimensional_simplices(self, dimension: int, nodes: dict, nodes_in_cubes: dict=None):
        if self.cubes_shape and nodes_in_cubes:
            return self._compute_n_dimensional_simplices_sectioned(dimension, nodes, nodes_in_cubes)
        else: 
            return self._compute_n_dimensional_simplices_trivial(dimension, nodes)
            
    def _num_cubes(self):
        assert(self.cubes_shape)
        result = 1
        for i in self.cubes_shape:
            result *= i
        return result

    def compute(self, nodes: dict):
        """Helper function to find edges of the overlapping clusters.
        Parameters
        ----------
        nodes:
            A dictionary with entires `{node id}:{list of ids in node}`

        Returns
        -------
        edges:
            A 1-skeleton of the nerve (intersecting  nodes)

        simplicies:
            Complete list of simplices


        """

        nodes_in_cubes = None
        # If there is information about the arrangement of the cubes, it is possible to use more efficient schemes for 
        # computing clusters. This requires the clusters to be accessed via their respective cubes.
        if self.cubes_shape:
            nodes_in_cubes = [[] for _ in range(self._num_cubes())]
            import re
            for node in nodes:
                # RegEx finds the number after "cube" in the string node, which is of format
                # e.g. "cube0_cluster0"
                cube = int(re.search("(?<=cube)\d*", node).group(0))
                nodes_in_cubes[cube].append(node)
        links = defaultdict(list)
        simplices = [[n] for n in nodes]

        # Sample candidates for each desired simplicial dimension        
        for simplicial_dimension in range (1, self.max_dimension+1):
            n_simplices = self._compute_n_dimensional_simplices(simplicial_dimension, nodes, nodes_in_cubes=nodes_in_cubes)
            simplices += n_simplices
            if(simplicial_dimension == 1):
                for edge in n_simplices:
                    links[edge[0]].append(edge[1])
        return links, simplices

