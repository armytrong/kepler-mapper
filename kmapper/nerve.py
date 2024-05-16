import itertools
from collections import defaultdict

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

    Warning: Not implemented yet.
    """
    def __init__(self, min_intersection=1, max_dimension=1):
        self.min_intersection = min_intersection
        self.max_dimension = max_dimension

    def __repr__(self):
        return "SimplicialNerve(min_intersection={})".format(self.min_intersection)
        
    def compute(self, nodes, links=None):
        result = defaultdict(list)
        simplices = [[n] for n in nodes]
        for simplicial_dimension in range (1, self.max_dimension+1):
            simplex_num_nodes = simplicial_dimension + 1
            candidates = itertools.combinations(nodes.keys(), simplex_num_nodes)
            for candidate in candidates:
                intersection = set(nodes[candidate[0]])
                for i in range(1,simplex_num_nodes):
                    intersection = intersection.intersection(nodes[candidate[i]])
                if(len(intersection) >= self.min_intersection):
                    if(simplicial_dimension == 1):
                        result[candidate[0]].append(candidate[1])

                    simplices.append(candidate)
        return result, simplices
