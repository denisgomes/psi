"""An object hierarchy data structure for primitives and complex primitives
used for culling and selection optimization.
"""

import sys
from itertools import izip

from psi.gui.bounds import AABB


sys.setrecursionlimit(1000)


class BVHNode(object):

    __slots__ = ("parent", "left", "right", "prims", "bound", "height")

    def __init__(self):
        self.parent = None
        self.left = None
        self.right = None
        self.prims = None   # list for leaf nodes
        self.bound = AABB()

    def is_leaf(self):
        return self.left is None and self.right is None


class BVHTree(object):
    """A bounding volume hierarchy with different build methods"""

    def __init__(self, min_prims=2, max_prims=5, bins=8):
        self._min_prims = min_prims     # merge criteria
        self._max_prims = max_prims     # split criteria
        self._bins = bins
        self._root = BVHNode()

    @property
    def root(self):
        return self._root

    @classmethod
    def build(cls, prims, mode=2, **kwargs):
        mode_map = {1: cls.median_build, 2: cls.binned_sah_build}
        return mode_map[mode](prims, **kwargs)

    @classmethod
    def median_build(cls, prims, min_prims=5, max_levels=5, **kwargs):
        """A recursive top down median split BVH tree builder that terminates
        if a minimum number of primitives is reached for a node or if the max
        level is reached.
        """
        self = cls(**kwargs)
        self._root = root = BVHNode()
        root_bound_merge = root.bound.merge
        for prim in prims:
            root_bound_merge(prim.bound)
        self._median_build(prims, root, level=0, min_prims=min_prims,
                           max_levels=max_levels)
        return self

    def _median_build(self, prims, node, level=0, min_prims=5, max_levels=5):
        if len(prims) <= min_prims or level == max_levels:  # terminate
            node.prims = []     # leaf node
            node.prims.extend(prims)
        else:
            size = node.bound.size()
            axis = size.index(max(size))  # largest axis
            value = node.bound.center()[axis]
            node.left = left_node = BVHNode()
            node.right = right_node = BVHNode()
            left_node.parent = right_node.parent = node
            left_prims = []
            right_prims = []

            # split and merge bounds
            for prim in prims:
                if prim.bound.center()[axis] < value:
                    left_prims.append(prim)
                    left_node.bound.merge(prim.bound)
                else:
                    right_prims.append(prim)
                    right_node.bound.merge(prim.bound)

            # next recursive call
            level += 1
            self._median_build(left_prims, left_node, level=level,
                               min_prims=min_prims, max_levels=max_levels)
            self._median_build(right_prims, right_node, level=level,
                               min_prims=min_prims, max_levels=max_levels)

    @classmethod
    def binned_sah_build(cls, prims, bins=16, min_prims=5, max_levels=5,
                         **kwargs):
        """Use the SAH (surface area heuristics) method to build the bounding
        volume hierarchy.  Based on the paper "On Fast Construction of
        SAH-based Bounding Volume Hierarchies" by Ingo Wald.

        1. Find longest axis
        2. Divide axis into bins with k cut planes.
        3. For each primitive
            a. Find number of prims in respective bin
            b. Grow bin AABB by merging prim AABBs
        4. Plane extraction
            a. Determine area and number of prims in left and right half
            b. Cost = Al*Nl + Ar*Nr
            c. Find lowest cost cutting plane using cost func
        5. Combine prims on the left and right of lowest cost plane
        6. Repeat until termination* criteria is reached

        * Termination
            a. A min number of prims is reached per node
            b. A maximum tree depth is reached
        """
        self = cls(**kwargs)
        self._root = root = BVHNode()
        root_bound_merge = root.bound.merge
        for prim in prims:
            root_bound_merge(prim.bound)
        self._binned_sah_build(prims, bins=bins, node=root,
                                min_prims=min_prims, max_levels=max_levels)
        return self

    def _binned_sah_build(self, prims, node, bins=16, min_prims=5,
                          max_levels=5, level=0):
        if len(prims) <= min_prims or level == max_levels:  # terminate
            node.prims = []     # leaf node
            node.prims.extend(prims)
        else:
            size = node.bound.size()
            value = max(size)
            axis = size.index(value)
            k = bins - 1
            eps = 0.1
            B = [AABB() for _ in range(bins)]   # bin bounds
            P = [list() for _ in range(bins)]   # bin prims
            k1 = bins*(1 - eps) / value
            k0 = node.bound.min[axis]
            for prim in prims:  # project prims into bins
                bid = int(k1 * (prim.bound.center()[axis]-k0))
                B[bid].merge(prim.bound)    # bin bounds
                P[bid].append(prim)     # bin prims

            # calculate cost equation terms for each plane
            nl, nr, al, ar = 0, 0, 0, 0
            NL, AL, NR, AR = [], [], [], []     # prim count and area
            rng = range(k)
            rev_rng = reversed(range(1, bins))
            for l, r in izip(rng, rev_rng):
                # count and merge left to right
                nl += len(P[l])
                al += B[l].area()
                NL.append(nl)
                AL.append(al)
                # right to left
                nr += len(P[r])
                ar += B[r].area()
                NR.append(nr)
                AR.append(ar)

            # calculate lowest cost plane - minimize cost func
            # minimize surface area to get tigher bounds
            plane = 0
            cost = NL[0]*AL[0] + NR[-1]*AR[-1]  # NR/AR - right to left
            for pl in range(1, k):
                new_cost = NL[pl]*AL[pl] + NR[-1-pl]*AR[-1-pl]
                if new_cost < cost:
                    cost = new_cost
                    plane = pl

            node.left = left_node = BVHNode()
            node.right = right_node = BVHNode()
            left_node.parent = right_node.parent = node
            left_prims = []
            right_prims = []

            # merge prims and bounds on left and right side of plane
            # go right to left since list.pop() is faster than list.pop(0)
            for i in reversed(range(k)):
                right_prims.extend(P.pop())
                right_node.bound.merge(B.pop())
                if i == plane:
                    break
            for prims, bounds in izip(P, B):
                left_prims.extend(prims)
                left_node.bound.merge(bounds)

            # next recursive call
            level += 1
            self._binned_sah_build(left_prims, left_node, bins=bins,
                                   min_prims=min_prims, max_levels=max_levels,
                                   level=level)
            self._binned_sah_build(right_prims, right_node, bins=bins,
                                   min_prims=min_prims, max_levels=max_levels,
                                   level=level)

    def prim_count(self, node=None):
        """Given a node, recursively determine the total number of primitives
        in the leaf nodes of node's subtree. The default is to return the
        total number of primitives in the tree.
        """
        if node is None:
            node = self._root

        if node.is_leaf():
            return len(node.prims)

        return self.prim_count(node.left) + self.prim_count(node.right)

    def leaf_count(self, node=None):
        """Given a node, recursively determine the total number of leaf nodes
        in the node's subtree. The default is to return the total number of
        leaf nodes in the tree.
        """
        if node is None:
            node = self._root

        if node.is_leaf():
            return 1

        return self.leaf_count(node.left) + self.leaf_count(node.right)

    def node_count(self, node=None):
        """Given a node, recursively determin the total number of nodes
        (including leaf nodes) in the node's subtree. The default is to
        return the total number of nodes in the tree.
        """
        if node is None:
            node = self._root

        if node.is_leaf():
            return 1

        return self.node_count(node.left) + 1 + self.node_count(node.right)

    def _prim_get_leaf(self, prim):
        """Determine the leaf node the primitive belongs to"""
        for leaf in self.iter_leaves():
            if prim in leaf.prims:
                return leaf

    def iter_prims(self, node=None):
        """Iterate through all the primtives in the tree.  Note that only the
        leaf nodes contain primitives.
        """
        if node is None:
            node = self._root

        if node.left and node.right:
            for left in self.iter_prims(node.left):
                yield left
            for right in self.iter_prims(node.right):
                yield right

        if node.is_leaf():
            for prim in node.prims:
                yield prim

    def iter_nodes(self, node=None):
        """Iterate through all nodes, including leaf nodes"""
        if node is None:
            node = self._root

        yield node

        if node.left and node.right:
            for left in self.iter_nodes(node.left):
                yield left

            for right in self.iter_nodes(node.right):
                yield right

    def iter_leaves(self, node=None):
        """Iterate through leaf nodes only"""
        if node is None:
            node = self._root

        if node.left and node.right:
            for left in self.iter_leaves(node.left):
                yield left

            for right in self.iter_leaves(node.right):
                yield right

        if node.is_leaf():
            yield node

    def get_picked_set(self, ray, node=None):
        """Selection culling"""
        if node is None:
            node = self.root

        for leaf in self.iter_leaves(node):

            if ray.hits_bound(leaf):    # ray-treenode intersection

                for prim in leaf.prims:
                    if ray.hits_prim(prim):
                        yield prim
