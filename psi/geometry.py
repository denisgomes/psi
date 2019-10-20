"""A generic centerline model of the piping system used to extract geometry
data.
"""


from psi.vector import Vec3


class Vertex:
    """A point in 3 dimensional space"""

    __slots__ = "coords"

    def __init__(self, x, y, z):
        self.coord = Vec3(x, y, z)

    @property
    def x(self):
        return self.coord.x

    @property
    def y(self):
        return self.coord.y

    @property
    def z(self):
        return self.coord.z


class Edge:
    """A line formed by joining two vertices.

    The direction of the edge is implicitly defined from v1 to v2.
    """

    __slots__ = ("v1", "v2")

    def __init__(self, v1, v2):
        self.v1 = v1
        self.v2 = v2

    def length(self):
        """Get the length of the edge"""
        r21 = self.v2.coord - self.v1.coord

        return r21.norm()


class Curve:
    pass


class WireRep:

    def __init__(self):
        self.vertices = {}
        self.edges = {}
        self.curves = {}

    def MV(self, x, y, z, vid=None):
        """Make a vertex"""
        v = Vertex(x, y, z)

        if vid is None:
            vid = 1
            while True:
                if vid in self.vertices:
                    vid += 1
                else:
                    break

        self.vertices[vid] = v

        return v

    def ME(self, v1, v2):
        """Make an edge by specifying two vertices"""
        e = Edge(v1, v2)

        return e

    def MC(self, edge):
        """Specific an edge to be a curve"""
        pass
