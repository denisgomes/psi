# Copyright (c) 2019 Denis Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
#    * Neither the name of Pipe Stress Infinity (PSI) nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""A WireRep data struture stores the centerline topology of the model"""

from psi.utils.euclid import Vector3
from psi.utils.bezier import bezier
from psi import units


@units.define(x="length", y="length", z="length")
class Vertex(object):
    def __init__(self, x=0, y=0, z=0):
        # self.co = [0, 0, 0]
        self.x = x
        self.y = y
        self.z = z
        self.edges = []
        self.data = {}

    @property
    def co(self):
        return self.x, self.y, self.z

    @co.setter
    def co(self, args):
        self.x, self.y, self.z = args


class Edge(object):

    def __init__(self):
        """Create an edge which goes from v1 to v2"""
        self.v1 = None
        self.v2 = None
        self.curve = None
        self.data = {}

    @property
    def length(self):
        return (Vector3(*self.v2.co)-Vector3(*self.v1.co)).magnitude()

    def shared_vert(self, other):
        """Return the shared vertex between self and other edges"""
        if self.v1 in (other.v1, other.v2):
            return self.v1
        elif self.v2 in (other.v1, other.v2):
            return self.v2

    def other_vert(self, v):
        if v is self.v1:
            return self.v2
        elif v is self.v2:
            return self.v1

    def has_vert(self, v):
        if self.v1 is v or self.v2 is v:
            return True
        else:
            return False


class Curve(object):
    """Bezier curves are used to represent a conic section such as a circle or
    an arc.  The curve representation of a circle is generated by stringing
    multiple curves together, where each curve is approximated by one or more
    edges.  When a curve is first created it is a simple edge.  By adding more
    control points, curves of varing complexity can be approximated.  An arc or
    bend for example can be approximated using a quadratic bezier curve with
    three control points or a cubic bezier curve with four control points.  For
    the quadratic case, the control points will be the bend near, elbow and far
    points.

    Bezier curves can also be rational in which case they are capable of
    representing a curve exactly.  A quadratic rational bezier curve can be
    used to represent a quarter circle using a weight value for the control
    point equal to cos(theta/2).  Rational bezier curves are like nurbs curves
    because a weight factor has to be provided for all the control points.
    When we are dealing with regular bezier curves, a quadratic approximation
    for a quater circle does not give accurate enough results for an arc
    length even if a large number of segments are used to subdivide the curve
    equation.  Therefore a cubic bezier curve is generally used to represent
    an arc with sufficient accurary and four cubic bezier curves are used to
    model a circle.

    Since accurary is very in CAD models, a rational bezier curve will be used
    as the mathematical model for a curve object.  As more and more edges are
    added a rational bezier curve will represent the underlying curve exactly.
    """

    _ndiv = 10  # number of approximating edges
    _max_circle_angle = 90  # max angle for arc before splitting

    def __init__(self):
        self.start_vert = None
        self.end_vert = None
        self.control_points = []    # list of lists
        self.weights = []   # weight for each control point
        self.edges = []      # approximating edges
        self.data = {}

    @property
    def length(self):
        """Returns the approximate arc-length of the curve"""
        length = 0
        for edge in self.edges:
            length += edge.length
        return length

    def approximate(self, ndiv):
        """Generate a list of points that approximate the exact bezier curve
        given by the control points and weights.  Ndiv is the number of edges
        that are used.
        """
        ts = [t/ndiv for t in range(ndiv+1)]
        bz = bezier(self._get_control_points(), self._get_weights())
        pnts = bz(ts)
        return pnts

    def _get_control_points(self):
        start_ctrl_pnt = [self.start_vert.co]   # start and end points added
        end_ctrl_pnt = [self.end_vert.co]       # to list
        return start_ctrl_pnt + self.control_points + end_ctrl_pnt

    def _get_weights(self):
        start_pnt_weight = [1]  # start and end weights are always 1
        end_pnt_weight = [1]
        return start_pnt_weight + self.weights + end_pnt_weight

    @property
    def vertices(self):
        verts = []
        for edge in self.edges:
            verts.extend([edge.v1, edge.v2])
        return verts


class Geometry(object):
    """A centerline wireframe representation of the piping model"""

    def __init__(self):
        self.vertices = {}
        self.edges = {}
        self.curves = {}

    def _new_id(self, entity):
        id = 1
        while id in entity:
            id += 1
        return id

    def _get_id(self, entity, entities):
        """Get the id of an entity"""
        for id, obj in entities.iteritems():
            if obj is entity:
                return id

    def MV(self, id=None, x=0, y=0, z=0):
        """Create a vertex"""
        if id is None:
            id = self._new_id(self.vertices)

        vert = self.vertices.get(id, None)
        assert vert is None, 'vertex exists'

        vert = Vertex(x, y, z)
        # vert.co[:] = x, y, z
        self.vertices[id] = vert
        return vert

    def KV(self, v):
        """Kill-Vertex"""
        assert not v.edges, "vertex attached to an edge, first delete the edge"
        v.data.clear()
        v.data = None
        id = self._get_id(v, self.vertices)
        del self.vertices[id]

    def ME(self, v1, v2):
        """Create a directed edge going from v1 to v2"""
        id = self._new_id(self.edges)
        edge = Edge()
        edge.v1, edge.v2 = v1, v2
        self.edges[id] = edge
        v1.edges.append(edge)   # add to disk cycle
        v2.edges.append(edge)
        return edge

    def KE(self, e):
        """Kill-Edge"""
        e.v1.edges.remove(e)
        e.v2.edges.remove(e)
        e.v1 = e.v2 = None
        e.data.clear()
        e.data = None
        id = self._get_id(e, self.edges)
        del self.edges[id]

    def SEMV(self, e, ref_v=None, length=None):
        """Split-Edge-Make-Vertex

        Split an edge and create a vertex in the process.  The default
        reference vertex is the first vertex of the edge. If the length
        parameter is not specified, the edge is split in half.
        """
        if ref_v:
            v1 = ref_v
            v2 = e.other_vert(v1)
        else:
            v1 = e.v1
            v2 = e.v2
        rv1 = Vector3(*v1.co)
        rv2 = Vector3(*v2.co)
        if not length:
            length = (rv2-rv1).magnitude() / 2
        rvn = (rv2-rv1).normalized() * length   # wrt to rv1
        rn = rvn + rv1  # rvn wrt global
        vn = self.MV(None, *rn[:])

        # delete edge but first save the data it contains
        data = e.data
        self.KE(e)

        # create the new edges and copy the data
        en = self.ME(v1, vn)    # create two new edges
        ek = self.ME(vn, v2)
        en.data = data.copy()   # copy edge data before deleting
        ek.data = data.copy()

        return en, ek, vn

    def JEKV(self, e1, e2):
        """Join-Edge-Kill-Vertex"""
        vs = e1.shared_vert(e2)
        v1 = e1.other_vert(vs)
        v2 = e2.other_vert(vs)
        self.KE(e1)
        self.KE(e2)
        self.KV(vs)
        # TODO: copy e1's data but note that e2's data might be different
        # consider loop direction to determine which e data to copy?
        en = self.ME(v1, v2)
        en.data = e1.data.copy()
        return en

    def MCu(self, e):
        """Make a curve given an edge"""
        id = self._new_id(self.curves)
        curve = Curve()
        curve.start_vert, curve.end_vert = e.v1, e.v2
        curve.edges.append(e)
        self.curves[id] = curve
        e.curve = curve
        e.data = curve.data
        return curve

    def KC(self, curve):
        """Kill-Curve"""
        del curve.control_points[:]
        del curve.weights[:]

        # delete all the approximating edges, do not delete the outer verts
        start_end_verts = curve.start_vert, curve.end_vert
        for edge in curve.edges:
            v1, v2 = edge.v1, edge.v2
            self.KE(edge)
            if v1 not in start_end_verts:
                self.KV(v1)
            if v2 not in start_end_verts:
                self.KV(v2)
        del curve.edges[:]
        curve.start_vert, curve.end_vert = None, None

        id = self._get_id(curve, self.curves)
        del self.curves[id]

        self.ME(*start_end_verts)   # create an edge in place of curve

    def ACPCu(self, curve, control_point, weight=1):
        """Add a control point to a curve"""
        curve.control_points.append(control_point)
        curve.weights.append(weight)

    def RCPCu(self, curve, control_point):
        """Remove a control point from a curve"""
        curve.control_points.remove(control_point)

    def UCP(self, curve, old_control_point, new_control_point):
        """Update a control point"""
        idx = curve.control_points.index(old_control_point)
        curve.control_points[idx][:] = new_control_point

    def UCu(self, curve, ndiv=25):
        """Update a curve.  Used after a control point is added, removed or
        mutated.
        """
        # delete all the approximating edges, do not delete the outer verts
        start_end_verts = curve.start_vert, curve.end_vert
        for edge in curve.edges:
            v1, v2 = edge.v1, edge.v2
            self.KE(edge)
            if v1 not in start_end_verts:
                self.KV(v1)
            if v2 not in start_end_verts:
                self.KV(v2)
        del curve.edges[:]

        # add edges and verts again
        if ndiv is None:
            ndiv = Curve._ndiv
        pnts = curve.approximate(ndiv)
        v1 = curve.start_vert
        for i in range(1, ndiv+1):
            if i == ndiv:
                v2 = curve.end_vert
            else:
                v2 = self.MV(None, *pnts[i])
            e = self.ME(v1, v2)
            e.curve = curve
            e.data = curve.data
            curve.edges.append(e)
            v1 = v2
