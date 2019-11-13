"""Piping finite elements.

To create a model, the user must first define a series of runs that make the
overall centerline of the model. After this process, additional components can
be added by giving a point and a direction in which to place the new item or by
"converting" a run to a different type of element.

Example
-------
First define the model centerline using Point and Run.

>>> Model("test")
>>> Point(10, 0, 0, 0)
>>> Bend(20, -2, radius=...)    # can be created by a point
>>> Bend(30, 0, 2, radius=...)
>>> Run(40, 0, 0, 2)
>>> Bend(50, 0, 0, 2, radius=...)
>>> Bend(60, 0, -2, radius=...)
>>> Run(70, 2)
>>> Point(80, -10)
>>> Run(90, 0, 2)
>>> Run(100, 0, 0.5)
>>> Point(90)
>>> Run(40)                     # defined twice and loop-closed (1) (2)

Now define the other components by point.

>>> elements(20, 30).split(25)          # split a bend - create run and bend
>>> Valve(25, length, mass, "mid")      # two elements, 3 nodes made
>>> Tee(40, "welding")                  # sif.Tee by point only
>>> Flange(10, length, mass, "down")    # flanged nozzle
>>> Flange(70, length, mass, "up")

There are three possible ways in which a point defined for a component can be
interpreted. If the point is defined as a midpoint, two elements are made, one
on each adjacent run. If the component is defined to be upstream, it will be
placed in the oposite direction to the flow, with respect to the element nodal
direction. If the component is defined to be downstream, it will be placed in
the same direction as the flow. In both the latter two cases, only one element
is created. When a component is defined at a boundary node, it is possible to
create another boundary node by choosing the direction that 'extends' the
model.

Components can also be created by picking a run element and 'converting' it to
a different type of element such as a valve or a flange.

>>> elements(50, 60).split(54, 0.5, ref=50)     # used node 50 as reference
>>> elements(54, 60).split(55, 0.5, ref=54)
>>> elements(54, 55).convert(Valve, mass=...)

Notes
-----
Only runs can be split and merged. To split a valve, it must first be changed
to a run and then split.


TODO
----
Dealing with element weight vs mass when it comes to units.

For Rigid elements the weight must be specified not the mass.

Bend and reducer having multiple run approximation.

How to easily create bends when defining geometry and how to update the
geometry when point coordinates or bend radius is updated.

How to covert different types of elements to runs and vice versa.

Unit conversion should be disabled before the analysis is performed.
"""

from math import cos, pi, sqrt

import numpy as np
import numpy.linalg as la

from psi.entity import Entity, EntityContainer
from psi.units import units
from psi.utils.euclid import Vector3


# TODO: Check app.points to determine if point exists, __new__


class Element(Entity):
    """Create an element in the active model."""

    nn = 2      # number of nodes
    ndof = 6    # nodal degrees of freedom

    def __new__(cls, name, *args, **kwargs):
        # call base class __new__
        inst = super(Element, cls).__new__(cls, name, *args, **kwargs)

        app = cls._app
        model = app.models.active_object
        assert model.active_point, "create or activate a point"

        # applies to all points including branch
        if name in model.geometry.vertices:
            raise NameError("name '%s' in use" % name)

        return inst

    def __init__(self, to_point, from_point, section, material):
        # element name is a tuple consisting of from_point and to_point
        super(Element, self).__init__((from_point, to_point))
        self.geometry = None
        self.section = section
        self.material = material

    def build(self, point):
        """Create and set element geometry"""
        raise NotImplementedError("abstract method")

    @property
    def from_point(self):
        return self.app.points(self.name[0])

    @from_point.setter
    def from_point(self, name):
        raise NotImplementedError("implement")

    @property
    def to_point(self):
        """The element name is used to store this parameter.

        In the case where a loop is closed, two elements will have the same
        name, but they have different from points.
        """
        return self.app.points(self.name[1])

    @to_point.setter
    def to_point(self, name):
        raise NotImplementedError("implement")

    @property
    def parent(self):
        return self.app.elements

    def next(self):
        """Return the next element in the downstream direction"""
        raise NotImplementedError("implement")

    def prev(self):
        """Return the previous element in the upstream direction"""
        raise NotImplementedError("implement")

    def dircos(self):
        """Element direction cosine matrix"""
        raise NotImplementedError("abstract method")

    def T(self):
        """Local dof to global DOF transformation matrix"""
        raise NotImplementedError("abstract method")

    def klocal(self):
        """Local stiffness matrix"""
        raise NotImplementedError("abstract method")

    def kglobal(self):
        """Global element stiffness matrix"""
        raise NotImplementedError("abstract method")

    def mlocal(self):
        """Local mass matrix"""
        raise NotImplementedError("abstract method")

    def mglobal(self):
        """Global element mass matrix"""
        raise NotImplementedError("abstract method")

    def flocal(self, load):
        """Local force matrix"""
        raise NotImplementedError("abstract method")

    def fglobal(self):
        """Global element force matrix"""
        raise NotImplementedError("abstract method")

    def __repr__(self):
        return "%s(%s, %s)" % (self.type, self.from_point.name,
                               self.to_point.name)


class Structural(Element):
    """Structural elements such as wide flanges and tee members"""
    pass


class Piping(Element):
    """Base class for piping elements such as runs and bends."""

    def __init__(self, to_point, from_point=None, section=None, material=None,
                 insulation=None, code=None):
        # active attributes
        model = self.app.models.active_object

        # these are all objects not names
        if from_point is None:
            from_point = model.active_point.name
        if section is None:
            section = model.active_section
        if material is None:
            material = model.active_material
        if insulation is None:
            insulation = model.active_insulation
        if code is None:
            code = model.active_code

        super(Piping, self).__init__(to_point, from_point, section, material)
        self.insulation = insulation
        self.code = code
        self.sifs = set()
        self.loads = set()
        self.loadcases = set()
        self.supports = set()
        model.active_point = to_point

    @classmethod
    def by_point(cls, from_point, to_point, section=None, material=None,
                 insulation=None, code=None):
        """Alternate constructor to create an element by points. This is the
        classical way to define a piping element mesh in which all the points
        can be defined and then all the elements.

        Note the method function signature defines an element direction which
        goes from the from_point to the to_point as opposed to the default way
        of element creation in which only the to_point is required and the
        from_point is assumed to be the active point.
        """
        inst = cls(to_point, from_point, section, material, insulation, code)

        return inst

    def convert(self, element_type, **kwargs):
        """Convert from one type of element to another"""
        raise NotImplementedError("abstract method")

    def mass(self):
        """The pipe weight and in the case of a valve or flange, the respective
        user defined weight.
        """
        mass = (self.material.rho.value * self.section.area *
                self.geometry.length)
        return mass

    def weight(self, accel):
        """Weight of the piping due to the acceleration"""
        return self.mass * accel

    def fluid_mass(self, fluden):
        di = self.section.od - 2*self.section.thke
        fluid_area = (pi/4) * (di**2)

        mass = fluden * fluid_area * self.geometry.length
        return mass

    def fluid_weight(self, fluden, accel):
        return self.fluid_mass(fluden) * accel

    def insulation_mass(self):
        dins = self.section.od + 2*self.insulation.thk
        do = self.section.od
        insulation_area = (pi/4)*(dins**2 - do**2)

        mass = (self.insulation.rho * insulation_area *
                self.geometry.length)
        return mass

    def insulation_weight(self, accel):
        return self.insulation_mass() * accel

    def total_mass(self, fluden):
        """Total weight of element including, pipe, insulation and contents."""
        return self.mass() + self.fluid_mass(fluden) + self.insulation_mass()

    def total_weight(self, fluden, accel):
        """Total weight of element including, pipe, insulation and contents."""
        return (self.weight(accel) + self.fluid_weight(fluden)
                + self.insulation_weight(accel))


@units.define(_dx="length", _dy="length", _dz="length")
class Run(Piping):
    """Define a pipe run by incremental offsets in the global x, y and z."""

    def __init__(self, point, dx, dy=0, dz=0, from_point=None, section=None,
                 material=None, insulation=None, code=None):
        super(Run, self).__init__(point, from_point, section, material,
                                  insulation, code)
        self._dx = dx
        self._dy = dy
        self._dz = dz
        self._split_edge = None     # for bends
        self.build(point)

    def sif(self):
        """Element stress intensification factor per the applicable code"""
        sifi = sifo = 1.0

        return (sifi, sifo)

    def kfac(self):
        """Element flexibility factor per the applicable code"""
        return 1.0

    def build(self, point):
        model = self.app.models.active_object
        elements = self.app.elements
        points = self.app.points

        from_point = self.from_point.name

        # create edge
        from_vert = model.geometry.vertices[from_point]
        to_vert = model.geometry.vertices.get(point, None)
        if to_vert is None:
            # make new vertex
            to_vec = Vector3(*from_vert.co) + Vector3(self.dx, self.dy,
                                                      self.dz)
            pnt = points.Point(point, *to_vec[:])
            # to_vert = model.geometry.MV(point, *to_vec[:])
            edge = model.geometry.ME(from_vert, pnt.vertex)
        else:
            # close the loop
            edge = model.geometry.ME(from_vert, to_vert)
            from_to_vec = Vector3(*to_vert.co[:]) - Vector3(*from_vert.co[:])
            self.dx, self.dy, self.dz = from_to_vec[:]

        # create run geometry
        self.geometry = edge

        # build bend element
        element = elements.objects.get(from_point, None)
        if element and isinstance(element, Bend):
            element.build(from_point)

    @property
    def dx(self):
        return self._dx

    @dx.setter
    def dx(self, value):
        """Change the vertex x-coordinate of the 'to' point"""
        xold = self._dx
        self._dx = value    # set x

        # TODO: Note that geometry cooridnates are managed
        e = self.geometry
        v1, v2 = e.v1, e.v2
        r1 = Vector3(*v1.co)
        r2 = Vector3(*v2.co)
        dx = self._dx       # managed

        x_new = dx - xold
        r3 = r2 - r1
        r3_new = r3 + Vector3(x_new, 0, 0)
        r2_new = r1 + r3_new

        # set value
        v2.co = r2_new[:]

    @property
    def dy(self):
        return self._dy

    @dy.setter
    def dy(self, value):
        """Change the vertex y-coordinate of the 'to' point"""
        yold = self._dy
        self._dy = value    # set y

        e = self.geometry
        v1, v2 = e.v1, e.v2
        r1 = Vector3(*v1.co)
        r2 = Vector3(*v2.co)
        dy = self._dy       # managed

        y_new = dy - yold
        r3 = r2 - r1
        r3_new = r3 + Vector3(0, y_new, 0)
        r2_new = r1 + r3_new

        # set value
        v2.co = r2_new[:]

    @property
    def dz(self):
        return self._dz

    @dz.setter
    def dz(self, value):
        """Change the vertex z-coordinate of the 'to' point"""
        zold = self._dz
        self._dz = value    # set y

        e = self.geometry
        v1, v2 = e.v1, e.v2
        r1 = Vector3(*v1.co)
        r2 = Vector3(*v2.co)
        dz = self._dz       # managed

        z_new = dz - zold
        r3 = r2 - r1
        r3_new = r3 + Vector3(0, 0, z_new)
        r2_new = r1 + r3_new

        # set value
        v2.co = r2_new[:]

    @property
    def length(self):
        dx, dy, dz = self.dx, self.dy, self.dz
        return sqrt(dx**2 + dy**2 + dz**2)

    @length.setter
    def length(self, value):
        """Set the length of the element."""
        e = self.geometry
        v1, v2 = e.v1, e.v2
        r1 = Vector3(*v1.co)
        r2 = Vector3(*v2.co)

        r3 = r2 - r1
        r3u = r3.normalized()   # unit vector
        r_new = r3u * value     # the new dx, dy, dz
        r2_new = r1 + r_new

        # set new value for v2
        dx, dy, dz = r2_new[:]
        self.dx = dx
        self.dy = dy
        self.dz = dz

    # def split(self, point, distance=None, reference=None):
    #     """Split the element by the distance from the fixed reference point.

    #     The default behaviour is to split the element in two and assign a new
    #     point. The original element is resized and a new element is created.
    #     All attributes of the original are copied to the newly created element.

    #     If the point already exists, a new point is defined and used instead.
    #     """
    #     self.parent.split(self, point, distance, reference)

    def dircos(self):
        """Return the direction cosine of the element in matrix form. The unit
        vector components of the local coordinate axes of the element determine
        the 3x3 transformation matrix.

        For straight elements, the local x axis is given by the direction from
        the 'from' point to the 'to' point. The local y is parallel to the
        global vertical direction and the local z is the cross product of the
        local x and y axes.

        An exceptional case is when the pipe runs in the global vertical
        direction. In this case the local x axis is the same as before, however
        it happens to be in the vertical direction. The local y is aligned with
        global x (this is arbitrary) and the local z is the cross product of
        local x with local y.
        """
        up = self.app.models.active_object.settings["core.vertical"]
        if up == "y":
            local_y = np.array([0., 1., 0.])
        elif up == "z":
            local_y = np.array([0., 0., 1.])

        from_vert = np.array(self.geometry.v1.co)
        to_vert = np.array(self.geometry.v2.co)
        local_x = to_vert - from_vert

        # check to see if local_x is parallel to global vertical
        if (1 - np.dot(local_x, local_y)) < 0.01:
            # yes parallel, vertical straight element
            local_y = np.array([1., 0., 0.])

        local_z = np.cross(local_x, local_y)

        # recalculate local y so that its orthogonal
        local_y = np.cross(local_x, local_z)

        return np.array([local_x/la.norm(local_x),
                         local_y/la.norm(local_y),
                         local_z/la.norm(local_z)]
                        )

    def T(self):
        """Local to global transformation matrix"""
        tf = np.empty((12, 12), dtype=np.float32)
        dc = self.dircos()

        tf[:3, :3] = tf[3:6, 3:6] = tf[6:9, 6:9] = tf[9:12, 9:12] = dc[:, :]

        return tf

    def klocal(self, temp, stiffness=None):
        """The local stiffness matrix of a straight piping element, i.e all
        elements with the expection of bends and reducers. The stiffness is a
        function of the temperature. If a stiffness parameter is provided it is
        directly applied to the stiffness matrix.
        """
        kmat = np.empty((12, 12), dtype=np.float32)

        a = self.geometry.length / 2    # for ease below

        E = self.material.eh[temp]
        if stiffness:
            E = stiffness   # overwrite

        nu = self.material.nu.value     # poisson's ratio
        G = E / (2 * (1 + nu))          # shear mod isotropic mats
        A = self.section.area
        J = self.section.ixx
        Iy = self.section.iyy
        Iz = self.section.izz

        kmat[0, 0] = A*E / 2*a
        kmat[0, 6] = kmat[6, 0] = -A*E / 2*a

        kmat[1, 1] = 3*E*Iz / 2*a**3
        kmat[1, 5] = kmat[5, 1] = 3*E*Iz / 2*a**2
        kmat[1, 7] = kmat[7, 1] = -3*E*Iz / 2*a**3
        kmat[1, 11] = kmat[11, 1] = 3*E*Iz / 2*a**2

        kmat[2, 2] = 3*E*Iy / 2*a**3
        kmat[2, 4] = kmat[4, 2] = -3*E*Iy / 2*a**2
        kmat[2, 8] = kmat[8, 2] = -3*E*Iy / 2*a**3
        kmat[2, 10] = kmat[10, 2] = -3*E*Iy / 2*a**2

        kmat[3, 3] = G*J / 2*a
        kmat[3, 9] = -G*J / 2*a

        kmat[4, 4] = 2*E*Iy / a
        kmat[4, 8] = kmat[8, 4] = 3*E*Iy / 2*a**2
        kmat[4, 10] = kmat[10, 4] = E*Iy / a

        kmat[5, 5] = 2*E*Iz / a
        kmat[5, 7] = kmat[7, 5] = -3*E*Iz / 2*a**2
        kmat[5, 11] = kmat[11, 5] = E*Iz / a

        kmat[6, 6] = A*E / 2*a

        kmat[7, 7] = 3*E*Iz / 2*a**3
        kmat[7, 11] = kmat[11, 7] = -3*E*Iz / 2*a**2

        kmat[8, 8] = 3*E*Iy / 2*a**3
        kmat[8, 10] = kmat[10, 8] = 3*E*Iy / 2*a**2

        kmat[9, 9] = G*J / 2*a

        kmat[10, 10] = 2*E*Iy / a

        kmat[11, 11] = 2*E*Iz / a

        return kmat

    def kglobal(self, temp, stiffness=None):
        """The element global stiffness matrix before assembly into the system
        matrix.

        Note that the transpose of T is taken because T is defined in row major
        format, in other words, the local axes components all appear on the
        same row. Taking the transpose, (i.e. inverse, sinces it's a rotation
        system.
        """
        T = self.T()    # build T once and reuse
        return T.transpose() * self.klocal(temp, stiffness) * T


@units.define(_radius="length")
class Bend(Run):
    """Define a pipe bend.

    The bend element is similar to the run in that it has a length to it.  A
    bend also has a radius which is used to create the underlying curve
    topology. A bend is modeled using a rational second order quadratic bezier
    curve.

    In order to build a bend, a valid radius must be provided. The bend and the
    adjacent run must be long enough to accomodate the bend, else the bend is
    not created. A bend cannot be defined where there are more than two
    elements defined, i.e. a tee intersection.

    The bend geometry is not created immediately upon definition but after the
    next element is defined. Initially, a run is created but the curve is built
    when the next adjacent run element is defined. The curve is also
    built/updated when the bend radius is set.

    A bend element is approximated using multiple runs element. The number of
    elements is a user defined parameter.
    """

    def __init__(self, point, dx, dy=0, dz=0, radius="long", tol=0.01,
                 from_point=None, section=None, material=None,
                 insulation=None, code=None):
        # calls build
        super(Bend, self).__init__(point, dx, dy, dz, from_point, section,
                                   material, insulation, code)
        self.runs = []          # internal runs
        self.radius = radius    # also calls build
        self.tol = tol

    def build(self, point):
        model = self.app.models.active_object

        if model.geometry.vertices.get(point, None) is None:
            # on initial super call
            super(Bend, self).build(point)
        else:
            # when radius is set
            assert self.radius != 0, "radius must be nonzero"

            # create a bend
            ctrlvert = model.geometry.vertices[point]

            try:
                # assuming disk cycle with two edges
                e1, e2 = ctrlvert.edges
            except ValueError:
                raise AssertionError("too many elements")

            # determine if lengths are long enough
            if (self.radius > e1.length+self.tol or
                    self.radius > e2.length+self.tol):
                raise ValueError("radius too large")

            # split edge to make room for bend
            verts = []
            for e in (e1, e2):
                elen = e.length
                if (self.radius <= elen+self.tol and
                        self.radius >= elen-self.tol):
                    # no need to split edge
                    v = e.other_vert(ctrlvert)
                    verts.append(v)
                elif self.radius < elen:
                    en, em, v = model.geometry.SEMV(e, ctrlvert,
                                                    self.radius)
                    # note en is closest to control point vert
                    model.geometry.KE(en)
                    verts.append(v)

            # create curve
            v1, v2 = verts
            rv1 = Vector3(*v1.co)
            rv2 = Vector3(*v2.co)
            ec = model.geometry.ME(v1, v2)  # on a diagonal
            curve = model.geometry.MCu(ec)
            model.geometry.ACPCu(curve, ctrlvert.co[:],
                                 cos(rv1.angle(rv2)/2))
            model.geometry.UCu(curve)

            self.geometry = curve

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        """Specify the bend radius of the element.

        If the next element is not defined, the curved element will not be
        created regardless of whether a valid bend radius is provided.
        """
        try:
            if isinstance(value, str):
                if value == "long":
                    value = 1.5 * self.section.od
                elif value == "short":
                    value = 1.0 * self.section.od
                elif value[-1] in ("d", "D"):
                    # "3D", "5d", "3.5D"
                    fac = float(value[:-1])
                    value = fac * self.section.od
        except:
            value = 0

        self._radius = value

        # update geometry - complain if radius is too big
        # TODO: add update code

        # if next element is defined
        model = self.app.models.active_object
        vert = model.geometry.vertices.get(self.to_point, None)
        if vert and len(vert.edges) == 2:
            self.build(self.to_point)

    def klocal(self, temp, stiffness=None):
        for run in self.runs:
            yield run.klocal(temp, stiffness)

    def kglobal(self, temp, stiffness=None):
        """Yield the global bend stiffness matrix for each element.

        The bend consists of multiple approximating elements and so the global
        matrix for the bend is the assembly of all the elements which make up
        the bend.
        """
        for run in self.runs:
            yield run.kglobal(temp, stiffness)

    def __repr__(self):
        return "%s %s" % (self.type, self.to_point)


class Reducer(Run):
    """Define a pipe reducer by giving a different section.

    The new section is activated automatically. A reducer is approximated using
    multiple runs with decreasing diameters.
    """

    def __init__(self, point, section2, dx, dy=0, dz=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Reducer, self).__init__(point, dx, dy, dz, from_point, section,
                                      material, insulation, code)
        section2.activate()     # automatic
        self.section2 = section2
        self.runs = []          # internal runs

    def build(self, point):
        """Similar to a bend a reducer consists of one curve which consists of
        multiple internal edges. Each edge has a different section associated
        with each with a smaller and smaller diameter to approximate the
        reducer geometry.
        """
        raise NotImplementedError("implement")

    def mass(self):
        """Mass of a reducer"""
        d1o = self.section.od
        d1i = self.section.od - 2*self.section.thke
        d2o = self.section2.od
        d2i = self.section2.od - 2*self.section2.thke

        vol = (1/12)*pi*self.geometry.length*(d1o**2+d1o*d2o+d2o**2 -
                                              d1i**2+d1i*d2i+d2i**2)
        return self.material.rho.value * vol

    def fluid_mass(self, fluden):
        """Fluid mass within the reducer"""
        d1i = self.section.od - 2*self.section.thke
        d2i = self.section2.od - 2*self.section2.thke

        fluvol = (1/12)*pi*self.geometry.length*(d1i**2+d1i*d2i+d2i**2)

        return fluden * fluvol

    def insulation_mass(self):
        """A conical shaped insulation profile is considered"""
        s, s2 = self.section, self.section2

        (s, s2) = ((s, s2) if (2*self.insulation.thk+s.od) >
                              (2*self.insulation.thk+s2.od) else
                              (s2, s))

        dt1 = s.od + 2*self.insulation.thk
        dt2 = s2.od + 2*self.insulation.thk
        dr1 = s.od
        dr2 = s2.od

        vt = (1/12) * (dt1**2 + dt1*dt2 + dt2**2)
        vr = (1/12) * (dr1**2 + dr1*dr2 + dr2**2)
        vi = vt - vr

        return s.inden*vi

    def klocal(self, temp, stiffness=None):
        for run in self.runs:
            yield run.klocal(temp, stiffness)

    def kglobal(self, temp, stiffness=None):
        """Yield the global reducer stiffness matrix for each element.

        The bend consists of multiple approximating elements and so the global
        matrix for the reducer is the assembly of all the elements which make
        up the bend.
        """
        for run in self.runs:
            yield run.kglobal(temp, stiffness)


@units.define(weight="force")
class Rigid(Run):
    """Rigid elements for which user defined weights must be provided.

    By default a weightless rigid element will have zero mass, however it will
    have mass due to insulation and fluid contents defined for the element. The
    user can zero out the additional mass if needed by setting the insulation
    and fluid attributes to None.
    """

    def __init__(self, point, dx, dy=0, dz=0, weight=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Rigid, self).__init__(point, dx, dy, dz, from_point, section,
                                    material, insulation, code)
        self.weight = weight

    def mass(self, accel):
        """Inverse where the mass is returned since the weight is a user
        defined parameter for rigid elements.
        """
        return self.weight / accel

    def klocal(self, temp):
        """The stiffness parameter for a rigid element is defined in the user
        specified options.
        """
        stiff = self.app.models.active_object.settings["core.rigid_stiffness"]

        return super(Rigid, self).stiffness(temp, stiff)


class Valve(Rigid):

    @classmethod
    def from_file(cls, point, dx, dy=0, dz=0, rating=150, valve_type="gate",
                  flange_type="weldneck", from_point=None, section=None,
                  material=None, insulation=None, code=None):
        """Select a valve from a data file"""
        raise NotImplementedError("implement")

    def __init__(self, point, dx, dy=0, dz=0, weight=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Valve, self).__init__(point, dx, dy, dz, weight, from_point,
                                    section, material, insulation, code)


class Flange(Rigid):

    @classmethod
    def from_file(cls, point, dx, dy=0, dz=0, rating=150,
                  flange_type="weldneck", from_point=None, section=None,
                  material=None, insulation=None, code=None):
        """Select a flange from a data file"""
        raise NotImplementedError("implement")

    def __init__(self, point, dx, dy=0, dz=0, weight=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Flange, self).__init__(point, dx, dy, dz, weight, from_point,
                                     section, material, insulation, code)


class ElementContainer(EntityContainer):

    def __init__(self):
        super(ElementContainer, self).__init__()
        self.Run = Run
        self.Bend = Bend
        self.Reducer = Reducer
        self.Rigid = Rigid
        self.Valve = Valve
        self.Flange = Flange

        # selected object
        self._active_objects = set()

    @property
    def active_objects(self):
        return self._active_objects

    def _iter_all(self, typ):
        """Helper generator method"""
        for obj in self._objects.values():
            if isinstance(obj, typ):
                yield obj

    @property
    def runs(self):
        """A generator containing all runs"""
        return self._iter_all(Run)

    @property
    def bends(self):
        """A generator containing all bends"""
        return self._iter_all(Bend)

    @property
    def reducers(self):
        """A generator containing all reducers"""
        return self._iter_all(Reducer)

    @property
    def rigids(self):
        """A generator containing all rigids"""
        return self._iter_all(Rigid)

    @property
    def valves(self):
        """A generator containing all valves"""
        return self._iter_all(Valve)

    @property
    def flanges(self):
        """A generator containing all flanges"""
        return self._iter_all(Flange)

    def iterate(self, element_type):
        """Iterate through the list of active elements"""
        for element in self._active_objects:
            yield element

    # def split(self, inst, point, distance, reference):
    #     raise NotImplementedError("implement")

    def __call__(self, from_point, to_point):
        """An element can be retrieved by from_point to to_point."""
        for val in self._objects.values():
            if (val.from_point.name == from_point
                    and val.to_point.name == to_point):
                return val

    def select(self, type_, item, *args):
        """Select a subset of elements"""

        if type_ == "all":
            self._active_objects.clear()
            self._active_objects.update(set(self._objects.values()))

        elif type_ == "none":
            self._active_objects.clear()

        elif type_ == "invert":
            selobjs = set(self._objects.values())
            selset = set(selobjs).difference(self._active_objects)
            if selset:
                self._active_objects.clear()
                self._active_objects.update(selset)

        elif type_ == "s":
            selset = set()
            if item == "type":
                for arg in args:
                    for elem in self._objects.values():
                        if isinstance(elem, arg):
                            selset.add(elem)
            elif item in "section":
                for arg in args:
                    for elem in self._objects.values():
                        if elem.section == arg:
                            selset.add(elem)
            elif item in "material":
                for arg in args:
                    for elem in self._objects.values():
                        if elem.material == arg:
                            selset.add(elem)

            self._active_objects.clear()
            self._active_objects.update(selset)
