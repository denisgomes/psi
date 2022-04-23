# Pipe Stress Infinity (PSI) - The pipe stress analysis and design software.
# Copyright (c) 2021 Denis Gomes <denisgomes@consultant.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Piping finite elements.

To create a model, the user must first define a series of runs that make the
overall centerline of the model. After this process, additional components can
be added by giving a point and a direction in which to place the new item or by
"converting" a run to a different type of element.

Example
-------
First define the model centerline using Point and Run.

.. code-block:: python

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

.. code-block:: python

    >>> elements(20, 30).split(25)      # split a bend - create run and bend
    >>> Valve(25, length, mass, "mid")      # two elements, 3 nodes made
    >>> Tee(40, "welding")                  # sif.Tee by point only
    >>> Flange(10, length, mass, "down")    # flanged nozzle
    >>> Flange(70, length, mass, "up")

There are three possible ways in which a point defined for a component can be
interpreted.

1. If the point is defined as a midpoint, two elements are made, one on each
   adjacent run.

2. If the component is defined to be upstream, it will be placed in the oposite
   direction to the flow, with respect to the element nodal direction.

3. If the component is defined to be downstream, it will be placed in the same
   direction as the flow. In both the latter two cases, only one element is
   created.

When a component is defined at a boundary node, it is possible to
create another boundary node by choosing the direction that 'extends' the
model.

Components can also be created by picking a run element and 'converting' it to
a different type of element such as a valve or a flange.

.. code-block:: python

    >>> elements(50, 60).split(54, 0.5, ref=50)     # used node 50 as reference
    >>> elements(54, 60).split(55, 0.5, ref=54)
    >>> elements(54, 55).convert(Valve, mass=...)

.. note::

    Only runs can be split and merged. To split a valve for example, it must
    first be changed to a run and then split.


The following items require extra attention:

* Dealing with element weight vs mass when it comes to units.
* For Rigid elements the weight must be specified not the mass.
* Bend and reducer having multiple run approximation.
* How to easily create bends when defining geometry and how to update the
  geometry when point coordinates or bend radius is updated.
* How to covert different types of elements to runs and vice versa.
* Unit conversion should be disabled before the analysis is performed.
"""

from math import pi, sqrt
from contextlib import redirect_stdout
import sys

import numpy as np
import numpy.linalg as la
from numpy.testing import assert_array_almost_equal
from numpy import sin, cos, arccos, arctan

from psi.entity import Entity, EntityContainer
from psi import units
from psi.utils.euclid import Vector3

from psi.sections import Pipe


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
    def length(self):
        x2, y2, z2 = self.to_point().xyz
        x1, y1, z1 = self.from_point().xyz

        return sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

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

    def select(self):
        """Select an element.

        The element is placed in the active element set.

        Example
        -------
        Active elements can then be assigned different properties such as loads
        in one step.

        .. code-block:: python

            ...
            >>> [elem.select() for elem in elements if elem.mat is mat1]
            >>> p1 = Weight("W1", 1)
            >>> p1.apply()  # assigns to all active elements
        """
        self.parent.select(self)

    def unselect(self):
        self.parent.unselect(self)

    def is_active(self):
        return self.parent.is_active(self)

    def angle(self, vector):
        """Returns the sin and cos of theta where theta is the angle between
        the element local x axis and an user vector.
        """
        a = (np.array(self.to_point.xyz, dtype=np.float64) -
             np.array(self.from_point.xyz, dtype=np.float64))
        b = np.array(vector, dtype=np.float64)

        # angle element makes with vector
        theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))

        sint = sin(theta)
        cost = cos(theta)

        return sint, cost

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
        """Convert a Run element to another and vice versa."""
        raise NotImplementedError("abstract method")

    def mass(self):
        """Pipe element mass"""
        mass = (self.material.rho.value * self.section.area *
                self.length)

        return mass

    def cladding_mass(self):
        """Internal cladding mass"""
        docld = self.section.id
        dicld = docld - 2*self.cladding.thk

        cld_area = (pi/4) * (docld**2-dicld**2)

        mass = (self.cladding.rho * cld_area * self.length)

        return mass

    def insulation_mass(self):
        """Insulation mass"""
        dins = self.section.od + 2*self.insulation.thk
        do = self.section.od
        insulation_area = (pi/4) * (dins**2-do**2)

        mass = (self.insulation.rho * insulation_area *
                self.length)

        return mass

    def refractory_mass(self):
        """Internal refractory mass"""
        doref = self.section.id
        diref = doref - 2*self.refractory.thk

        ref_area = (pi/4) * (doref**2-diref**2)

        mass = (self.refractory.rho * ref_area * self.length)

        return mass

    def delete(self):
        """Deleting an element other than a run has the effect of converting to
        a run, so the points and underlying vertex data is unaffected.

        Deleting a run....

        Deleting a point has the effect of converting the element to a run and
        merging adjacent runs together.
        """
        raise NotImplementedError("abstract method")


@units.define(_dx="length", _dy="length", _dz="length")
class Run(Piping):
    """Define a pipe run by incremental offsets in the global x, y and z
    directions.
    """

    def __init__(self, point, dx, dy=0, dz=0, from_point=None, section=None,
                 material=None, insulation=None, code=None):
        super(Run, self).__init__(point, from_point, section, material,
                                  insulation, code)
        self._dx = dx
        self._dy = dy
        self._dz = dz
        self._split_edge = None     # for bends
        self.build(point)

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

    @property
    def dx(self):
        """Get and set the vertex x-coordinate of the 'to' point."""
        return self._dx

    @dx.setter
    def dx(self, value):
        xold = self._dx
        self._dx = value    # set x

        # TODO: Note that geometry coordinates are managed
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
        """Get and set the vertex y-coordinate of the 'to' point."""
        return self._dy

    @dy.setter
    def dy(self, value):
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
        """Get and set the vertex z-coordinate of the 'to' point."""
        return self._dz

    @dz.setter
    def dz(self, value):
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
        """Get and set the length of the element."""
        dx, dy, dz = self.dx, self.dy, self.dz
        return sqrt(dx**2 + dy**2 + dz**2)

    @length.setter
    def length(self, value):
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
        up = self.app.models.active_object.settings.vertical
        # note, local y is being set to a global direction
        if up == "y":
            local_y = np.array([0., 1., 0.], dtype=np.float64)
        elif up == "z":
            local_y = np.array([0., 0., 1.], dtype=np.float64)

        from_vert = np.array(self.geometry.v1.co)
        to_vert = np.array(self.geometry.v2.co)
        local_x = to_vert - from_vert

        if self.is_vertical:    # set to global x
            local_y = np.array([1., 0., 0.], dtype=np.float64)

        local_z = np.cross(local_x, local_y)

        # recalculate local y so that it's orthogonal
        local_y = np.cross(local_z, local_x)

        dc = np.array([local_x/la.norm(local_x),
                       local_y/la.norm(local_y),
                       local_z/la.norm(local_z)], dtype=np.float64)

        return dc

    @property
    def is_vertical(self):
        up = self.app.models.active_object.settings.vertical
        # note, local y is being set to a global direction
        if up == "y":
            vertical = np.array([0., 1., 0.], dtype=np.float64)
        elif up == "z":
            vertical = np.array([0., 0., 1.], dtype=np.float64)

        from_vert = np.array(self.geometry.v1.co)
        to_vert = np.array(self.geometry.v2.co)
        local_x = to_vert - from_vert

        try:
            assert_array_almost_equal(np.abs(local_x) / la.norm(local_x),
                                      np.abs(vertical) / la.norm(vertical),
                                      decimal=5)
            return True
        except AssertionError:
            return False

    @property
    def is_horizontal(self):
        up = self.app.models.active_object.settings.vertical

        from_vert = np.array(self.geometry.v1.co)
        to_vert = np.array(self.geometry.v2.co)
        local_x = to_vert - from_vert

        if up == "y":
            vertical = np.array([0., 1., 0.], dtype=np.float64)
        elif up == "z":
            vertical = np.array([0., 0., 1.], dtype=np.float64)

        if np.dot(local_x, vertical) == 0:
            return True

        return False

    def T(self):
        """Local to global transformation matrix"""
        tf = np.zeros((12, 12), dtype=np.float64)
        dc = self.dircos()

        tf[:3, :3] = tf[3:6, 3:6] = tf[6:9, 6:9] = tf[9:12, 9:12] = dc[:, :]

        return tf

    def klocal(self, temp, sfac=1.0):
        """The local stiffness matrix of a straight piping element, i.e all
        elements with the expection of bends and reducers are derived from a
        run. The stiffness is a function of the temperature.

        sfac is used for Rigid types only so that the thickness for these types
        can be increased by a factor of 10.

        Note, the code flexibility factor is applied to the element for the
        bending directions. The stiffness is divided by the flexibility factor,
        effectively making the element more flexible in the transverse bending
        directions. The torsional stiffness is not altered.

        Stiffness matrix from 'Theory of Matrix Structural Analysis' by J.S.
        Przemieniecki.
        """
        kmat = np.zeros((12, 12), dtype=np.float64)

        L = self.length
        E = self.material.ymod[temp]
        nu = self.material.nu.value     # poisson's ratio
        G = E / (2 * (1 + nu))          # shear mod, isotropic mats

        # used for rigids
        self.section.thk *= sfac

        J = self.section.ixx
        Iy = self.section.iyy
        Iz = self.section.izz

        A = self.section.area
        Ay = 0  # shear area y
        Az = 0  # shear area z
        phi_y = 0
        phi_z = 0
        if self.app.models.active_object.settings.timoshenko:
            if self.section.is_thin_wall:
                Ay = Az = (1/2) * A     # pipe shear shape factors
            else:
                Ay = Az = (27/32) * A   # heavy wall shape factors

            phi_y = (12*E*Iz) / (G*Ay*L**2)
            phi_z = (12*E*Iy) / (G*Az*L**2)

        # flex factor
        kfac = self.code.kfac(self)

        kmat[0, 0] = (E*A) / L
        kmat[0, 6] = kmat[6, 0] = (-E*A) / L

        kmat[1, 1] = (12*E*Iz / (L**3*(1+phi_y))) / kfac
        kmat[1, 5] = kmat[5, 1] = (6*E*Iz / (L**2*(1+phi_y))) / kfac
        kmat[1, 7] = kmat[7, 1] = (-12*E*Iz / (L**3*(1+phi_y))) / kfac
        kmat[1, 11] = kmat[11, 1] = (6*E*Iz / (L**2*(1+phi_y))) / kfac

        kmat[2, 2] = (12*E*Iy / (L**3*(1+phi_z))) / kfac
        kmat[2, 4] = kmat[4, 2] = (-6*E*Iy / (L**2*(1+phi_z))) / kfac
        kmat[2, 8] = kmat[8, 2] = (-12*E*Iy / (L**3*(1+phi_z))) / kfac
        kmat[2, 10] = kmat[10, 2] = (-6*E*Iy / (L**2*(1+phi_z))) / kfac

        kmat[3, 3] = (G*J) / L
        kmat[9, 3] = kmat[3, 9] = (-G*J) / L

        kmat[4, 4] = ((4+phi_z)*E*Iy / (L*(1+phi_z))) / kfac
        kmat[4, 8] = kmat[8, 4] = (6*E*Iy / (L**2*(1+phi_z))) / kfac
        kmat[4, 10] = kmat[10, 4] = ((2-phi_z)*E*Iy / (L*(1+phi_z))) / kfac

        kmat[5, 5] = ((4+phi_y)*E*Iz / (L*(1+phi_y))) / kfac
        kmat[5, 7] = kmat[7, 5] = (-6*E*Iz / (L**2*(1+phi_y))) / kfac
        kmat[5, 11] = kmat[11, 5] = ((2-phi_y)*E*Iz / (L*(1+phi_y))) / kfac

        kmat[6, 6] = (E*A) / L

        kmat[7, 7] = (12*E*Iz / (L**3*(1+phi_y))) / kfac
        kmat[7, 11] = kmat[11, 7] = (-6*E*Iz / (L**2*(1+phi_y))) / kfac

        kmat[8, 8] = (12*E*Iy / (L**3*(1+phi_z))) / kfac
        kmat[8, 10] = kmat[10, 8] = (6*E*Iy / (L**2*(1+phi_z))) / kfac

        kmat[9, 9] = (G*J) / L

        kmat[10, 10] = ((4+phi_z)*E*Iy / (L*(1+phi_z))) / kfac

        kmat[11, 11] = ((4+phi_y)*E*Iz / (L*(1+phi_y))) / kfac

        # with redirect_stdout(sys.__stdout__):
        #     print(kmat)

        return kmat

    def kglobal(self, temp, sfac=1.0):
        """The element global stiffness matrix before assembly into the system
        matrix.

        Note that the transpose of T is taken because T is defined in row major
        format, in other words, the local axes components all appear on the
        same row. Taking the transpose, (i.e. inverse, sinces it's a rotation
        system.
        """
        T = self.T()    # build T once and reuse

        return T.transpose() @ self.klocal(temp, sfac) @ T

    def mlocal(self):
        """Local mass matrix of run element including both the linear and
        rotatory inertias. It also includes the coupling betweeen the degrees
        of freedom not accounted for with the lumped mass approach, which
        has be shown to less accurate.

        From 'Theory of Matrix Structural Analysis' by J.S. Przemieniecki.
        """
        mmat = np.zeros((12, 12), dtype=np.float64)

        L = self.length
        rho = self.material.rho.value
        A = self.section.area

        J = self.section.ixx
        Iy = self.section.iyy
        Iz = self.section.izz

        # remove rotational inertia
        # J = 0
        # Iy = 0
        # Iz = 0

        mmat[0, 0] = mmat[6, 6] = 1/3
        mmat[1, 1] = mmat[7, 7] = 13/35 + 6*Iz/(5*A*L**2)
        mmat[2, 2] = mmat[8, 8] = 13/35 + 6*Iy/(5*A*L**2)
        mmat[3, 3] = mmat[9, 9] = J/(3*A)
        mmat[4, 4] = mmat[10, 10] = L**2/105 + 2*Iy/(15*A)
        mmat[5, 5] = mmat[11, 11] = L**2/105 + 2*Iz/(15*A)

        mmat[4, 2] = -(11*L)/210 - Iy/(10*A*L)
        mmat[7, 5] = 13*L/420 - Iz/(10*A*L)
        mmat[10, 8] = (11*L)/210 + Iy/(10*A*L)

        mmat[5, 1] = (11*L)/210 + Iz/(10*A*L)
        mmat[8, 4] = -13*L/420 + Iy/(10*A*L)
        mmat[11, 7] = -(11*L)/210 - Iz/(10*A*L)

        mmat[7, 1] = 9/70 - 6*Iz/(5*A*L**2)
        mmat[8, 2] = 9/70 - 6*Iy/(5*A*L**2)
        mmat[9, 3] = J/(6*A)
        mmat[10, 4] = -L**2/140 - Iy/(30*A)
        mmat[11, 5] = -L**2/140 - Iz/(30*A)

        mmat[10, 2] = 13*L/420 - Iy/(10*A*L)

        mmat[11, 1] = -13*L/420 + Iz/(10*A*L)

        return (rho*A*L) * mmat

    def mglobal(self):
        """Global mass matrix of run"""

        T = self.T()    # build T once and reuse

        return T.transpose() @ self.mlocal() @ T


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
    elements is a user defined parameter. The underlying curve data for the
    geometry is used to create additional runs and solved.  The temporary runs
    are then deleted so that all the model data is not affected overall.

    The internal run approximations and results are not accessible by the user.
    """

    def __init__(self, point, dx, dy=0, dz=0, radius="long", flange=0,
                 tol=0.01, from_point=None, section=None, material=None,
                 insulation=None, code=None):
        # calls build
        super(Bend, self).__init__(point, dx, dy, dz, from_point, section,
                                   material, insulation, code)
        self.radius = radius    # also calls build
        self.flange = flange    # 0=None, 1=single side, 2=both sides
        self.tol = tol

        self._runs = []         # internal runs

    @property
    def runs(self):
        return self._runs

    def build(self, point):
        """This build should happen during solver pre-processor step to create
        the approximating run elements and populating the runs list.
        """
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

    @property
    def near_point(self):
        raise NotImplementedError("implement")

    @property
    def far_point(self):
        raise NotImplementedError("implement")

    @property
    def mid_point(self):
        raise NotImplementedError("implement")

    def klocal(self, temp, sfac=1.0):
        # poluate the run list with temporary elements first

        for run in self._runs:
            yield run.klocal(temp, sfac)

    def kglobal(self, temp, sfac=1.0):
        """Yield the global bend stiffness matrix for each element.

        The bend consists of multiple approximating elements and so the global
        matrix for the bend is the assembly of all the elements which make up
        the bend.
        """
        for run in self.runs:
            yield run.kglobal(temp, sfac)

    def __repr__(self):
        return "%s %s" % (self.type, self.to_point)


class Reducer(Run):
    """Define a pipe reducer by giving a different section.

    The new section is activated automatically. A reducer is approximated by
    multiple runs with decreasing diameters.

    A reducer element is approximated using multiple run elements. The number
    of elements is a user defined parameter. The diameter of each element is
    stepped down from the larger to smaller pipe diameter. Similar to a bend,
    a reducer is also approximated using a curve, but in this case a straight
    curve. The straight curve can be easily subdivided by calling the euler
    operations on the curve.

    The code flexibility factor is applied to all approximating elements of the
    reducer. The internal points, elements and sections are not accessible by
    the end user and do not live in the model points or elements list.
    """

    def __init__(self, point, dx, dy=0, dz=0, section2=None, is_concentric=True,
                 from_point=None, section=None, material=None, insulation=None,
                 code=None):
        self.section2 = section2
        self.is_concentric = is_concentric
        self._runs = []             # internal run approximations
        self._sections = []         # internal section approximations

        super(Reducer, self).__init__(point, dx, dy, dz, from_point, section,
                                      material, insulation, code)

    def build(self, point):
        """Each subrun has a different section associated with it with each run
        getting smaller and smaller to approximate the reducer geometry.
        """
        super(Reducer, self).build(point)
        edge = self.geometry

        if self.section2 is None:
            self.section2 = self.section

        v1, v2 = Vector3(*edge.v1.co), Vector3(*edge.v2.co)
        v = v2 - v1
        v.normalize()   # unit vector

        L = edge.length         # length of reducer
        D = self.section.od     # diameter at from point
        d = self.section2.od    # diameter at to point
        T = self.section.thk    # thickness at from point
        t = self.section2.thk   # thickness at to point

        mat = self.material
        insul = self.insulation
        code = self.code

        N = 9   # number of division

        for i in range(0, N):
            x_i = i * (L/N)         # element start length
            x_j = (i+1) * (L/N)     # element end length

            # diameters
            de_i = (d-D)/L * x_i + D
            de_j = (d-D)/L * x_j + D
            de = (de_i+de_j) / 2

            # thickness
            te_i = (t-T)/L * x_i + T
            te_j = (t-T)/L * x_j + T
            te = (te_i+te_j) / 2

            # with redirect_stdout(sys.__stdout__):
            #     print(de, te)

            # create sections specific to the reducer
            # names and fix the runs, also section corrosion allowance
            sec = Pipe(str(self.name)+str(i), de, te)
            self._sections.append(sec)

            # create all new points and runs specific to reducer
            ve = v1 + x_j*v
            run = Run(str(self.name)+str(i), dx=ve.x, dy=ve.y, dz=ve.z,
                      from_point=from_point, section=sec, material=mat,
                      insulation=insul, code=code)
            self._runs.append(run)

    def mass(self):
        """Mass of a reducer.

        Combined mass of all the internal run elements
        """
        d1o = self.section.od
        d1i = self.section.od - 2*self.section.thke
        d2o = self.section2.od
        d2i = self.section2.od - 2*self.section2.thke

        vol = (1/12)*pi*self.length*(d1o**2+d1o*d2o+d2o**2 -
                                              d1i**2+d1i*d2i+d2i**2)
        return self.material.rho.value * vol

    def fluid_mass(self, fluden):
        """Fluid mass within the reducer.

        Combined fluid mass of all the internal run elements
        """
        d1i = self.section.od - 2*self.section.thke
        d2i = self.section2.od - 2*self.section2.thke

        fluvol = (1/12)*pi*self.length*(d1i**2+d1i*d2i+d2i**2)

        return fluden * fluvol

    def insulation_mass(self):
        """A conical shaped insulation profile is considered.

        Combined insulation mass of all the internal run elements
        """
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

    def klocal(self, temp, sfac=1.0):
        """The list of internal run elements are condensed to give a local
        12x12 stiffness matrix. The degrees of freedom are reduced using the
        formulation from Bathe, pg 718.

        Kaa_bar = Kaa - Kac*Kcc^-1*Kca  # condensed stiffness matrix
        Ra_bar = Ra - Kac*Kcc^-1*Rc     # condensed load vector

        See the static condensation pdf in the reference folder for the
        rest of the equations.

        The procedure involves rearranging the internal degrees of freedom to
        the top and to the left and then applying the equations above to
        compute the condensed 12x12 matrix.

        The force vector should also be condensed accordingly given the runs.
        """
        for run in self.runs:
            yield run.klocal(temp, sfac)

    @property
    def alpha(self):
        # sloped portion assumed 60% of length
        D1 = self.section.od
        D2 = self.section2.od
        L = 0.6 * self.length

        alp = math.degrees(atan(0.5*(D1-D2)/L))
        return alp


@units.define(weight="force")
class Rigid(Run):
    """Rigid elements for which user defined weights must be provided.

    By default a weightless rigid element with zero mass is assumed, however it
    will have mass due to insulation and fluid contents defined for derived
    classs such as a Valve. The user can zero out the additional mass if
    needed by setting the insulation and fluid attributes to None.

    .. note::

        The wall thickness of a rigid element is 10 times the thickness of a
        run.
    """

    def __init__(self, point, dx, dy=0, dz=0, weight=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Rigid, self).__init__(point, dx, dy, dz, from_point, section,
                                    material, insulation, code)
        self.weight = weight
        self.include_thermal = True     # include temperature
        self.include_fluid = True       # include contents

    def mass(self, accel):
        """The mass is returned since the weight is a user defined parameter
        for rigid elements.
        """
        return self.weight / accel

    def klocal(self, temp, sfac=10.0):
        """The thickness of rigid elements are multiplied by a factor of 10."""
        return super(Rigid, self).klocal(temp, sfac)

    def kglobal(self, temp, sfac=10.0):
        """The thickness of rigid elements are multiplied by a factor of 10."""
        return super(Rigid, self).kglobal(temp, sfac)


class Valve(Rigid):
    """A valve is modeled as a rigid element with a user defined weight.

    User defined weights can be loaded from a data file. Insulation weight is
    multiplied by 1.75 to account for additional material.
    """

    def __init__(self, point, dx, dy=0, dz=0, weight=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Valve, self).__init__(point, dx, dy, dz, weight, from_point,
                                    section, material, insulation, code)

    @classmethod
    def from_file(cls, point, dx, dy=0, dz=0, rating=150, valve_type="gate",
                  flange_type="weldneck", from_point=None, section=None,
                  material=None, insulation=None, code=None):
        """Select a valve from a data file"""
        raise NotImplementedError("implement")

    def insulation_mass(self):
        """Valves have more insulation to cover surface area."""
        return 1.75 * super(Valve, self).insulation_mass()


class Flange(Rigid):
    """A flange piping component"""

    def __init__(self, point, dx, dy=0, dz=0, weight=0, from_point=None,
                 section=None, material=None, insulation=None, code=None):
        super(Flange, self).__init__(point, dx, dy, dz, weight, from_point,
                                     section, material, insulation, code)

    @classmethod
    def from_file(cls, point, dx, dy=0, dz=0, rating=150,
                  flange_type="weldneck", from_point=None, section=None,
                  material=None, insulation=None, code=None):
        """Select a flange from a data file"""
        raise NotImplementedError("implement")

    def insulation_mass(self):
        """Flanges have more insulation to cover surface area."""
        return 1.75 * super(Flange, self).insulation_mass()


class Bellow(Run):
    """An expansion joint"""

    pass


class ElementContainer(EntityContainer):

    def __init__(self):
        super(ElementContainer, self).__init__()
        self.Run = Run
        self.Bend = Bend
        self.Reducer = Reducer
        self.Rigid = Rigid
        self.Valve = Valve
        self.Flange = Flange
        self.Bellow = Bellow

        # selected set of elements
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

    def select(self, inst=None):
        """Add the element instance to the active set"""
        if inst:
            self._active_objects.add(inst)
        else:
            # select all elements
            self._active_objects.update(self._objects.values())

    def unselect(self, inst=None):
        """Remove the element instance from the active set"""
        if inst:
            self._active_objects.remove(inst)
        else:
            # unselect all
            self._active_objects.clear()

    def is_active(self, inst):
        """Check to see if the element is in the active set"""
        return inst in self._active_objects
