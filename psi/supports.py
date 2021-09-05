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

"""Implementation of different types of theoritical pipe supports.

Supports are implemented using the penalty approach where the global system
stiffness and force matrices are modified by the support stiffness value and
displacement respectively.

Axis aligned supports directly scale the diagonal elements of the global
stiffness matrix. Inclined or skewed supports are implemented using constraint
equations and modify more than the diagonal terms as the displacements are
coupled via the direction cosines.

For support displacements, the support displacement vector is multiplied by the
support stiffness and added to the corresponding force vector. Refer to
"Introduction to Finite Elements in Engineering" by Chandrupatla and Belegundu
for additional details.

A stiffness value of 1000*K, where K is the largest stiffness in the global
stiffness matrix has shown to produce good results based on textbook examples.
Reasonable default values for translation and rotation stiffness are specified
for each support. The user can change the default values via model settings.

X, Y and Z supports are inherited from Inclined supports and can take several
different forms. They can be snubbers (only active in occasional load cases),
single or bi-directional, translational or rotational and/or define friction
and gaps.
"""

import warnings
import sys
from contextlib import redirect_stdout

import numpy as np

from psi.entity import Entity, EntityContainer
from psi import units
from psi.units import DEFAULT_UNITS


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class Support(Entity):

    def __init__(self, name, point):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        translation_stiffnesss : float
            Stiffness in the translational directions.

            .. note::
               The default value is based on english units.

        rotation_stiffness : float
            Stiffness in the rotational directions.

            .. note::
               The default value is based on english units.
        """
        super(Support, self).__init__(name)
        self.point = point
        self.elememt = None

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness

    def apply(self, element):
        """Apply the support to the element.

        Parameters
        ----------
        element : Element
            An element object.

        Supports are applied to an element at a node on the element.

        Example
        -------
        Create a support at node 20 of run element 20.

        .. code-block:: python

            >>> run20 = elements(10, 20)
            >>> anc = Anchor("anc20", 20)
            >>> anc.apply([run20])
        """
        self.parent.apply([self], element)

    @property
    def parent(self):
        return self.app.supports

    def kglobal(self, element):
        """The support stiffness matrix used to modify the system stiffness
        matrix.

        .. note::

            kglobal is a matrix. For axis aligned supports only the diagonal
            terms consist of the penalty terms. This is not the case for
            skewed supports which modify more than just the diagonal elements.
        """
        k = np.zeros((12, 12), dtype=np.float64)

        c = self.cglobal(element)   # penalty vector
        di = np.diag_indices(6)

        if self.point == element.from_point.name:
            k[:6, :6][di] = c[:6, 0]

        elif self.point == element.to_point.name:
            k[6:12, 6:12][di] = c[6:12, 0]

        return k

    def fglobal(self, element):
        """The support force required to produce a prescribed displacement at a
        node.
        """
        pass

    def cglobal(self, element):
        """The support penalty vector consisting of translation and rotation
        terms. The first 6 DOFs apply to node i and last 6 to node j. The
        first 3 DOFs are translational DOFs for node i. DOFs 3 to 6 are the
        rotational DOFs for node i.

        .. note::

            The cglobal vector appears in the kglobal diagonal terms when the
            support is axis aligned.
        """
        raise NotImplementedError("implement")

    def dglobal(self, element):
        """Nodal displacement vector used for penalty method.

        By default all supports except for the Displacement support are assumed
        to have 0 imposed displacement.
        """
        a = np.zeros((12, 1), dtype=np.float64)

        return a


class Anchor(Support):
    """Support with all 6 degrees of freedom at a node fixed."""

    def __init__(self, name, point):
        """Create an anchor support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        translation_stiffnesss : float
            Stiffness in the translational direction.

        rotation_stiffness : float
            Stiffness in the rotational direction.
        """
        super(Anchor, self).__init__(name, point)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            c[:3, 0] = [self.translation_stiffness] * 3
            c[3:6, 0] = [self.rotation_stiffness] * 3

        elif self.point == element.to_point.name:
            c[6:9, 0] = [self.translation_stiffness] * 3
            c[9:12, 0] = [self.rotation_stiffness] * 3

        return c


@units.define(gap="length")
class RigidSupport(Support):    # Abstract Base Class
    """A generic rigid support with different configurations.

        1. Translational or rotational (X, Y, Z, RX, RY, RZ)
        2. Single of double-acting, (+/-) (X, Y, Z, RX, RY, RZ)
        3. Gap and friction,
        4. Snubber (+/-) (XSNB, YSNB, ZSNB).

    The default support is a rigid double-acting translational support with no
    friction or gap, in other words a linear support.

    .. note::

        Snubbers are translational linear supports, so gaps and friction are
        not valid. They are active only for occasional cases.

    .. note::

        Rotational supports do not have friction.

    X, Y, Z - support for friction

    +X, +Y, +Z, +RX, +RY, +RZ - support for gaps
    +X, +Y, +Z - support for friction

    -X, -Y, -Z, -RX, -RY, -RZ - support for gaps
    -x, -Y, -Z - support for friction

    (+/-) SNB - linear support
    """

    def __init__(self, name, point, direction=None, gap=0.0, mu=0.0,
                 is_rotational=False, is_snubber=False):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        direction : str
            Support direction. Default is None. "+" and "-" is used to
            specify a directional support (non-linear).

        mu : float
            Support friction (non-linear).

        is_rotational : bool
            True if the support is a rotational restraint.

        is_snubber : bool
            True if the support is snubber.

        translation_stiffnesss : float
            Stiffness in the translational direction.

        rotation_stiffness : float
            Stiffness in the rotational direction.

        gap : float
            Support gap (non-linear).
        """
        super(RigidSupport, self).__init__(name, point)
        self.direction = direction
        self.mu = mu
        self._is_rotational = is_rotational
        self._is_snubber = is_snubber

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness
            self.gap = gap

    @property
    def is_rotational(self):
        return self._is_rotational

    @is_rotational.setter
    def is_rotational(self, value):
        if value is True:
            self._is_rotational = True
            self._is_snubber = False
            self.mu = 0.0

    @property
    def is_snubber(self):
        return self._is_snubber

    @is_snubber.setter
    def is_snubber(self, value):
        if value is True:
            self._is_snubber = True
            self._is_rotational = False
            self.mu = 0.0
            self.gap = 0.0

    @property
    def is_nonlinear(self):
        return self.has_friction or self.has_gap or not self.is_snubber

    @property
    def has_friction(self):
        return self.mu > 0.0

    @property
    def has_gap(self):
        return self.gap > 0.0 or str(self.direction) in "+-"

    @property
    def type(self):
        supdir = self.direction if self.direction else ""

        if self.is_snubber:
            _type = supdir + super(RigidSupport, self).type + "SNB"
        elif self.is_rotational:
            _type = supdir + "R" + super(RigidSupport, self).type
        else:
            _type = super(RigidSupport, self).type

        return _type


class Inclined(RigidSupport):
    """A skewed roller support not aligned with a global axis.

    Parameters
    ----------
    dircos : tuple
        Support direction cosine. Also pointing in the positive support
        direction.

        .. note::
            The support direction gives the sense of the support, "+" or "-".
    """

    def __init__(self, name, point, direction=None, gap=0.0, mu=0.0,
                 is_rotational=False, is_snubber=False, dircos=None):

        super(Inclined, self).__init__(name, point, direction, gap, mu,
                                       is_rotational, is_snubber)
        self.dircos = dircos

    def kglobal(self, element):
        k = np.zeros((12, 12), dtype=np.float64)
        B1, B2, B3 = self.dircos
        B = np.array([[B1], [B2], [B3]], dtype=np.float64) @  \
            np.array([[B1, B2, B3]], dtype=np.float64)

        if self.is_rotational is False:
            k_trans = self.translation_stiffness * B

            if self.point == element.from_point.name:
                k[:3, :3] = k_trans[:, :]
            elif self.point == element.to_point.name:
                k[6:9, 6:9] = k_trans[:, :]

        else:
            k_rot = self.rotation_stiffness * B

            if self.point == element.from_point.name:
                k[3:6, 3:6] = k_rot[:, :]
            elif self.point == element.to_point.name:
                k[9:12, 9:12] = k_rot[:, :]

        return k

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        # does c need to be multiplied by dircos?
        # does d need to be multiplied by dircos?

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[:3, 0] = [self.translation_stiffness if dc != 0 else 0
                            for dc in self.dircos]
            elif self.point == element.to_point.name:
                c[6:9, 0] = [self.translation_stiffness if dc != 0 else 0
                             for dc in self.dircos]
        else:
            if self.point == element.from_point.name:
                c[3:6, 0] = [self.rotation_stiffness if dc != 0 else 0
                             for dc in self.dircos]
            elif self.point == element.to_point.name:
                c[9:12, 0] = [self.rotation_stiffness if dc != 0 else 0
                              for dc in self.dircos]

        return c

    @property
    def type(self):
        supdir = self.direction if self.direction else ""

        if self.is_snubber:
            _type = supdir + "I" + "SNB"
        elif self.is_rotational:
            _type = supdir + "R" + "I"
        else:
            _type = "I"

        # add on the dircos
        _type += repr(self.dircos)

        return _type


class X(RigidSupport):
    """Support aligned with the global x direction."""

    @property
    def dircos(self):
        return (1, 0, 0)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[0, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[6, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[3, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[9, 0] = self.rotation_stiffness

        return c


class Y(RigidSupport):
    """Support aligned with the global y direction."""

    @property
    def dircos(self):
        return (0, 1, 0)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[1, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[7, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[4, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[10, 0] = self.rotation_stiffness

        return c


class Z(RigidSupport):
    """Support aligned with the global z direction."""

    @property
    def dircos(self):
        return (0, 0, 1)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[11, 0] = self.rotation_stiffness

        return c


class LineStop(RigidSupport):
    """Support aligned with the axial direction of the pipe.

    LineStop supports are used to redirect thermal movement. They are commonly
    used for rack piping with expansion loops.
    """

    _dircos = None

    @property
    def dircos(self):
        return self._dircos

    def apply(self, element):
        super(LineStop, self).apply(element)

        # dircos depends on element direction
        r = (np.array(self.element.to_point.xyz, dtype=np.float64) -
             np.array(self.element.from_point.xyz, dtype=np.float64))
        self._dircos = tuple(r / np.linalg.norm(r))

    def clocal(self, element):
        """The local element x direction is the inline/axial direction"""
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[0, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[6, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[3, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[9, 0] = self.rotation_stiffness

        return c

    def cglobal(self, element):
        """Transform such that support is always inline with respect to the
        element coordinate system.
        """
        T = element.T()

        return T.transpose() @ self.clocal(element)


class Guide(RigidSupport):
    """Support perpendicular to the pipe run direction.

    An exceptional case is a guided riser support which restricts movement in
    the horizontal plane.
    """

    _dir1 = True
    _dircos1 = None
    _dircos2 = None

    @property
    def toggle(self):
        if self._dir1:
            self._dir1 = False
        else:
            self._dir1 = True

    @property
    def dircos(self):
        if self._dir1:
            return self._dircos1
        else:
            return self._dircos2

    def apply(self, element):
        super(Guide, self).apply(element)

        dc = self.element.dircos()  # element dircos
        self._dircos1 = tuple(dc[:, 1])
        self._dircos2 = tuple(dc[:, 2])

    def clocal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                if self._dir1:
                    c[1, 0] = self.translation_stiffness
                else:
                    c[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                if self._dir1:
                    c[7, 0] = self.translation_stiffness
                else:
                    c[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                if self._dir1:
                    c[4, 0] = self.rotation_stiffness
                else:
                    c[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                if self._dir1:
                    c[10, 0] = self.rotation_stiffness
                else:
                    c[11, 0] = self.rotation_stiffness

        return c

    def cglobal(self, element):
        T = element.T()

        return T.transpose() @ self.clocal(element)


@units.define(_dx="length", _dy="length", _dz="length",
              _mx="rotation", _my="rotation", _mz="rotation")
class Displacement(Support):
    """A displacement support.

    Displacements are applied to a stiffness matrix similar to how supports
    are. Supports are in essence a special case with 0 movement in the
    direction of stiffness. Using the penalty approach, the stiffness and force
    terms in the global system matrix are modified.

    Support displacements are associated to an operating case and typically
    used with a thermal case to model equipment nozzle movements.

    .. note::

        If a displacement is not explicitly defined for a particular direction,
        (i.e. None) the pipe is free to move in that direction.
    """

    def __init__(self, name, opercase, point, dx=None, dy=None, dz=None,
                 rx=None, ry=None, rz=None):
        """Create a displacement support instance."""
        Support.__init__(self, name, point)
        self.opercase = opercase
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rx = rx
        self.ry = ry
        self.rz = rz

    @property
    def dx(self):
        return (None if self._dx is np.nan else self._dx)

    @dx.setter
    def dx(self, value):
        self._dx = np.nan if value is None else value

    @property
    def dy(self):
        return (None if self._dy is np.nan else self._dy)

    @dy.setter
    def dy(self, value):
        self._dy = np.nan if value is None else value

    @property
    def dz(self):
        return (None if self._dz is np.nan else self._dz)

    @dz.setter
    def dz(self, value):
        self._dz = np.nan if value is None else value

    @property
    def rx(self):
        return (None if self._rx is np.nan else self._rx)

    @rx.setter
    def rx(self, value):
        self._rx = np.nan if value is None else value

    @property
    def ry(self):
        return (None if self._ry is np.nan else self._ry)

    @ry.setter
    def ry(self, value):
        self._ry = np.nan if value is None else value

    @property
    def rz(self):
        return (None if self._rz is np.nan else self._rz)

    @rz.setter
    def rz(self, value):
        self._rz = np.nan if value is None else value

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)
        disp = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            c[:3, 0] = [self.translation_stiffness] * 3
            c[3:6, 0] = [self.rotation_stiffness] * 3

            disp[:6, 0] = [self._dx, self._dy, self._dz,
                           self._rx, self._ry, self._rz]

        elif self.point == element.to_point.name:
            c[6:9, 0] = [self.translation_stiffness] * 3
            c[9:12, 0] = [self.rotation_stiffness] * 3

            disp[6:12, 0] = [self._dx, self._dy, self._dz,
                             self._rx, self._ry, self._rz]

        # zero out stiffness with corresponding 'nan' displacement, i.e
        # a free DOF
        c[np.isnan(disp)] = 0

        return c

    def dglobal(self, element):
        """Nodal displacements are multiplied by the support stiffness and
        added to the force vector.
        """
        a = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            a[:6, 0] = [self._dx, self._dy, self._dz,
                        self._rx, self._ry, self._rz]

        elif self.point == element.to_point.name:
            a[6:12, 0] = [self._dx, self._dy, self._dz,
                          self._rx, self._ry, self._rz]

        # zero out 'nan' displacement
        a[np.isnan(a)] = 0

        return a


@units.define(spring_rate="translation_stiffness", cold_load="force")
class Spring(Support):

    def __init__(self, name, point, spring_rate, cold_load, variability=25,
                 is_constant=False):
        warnings.warn("do not use, implementation incomplete")

        super(Spring, self).__init__(name, point)
        self.spring_rate = spring_rate
        self.cold_load = cold_load
        self._variability = variability     # 25% per MSS-SP-58
        self._is_constant = is_constant     # true if constant effort

    @classmethod
    def from_table(cls, name, point, hot_load, movement, vendor="anvil"):
        """Pick a spring from a vendor hanger table based on the hot load
        and movement.
        """
        pass

    @property
    def variability(self):
        return self._variability

    @variability.setter
    def variability(self, variability):
        """The change in load from the cold load to the hot load"""
        self._variability = variability
        if variability == 0:
            self._is_constant = True

    @property
    def is_constant(self):
        return self._is_constant

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        vert = self.app.models.active_object.settings.vertical

        if vert == "y":
            if self.point == element.from_point.name:
                c[1, 0] = self.spring_rate
            elif self.point == element.to_point.name:
                c[7, 0] = self.spring_rate

        elif vert == "z":
            if self.point == element.from_point.name:
                c[2, 0] = self.spring_rate
            elif self.point == element.to_point.name:
                c[8, 0] = self.spring_rate

        return c


class SupportContainer(EntityContainer):

    def __init__(self):
        super(SupportContainer, self).__init__()
        self.Anchor = Anchor
        self.Inclined = Inclined
        self.X = X
        self.Y = Y
        self.Z = Z
        self.LineStop = LineStop
        self.Guide = Guide
        self.Spring = Spring
        self.Displacement = Displacement

    def apply(self, supports=[], elements=[]):
        """Apply supports to elements.

        A reference of the support is attached to each element, a one to one
        assignment.

        Parameters
        ----------
        supports : list
            A list of supports

        elements : list
            A list of elements. If elements is None, loads are applied to all
            elements.
        """
        for support, element in zip(supports, elements):
            support.element = element
            element.supports.add(support)
