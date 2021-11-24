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
for additional details. The user specified global support displacement is
projected in the support direction to determine the displacement components. If
a compoment is 0, the support will not be active along that axis, and the pipe
will be able to move in that direction.

A stiffness value of 1000*K, where K is the largest stiffness in the global
stiffness matrix has shown to produce good results based on textbook examples.
Reasonable default values for translation and rotation stiffness are specified
for each support. The user can change the default values via model settings.

X, Y and Z supports are inherited from RigidSupport and can take several
different forms. They can be snubbers (only active in occasional load cases),
single or bi-directional, translational or rotational and/or define friction
and gaps.
"""

import warnings
import sys
from contextlib import redirect_stdout

import numpy as np

from psi import units
from psi.entity import Entity, EntityContainer
from psi.loads import Displacement
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

    def cglobal(self, element):
        """The support penalty vector consisting of translation and rotation
        terms. The first 6 DOFs apply to node i and last 6 to node j. The first
        3 DOFs are translational DOFs for node i. DOFs 3 to 6 are the
        rotational DOFs for node i.

        .. note::

            The cglobal vector appears in the kglobal diagonal terms when the
            support is axis aligned.
        """
        raise NotImplementedError("implement")

    def dglobal(self, element, loadcase):
        """Nodal displacement vector used for penalty method.

        By default all supports except for a support with displacement load is
        assumed to have 0 imposed displacement.

        The displacement is given as a vector given with dx, dy, dz, rx, ry and
        rz at the support node.
        """
        a = np.zeros((12, 1), dtype=np.float64)

        disps = [disp for disp in element.loads if isinstance(disp,
                 Displacement)]
        for disp in disps:
            if disp in loadcase and self.element.from_point == disp.point:
                a[:6, 0] += disp.to_list()

            elif disp in loadcase and self.element.to_point == disp.point:
                a[6:12, 0] += disp.to_list()

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

    SNB - linear support
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
        assert isinstance(value, bool), "must be boolean"

        self.is_rotational = value
        if value is True:
            self._is_snubber = False
            self.mu = 0.0

    @property
    def is_snubber(self):
        return self._is_snubber

    @is_snubber.setter
    def is_snubber(self, value):
        assert isinstance(value, bool), "must be boolean"

        self._is_snubber = value
        if value is True:
            self._is_rotational = False
            self.direction = None
            self.mu = 0.0
            self.gap = 0.0

    @property
    def is_nonlinear(self):
        return (self.has_friction or self.has_gap) and not self.is_snubber

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

    @property
    def dircos(self):
        # dircos depends on element direction
        r = (np.array(self.element.to_point.xyz, dtype=np.float64) -
             np.array(self.element.from_point.xyz, dtype=np.float64))
        dc = tuple(r / np.linalg.norm(r))

        return dc

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

    The support direction can be toggled between the two local element lateral
    directions of the pipe.

    .. note::
        If a pipe is horizontal, there is only one lateral direction. Any other
        configuration will result in two lateral directions (perpendicular to
        the pipe run direction), including for skewed piping.
    """

    _dir1 = True

    def flip(self):
        """Select alternate support lateral direction.

        .. note::
            This method does not affect the support direction for a support on
            horizontal piping.
        """
        if self._dir1:
            self._dir1 = False
        else:
            self._dir1 = True

    @property
    def dircos(self):
        vert = self.app.models.active_object.settings.vertical
        dc = self.element.dircos()

        dc1 = tuple(dc[1, :])
        dc2 = tuple(dc[2, :])

        local_x = dc[0, :]
        if vert == "y":
            up = np.array([[0], [1], [0]], dtype=np.float32)

            if np.dot(local_x, up) == 0:
                dc1 = dc2 = tuple(dc[2, :])

        elif vert == "z":
            up = np.array([[0], [0], [1]], dtype=np.float32)

            if np.dot(local_x, up) == 0:
                dc1 = dc2 = tuple(dc[1, :])

        if self._dir1:
            return dc1
        else:
            return dc2

    def clocal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                if element.is_horizontal:
                    # element local y is always aligned with global y
                    # for horizontal elements, thus lateral in local z
                    c[2, 0] = self.translation_stiffness
                elif self._dir1:
                    c[1, 0] = self.translation_stiffness
                else:
                    c[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                if element.is_horizontal:
                    c[8, 0] = self.translation_stiffness
                elif self._dir1:
                    c[7, 0] = self.translation_stiffness
                else:
                    c[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                if element.is_horizontal:
                    c[5, 0] = self.rotation_stiffness
                elif self._dir1:
                    c[4, 0] = self.rotation_stiffness
                else:
                    c[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                if element.is_horizontal:
                    c[11, 0] = self.rotation_stiffness
                elif self._dir1:
                    c[10, 0] = self.rotation_stiffness
                else:
                    c[11, 0] = self.rotation_stiffness

        return c

    def cglobal(self, element):
        T = element.T()

        return T.transpose() @ self.clocal(element)


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
