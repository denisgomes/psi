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

"""Different types of pipe supports.

Supports are implemented using the penalty approach where the global system
stiffness and force matrices are modified by the support stiffness value. A
stiffness value of 1000*K, where K is the largest stiffness in the global
stiffness matrix has shown to produce good results.
"""

import warnings

import numpy as np

from psi.entity import Entity, EntityContainer
from psi.loads import Load
from psi import units
from psi.units import DEFAULT_UNITS


class Support(Entity):

    def __init__(self, name, point):
        super(Support, self).__init__(name)
        self.point = point

    def apply(self, elements=None):
        """Apply the support to the elements.

        Parameters
        ----------
        elements : list
            A list of elements. If elements is None, support is applied to the
            active elements.
        """
        self.parent.apply([self], elements)

    @property
    def parent(self):
        return self.app.supports

    def kglobal(self, element):
        raise NotImplementedError("implement")

    def dglobal(self, element):
        """Nodal displacements used for penalty method.

        By default all supports are assumed to have 0 displacements.
        """
        a = np.zeros((12, 1), dtype=np.float64)

        return a


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
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
            Stiffness in the translational directions.

            .. note::
               The default value is based on imperial units.

        rotation_stiffness : float
            Stiffness in the rotational directions.

            .. note::
               The default value is based on imperial units.
        """
        super(Anchor, self).__init__(name, point)

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            k[:3, 0] = [self.translation_stiffness] * 3
            k[3:6, 0] = [self.rotation_stiffness] * 3

        elif self.point == element.to_point.name:
            k[6:9, 0] = [self.translation_stiffness] * 3
            k[9:12, 0] = [self.rotation_stiffness] * 3

        return k


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness",
              gap="length")
class RigidSupport(Support):
    """Rigid supports can be translational or rotational. They can also be
    double-acting or directional.

    The default support is a rigid double-acting translational support with no
    gap.
    """

    def __init__(self, name, point, friction=None, is_rotational=False,
                 is_nonlinear=False, is_snubber=False, gap=0):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        friction : float
            Support friction (non-linear).

        is_rotational : bool
            True if the support is a rotational restraint.

        is_nonlinear : bool
            True if the support is a nonlinear restraint (i.e can liftoff)

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
        self.friction = friction
        self.is_rotational = is_rotational
        self.is_nonlinear = is_nonlinear
        self.is_snubber = is_snubber

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness
            self.gap = gap


class GlobalX(RigidSupport):
    """Support aligned with the global x direction."""

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                k[0, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                k[6, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                k[3, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                k[9, 0] = self.rotation_stiffness

        return k


class GlobalY(RigidSupport):
    """Support aligned with the global y direction."""

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                k[1, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                k[7, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                k[4, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                k[10, 0] = self.rotation_stiffness

        return k


class GlobalZ(RigidSupport):
    """Support aligned with the global z direction."""

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                k[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                k[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                k[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                k[11, 0] = self.rotation_stiffness

        return k


class LineStop(RigidSupport):
    """Support aligned with the axial direction of the pipe.

    LineStop supports are used to redirect thermal movement. They are commonly
    used for rack piping with expansion loops.
    """

    def klocal(self, element):
        """The local element x direction is the inline/axial direction"""
        k = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                k[0, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                k[6, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                k[3, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                k[9, 0] = self.rotation_stiffness

        return k

    def kglobal(self, element):
        """Transform such that support is always inline with respect to the
        element coordinate system.
        """
        T = element.T()

        return T.transpose() @ self.klocal(element)


class Guide(RigidSupport):
    """Support perpendicular to the pipe direction.

    An exceptional case is a guided riser support which restricts movement in
    the horizontal plane.
    """

    def klocal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                k[1, 0] = self.translation_stiffness
                k[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                k[7, 0] = self.translation_stiffness
                k[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                k[4, 0] = self.rotation_stiffness
                k[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                k[10, 0] = self.rotation_stiffness
                k[11, 0] = self.rotation_stiffness

        return k

    def kglobal(self, element):
        T = element.T()

        return T.transpose() @ self.klocal(element)


@units.define(dx="length", dy="length", dz="length",
              mx="rotation", my="rotation", mz="rotation",
              translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class Displacement(Support, Load):
    """A generic global displacement vector.

    Displacements are applied similar to how supports are. Supports are in
    essence a special case with 0 movement in the direction of stiffness.
    Using the penalty approach, the stiffness and force terms in the global
    system matrix are modified.

    Displacements are associated to an operating case and typically used with a
    thermal case.

    TODO:

    Add a way to specify a 'free' nonzero displacement; dx, dy etc for
    example must be a number value so it is set to 0 by default otherwise
    free.

    Move to loads.py so that displacement is opercase dependent and update
    solver.py such that the global stiffness matrix is reset before solving for
    each loadcase.
    """

    def __init__(self, name, opercase, point, dx=0, dy=0, dz=0,
                 rx=0, ry=0, rz=0):
        """Create a displacement support instance."""
        Support.__init__(self, name, point)
        Load.__init__(self, name, opercase)
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rx = rx
        self.ry = ry
        self.rz = rz

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            k[:3, 0] = [self.translation_stiffness] * 3
            k[3:6, 0] = [self.rotation_stiffness] * 3

        elif self.point == element.to_point.name:
            k[6:9, 0] = [self.translation_stiffness] * 3
            k[9:12, 0] = [self.rotation_stiffness] * 3

        return k

    def dglobal(self, element):
        """Nodal displacements used for penalty method"""
        a = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            a[:6, 0] = [self.dx, self.dy, self.dz,
                        self.rx, self.ry, self.rz]

        elif self.point == element.to_point.name:
            a[6:12, 0] = [self.dx, self.dy, self.dz,
                          self.rx, self.ry, self.rz]

        return a

    def fglobal(self, element):
        return self.kglobal(element) * self.dglobal()


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

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        vert = self.app.models.active_object.settings.vertical

        if vert == "y":
            if self.point == element.from_point.name:
                k[1, 0] = self.spring_rate
            elif self.point == element.to_point.name:
                k[7, 0] = self.spring_rate

        elif vert == "z":
            if self.point == element.from_point.name:
                k[2, 0] = self.spring_rate
            elif self.point == element.to_point.name:
                k[8, 0] = self.spring_rate

        return k


class Incline(RigidSupport):
    """Skewed support implemented using direction cosine."""

    def __init__(self, name, point, friction=None, is_rotational=False,
                 is_nonlinear=False, is_snubber=False, gap=0,
                 dircos=(1, 0, 0)):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        friction : float
            Support friction (non-linear).

        is_rotational : bool
            True if the support is a rotational restraint.

        is_nonlinear : bool
            True if the support is a nonlinear restraint (i.e can liftoff)

        is_snubber : bool
            True if the support is snubber.

        translation_stiffnesss : float
            Stiffness in the translational direction.

        rotation_stiffness : float
            Stiffness in the rotational direction.

        gap : float
            Support gap (non-linear).

        dircos : tuple
            Support direction cosine vector. Defaults to global x support.
        """
        super(Incline, self).__init__(name, point, friction=friction,
                                      is_rotational=is_rotational,
                                      is_nonlinear=is_nonlinear,
                                      is_snubber=is_snubber, gap=gap)
        self.dircos = dircos

    def kglobal(self, element):
        k = np.zeros((12, 1), dtype=np.float64)

        cosx, cosy, cosz = self.dircos

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                k[0, 0] = self.translation_stiffness * cosx
                k[1, 0] = self.translation_stiffness * cosy
                k[2, 0] = self.translation_stiffness * cosz
            elif self.point == element.to_point.name:
                k[6, 0] = self.translation_stiffness * cosx
                k[7, 0] = self.translation_stiffness * cosy
                k[8, 0] = self.translation_stiffness * cosz

        else:
            if self.point == element.from_point.name:
                k[3, 0] = self.rotation_stiffness * cosx
                k[4, 0] = self.rotation_stiffness * cosy
                k[5, 0] = self.rotation_stiffness * cosz
            elif self.point == element.to_point.name:
                k[9, 0] = self.rotation_stiffness * cosx
                k[10, 0] = self.rotation_stiffness * cosy
                k[11, 0] = self.rotation_stiffness * cosz

        return k


class SupportContainer(EntityContainer):

    def __init__(self):
        super(SupportContainer, self).__init__()
        self.Anchor = Anchor
        self.GlobalX = GlobalX
        self.GlobalY = GlobalY
        self.GlobalZ = GlobalZ
        self.LineStop = LineStop
        self.Guide = Guide
        self.Spring = Spring
        self.Incline = Incline
        self.Displacement = Displacement

    def apply(self, supports=[], elements=None):
        """Apply supports to elements.

        A copy of the support is attached to each element, where each copy
        shares the same name.

        Parameters
        ----------
        supports : list
            A list of supports
        elements : list
            A list of elements. If elements is None, loads are applied to all
            elements.
        """
        if elements is None:
            elements = []

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.supports.update(supports)
