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

"""Different types of pipe supports.

Supports are implemented using the penalty approach where the global system
stiffness and force matrices are modified by the support stiffness value. A
stiffness value of 1000*K, where K is the largest stiffness in the global
stiffness matrix has shown to produce good results.
"""

import numpy as np

from psi.entity import Entity, EntityContainer
from psi import units


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

    def __init__(self, name, point, translation_stiffness=1e12,
                 rotation_stiffness=1e12):
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

        with units.Units(user_units="english"):
            self.translation_stiffness = translation_stiffness
            self.rotation_stiffness = rotation_stiffness

    def kglobal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            f[:3, 0] = [self.translation_stiffness] * 3
            f[3:6, 0] = [self.rotation_stiffness] * 3

        elif self.point == element.to_point.name:
            f[6:9, 0] = [self.translation_stiffness] * 3
            f[9:12, 0] = [self.rotation_stiffness] * 3

        return f


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness",
              gap="length")
class RigidSupport(Support):
    """Rigid supports can be translational or rotational. They can also be
    double-acting or directional.

    The default support is a rigid double-acting translational support with no
    gap.
    """

    def __init__(self, name, point, dircos=None, friction=None,
                 is_rotational=False, is_directional=False, is_snubber=False,
                 translation_stiffness=1e12, rotation_stiffness=1e12,
                 gap=0):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        dircos : tuple
            Support direction cosine.

        friction : float
            Support friction (non-linear).

        is_rotational : bool
            True if the support is a rotational restraint.

        is_directional : bool
            True if the support is a translational restraint.

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
        self._dircos = dircos   # direction cosine
        self.friction = friction
        self.is_rotational = is_rotational
        self.is_directional = is_directional
        self.is_snubber = is_snubber

        with units.Units(user_units="english"):
            self.translation_stiffness = translation_stiffness
            self.rotation_stiffness = rotation_stiffness
            self.gap = gap

    @property
    def stiffness(self):
        if self.is_rotational is False:
            return self.translation_stiffness
        else:
            return self.rotation_stiffness

    @stiffness.setter
    def stiffness(self, value):
        if self.is_rotational is False:
            self.translation_stiffness = value
        else:
            self.rotation_stiffness = value

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, cosx, cosy, cosz):
        self._direction = (cosx, cosy, cosz)


class GlobalX(RigidSupport):

    def kglobal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            f[0, 0] = self.translation_stiffness

        elif self.point == element.to_point.name:
            f[6, 0] = self.translation_stiffness

        return f


class GlobalY(RigidSupport):

    def kglobal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            f[1, 0] = self.translation_stiffness

        elif self.point == element.to_point.name:
            f[7, 0] = self.translation_stiffness

        return f


class GlobalZ(RigidSupport):

    def kglobal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            f[2, 0] = self.translation_stiffness

        elif self.point == element.to_point.name:
            f[8, 0] = self.translation_stiffness

        return f


@units.define(stiffness="translation_stiffness", cold_load="force")
class Spring(Support):

    def __init__(self, name, point, stiffness, cold_load, variability=25,
                 is_constant=False):
        super(Spring, self).__init__(name, point)
        self.stiffness = stiffness
        self.cold_load = cold_load
        self._variability = variability     # 25% allowable variation
        self._is_constant = is_constant     # true if constant effort

    @classmethod
    def from_table(cls, name, point, hot_load, movement, vendor="Lisega"):
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
        f = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            f[1, 0] = self.stiffness

        elif self.point == element.to_point.name:
            f[7, 0] = self.stiffness

        return f


@units.define(dx="length", dy="length", dz="length",
              mx="rotation", my="rotation", mz="rotation",
              translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class Displacement(Support):
    """A generic global displacement vector.

    Displacements are applied similar to how supports are. Supports are in
    essence a special case with 0 movement in the direction of stiffness.
    Using the penalty approach, the stiffness and force terms in the global
    system matrix are modified.

    Displacements are associated to an operating case and typically used with a
    thermal case.

    TODO: how to specify a 'free' nonzero displacement; dx, dy etc for example
    must be a number value so it is set to 0 by default
    """

    def __init__(self, name, opercase, point,
                 dx=0, dy=0, dz=0, rx=0, ry=0, rz=0,
                 translation_stiffness=1e12, rotation_stiffness=1e12):
        """Create a displacement support instance."""
        super(Displacement, self).__init__(name, point)
        self.opercase = opercase
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.translation_stiffness = translation_stiffness
        self.rotation_stiffness = rotation_stiffness

    def kglobal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            f[:3, 0] = [self.translation_stiffness] * 3
            f[3:6, 0] = [self.rotation_stiffness] * 3

        elif self.point == element.to_point.name:
            f[6:9, 0] = [self.translation_stiffness] * 3
            f[9:12, 0] = [self.rotation_stiffness] * 3

        return f

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


class SupportContainer(EntityContainer):

    def __init__(self):
        super(SupportContainer, self).__init__()
        self.Anchor = Anchor
        self.GlobalX = GlobalX
        self.GlobalY = GlobalY
        self.GlobalZ = GlobalZ
        self.Spring = Spring
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
