"""Different types of pipe supports"""

import numpy as np

from psi.entity import Entity, EntityContainer
from psi import units
from psi.settings import options


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


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class Anchor(Support):
    """Suppport with all 6 degrees of freedom at a node fixed"""

    def __init__(self, name, point, translation_stiffness=1e12,
                 rotation_stiffness=1e12):
        super(Anchor, self).__init__(name, point)
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


@units.define(gap="length", translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class RigidSupport(Support):
    """Rigid supports can be translational or rotational. They can also be
    double-acting or directional. The default support is a rigid double-acting
    translational support with no gap.
    """

    def __init__(self, name, point, dircos=None, gap=0, friction=None,
                 is_rotational=False, is_directional=False, is_snubber=False,
                 translation_stiffness=1e12, rotation_stiffness=1e12):
        super(RigidSupport, self).__init__(name, point)
        self._dircos = dircos   # direction cosine
        self.gap = gap
        self.translation_stiffness = translation_stiffness
        self.rotation_stiffness = rotation_stiffness
        self.friction = friction
        self.is_rotational = is_rotational
        self.is_directional = is_directional
        self.is_snubber = is_snubber

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


class SupportContainer(EntityContainer):

    def __init__(self):
        super(SupportContainer, self).__init__()
        self.Anchor = Anchor
        self.GlobalX = GlobalX
        self.GlobalY = GlobalY
        self.GlobalZ = GlobalZ
        self.Spring = Spring

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
