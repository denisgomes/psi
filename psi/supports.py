"""Different types of pipe supports"""


from psi.entity import Entity, EntityContainer
from psi.units import units


class Support(Entity):

    def __init__(self, name, point):
        super(Support, self).__init__(name)
        self.point = point

    @property
    def parent(self):
        return self.app.supports


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class Anchor(Support):

    def __init__(self, name, point, translation_stiffness=1e12,
                 rotation_stiffness=1e12):
        super(Anchor, self).__init__(name, point)
        self.translation_stiffness = translation_stiffness
        self.rotation_stiffness = rotation_stiffness


@units.define(gap="length", translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class RigidSupport(Support):
    """Rigid supports can be translational or rotational. They can also be
    double-acting or directional. The default support is a rigid double-acting
    translational support with no gap.
    """

    def __init__(self, name, point, direction=None, gap=0., friction=None,
                 is_rotational=False, is_directional=False, is_snubber=False,
                 translation_stiffness=1e12, rotation_stiffness=1e12):
        super(RigidSupport, self).__init__(name, point)
        self._direction = direction     # direction cosine
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


class Vertical(RigidSupport):
    pass


class Lateral(RigidSupport):
    pass


class Axial(RigidSupport):
    pass


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


class SupportContainer(EntityContainer):

    def __init__(self):
        super(SupportContainer, self).__init__()
        self.Anchor = Anchor
        self.Vertical = Vertical
        self.Lateral = Lateral
        self.Axial = Axial
        self.Spring = Spring
