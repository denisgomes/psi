"""Different results extracted from the loadcases solved"""

import psi
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi import units


class Result(Entity, ActiveEntityMixin):
    pass


class Stress(Result):
    """Code stress output results"""
    pass


class Movement(Result):
    """Nodal displacement results"""
    pass


class Reaction(Result):
    """Support reaction results"""
    pass


class ResultContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(ResultContainer, self).__init__()
        self.Stress = Stress
        self.Movement = Movement
        self.Reaction = Reaction
