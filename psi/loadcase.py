"""Each loadcase consists of one or many different types of loads applied to
the piping system at once.  For instance, the sustained case for a system
without springs consists of the weight(W) case along with a pressure case(P1).
Both these loads are first applied to the underlying finite element model and
the solution is generated for that particular case.  This is a primary load
case since it is composed of primary loads.

A load combination case consists of one or more primary load cases combined
using a specific combination method.  The solution vectors from the primary
load cases are added, subtracted, etc., in order to produce these kind of
secondary load cases.

The load case may consist of primary loads such as Weight, Thermal, etc., or it
may be a be a combination of two or more secondary (ie. result cases) combined
together.

The results for the load case are stored internally in the load case object.
"""

from psi.entity import Entity, EntityContainer
from psi.utils.orderedset import OrderedSet


class BaseCase(Entity):

    def __init__(self, name, loads, stype="SUS"):
        super(BaseCase, self).__init__(name)

        self.cases = OrderedSet()   # set of loads
        self.cases.update(loads)

        self.stype = stype  # HRG, HYD, SUS, EXP, OCC, OPE, FAT
        self.results = {}   # node -> value

    @property
    def parent(self):
        return self.app.loadcases


class LoadCase(BaseCase):
    """A set of primary load cases consisting of different types of loads.

    Note that the loads for a given case must be unique, in other words, two
    sustained loads cannot be defined for one case.
    """

    def __int__(self, name, loads, stype="SUS"):
        super(LoadCase, self).__init__(name, loads, stype)

    def solve(self):
        """Solve the particular load case"""
        pass


class LoadComb(BaseCase):

    def __init__(self, name, loadcases, stype, combmeth):
        super(LoadCase, self).__init__(name, loadcases, stype, combmeth)

    def combine(self):
        pass


class LoadCaseContainer(EntityContainer):

    def __init__(self):
        super(LoadCaseContainer, self).__init__()
        self.LoadCase = LoadCase
        self.LoadComb = LoadComb

    def defaults(self):
        """A set of default cases automatically created for the user"""
        pass

    def solve(self):
        """Solve all primary cases and combination cases"""
        pass


class Stress(object):
    """Stress results"""
    pass


class Displacements(object):
    """Displacements results"""
    pass


class Restraint(object):
    """Restraint results"""
    pass
