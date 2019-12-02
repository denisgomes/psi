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

A loadcase may consist of primary loads such as Weight, Thermal, etc., or it
may be a be a combination of two or more secondary (ie. result cases) combined
together.

The results for the loadcase is stored internally in the load case object.
"""

from psi.entity import Entity, EntityContainer
from psi.utils.orderedset import OrderedSet


class BaseCase(Entity):

    def __init__(self, name, stype="SUS"):
        super(BaseCase, self).__init__(name)
        self.stype = stype  # HRG, HYD, SUS, EXP, OCC, OPE, FAT
        self.results = {}   # node -> value

    @property
    def parent(self):
        return self.app.loadcases


class LoadCase(BaseCase):
    """A set of primary load cases consisting of different types of loads.

    Note that the loads for a given case must be unique, in other words, the
    same load cannot be specified twice.
    """
    def __init__(self, name, stype="SUS", loads=[]):
        super(LoadCase, self).__init__(name, stype)
        self.loads = OrderedSet()

        for load in loads:
            self.loads.add(load)

    def solve(self):
        """Solve the particular load case"""
        pass


class LoadComb(BaseCase):

    def __init__(self, name, stype="SUS", loadcases=[]):
        super(LoadComb, self).__init__(name, stype)

        self.loadcases = OrderedSet()

        for loadcase in loadcases:
            self.loadcases.add(loadcase)

    def __add__(self, other):
        if isinstance(other, LoadComb):
            # add results here
            pass

    def __sub__(self, other):
        pass

    def __mult__(self, other):
        pass

    def __rmul__(self, other):
        pass

    def __truediv__(self, other):
        pass

    __floordiv__ = __truediv__

    def __pow__(self, other, value):
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
