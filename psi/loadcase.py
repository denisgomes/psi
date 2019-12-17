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

The displacement, support reactions and internal force results for the loadcase
are stored internally in the load case object. Note that properties with units
for different values are stored separately and combined on the fly.
"""

from itertools import zip_longest, product
from collections import defaultdict
from contextlib import redirect_stdout
import sys

import numpy as np

from psi.entity import Entity, EntityContainer
from psi.utils.orderedset import OrderedSet
from psi import units


@units.define(_values="length")
class Translation(object):
    """Translation components of the displacement row vector"""

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """A numpy array of translations at each node.

        Parameters
        ----------
        data : numpy.array row vector
        """
        x = data[::6]
        y = data[1::6]
        z = data[2::6]
        xyz = np.array(list(zip(x, y, z)), dtype=np.float64)

        # column vector of translation
        self._values = xyz.flatten().reshape((-1, 1))


@units.define(_values="rotation")
class Rotation(object):
    """Rotation components of the displacement row vector"""

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """A numpy array of rotations at each node.

        Parameters
        ----------
        data : numpy.array row vector
        """
        rx = data[3::6]
        ry = data[4::6]
        rz = data[5::6]
        rxyz = np.array(list(zip(rx, ry, rz)), dtype=np.float64)

        self._values = rxyz.flatten().reshape((-1, 1))


class Movements(object):

    def __init__(self):
        self._translation = Translation()
        self._rotation = Rotation()

    def __call__(self, point):
        pass

    @property
    def results(self):
        x = self._translation.results[::3]
        y = self._translation.results[1::3]
        z = self._translation.results[2::3]
        rx = self._rotation.results[::3]
        ry = self._rotation.results[1::3]
        rz = self._rotation.results[2::3]

        disp = np.array(list(zip(x, y, z, rx, ry, rz)), dtype=np.float64)
        values = disp.flatten().reshape((-1, 1))

        return values

    @results.setter
    def results(self, data):
        self._translation.results = data
        self._rotation.results = data


@units.define(_values="force")
class Force(object):

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """A numpy array of forces at each node.

        Parameters
        ----------
        values : numpy.array
        """
        fx = data[::6]
        fy = data[1::6]
        fz = data[2::6]
        fxyz = np.array(list(zip(fx, fy, fz)), dtype=np.float64)

        # column vector of translation
        self._values = fxyz.flatten().reshape((-1, 1))


@units.define(_values="moment_output")
class Moment(object):

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """A numpy array of moments at each node.

        Parameters
        ----------
        values : numpy.array
        """
        mx = data[3::6]
        my = data[4::6]
        mz = data[5::6]
        mxyz = np.array(list(zip(mx, my, mz)), dtype=np.float64)

        self._values = mxyz.flatten().reshape((-1, 1))


class Reactions(object):
    """Support reactions loads"""

    def __init__(self):
        self._force = Force()
        self._moment = Moment()

    @property
    def results(self):
        fx = self._force.results[::3]
        fy = self._force.results[1::3]
        fz = self._force.results[2::3]
        mx = self._moment.results[::3]
        my = self._moment.results[1::3]
        mz = self._moment.results[2::3]

        disp = np.array(list(zip(fx, fy, fz, mx, my, mz)), dtype=np.float64)
        values = disp.flatten().reshape((-1, 1))

        return values

    @results.setter
    def results(self, data):
        self._force.results = data
        self._moment.results = data


class Forces(object):
    """Member internal forces and moments"""

    def __init__(self):
        self._force = Force()
        self._moment = Moment()

    @property
    def results(self):
        fx = self._force.results[::3]
        fy = self._force.results[1::3]
        fz = self._force.results[2::3]
        mx = self._moment.results[::3]
        my = self._moment.results[1::3]
        mz = self._moment.results[2::3]

        disp = np.array(list(zip(fx, fy, fz, mx, my, mz)), dtype=np.float64)
        values = disp.flatten().reshape((-1, 1))

        return values

    @results.setter
    def results(self, data):
        self._force.results = data
        self._moment.results = data


class BaseCase(Entity):

    def __init__(self, name, stype="SUS"):
        super(BaseCase, self).__init__(name)
        self.stype = stype  # HRG, HYD, SUS, EXP, OCC, OPE, FAT

    @property
    def points(self):
        return list(self.app.points)

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

        # results objects
        self._movements = Movements()
        self._reactions = Reactions()
        self._forces = Forces()

    @property
    def movements(self):
        return self._movements

    @property
    def reactions(self):
        return self._reactions

    @property
    def forces(self):
        return self._forces


class LoadComb(BaseCase):
    """Add loadcases using different combination methods.

    Note: Combinations pull stored data from loadcases on the fly and do the
    necessary operations.

    Parameters
    ----------
    name : str
        Unique name for load combination object.

    stype : str
        Type of code stress. Defaults to sustained stress.

        HGR - hanger load case
        HYD - hydro load case
        SUS - sustained stress case
        EXP - thermal expansion stress.case
        OCC - occasional stress case
        OPE - operating stress case
        FAT - Fatigue stress case

    method : str
        Result combination method.

        Algebriac - Disp/force results added vectorially. Stresses are derived
        from the force results.

        Scalar - Disp/force/ results added vectorially similar to the algebraic
        method. Stresses are added together.

        SRSS - Square root of the sum squared.

        Abs - Absolute summation.

        Signmax - Signed max.

        Signmin - Signed min.

    loadcases : list of loadcases.
        List of load cases.

    factors : list of numbers.
        A list of factors corresponding to each loadcase. If a factor is not
        given, a default value of 1 is used.
    """

    def __init__(self, name, stype="ope", method="algebraic", loadcases=[],
                 factors=[]):
        super(LoadComb, self).__init__(name, stype)
        self.method = method
        self.loadcases = OrderedSet()
        self.factors = factors

        for loadcase in loadcases:
            if isinstance(loadcase, LoadCase):
                self.loadcases.add(loadcase)

    @property
    def movements(self):
        movements = Movements()
        movements.results = np.zeros(len(self.points) * 6, dtype=np.float64)

        for factor, loadcase in zip_longest(self.factors, self.loadcases,
                                            fillvalue=1):
            if self.method in ("algebraic", "scalar"):
                movements.results += (factor * loadcase.movements.results)
            elif self.method == "srss":
                # note: sign of factor has no effect, always positive
                movements.results += (factor * loadcase.movements.results)**2
                movements.results = np.sqrt(movements.results)
            elif self.method == "abs":
                movements.results += (factor * np.abs(loadcase.movements.results))
            elif self.method == "signmax":
                # overwrite every time through the loop
                movements.results = np.maximum(movements.results,
                                               loadcase.movements.results)
            elif self.method == "signmin":
                movements.results = np.minimum(movements.results,
                                               loadcase.movements.results)

        return movements

    @property
    def reactions(self):
        reactions = Reactions()
        reactions.results = np.zeros(len(self.points) * 6, dtype=np.float64)

        for factor, loadcase in zip_longest(self.factors, self.loadcases,
                                            fillvalue=1):
            if self.method in ("algebraic", "scalar"):
                reactions.results += (factor * loadcase.reactions.results)
            elif self.method == "srss":
                # note: sign of factor has no effect, always positive
                reactions.results += (factor * loadcase.reactions.results)**2
                reactions.results = np.sqrt(reactions.results)
            elif self.method == "abs":
                reactions.results += (factor * np.abs(loadcase.reactions.results))
            elif self.method == "signmax":
                # overwrite every time through the loop
                reactions.results = np.maximum(reactions.results,
                                               loadcase.movements.results)
            elif self.method == "signmin":
                reactions.results = np.minimum(reactions.results,
                                               loadcase.movements.results)

        return reactions

    @property
    def forces(self):
        forces = Forces()
        forces.results = np.zeros(len(self.points) * 6, dtype=np.float64)

        for factor, loadcase in zip_longest(self.factors, self.loadcases,
                                            fillvalue=1):
            if self.method in ("algebraic", "scalar"):
                forces.results += (factor * loadcase.forces.results)
            elif self.method == "srss":
                # note: sign of factor has no effect, always positive
                forces.results += (factor * loadcase.forces.results)**2
                forces.results = np.sqrt(forces.results)
            elif self.method == "abs":
                forces.results += (factor * np.abs(loadcase.forces.results))
            elif self.method == "signmax":
                # overwrite every time through the loop
                forces.results = np.maximum(forces.results,
                                            loadcase.movements.results)
            elif self.method == "signmin":
                forces.results = np.minimum(forces.results,
                                            loadcase.movements.results)

        return forces

    @property
    def stresses(self):
        raise NotImplementedError("implement")


class LoadCaseContainer(EntityContainer):

    def __init__(self):
        super(LoadCaseContainer, self).__init__()
        self.LoadCase = LoadCase
        self.LoadComb = LoadComb

    def defaults(self):
        """A set of cases automatically created based on the thermal cases.
        Each temperature case is an operating case. Different combinations of
        weight and pressure loads associated with each thermal case is
        considered. In addition, loads such as wind, seismic and force loads
        (i.e. SRV loads) are applied to each of the operating cases.
        """
        loadtypes = defaultdict(list)
        for load in self.app.loads:
            loadtypes[load.type].append(load)

        operating = []
        sustained = []
        thermal = []
        for thermal in loadtypes["Thermal"]:
            operating.append(list(product(loadtypes["Weight"],
                                          loadtypes["Pressure"],
                                          # loadtypes["Fluid"],
                                          [thermal],
                                          )))
            sustained.append(list(product(loadtypes["Weight"],
                                          loadtypes["Pressure"],
                                          #loadtypes["Fluid"],
                                          )))

        # operating, sustained, occassional loads - primary loadcases
        count = 1
        for oper in operating:
            LoadCase("L%s" % count, "ope", loads=oper[0])
            count += 1

        for sus in sustained:
            LoadCase("L%s" % count, "sus", loads=sus[0])
            count += 1

        # thermal and occasional stresses - secondary loadcombs

        with redirect_stdout(sys.__stdout__):
            for case in operating:
                print(list(case))
