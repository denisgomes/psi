# Pipe Stress Infinity (PSI) - The pipe stress design and analysis software.
# Copyright (c) 2019 Denis Gomes

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
from psi.point import Point
from psi.utils.orderedset import OrderedSet
from psi import units


@units.define(_values="length")
class Translation:
    """Translation components of the displacement row vector.

    .. note::
       The units for the row vector are converted on the fly to base units and
       back to user units.
    """

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """A numpy array of translations at each node.

        Parameters
        ----------
        data : numpy.array row vector

        .. note::
            Data is passed in as a row vector and then converted to a column
            vector.
        """
        dx = data[::6]
        dy = data[1::6]
        dz = data[2::6]
        xyz = np.array(list(zip(dx, dy, dz)), dtype=np.float64)

        # column vector of translation
        self._values = xyz.flatten().reshape((-1, 1))


@units.define(_values="rotation")
class Rotation:
    """Rotation components of the displacement row vector.

    .. note::
       The units for the row vector are converted on the fly to base units and
       back to user units.
    """

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """A numpy array of rotations at each node.

        Parameters
        ----------
        data : numpy.array row vector

        .. note::
            Data is passed in as a row vector and then converted to a column
            vector.
        """
        rx = data[3::6]
        ry = data[4::6]
        rz = data[5::6]
        rxyz = np.array(list(zip(rx, ry, rz)), dtype=np.float64)

        self._values = rxyz.flatten().reshape((-1, 1))


class Movements:
    """Both the translation and rotation components of the displacement row
    vector.

    .. note::
       The translation and rotation components are seperately stored because
       they have different units.

    .. code-block: python

        >>> r1.movements[pt10]          # returns all comps
        >>> r1.movements[pt10][0:3]     # returns dx, dy, dz
        >>> r1.movements[pt10][3:6]     # returns rx, ry, rz

    Slices with strides can also be used such as:

    .. code-block: python

        >>> r1.movements[pt10][::2]     # returns dx, dz, ry

    Data index and corresponding nodal displacement.

    * 0 = dx
    * 1 = dy
    * 2 = dz
    * 3 = rx
    * 4 = ry
    * 5 = rz
    """

    def __init__(self, app):
        self._app = app
        self._translation = Translation()
        self._rotation = Rotation()

    def __getitem__(self, item):
        """Get nodal results."""
        ndof = 6    # degrees of freedom per node
        if isinstance(item, Point):
            try:
                model = self._app.models.active_object
                idxi = list(model.points).index(item)
                niqi, niqj = idxi*ndof, idxi*ndof + ndof

                # note that .results is a column vector
                return self.results[niqi:niqj, 0]
            except:
                raise ValueError("node is not in results array")

        else:
            raise ValueError("value not found")

    @property
    def results(self):
        dx = self._translation.results[::3]
        dy = self._translation.results[1::3]
        dz = self._translation.results[2::3]
        rx = self._rotation.results[::3]
        ry = self._rotation.results[1::3]
        rz = self._rotation.results[2::3]

        disp = np.array(list(zip(dx, dy, dz, rx, ry, rz)), dtype=np.float64)
        values = disp.flatten().reshape((-1, 1))

        return values

    @results.setter
    def results(self, data):
        self._translation.results = data
        self._rotation.results = data


@units.define(_values="force")
class Force:
    """Translation components of the force row vector.

    .. note::
       The units for the row vector are converted on the fly to base units and
       back to user units.
    """

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
class Moment:
    """Rotation components of the force row vector.

    .. note::
       The units for the row vector are converted on the fly to base units and
       back to user units.
    """

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


class Reactions:
    """Support reaction loads.

    .. note::
       The translation and rotation components are separated because they have
       different units.

    .. code-block: python

        >>> r1.reactions[pt10]          # returns all comps
        >>> r1.reactions[pt10][0:3]     # returns dx, dy, dz
        >>> r1.reactions[pt10][3:6]     # returns rx, ry, rz

    Slices with strides can also be used such as:

    .. code-block: python

        >>> r1.reactions[pt10][::2]     # returns dx, dz, ry

    Data index and corresponding nodal force/moment.

    * 0 = fx
    * 1 = fy
    * 2 = fz
    * 3 = mx
    * 4 = my
    * 5 = mz
    """

    def __init__(self, app):
        self._app = app
        self._force = Force()
        self._moment = Moment()

    def __getitem__(self, item):
        """Get nodal results."""
        ndof = 6    # degrees of freedom per node
        if isinstance(item, Point):
            try:
                model = self._app.models.active_object
                idxi = list(model.points).index(item)
                niqi, niqj = idxi*ndof, idxi*ndof + ndof

                # note that .results is a column vector
                return self.results[niqi:niqj, 0]
            except:
                raise ValueError("node is not in results array")

        else:
            raise ValueError("value not found")

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


class Forces:
    """Member internal forces and moments.

    .. note::
       The translation and rotation components are separated because they have
       different units.

    .. code-block: python

        >>> r1.forces[pt10]          # returns all comps
        >>> r1.forces[pt10][0:3]     # returns dx, dy, dz
        >>> r1.forces[pt10][3:6]     # returns rx, ry, rz

    Slices with strides can also be used such as:

    .. code-block: python

        >>> r1.forces[pt10][::2]     # returns dx, dz, ry

    Data index and corresponding nodal force/moment.

    * 0 = fx
    * 1 = fy
    * 2 = fz
    * 3 = mx
    * 4 = my
    * 5 = mz
    """

    def __init__(self, app):
        self._app = app
        self._force = Force()
        self._moment = Moment()

    def __getitem__(self, item):
        """Get nodal results."""
        ndof = 6    # degrees of freedom per node
        if isinstance(item, Point):
            try:
                model = self._app.models.active_object
                idxi = list(model.points).index(item)
                niqi, niqj = idxi*ndof, idxi*ndof + ndof

                # note that .results is a column vector
                return self.results[niqi:niqj, 0]
            except:
                raise ValueError("node is not in results array")

        else:
            raise ValueError("value not found")

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


@units.define(_values="stress")
class Stress:
    """Element hoop stress due to pressure"""

    @property
    def results(self):
        return self._values

    @results.setter
    def results(self, data):
        """Convert row to column vector"""
        self._values = data.flatten().reshape((-1, 1))


@units.define(_values="stress")
class Shoop(Stress):
    """Element hoop stress due to pressure"""
    pass


@units.define(_values="stress")
class Slp(Stress):
    """Element longitudinal stress due to pressure"""
    pass


@units.define(_values="stress")
class Slb(Stress):
    """Element longitudinal stress due to bending"""
    pass


@units.define(_values="stress")
class Sl(Stress):
    """Total element longitudinal stress"""
    pass


@units.define(_values="stress")
class Stor(Stress):
    """Element torsional stress"""
    pass


@units.define(_values="stress")
class Sax(Stress):
    """Element axial (F/A) stress"""
    pass


@units.define(results="stress")
class Sallow:
    """Element code stress at a node"""
    pass


class Stresses:
    """A table of nodal stresses.

    Each row of the table corresponds to a different node.
    """

    def __init__(self, app):
        self._app = app
        self._shoop = Shoop()
        self._sax = Sax()
        self._stor = Stor()
        self._slp = Slp()
        self._slb = Slb()
        self._sl = Sl()
        self._sifi = None
        self._sifo = None
        self._sallow = Sallow()
        self._sratio = None
        self._scode = None

    def __getitem__(self, item):
        """Get nodal results."""
        ndof = 6    # degrees of freedom per node
        if isinstance(item, Point):
            try:
                model = self._app.models.active_object
                idxi = list(model.points).index(item)
                niqi, niqj = idxi*ndof, idxi*ndof + ndof

                # note that .results is a column vector
                return self.results[niqi:niqj, 0]
            except:
                raise ValueError("node is not in results array")

        else:
            raise ValueError("value not found")

    @property
    def results(self):
        shoop = self._shoop.results
        sax = self._sax.results
        stor = self._stor.results
        slp = self._slp.results
        slb = self._slb.results
        sl = self._sl.results
        sifi = self._sifi
        sifo = self._sifo
        sallow = self._sallow.results
        sratio = self._sratio
        scodes = self._scodes

        data = zip(shoop, sax, stor, slp, slb, sl, sifi, sifo, sallow, sratio,
                   scodes)

        values = list(data)

        return values

    @results.setter
    def results(self, data):
        results, codes = data

        self._shoop.results = results[:, 0]
        self._sax.results = results[:, 1]
        self._stor.results = results[:, 2]
        self._slp.results = results[:, 3]
        self._slb.results = results[:, 4]
        self._sl.results = results[:, 5]
        self._sifi = results[:, 6]
        self._sifo = results[:, 7]
        self._sallow.results = results[:, 8]
        self._sratio = results[:, 9]
        self._scodes = codes


class BaseCase(Entity):

    def __init__(self, name, stype="sus"):
        super(BaseCase, self).__init__(name)
        self._stype = stype  # HRG, HYD, SUS, EXP, OCC, OPE, FAT

    @property
    def stype(self):
        return self._stype

    @property
    def points(self):
        return list(self.app.points)

    @property
    def parent(self):
        return self.app.loadcases


class LoadCase(BaseCase):
    """A set of primary load cases consisting of different types of loads.

    .. note::
        The loads for a given case must be unique. In other words, the same
        load cannot be specified twice. An equality check is performed based
        on the load type, name and operating case it belongs to.
    """

    def __init__(self, name, stype="sus", loadtypes=[], opercases=[]):
        super(LoadCase, self).__init__(name, stype)
        self._loads = OrderedSet()
        self._loadtypes = loadtypes
        self._opercases = opercases

        for load in zip(loadtypes, opercases):
            self._loads.add(load)

        # results objects
        self._movements = Movements(self.app)
        self._reactions = Reactions(self.app)
        self._forces = Forces(self.app)
        self._stresses = Stresses(self.app)

    @property
    def loadtypes(self):
        return tuple(self._loadtypes)

    @property
    def opercases(self):
        return tuple(self._opercases)

    @property
    def movements(self):
        """Return the nodal displacement results array."""
        return self._movements

    @property
    def reactions(self):
        """Return the nodal reaction results array."""
        return self._reactions

    @property
    def forces(self):
        """Return the force reaction results array."""
        return self._forces

    @property
    def stresses(self):
        """Return the nodal stress results array."""
        return self._stresses

    @property
    def label(self):
        """Return the loadcase label used in the reports."""
        lbl = []
        for loadtype, opercase in zip(self._loadtypes, self._opercases):
            lbl.append("%s[%s]" % (loadtype.label, opercase))

        return " + ".join(lbl)


class LoadComb(BaseCase):
    """Combine primary loadcases using different combination methods.

    .. note::
        Combinations pull stored data from loadcases on the fly and do the
        necessary combination operations.

    .. attention::
        A loadcase and a loadcomb are derived from the same basecase and so
        they have the same namespace when it comes to name.
    """

    def __init__(self, name, stype="ope", method="algebraic", loadcases=[],
                 factors=[]):
        """Create a loadcomb instance.

        Parameters
        ----------
        name : str
            Unique name for load combination object.

        stype : str
            Type of code stress. Defaults to sustained stress.

                * HGR - Hanger load case
                * HYD - Hydro load case
                * SUS - Sustained stress case
                * EXP - Thermal expansion stress.case
                * OCC - Occasional stress case
                * OPE - Operating stress case
                * FAT - Fatigue stress case

        method : str
            Result combination method.

                * Algebriac - Disp/force results added vectorially. Stresses
                  are derived from the vectorially added force results.
                * Scalar - Disp/force results added vectorially similar to the
                  algebraic method. Stresses are added together.
                * SRSS - Square root of the sum squared. Direction independant.
                * Abs - Absolute summation.
                * Signmax - Signed max.
                * Signmin - Signed min.

        loadcases : list of loadcases.
            List of load cases.

        factors : list of numbers.
            A list of factors corresponding to each loadcase.

            .. note::
               If a factor is not given, a default value of 1 is used. Also,
               the number of factors must match the number of loadcases.
        """
        super(LoadComb, self).__init__(name, stype)
        self._method = method
        self._loadcases = OrderedSet()
        self._factors = factors

        for loadcase in loadcases:
            if isinstance(loadcase, LoadCase):
                self._loadcases.add(loadcase)

        # results objects
        self._movements = Movements(self.app)
        self._reactions = Reactions(self.app)
        self._forces = Forces(self.app)
        self._stresses = Stresses(self.app)

    @property
    def method(self):
        return self._method

    @property
    def loadcases(self):
        return tuple(self._loadcases)

    @property
    def factors(self):
        return tuple(self._factors)

    @property
    def movements(self):
        """Return the combined nodal displacement array."""
        self._movements.results = np.zeros(len(self.points) * 6, dtype=np.float64)

        for factor, loadcase in zip_longest(self._factors, self._loadcases,
                                            fillvalue=1):
            if self._method in ("algebraic", "scalar"):
                self._movements.results += (factor * loadcase.movements.results)
            elif self._method == "srss":
                # note: sign of factor has no effect, always positive
                self._movements.results += (factor * loadcase.movements.results)**2
            elif self._method == "abs":
                self._movements.results += (factor * np.abs(loadcase.movements.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._movements.results = np.maximum(self._movements.results,
                                                     loadcase.movements.results)
            elif self._method == "signmin":
                self._movements.results = np.minimum(self._movements.results,
                                                     loadcase.movements.results)

        if self._method == "srss":
            # this should be done once
            self._movements.results = np.sqrt(self._movements.results)

        return self._movements

    @property
    def reactions(self):
        """Return the combined nodal reaction array."""
        self._reactions.results = np.zeros(len(self.points) * 6, dtype=np.float64)

        for factor, loadcase in zip_longest(self._factors, self._loadcases,
                                            fillvalue=1):
            if self._method in ("algebraic", "scalar"):
                self._reactions.results += (factor * loadcase.reactions.results)
            elif self._method == "srss":
                # note: sign of factor has no effect, always positive
                self._reactions.results += (factor * loadcase.reactions.results)**2
            elif self._method == "abs":
                self._reactions.results += (factor * np.abs(loadcase.reactions.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._reactions.results = np.maximum(self._reactions.results,
                                                     loadcase.movements.results)
            elif self._method == "signmin":
                self._reactions.results = np.minimum(self._reactions.results,
                                                     loadcase.movements.results)

        if self._method == "srss":
            self._reactions.results = np.sqrt(self._reactions.results)

        return self._reactions

    @property
    def forces(self):
        """Return the combined nodal forces array."""
        self._forces.results = np.zeros(len(self.points) * 6, dtype=np.float64)

        for factor, loadcase in zip_longest(self._factors, self._loadcases,
                                            fillvalue=1):
            if self._method in ("algebraic", "scalar"):
                self._forces.results += (factor * loadcase.forces.results)
            elif self._method == "srss":
                # note: sign of factor has no effect, always positive
                self._forces.results += (factor * loadcase.forces.results)**2
            elif self._method == "abs":
                self._forces.results += (factor * np.abs(loadcase.forces.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._forces.results = np.maximum(self._forces.results,
                                                  loadcase.movements.results)
            elif self._method == "signmin":
                self._forces.results = np.minimum(self._forces.results,
                                                  loadcase.movements.results)

        if self._method == "srss":
            self._forces.results = np.sqrt(self._forces.results)

        return self._forces

    @property
    def stresses(self):
        """Return the combined nodal stresses array.

        .. note:: Stresses are calculated based on the stress type and the
            element code.

        The individual loadcases that make up the load combination can each
        have a different allowable stress. The allowable is determined by the
        stress type attribute of the load combination.  The following logic is
        applied to determine the allowable:

            1. All loadcases with the same stress type as that of the load
            combination is considered and the one with the smallest allowable
            is used as the combined allowable stress.

            2. If the loadcases are not of the same stress type as the load
            combination, again the smallest allowable of the all the loadcases
            is taken.
        """
        self._stresses.results = np.zeros((len(self.points), 10), dtype=np.float64)

        # copy sifi, sifo, as they are unchanged between loadcases
        for factor, loadcase in zip_longest(self._factors, self._loadcases,
                                            fillvalue=1):
            if self._method == "algebraic":
                # add internal forces algebraically first then calculate code
                # stress, this is done by the solver in solver.py
                return self._stresses
            elif self._method == "scaler":
                self._stresses.results[:, :6] += (factor * loadcase.stresses.results[:, :6])
            elif self._method == "srss":
                # note: sign of factor has no effect, always positive
                self._stresses.results[:, :6] += (factor * loadcase.stresses.results[:, :6])**2
            elif self._method == "abs":
                self._stresses.results[:, :6] += (factor * np.abs(loadcase.stresses.results[:, :6]))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._stresses.results[:, :6] = np.maximum(self._stresses.results[:, :6],
                                                           loadcase.stresses.results[:, :6])
            elif self._method == "signmin":
                self._stresses.results[:, :6] = np.minimum(self._stresses.results[:, :6],
                                                           loadcase.stresses.results[:, :6])

        # final IR calculated based on combined stress and the allowable stress
        self._stresses.results[:, -1] = self._stresses.results[:, 5] / self._stresses.results[:, -2]

        if self._method == "srss":
            self._stresses.results[:, :6] = np.sqrt(self._stresses.results[:, :6])

        return self._stresses

    @property
    def label(self):
        """Return the loadcomb label used in the reports."""
        lbl = []
        for loadcase in self._loadcases:
            lbl.append("%s" % loadcase.name)

        return " + ".join(lbl)


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
