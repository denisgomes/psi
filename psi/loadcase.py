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
class Translation(object):
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
class Rotation(object):
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


class Movements(object):
    """Both the translation and rotation components of the displacement row
    vector.

    .. note::
       The translation and rotation components are separated because they have
       different units.
    """

    def __init__(self, app):
        self._app = app
        self._translation = Translation()
        self._rotation = Rotation()

    def __getitem__(self, item):
        """Get nodal results.

        .. code-block: python

            >>> r1.movements[pt10]          # returns all comps
            >>> r1.movements[pt10][0:3]     # returns dx, dy, dz
            >>> r1.movements[pt10][3:6]     # returns rx, ry, rz

        Slices with strides can also be used such as:

        .. code-block: python

            >>> r1.movements[pt10][::2]     # returns dx, dz, ry
        """
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
class Force(object):
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
class Moment(object):
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


class Reactions(object):
    """Support reaction loads.

    .. note::
       The translation and rotation components are separated because they have
       different units.
    """

    def __init__(self, app):
        self._app = app
        self._force = Force()
        self._moment = Moment()

    def __getitem__(self, item):
        """Get nodal results.

        .. code-block: python

            >>> r1.reactions[pt10]          # returns all comps
            >>> r1.reactions[pt10][0:3]     # returns dx, dy, dz
            >>> r1.reactions[pt10][3:6]     # returns rx, ry, rz

        Slices with strides can also be used such as:

        .. code-block: python

            >>> r1.reactions[pt10][::2]     # returns dx, dz, ry
        """
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


class Forces(object):
    """Member internal forces and moments.

    .. note::
       The translation and rotation components are separated because they have
       different units.
    """

    def __init__(self, app):
        self._app = app
        self._force = Force()
        self._moment = Moment()

    def __getitem__(self, item):
        """Get nodal results.

        .. code-block: python

            >>> r1.forces[pt10]          # returns all comps
            >>> r1.forces[pt10][0:3]     # returns dx, dy, dz
            >>> r1.forces[pt10][3:6]     # returns rx, ry, rz

        Slices with strides can also be used such as:

        .. code-block: python

            >>> r1.forces[pt10][::2]     # returns dx, dz, ry
        """
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


class BaseCase(Entity):

    def __init__(self, name, stype="sus"):
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

    .. note::
        The loads for a given case must be unique. In other words, the same
        load cannot be specified twice. A equality check is performed based
        on the load type, name and operating case it belongs to.
    """

    def __init__(self, name, stype="sus", loadtypes=[], opercases=[]):
        super(LoadCase, self).__init__(name, stype)
        self.loads = OrderedSet()

        for load in zip(loadtypes, opercases):
            self.loads.add(load)

        # results objects
        self._movements = Movements(self.app)
        self._reactions = Reactions(self.app)
        self._forces = Forces(self.app)

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
        raise NotImplementedError("implement")

    @property
    def label(self):
        """Return the loadcase label used in the reports."""
        lbl = []
        for loadtype, opercase in self.loads:
            lbl.append("%s[%s]" % (loadtype.label, opercase))

        return " + ".join(lbl)


class LoadComb(BaseCase):
    """Combine loadcases using different combination methods.

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
                  are derived from the force results.
                * Scalar - Disp/force/ results added vectorially similar to the
                  algebraic method. Stresses are added together.
                * SRSS - Square root of the sum squared.
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
        self.method = method
        self.loadcases = OrderedSet()
        self.factors = factors

        for loadcase in loadcases:
            if isinstance(loadcase, LoadCase):
                self.loadcases.add(loadcase)

    @property
    def movements(self):
        """Return the combined nodal displacement array."""
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
        """Return the combined nodal reaction array."""
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
        """Return the combined nodal forces array."""
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
        """Return the combined nodal stresses array.

        .. note:: Stresses are calculated based on the stress type and the
           element code.
        """
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
