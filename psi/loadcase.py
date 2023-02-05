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

from collections import namedtuple
from itertools import zip_longest, product
from collections import defaultdict
from contextlib import redirect_stdout
import sys

import numpy as np

from psi.entity import Entity, EntityContainer
from psi.point import Point
from psi.elements import Element
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
       The translation and rotation components are separately stored because
       they have different units.

    .. code-block: python

        >>> r1.movements[pt10]          # returns all comps
        >>> r1.movements[pt10][0:3]     # returns dx, dy, dz
        >>> r1.movements[pt10][3:6]     # returns rx, ry, rz
        >>> r1.movements[pt10].dx       # return dx
        >>> r1.movements[pt10].rz       # return rz

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

    Result = namedtuple("Result", "dx dy dz rx ry rz")

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
                return Movements.Result(*self.results[niqi:niqj, 0])
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

    def dx(self, func=None):
        """Perform an operation on the translation array - x components"""
        if func:
            return func(self._translation.results[::3])

        return self._translation.results[::3]

    def dy(self, func=None):
        """Perform an operation on the translation array - y components"""
        if func:
            return func(self._translation.results[1::3])

        return self._translation.results[1::3]

    def dz(self, func=None):
        """Perform an operation on the translation array - z components"""
        if func:
            return func(self._translation.results[2::3])

        return self._translation.results[2::3]

    def rx(self, func=None):
        """Perform an operation on the rotation array - rx components"""
        if func:
            return func(self._rotation.results[::3])

        return self._rotation.results[::3]

    def ry(self, func=None):
        """Perform an operation on the rotation array - ry components"""
        if func:
            return func(self._rotation.results[1::3])

        return self._rotation.results[1::3]

    def rz(self, func=None):
        """Perform an operation on the rotation array - rz components"""
        if func:
            return func(self._rotation.results[2::3])

        return self._rotation.results[2::3]


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
        >>> r1.reactions[pt10][0:3]     # returns fx, fy, fz
        >>> r1.reactions[pt10][3:6]     # returns mx, my, mz
        >>> r1.reactions[pt10].fx       # returns fx
        >>> r1.reactions[pt10].mz       # returns mz

    Slices with strides can also be used such as:

    .. code-block: python

        >>> r1.reactions[pt10][::2]     # returns fx, fz, my

    Data index and corresponding nodal force/moment.

    * 0 = fx
    * 1 = fy
    * 2 = fz
    * 3 = mx
    * 4 = my
    * 5 = mz
    """

    Result = namedtuple("Result", "fx fy fz mx my mz")

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
                return Reactions.Result(*self.results[niqi:niqj, 0])
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

    def fx(self, func=None):
        if func:
            return func(self._force.results[::3])

        return self._force.results[::3]

    def fy(self, func=None):
        if func:
            return func(self._force.results[1::3])

        return self._force.results[1::3]

    def fz(self, func=None):
        if func:
            return func(self._force.results[2::3])

        return self._force.results[2::3]

    def mx(self, func=None):
        if func:
            return func(self._moment.results[::3])

        return self._moment.results[::3]

    def my(self, func=None):
        if func:
            return func(self._moment.results[1::3])

        return self._moment.results[1::3]

    def mz(self, func=None):
        if func:
            return func(self._moment.results[2::3])

        return self._moment.results[2::3]


class Forces:
    """Member internal forces and moments.

    .. note::
       The translation and rotation components are separated because they have
       different units.

    .. code-block: python

        >>> r1.forces[pt10]         # returns all comps
        >>> r1.forces[pt10][0:3]    # returns fx, fy, fz
        >>> r1.forces[pt10][3:6]    # returns mx, my, mz
        >>> r1.forces[pt10].fx      # return fx
        >>> r1.forces[pt10].mz      # returns mz

    Slices with strides can also be used such as:

    .. code-block: python

        >>> r1.forces[pt10][::2]     # returns fx, fz, my

    Data index and corresponding nodal force/moment.

    * 0 = fx
    * 1 = fy
    * 2 = fz
    * 3 = mx
    * 4 = my
    * 5 = mz
    """

    Result = namedtuple("Result", "fx fy fz mx my mz")

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
                return Forces.Result(*self.results[niqi:niqj, 0])
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

    def fx(self, func=None):
        if func:
            return func(self._force.results[::3])

        return self._force.results[::3]

    def fy(self, func=None):
        if func:
            return func(self._force.results[1::3])

        return self._force.results[1::3]

    def fz(self, func=None):
        if func:
            return func(self._force.results[2::3])

        return self._force.results[2::3]

    def mx(self, func=None):
        if func:
            return func(self._moment.results[::3])

        return self._moment.results[::3]

    def my(self, func=None):
        if func:
            return func(self._moment.results[1::3])

        return self._moment.results[1::3]

    def mz(self, func=None):
        if func:
            return func(self._moment.results[2::3])

        return self._moment.results[2::3]


class MemberForces:
    """Member internal forces and moments with respect to the element local
    coordinate system.

    .. note::
       The translation and rotation components are separated because they have
       different units.

    .. code-block: python

        >>> r1.mforces[elem10]         # returns all comps
        >>> r1.mforces[elem10][0:3]    # returns fx, fy, fz node i
        >>> r1.mforces[elem10][3:6]    # returns mx, my, mz node i
        >>> r1.mforces[elem10].fxi     # return fx node i
        >>> r1.mforces[elem10].mzi     # returns mz node i

    Slices with strides can also be used such as:

    .. code-block: python

        >>> r1.mforces[elem10][::2]      # returns fxi, fzi, myi, ...

    Data index and corresponding nodal force/moment.

    * 0 = fxi
    * 1 = fyi
    * 2 = fzi
    * 3 = mxi
    * 4 = myi
    * 5 = mzi
    * 6 = fxj
    * 7 = fyj
    * 8 = fzj
    * 9 = mxj
    * 10 = myj
    * 11 = mzj
    """

    Result = namedtuple("Result", "fxi fyi fzi mxi myi mzi \
                                   fxj fyj fzj mxj myj mzj")

    def __init__(self, app):
        self._app = app
        self._forcei = Force()
        self._forcej = Force()
        self._momenti = Moment()
        self._momentj = Moment()

    def __getitem__(self, item):
        """Get element results."""
        edof = 12   # degrees of freedom per element
        if isinstance(item, Element):
            try:
                model = self._app.models.active_object
                idxi = list(model.elements).index(item)
                niqi, niqj = idxi*edof, idxi*edof + edof

                # note that .results is a column vector
                return MemberForces.Result(*self.results[niqi:niqj, 0])
            except:
                raise ValueError("element is not in results array")

        else:
            raise ValueError("value not found")

    @property
    def results(self):
        fxi = self._forcei.results[::3]
        fyi = self._forcei.results[1::3]
        fzi = self._forcei.results[2::3]
        mxi = self._momenti.results[::3]
        myi = self._momenti.results[1::3]
        mzi = self._momenti.results[2::3]

        fxj = self._forcej.results[::3]
        fyj = self._forcej.results[1::3]
        fzj = self._forcej.results[2::3]
        mxj = self._momentj.results[::3]
        myj = self._momentj.results[1::3]
        mzj = self._momentj.results[2::3]

        disp = np.array(list(zip(fxi, fyi, fzi, mxi, myi, mzi,
                                 fxj, fyj, fzj, mxj, myj, mzj)),
                                 dtype=np.float64)
        values = disp.flatten().reshape((-1, 1))

        return values

    @results.setter
    def results(self, data):
        dataij = data.reshape((-1, 6))
        datai = dataij[::2].flatten()       # elements node i
        dataj = dataij[1::2].flatten()      # elements node j
        self._forcei.results = datai
        self._momenti.results = datai
        self._forcej.results = dataj
        self._momentj.results = dataj

    def fxi(self, func=None):
        if func:
            return func(self._forcei.results[::3])

        return self._forcei.results[::3]

    def fyi(self, func=None):
        if func:
            return func(self._forcei.results[1::3])

        return self._forcei.results[1::3]

    def fzi(self, func=None):
        if func:
            return func(self._forcei.results[2::3])

        return self._forcei.results[2::3]

    def mxi(self, func=None):
        if func:
            return func(self._momenti.results[::3])

        return self._momenti.results[::3]

    def myi(self, func=None):
        if func:
            return func(self._momenti.results[1::3])

        return self._momenti.results[1::3]

    def mzi(self, func=None):
        if func:
            return func(self._momenti.results[2::3])

        return self._momenti.results[2::3]

    def fxj(self, func=None):
        if func:
            return func(self._forcej.results[::3])

        return self._forcej.results[::3]

    def fyj(self, func=None):
        if func:
            return func(self._forcej.results[1::3])

        return self._forcej.results[1::3]

    def fzj(self, func=None):
        if func:
            return func(self._forcej.results[2::3])

        return self._forcej.results[2::3]

    def mxj(self, func=None):
        if func:
            return func(self._momentj.results[::3])

        return self._momentj.results[::3]

    def myj(self, func=None):
        if func:
            return func(self._momentj.results[1::3])

        return self._momentj.results[1::3]

    def mzj(self, func=None):
        if func:
            return func(self._momentj.results[2::3])

        return self._momentj.results[2::3]


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


@units.define(_values="stress")
class S1(Stress):
    """Maximum principal stress"""
    pass


@units.define(_values="stress")
class S2(Stress):
    """Minimum principal stress"""
    pass


@units.define(_values="stress")
class MaxShear(Stress):
    """Maximum shear stress"""
    pass


@units.define(_values="stress")
class SInt(Stress):
    """Stress Intensity"""
    pass


@units.define(_values="stress")
class SVon(Stress):
    """Von Mises Stress"""
    pass


class Stresses:
    """A table of nodal stress results.

    .. note::
        Each row of the table corresponds to a different node.
    """

    Result = namedtuple("Result", ["shoop", "sax", "stor", "slp", "slb", "sl",
                                   "sifi", "sifo", "sallow", "sratio",
                                   "s1", "s2", "mshear", "sint", "svon",
                                   "scode"])

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
        self._s1 = S1()
        self._s2 = S2()
        self._mshear = MaxShear()
        self._sint = SInt()
        self._svon = SVon()
        self._scode = None

    def __getitem__(self, item):
        """Get nodal results."""
        if isinstance(item, Point):
            try:
                model = self._app.models.active_object
                idxi = list(model.points).index(item)

                # note that .results is a table, each row corresponds to node
                # the stress terms are in arrays because of units
                (shoop, sax, stor, slp, slb, sl, sifi, sifo, sallow,
                 sratio, s1, s2, mshear, sint, svon, scode) = self.results[idxi]
                return Stresses.Result(shoop[0], sax[0], stor[0], slp[0],
                                       slb[0], sl[0], sifi, sifo,
                                       sallow, sratio, s1[0], s2[0],
                                       mshear[0], sint[0], svon[0], scode)
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
        s1 = self._s1.results
        s2 = self._s2.results
        mshear = self._mshear.results
        sint = self._sint.results
        svon = self._svon.results
        scodes = self._scodes

        data = zip(shoop, sax, stor, slp, slb, sl, sifi, sifo, sallow, sratio,
                   s1, s2, mshear, sint, svon, scodes)
        values = list(data)     # mixed datatype

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
        self._s1.results = results[:, 10]
        self._s2.results = results[:, 11]
        self._mshear.results = results[:, 12]
        self._sint.results = results[:, 13]
        self._svon.results = results[:, 14]
        self._scodes = codes

    def shoop(self, func=None):
        if func:
            return func(self._shoop.results)

        return self._shoop.results

    def sax(self, func=None):
        if func:
            return func(self._sax.results)

        return self._sax.results

    def stor(self, func=None):
        if func:
            return func(self._stor.results)

        return self._stor.results

    def slp(self, func=None):
        if func:
            return func(self._slp.results)

        return self._slp.results

    def slb(self, func=None):
        if func:
            return func(self._slb.results)

        return self._slb.results

    def sl(self, func=None):
        if func:
            return func(self._sl.results)

        return self._sl.results

    def sifi(self, func=None):
        if func:
            return func(self._sifi.results)

        return self._sifi.results

    def sifo(self, func=None):
        if func:
            return func(self._sifo.results)

        return self._sifo.results

    def sratio(self, func=None):
        if func:
            return func(self._sratio.results)

        return self._sratio.results

    def sig1(self, func=None):
        if func:
            return func(self._s1.results)

        return self._s1.results

    def sig2(self, func=None):
        if func:
            return func(self._s2.results)

        return self._s2.results

    def mshear(self, func=None):
        if func:
            return func(self._mshear.results)

        return self._mshear.results

    def sint(self, func=None):
        if func:
            return func(self._sint.results)

        return self._sint.results

    def svon(self, func=None):
        if func:
            return func(self._svon.results)

        return self._svon.results


class BaseCase(Entity):

    stypes = ["hgr", "hyd", "sus", "exp", "occ", "ope", "fat"]

    def __init__(self, name, stype="sus"):
        super(BaseCase, self).__init__(name)
        self._stype = stype  # HRG, HYD, SUS, EXP, OCC, OPE, FAT

    @property
    def stype(self):
        return self._stype

    @stype.setter
    def stype(self, value):
        assert value in BaseCase.stypes, "invalid stress type"
        self._stype = value

    @property
    def points(self):
        return list(self.app.points)

    @property
    def elements(self):
        return list(self.app.elements)

    @property
    def supports(self):
        return list(self.app.supports)

    @property
    def parent(self):
        return self.app.loadcases


class LoadCase(BaseCase):
    """A set of primary load cases consisting of different types of loads and
    the operating case the load belongs to.

    Example
    -------
    Create a deadweight loadcase.

    .. code-block:: python

        >>> w1 = Weight('w1', 1)
        >>> p1 = Pressure('p1', 1)
        ...
        >>> lc1 = LoadCase('lc1', 'sus', [Weight, Pressure], [1, 1])

    .. note::
        The weight and pressure load have been defined for operating case 1.

    .. attention::
        The loads for a given case must be unique. In other words, the same
        load cannot be specified twice. An equality check is performed based
        on the load type, name and operating case it belongs to.
    """

    def __init__(self, name, stype="sus", loadtypes=[], opercases=[]):
        """Create a load case instance.

        Parameters
        ----------
        name : str
            Unique name for load case object.

        stype : str
            Type of code stress. Defaults to sustained stress.

                * HGR - Hanger load case
                * HYD - Hydro load case
                * SUS - Sustained stress case
                * EXP - Thermal expansion stress.case
                * OCC - Occasional stress case
                * OPE - Operating stress case
                * FAT - Fatigue stress case

        loadtypes : list of Load classes

        opercases : list of corresponding operating case numbers
        """
        super(LoadCase, self).__init__(name, stype)
        # self._loads = OrderedSet()
        self._loadtypes = loadtypes
        self._opercases = opercases

        # for load in zip(loadtypes, opercases):
        #     self._loads.add(load)

        # results objects
        self._movements = Movements(self.app)
        self._reactions = Reactions(self.app)
        self._forces = Forces(self.app)
        self._mforces = MemberForces(self.app)
        self._stresses = Stresses(self.app)

    def __contains__(self, load):
        """Check load in loadcase"""
        return (load.__class__ in self.loadtypes and
                load.opercase in self.opercases)

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
    def mforces(self):
        return self._mforces

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
                * Max - Largest value. Sign tacked on after taking abs.
                * Min - Smallest value. Sign tacked on after taking abs.
                * Signmax - Signed max. Sign accounted for.
                * Signmin - Signed min. Sign accounted for.

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
        self._mforces = MemberForces(self.app)
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
    def loadtypes(self):
        ltypes = []
        for loadcase in self.loadcases:
            ltypes.extend(loadcase.loadtypes)

        return tuple(ltypes)

    @property
    def opercases(self):
        opers = []
        for loadcase in self.loadcases:
            opers.extend(loadcase.opercases)

        return tuple(opers)

    @staticmethod
    def _maxfunc(a, b):
        """Return the element with the largest value not accounting for
        sign.
        """
        return a if abs(a) > abs(b) else b

    @staticmethod
    def _minfunc(a, b):
        """Return the element with the smallest value not accounting for
        sign.
        """
        return a if abs(a) < abs(b) else b

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
            elif self._method == "max":
                # overwrite every time through the loop
                self._movements.results = np.array(map(LoadComb._maxfunc,
                                                       self._movements.results,
                                                       factor * loadcase.movements.results))
            elif self._method == "min":
                self._movements.results = np.array(map(LoadComb._minfunc,
                                                       self._movements.results,
                                                       factor * loadcase.movements.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._movements.results = np.maximum(self._movements.results,
                                                     factor * loadcase.movements.results)
            elif self._method == "signmin":
                self._movements.results = np.minimum(self._movements.results,
                                                     factor * loadcase.movements.results)

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
            elif self._method == "max":
                # overwrite every time through the loop
                self._reactions.results = np.array(map(LoadComb._maxfunc,
                                                       self._reactions.results,
                                                       factor * loadcase.reactions.results))
            elif self._method == "min":
                self._reactions.results = np.array(map(LoadComb._minfunc,
                                                       self._reactions.results,
                                                       factor * loadcase.reactions.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._reactions.results = np.maximum(self._reactions.results,
                                                     factor * loadcase.reactions.results)
            elif self._method == "signmin":
                self._reactions.results = np.minimum(self._reactions.results,
                                                     factor * loadcase.reactions.results)

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
            elif self._method == "max":
                # overwrite every time through the loop
                self._forces.results = np.array(map(LoadComb._maxfunc,
                                                    self._forces.results,
                                                    factor * loadcase.forces.results))
            elif self._method == "min":
                self._forces.results = np.array(map(LoadComb._minfunc,
                                                    self._forces.results,
                                                    factor * loadcase.forces.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._forces.results = np.maximum(self._forces.results,
                                                  factor * loadcase.forces.results)
            elif self._method == "signmin":
                self._forces.results = np.minimum(self._forces.results,
                                                  factor * loadcase.forces.results)

        if self._method == "srss":
            self._forces.results = np.sqrt(self._forces.results)

        return self._forces

    @property
    def mforces(self):
        """Return the combined member forces array."""
        self._mforces.results = np.zeros(len(self.elements) * 12, dtype=np.float64)

        for factor, loadcase in zip_longest(self._factors, self._loadcases,
                                            fillvalue=1):
            if self._method in ("algebraic", "scalar"):
                self._mforces.results += (factor * loadcase.mforces.results)
            elif self._method == "srss":
                # note: sign of factor has no effect, always positive
                self._mforces.results += (factor * loadcase.mforces.results)**2
            elif self._method == "abs":
                self._mforces.results += (factor * np.abs(loadcase.mforces.results))
            elif self._method == "max":
                # overwrite every time through the loop
                self._mforces.results = np.array(map(LoadComb._maxfunc,
                                                    self._mforces.results,
                                                    factor * loadcase.mforces.results))
            elif self._method == "min":
                self._mforces.results = np.array(map(LoadComb._minfunc,
                                                    self._mforces.results,
                                                    factor * loadcase.mforces.results))
            elif self._method == "signmax":
                # overwrite every time through the loop
                self._mforces.results = np.maximum(self._mforces.results,
                                                  factor * loadcase.mforces.results)
            elif self._method == "signmin":
                self._mforces.results = np.minimum(self._mforces.results,
                                                  factor * loadcase.mforces.results)

        if self._method == "srss":
            self._mforces.results = np.sqrt(self._mforces.results)

        return self._mforces

    @property
    def stresses(self):
        """Return the combined nodal stresses array.

        .. note:: Stresses are calculated based on the stress type and the
            element code.

        The individual loadcases that make up the load combination can each
        have a different allowable stress. The allowable for the combination is
        determined by the stress type attribute of the load combination.  The
        following logic is applied to determine the allowable:

            1. All loadcases with the same stress type as that of the load
            combination is considered and the one with the smallest allowable
            is used as the combined allowable stress.

            2. If the loadcases are not of the same stress type as the load
            combination, again the smallest allowable of the all the loadcases
            is taken.
        """
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
