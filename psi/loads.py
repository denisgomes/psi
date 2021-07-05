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

"""Applying loads to model nodes and elements.


Example
-------
Suppose you have a list of elements 'element_list' that you want to apply a
single or multiple loads to.

To create a Thermal load for operating case 1, type:

.. code-block:: python

    >>> t1 = Thermal('T1', 1, 500)
    >>> t1.apply(element_list)

Similarly, to add a deadweight load to all active elements corresponding to
operating case 1, type:

.. code-block:: python

    >>> w1 = Weight('W1', 1)
    >>> w1.apply()

For multiple loads use the LoadContainer apply method:

.. code-block:: python

    >>> loads.apply([L1, L2, ..., LN], element_list)

.. warning::

    An element can have multiple loads of the same type however each load
    should be of a different operating case. This is not currently enforced by
    the program so an element can have two Weight loads defined for the same
    operating case in which case the weights will be added together.
"""

from __future__ import division

import csv
from math import pi
# import sys
# from contextlib import redirect_stdout

import numpy as np

import psi
from psi.entity import Entity, EntityContainer
from psi.elements import Rigid
from psi.sections import Pipe
from psi import units
from psi.units import DEFAULT_UNITS


class Load(Entity):

    def __init__(self, name, opercase):
        super(Load, self).__init__(name)
        self.opercase = opercase

    @property
    def parent(self):
        """Returns the LoadContainer instance."""
        return self.app.loads

    def apply(self, elements=None):
        """Apply the load to the elements.

        Parameters
        ----------
        elements : list
            A list of elements. If elements is None, load is applied to the
            active elements.
        """
        self.parent.apply([self], elements)

    def flocal(self, element):
        """Return the element local force vector for the load."""
        raise NotImplementedError("implement")

    def fglobal(self, element):
        """Return the element global force vector for the load."""
        T = element.T()    # build T once and reuse

        return T.transpose() @ self.flocal(element)


@units.define(gfac="g_load")
class Weight(Load):
    """The element deadweight load.

    Includes pipe, insulation, cladding and refractory weight. The weight load
    is applied as an uniform load across the length of the entire element. The
    direction of the load vector is always down with respect to the vertical.
    """

    label = "W"

    def __init__(self, name, opercase, gfac=1.0):
        """The weight load for each element.

        Parameters
        ----------
        name : str
            Unique name for load.

        opercase : int
            Operating case the load belongs to.

        gfac : float
            The gfac is set to 1 by default for earth gravity.
        """
        super(Weight, self).__init__(name, opercase)
        self.gfac = gfac

    def pipe(self, element):
        """Pipe component weight.

        Parameters
        ----------
        element : Element
            A piping element.

        Returns
        -------
        The weight of the element.
        """
        if isinstance(element, Rigid):
            return element.mass(self.gfac) * self.gfac

        return element.mass() * self.gfac

    def cladding(self, element):
        """Weight of cladding based on thickness"""
        try:
            return element.cladding_mass() * self.gfac
        except:
            return 0    # implement

    def insulation(self, element):
        """Insulation weight"""
        try:
            return element.insulation_mass() * self.gfac
        except:
            return 0

    def refractory(self, element):
        """Weight of refractory based on thickness"""
        try:
            return element.refractory_mass() * self.gfac
        except:
            return 0    # implement

    def total(self, element):
        """Total weight of piping element including pipe, cladding, insulation
        and refractory.

        If the weight is set to zero for rigid elements all terms are zero.
        """
        if isinstance(element, Rigid) and element.weight == 0:
            return 0

        return (self.pipe(element) + self.insulation(element) +
                self.cladding(element) + self.refractory(element))

    def flocal(self, element):
        """The global force vector associated with the body load is transformed
        to local coordinates and the typical textbook force and moment formulas
        for a fixed beam with uniform loading is used to create the local force
        vector. The returned local force vector is then transformed back to
        global coordinates by multiplying by the transformation matrix.
        """
        L = element.length

        # the vertical direction
        vert = self.app.models.active_object.settings.vertical

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        w = self.total(element) / L
        if vert == "y":
            # w transformed to local coord
            wl = element.dircos() @ np.array([[0], [-w], [0]])
        elif vert == "z":
            # w transformed to local coord
            wl = element.dircos() @ np.array([[0], [0], [-w]])

        wxl, wyl, wzl = wl
        f[:, 0] = [wxl*L/2, wyl*L/2, wzl*L/2, 0, -wzl*L**2/12, wyl*L**2/12,
                   wxl*L/2, wyl*L/2, wzl*L/2, 0, wzl*L**2/12, -wyl*L**2/12]

        return f


@units.define(pres="pressure")
class Pressure(Load):
    """The element pressure load.

    The calculated sustained stress in most piping codes uses the internal
    pressure to calculate a primary longitudinal stress typically equal to
    P*D/4*t which is added to the primary bending stress.

    The pressure on capped ends due to the system being closed also has the
    effect of pulling on the pipe axially and imparting a tensile thrust load.

    .. note::

        If the pressure stress per applicable piping code is accounted for,
        the pressure due to capped ends need not be considered as doing so
        will in effect double the pressure stress.

    The pressure can have a stiffening effect on the piping system called the
    Bourdon effect. For large D/t ratio pipe bends with large system pressures,
    the effects of ovalization can be made worse due to this effect. When a
    straight pipe is pressurized it wants to shrink axially due to the radial
    growth whereas a pipe bend wants to open up.

    A pressure dependent hoop stress can also be calculated and typically used
    for pipe wall sizing based on code requirements. The hoop stress is based
    on the highest pressure the system is expected to have. The final wall
    thickness is determined by adding the corrosion allowable, mill tolerance
    and also accounting for reductions due to threads. The calculated value is
    then added to the minimum thickess calculation and rounded up.
    """

    label = "P"

    def __init__(self, name, opercase, pres=0):
        """The pressure load for each element.

        Parameters
        ----------
        name : str
            Unique name for load.

        opercase : int
            Operating case the load belongs to.

        pres : float
            The pressure is set to 0 by default.
        """
        super(Pressure, self).__init__(name, opercase)
        self.pres = pres

    def thrust(self, element):
        """Force due to pressure thrust.

        Force causing axial elongation due to the system being pressurized and
        closed.

        .. note::

            The thrust force is based on the inner diameter of the pipe and the
            pipe maximum internal pressure of the system.
        """
        pipe = element.section

        id_ = pipe.od - 2*pipe.thk  # nominal inner diameter

        return (self.pres*np.pi*id_**2) / 4

    def bourdon(self, element):
        """The force due to Bourdon effect.

        The axial shortening of piping due to the internal pressure. This is a
        byproduct of the poisson's ratio effect. This effect is most apparent
        in underground and unrestrained piping.

        For steel piping the bourdon effect is negligible however it is more
        apparent for plastic piping.
        """
        pipe = element.section
        mat = element.material

        id_ = pipe.od - 2*pipe.thk
        nu = mat.nu.value

        return 2 * nu * self.pres * (np.pi*id_**2)/4

    def flocal(self, element):
        """If pressure thrust or bourdon effects are activated, a force is
        applied in the axial (local x) direction of the piping.

        The thrust force will put the element in tension and cause it to
        elongate. On the other hand, the bourdon pressure has a compressive
        shortening effect on the pipe length.

        The local force vector shows the contribution of both forces on an
        element.
        """
        f = np.zeros((12, 1), dtype=np.float64)

        # pressure thrust force - elongation
        ft = 0
        if self.app.models.active_object.settings.pressure_thrust:
            ft = self.thrust(element)

        # bourdon effect - shortening
        fb = 0
        if self.app.models.active_object.settings.bourdon_effect:
            fb = self.bourdon(element)

        f[:, 0] = [(-ft+fb)/2, 0, 0, 0, 0, 0,
                   (ft-fb)/2, 0, 0, 0, 0, 0]

        return f


@units.define(pres="pressure")
class Hydro(Pressure):
    """Hydro test pressure.

    Test pressure is typically 1.5 times the design pressure. Hydro testing is
    typically performed in the cold installed configuration to ensure there are
    no leaks.

    During testing, spring cans are locked using travel stops so that they
    behave as full Y supports. Spring can bodies are designed to support
    additional deadweight loads imposed during testing. Extra precaution should
    be taken to ensure these loads are acceptable for very large piping.

    Pneumatic testing can also be used along with RT (x-ray) to avoid
    overloading a system. The individual spool pieces can also be tested
    on the factory floor by capping both ends with a blind flange and then
    pumping it with water.

    .. warning::

        The hydro load only accounts for the sustained pressure stress due to
        test pressure. A fluid weight load should also be added in conjunction
        to account for the mechanical loading. See code below:

        .. code-block::

            >>> hp = Hydro('HP', 1, 200)
            >>> fl = Fluid.from_file('F1', 1, "water")
            >>> loads.apply([hp, fl])   # active elements
            ...
            >>> lc1 = LoadCase('L1', 'sus', [Hydro, Fluid], [1, 1])
    """

    label = "HP"

    def __init__(self, name, opercase, pres=0):
        super(Hydro, self).__init__(name, opercase)
        self.pres = pres


@units.define(temp="temperature", tref="temperature")
class Thermal(Load):
    """Thermal expansion load."""

    label = "T"

    def __init__(self, name, opercase, temp, tref):
        """The thermal load for each element.

        Parameters
        ----------
        name : str
            Unique name for load.

        opercase : int
            Operating case the load belongs to.

        temp : float
            The temperature of the element(s).

        tref : float
            The reference temperature used to calculate delta T.
        """
        super(Thermal, self).__init__(name, opercase)
        self.temp = temp
        self.tref = tref

    def flocal(self, element):
        # NOTE: the value of E determines the elongation,
        # verify if this should be expected behavior
        delT = self.temp - self.tref

        alp = element.material.alp[self.temp]
        E = element.material.ymod[self.tref]
        A = element.section.area

        # thermal axial force causing displacement
        # larger E produces larger force and displacement
        fa = (E*A) * (alp*delT)

        # thermal effect turned off for rigids
        if isinstance(element, Rigid) and not element.include_thermal:
            fa = 0

        f = np.zeros((12, 1), dtype=np.float64)

        f[:, 0] = [-fa, 0, 0, 0, 0, 0,
                   fa, 0, 0, 0, 0, 0]

        return f


@units.define(rho="fluid_density", gfac="g_load")
class Fluid(Load):
    """Contents load"""

    label = "FL"

    def __init__(self, name, opercase, rho, gfac=1.0):
        super(Fluid, self).__init__(name, opercase)
        self.rho = rho
        self.gfac = gfac

    def mass(self, element):
        """Mass of the fluid"""
        di = element.section.od - 2*element.section.thke
        fluid_area = (pi/4) * (di**2)
        mass = self.rho * fluid_area * element.geometry.length

        return mass

    def weight(self, element):
        """Weight of the fluid"""
        return self.mass(element) * self.gfac

    @classmethod
    def from_file(cls, name, opercase, fluid="water", gfac=1.0, fname=None,
                  default_units="english"):
        if fname is None:
            fname = psi.FLUID_DATA_FILE

        with open(fname, "r") as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if row["fluid"] == fluid:
                    rho = float(row["rho"])

                    with units.Units(user_units=default_units):
                        return cls(name, opercase, rho, gfac)
            else:
                return None

    def flocal(self, element):
        L = element.length

        # the vertical direction
        vert = self.app.models.active_object.settings.vertical

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        # note that for riser piping the fluid weigth in the column can be
        # larger because the water does not stick to the pipe
        w = self.weight(element) / L

        # fluid weight effect turned off for rigids
        if isinstance(element, Rigid) and not element.include_fluid:
            w = 0

        if vert == "y":
            # w transformed to local coord
            wl = element.dircos() @ np.array([[0], [-w], [0]])
        elif vert == "z":
            # w transformed to local coord
            wl = element.dircos() @ np.array([[0], [0], [-w]])

        wxl, wyl, wzl = wl
        f[:, 0] = [wxl*L/2, wyl*L/2, wzl*L/2, 0, -wzl*L**2/12, wyl*L**2/12,
                   wxl*L/2, wyl*L/2, wzl*L/2, 0, wzl*L**2/12, -wyl*L**2/12]

        return f


@units.define(fx="force", fy="force", fz="force", mx="moment_input",
              my="moment_input", mz="moment_input")
class Force(Load):
    """A generic global force vector"""

    label = "F"

    def __init__(self, name, opercase, point,
                 fx=0, fy=0, fz=0, mx=0, my=0, mz=0):
        super(Force, self).__init__(name, opercase)
        self.point = point
        self.fx = fx
        self.fy = fy
        self.fz = fz
        self.mx = mx
        self.my = my
        self.mz = mz

    def fglobal(self, element):
        """The forces are applied in the global directions"""
        f = np.zeros((12, 1), dtype=np.float64)

        fx = self.fx
        fy = self.fy
        fz = self.fz
        mx = self.mx
        my = self.my
        mz = self.mz

        if self.point == element.from_point.name:
            f[0:6, 0] = [fx, fy, fz, mx, my, mz]

        elif self.point == element.to_point.name:
            f[6:12, 0] = [fx, fy, fz, mx, my, mz]

        return f


class Seismic(Weight):
    """Three directional seismic loading applied as a gravity (g) factor in
    global coordinates.

    Parameters
    ----------
    name : str
        Unique name for load.

    opercase : int
        Operating case the load belongs to.

    gx : float
        Seismic g load factor in the global x direction. Defaults to 0.

    gy : float
        Seismic g load factor in the global y direction. Defaults to 0.

    gz : float
        Seismic g load factor in the global z direction. Defaults to 0.

    gfac: float
        Gravity factor with 1 being earth gravity. Defaults to 1.
    """

    label = "S"

    def __init__(self, name, opercase, gx=0.0, gy=0.0, gz=0.0, gfac=1.0):
        super(Seismic, self).__init__(name, opercase, gfac)
        self.gx = gx
        self.gy = gy
        self.gz = gz

    @classmethod
    def from_ASCE716(cls, name, **kwargs):
        raise NotImplementedError("implement")

    def flocal(self, element):
        L = element.length

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        # gx, gy, and gz are fractions of earth gravity
        wx = (self.total(element)*self.gx) / L
        wy = (self.total(element)*self.gy) / L
        wz = (self.total(element)*self.gz) / L

        # w transformed to local coord
        wl = element.dircos() @ np.array([[wx], [wy], [wz]])

        wxl, wyl, wzl = wl
        f[:, 0] = [wxl*L/2, wyl*L/2, wzl*L/2, 0, -wzl*L**2/12, wyl*L**2/12,
                   wxl*L/2, wyl*L/2, wzl*L/2, 0, wzl*L**2/12, -wyl*L**2/12]

        return f


@units.define(gelev="elevation")
class Wind(Load):
    """Wind force applied as uniform loading.

    The pressure due to wind is applied as a uniform force. It is a function
    of the pipe elevation. A pressure profile versus elevation table is used
    to determine the pressure at node i and j of a piping element. Then the
    computed average of the two pressures is applied over the entire element
    projected length as an uniform load.

    .. note::
        For more accurate results make sure to create a node anywhere the
        piping system crosses the ground elevation.

        If the pipe elevation at a node is less than the ground elevation the
        pressure contribution from that node is 0. If both from and to point
        elevations are less than ground, both pressures are 0 and thus the
        average pressure is also 0.

    Parameters
    ----------
    name : str
        Unique name for load.

    opercase : int
        Operating case the load belongs to.

    profile : list of tuples
        Wind pressure profile.

        List of elevation versus pressure tuples given in the format [(elev1,
        pres1), (elev2, pres2), ... (elevn, presn)] with respect to the ground
        elevation.

        .. note::
            The first elevation (elev1) corresponds to the ground elevation and
            must be set to zero. In other words, the 0 reference of the profile
            is located at ground elevation. Use the Wind.gelev attribute to set
            the global vertical position of ground.

    shape : float
        The element wind shape factor. Set to a typical default value of 0.7.

    dirvec : tuple
        The wind vector direction given in global coordinates (x, y, z).

    gelev : float
        The ground elevation with respect to global coordinates.  Elevation
        below which the wind pressure is zero.  Default is set to 0.

    is_projected : bool
        If true, the load is applied over the projected length of the element
        in the respective global direction. Default is set to True.
    """

    label = "WIN"

    @units.define(_keys="elevation", _values="wind_load")
    class Profile(object):
        """Tabular wind profile data given in list of (key, value) pairs."""

        def __init__(self, data):
            self._keys = []
            self._values = []
            self.table = data

        @property
        def _data(self):
            return list(zip(self._keys, self._values))

        @_data.setter
        def _data(self, data):
            self._keys, self._values = list(zip(*data))

        @property
        def table(self):
            """Return tabular list of (key, value) pairs"""
            return self._data

        @table.setter
        def table(self, data):
            """Set tabular list of (key, value) pairs"""
            self._data = sorted(data, key=lambda val: val[0])

        def __getitem__(self, item):
            """Key based value interpolation"""
            tvals = self._data  # keys vs. vals

            if len(tvals) == 1:
                # only one elevation entered
                return tvals[0][1]
            elif len(tvals) > 1:
                # interpolate
                sorted(tvals, key=lambda p: p[0])   # by key
                keys, vals = zip(*tvals)

                # user not expected to enter key with 6 digits
                # of precision
                keys = [round(key, 6) for key in keys]

                # pick top or bottom limit
                # if item is out of bounds
                minkey = min(keys)
                maxkey = max(keys)
                if item < minkey:
                    return vals[0]
                elif item > maxkey:
                    return vals[-1]

                # find bounding indices and interpolate
                for idx, key in enumerate(keys):
                    if key >= item:
                        uidx = idx
                        lidx = uidx - 1
                        break

                if key == item:    # exact match
                    return vals[uidx]
                else:
                    x1, y1 = keys[lidx], vals[lidx]
                    x2, y2 = keys[uidx], vals[uidx]

                    try:
                        return y1 - (y2-y1)/(x2-x1) * (x1-item)
                    except TypeError:   # None value
                        raise ValueError("pressure interpolation failed")
            else:
                raise ValueError("pressure undetermined")

        def clear(self):
            self._keys = []
            self._values = []

        def __repr__(self):
            if len(self._data) == 1:
                return str(self._data[0][1])

            return str(self.table)

    def __init__(self, name, opercase, profile=[], dirvec=(1, 0, 0),
                 shape=0.7, gelev=0, is_projected=True):
        super(Wind, self).__init__(name, opercase)
        self.profile = Wind.Profile(profile)
        self.dirvec = dirvec
        self.shape = shape
        self.gelev = gelev
        self.is_projected = is_projected

    @classmethod
    def from_ASCE7(cls, name, opercase, **kwargs):
        """Generate a profile and define shape factor based on ASCE criteria.
        """
        raise NotImplementedError

    def flocal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        L = element.length

        # wind unit direction
        d = np.array(self.dirvec, dtype=np.float64)
        du = d / np.linalg.norm(d)  # wind unit vector

        # take average pressure from node i and j
        vert = self.app.models.active_object.settings.vertical
        if vert == "y":
            elev1 = element.from_point.y
            elev2 = element.to_point.y
        elif vert == "z":
            elev1 = element.from_point.z
            elev2 = element.to_point.z
        # round to avoid difference from unit conversion
        # when doing floating point comparison
        elev1, elev2 = round(elev1, 6), round(elev2, 6)
        gelev = round(self.gelev, 6)
        # profile pressure taken with respect to ground elevation
        pres1 = 0 if elev1 <= gelev else self.profile[elev1-gelev]
        pres2 = 0 if elev2 <= gelev else self.profile[elev2-gelev]
        pavg = (pres1+pres2) / 2    # magnitude

        # with redirect_stdout(sys.__stdout__):
        #     print(elev1, self.gelev, pres1)
        #     print(elev2, self.gelev, pres2)

        # pipe and insulation diameter exposed to wind
        if type(element.section) is Pipe:
            dp = element.section.od
            if element.insulation:
                dp += (2*element.insulation.thk)

        if self.is_projected:
            sint, _ = element.angle(du)
            Ap = dp*L*sint
        else:   # normal direction
            Ap = dp*L

        # vector force due to wind
        fw = (Ap*self.shape*pavg) * du

        # force applied over element length
        fwx, fwy, fwz = fw / L

        # fw transformed to local coord
        wl = element.dircos() @ np.array([[fwx], [fwy], [fwz]])

        wxl, wyl, wzl = wl
        f[:, 0] = [wxl*L/2, wyl*L/2, wzl*L/2, 0, -wzl*L**2/12, wyl*L**2/12,
                   wxl*L/2, wyl*L/2, wzl*L/2, 0, wzl*L**2/12, -wyl*L**2/12]

        return f


class LoadContainer(EntityContainer):

    def __init__(self):
        super(LoadContainer, self).__init__()
        self.Weight = Weight
        self.Pressure = Pressure
        self.Hydro = Hydro
        self.Thermal = Thermal
        self.Fluid = Fluid
        self.Force = Force
        self.Seismic = Seismic
        self.Wind = Wind

    def apply(self, loads=[], elements=None):
        """Apply loads to elements.

        A reference for each load is assigned to each element.

        .. note::
            One pipe element can be assigned multiple loads of the same type.
            Each load must be associated with a different operating case.

        Parameters
        ----------
        loads : list
            A list of loads

        elements : list
            A list of elements. If elements is None, loads are applied to all
            active elements.
        """
        if elements is None:
            elements = []

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.loads.update(loads)
