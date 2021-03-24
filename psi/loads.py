# Pipe Stress Infinity (PSI) - The pipe stress analysis and design software.
# Copyright (c) 2020 Denis Gomes

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

    >>> t1 = Thermal('t1', 1, 500)
    >>> t1.apply(element_list)

For multiple loads use the LoadContainer apply method:

.. code-block:: python

    >>> loads.apply([L1, L2, ..., LN], element_list)

.. note::
    An element can have multiple loads of the same type however each load must
    be of a different operating case. This applies to loads that are operating
    case dependent such as Thermal and Pressure. Forces on the other hand do
    not have to be opercase dependent. This is not enforced by the program at
    the moment so an element can have two Weight loads defined for the same
    operating case, in which the effect of the weights will added together.
"""

from __future__ import division

import csv
from math import pi
import sys
from contextlib import redirect_stdout

import numpy as np
from numpy import sin, cos, arccos

import psi
from psi.entity import Entity, EntityContainer
from psi.elements import Rigid
from psi import units


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

    Includes pipe, insulation, cladding and refractory. The weight load is
    applied as an uniform load across the length of the entire element. The
    direction of the load vector is always down with respect to the vertical.

    The load vector is defined in local coordinates such that when it is
    transformed to global coordinates, the two components of the local gravity
    vector will resolve to the one component global gravity force vector
    parallel to the vertical down direction.
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
            return 0.0

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
        L = element.length

        # the vertical direction
        vert = self.app.models.active_object.settings.vertical

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        w = self.total(element) / L

        if vert == "y":
            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 1., 0.], dtype=np.float64)

            # angle element makes with vertical
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            f[:, 0] = [-w*L*cost/2, -w*L*sint/2, 0, 0, 0, -w*L**2*sint/12,
                       -w*L*cost/2, -w*L*sint/2, 0, 0, 0, w*L**2*sint/12]

        elif vert == "z":
            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 0., 1.], dtype=np.float64)

            # angle element makes with vertical
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            f[:, 0] = [-w*L*cost/2, -w*L*sint/2, 0, 0, 0, -w*L**2*sint/12,
                       -w*L*cost/2, -w*L*sint/2, 0, 0, 0, w*L**2*sint/12]

        return f


@units.define(pres="pressure")
class Pressure(Load):
    """Internal or external pressure.

    The pressure can have a stiffening effect on the piping system called the
    Bourdon effect. For large D/t ratio pipe bends with large system pressures,
    the effects of ovalization can be made worse due to this effect. When a
    straight pipe is pressurized it wants to shrink axially due to the radial
    growth whereas a pipe bend wants to open up.

    The pressure on capped ends also has the effect of pulling on the pipe
    axially. The sustained stress in most piping codes uses the internal
    pressure to calculate a longitudinal stress which is added to the bending
    stress. For B31.1 piping for example, the longitudinal stress due to
    pressure is equal to P*D/4*t. If pressure stress per code is accounted for
    the pressure due to capped ends need not be considered as doing so will in
    effect double the pressure stress.

    A hoop stress can also be calculated to determine if the pipe wall is sized
    properly based on code requirements.
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
        """
        pipe = element.section

        od = pipe.od
        id_ = pipe.od - 2*pipe.thk  # inner diameter

        return self.pres * np.pi * (od - id_)**2 / 4

    def bourdon(self, element):
        """The force due to the Bourdon effect.

        The axial shortening of piping due to the internal pressure. This is a
        byproduct of the poisson's ratio effect. This effect is present for
        underground and unrestrained piping.
        """
        pipe = element.section
        mat = element.material

        od = pipe.od
        id_ = pipe.od - 2*pipe.thk
        nu = mat.nu.value
        area = pipe.area

        return 2 * nu * self.pres * area * (id_**2 / (od**2-id_**2))

    def flocal(self, element):
        """If pressure thrust or bourdon effects are activated, a force is
        applied in the axial direction of the piping similar to how thermal
        loads are applied.
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

        f[:, 0] = [-ft+fb, 0, 0, 0, 0, 0,
                   ft-fb, 0, 0, 0, 0, 0]

        return f


@units.define(pres="pressure")
class Hydro(Pressure):
    """Hydro test pressure.

    Test pressure is typically 1.5 times the design pressure. Hydro testing is
    performed to ensure there are no leaks. Pneumatic testing can also be used
    along with RT (x-ray).

    During hydrotesting, spring cans are locked using travel stops so that they
    behave as full Y supports. Springs can bodies are designed to sustain any
    additional deadweight load imposed during testing. Extra precaution should
    be taken to ensure the loads are acceptable for very large piping.
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

        # with redirect_stdout(sys.__stdout__):
        #     print(E)
        #     print(delT)
        #     print(alp)

        f = np.zeros((12, 1), dtype=np.float64)

        f[:, 0] = [-fa, 0, 0, 0, 0, 0,
                   fa, 0, 0, 0, 0, 0]

        return f


@units.define(rho="density", gfac="g_load")
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
    def from_file(cls, name, fluid, fname=None):
        if fname is None:
            fname = psi.FLUID_DATA_FILE

        with open(fname, "r") as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if row["fluid"] == fluid:
                    rho = float(row["rho"])
                    return cls(name, rho)
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
            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 1., 0.], dtype=np.float64)

            # angle element makes with vertical
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            f[:, 0] = [-w*L*cost/2, -w*L*sint/2, 0, 0, 0, -w*L**2*sint/12,
                       -w*L*cost/2, -w*L*sint/2, 0, 0, 0, w*L**2*sint/12]

        elif vert == "z":
            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 0., 1.], dtype=np.float64)

            # angle element makes with vertical
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            f[:, 0] = [-w*L*cost/2, -w*L*sint/2, 0, 0, 0, -w*L**2*sint/12,
                       -w*L*cost/2, -w*L*sint/2, 0, 0, 0, w*L**2*sint/12]

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


@units.define(ux="uniform_load", uy="uniform_load", uz="uniform_load")
class Uniform(Load):
    """Generic uniform load with respect to the global coordinate system.

    Parameters
    ----------
    name : str
        Unique name for load.

    opercase : int
        Operating case the load belongs to.

    ux : float
        Gravity load factor in the global x direction.

    uy : float
        Gravity load factor in the global y direction.

    uz : float
        Gravity load factor in the global z direction.

    projected : bool
        The load is applied over the projected length of the piping element in
        the respective global axes. Default is set to True.
    """

    label = "U"

    def __init__(self, name, opercase, ux=0, uy=0, uz=0, projected=True):
        super(Uniform, self).__init__(name, opercase)
        self.ux = ux
        self.uy = uy
        self.uz = uz
        self.projected = projected

    def flocal(self, element):
        L = element.length

        # combined force vector
        f = np.zeros((12, 1), dtype=np.float64)

        # directional force vectors
        fx = np.zeros((12, 1), dtype=np.float64)
        fy = np.zeros((12, 1), dtype=np.float64)
        fz = np.zeros((12, 1), dtype=np.float64)

        if self.ux:
            wx = self.ux

            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([1., 0., 0.], dtype=np.float64)

            # angle element makes with global x
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            # projected length in global x direction
            if self.projected:
                L = L*sint

            fx[:, 0] = [-wx*L*cost/2, -wx*L*sint/2, 0, 0, 0, -wx*L**2*sint/12,
                        -wx*L*cost/2, -wx*L*sint/2, 0, 0, 0, wx*L**2*sint/12]

        if self.uy:
            wy = self.uy

            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 1., 0.], dtype=np.float64)

            # angle element makes with global y
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            # projected length in global x direction
            if self.projected:
                L = L*sint

            fy[:, 0] = [-wy*L*cost/2, -wy*L*sint/2, 0, 0, 0, -wy*L**2*sint/12,
                        -wy*L*cost/2, -wy*L*sint/2, 0, 0, 0, wy*L**2*sint/12]

        if self.uz:
            wz = self.uz

            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 0., 1.], dtype=np.float64)

            # angle element makes with global y
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            # projected length in global z direction
            if self.projected:
                L = L*sint

            fz[:, 0] = [-wz*L*cost/2, -wz*L*sint/2, 0, 0, 0, -wz*L**2*sint/12,
                        -wz*L*cost/2, -wz*L*sint/2, 0, 0, 0, wz*L**2*sint/12]

        # total uniform sum
        f = fx + fy + fz

        return f


@units.define(gx="g_load", gy="g_load", gz="g_load")
class Seismic(Load):
    """Three directional seismic loading applied as a gravity (g) factor in
    global coordinates.

    Parameters
    ----------
    name : str
        Unique name for load.

    opercase : int
        Operating case the load belongs to.

    gx : float
        Gravity load factor in the global x direction.

    gy : float
        Gravity load factor in the global y direction.

    gz : float
        Gravity load factor in the global z direction.
    """

    label = "S"

    def __init__(self, name, opercase, gx=0.0, gy=0.0, gz=0.0):
        super(Seismic, self).__init__(name, opercase)
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

        # directional force vectors
        fx = np.zeros((12, 1), dtype=np.float64)
        fy = np.zeros((12, 1), dtype=np.float64)
        fz = np.zeros((12, 1), dtype=np.float64)

        # the gfac is 1 by default for earth gravity
        wx = (self.total(element)*self.gx) / L
        wy = (self.total(element)*self.gy) / L
        wz = (self.total(element)*self.gz) / L

        if self.ux:
            wx = self.ux

            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([1., 0., 0.], dtype=np.float64)

            # angle element makes with global x
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            fx[:, 0] = [-wx*L*cost/2, -wx*L*sint/2, 0, 0, 0, -wx*L**2*sint/12,
                        -wx*L*cost/2, -wx*L*sint/2, 0, 0, 0, wx*L**2*sint/12]

        if self.uy:
            wy = self.uy

            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 1., 0.], dtype=np.float64)

            # angle element makes with global y
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            fy[:, 0] = [-wy*L*cost/2, -wy*L*sint/2, 0, 0, 0, -wy*L**2*sint/12,
                        -wy*L*cost/2, -wy*L*sint/2, 0, 0, 0, wy*L**2*sint/12]

        if self.uz:
            wz = self.uz

            a = (np.array(element.to_point.xyz, dtype=np.float64) -
                 np.array(element.from_point.xyz, dtype=np.float64))
            b = np.array([0., 0., 1.], dtype=np.float64)

            # angle element makes with global y
            theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
            sint = sin(theta)
            cost = cos(theta)

            fz[:, 0] = [-wz*L*cost/2, -wz*L*sint/2, 0, 0, 0, -wz*L**2*sint/12,
                        -wz*L*cost/2, -wz*L*sint/2, 0, 0, 0, wz*L**2*sint/12]

        # total uniform sum
        f = fx + fy + fz

        return f


class Wind(Load):
    """Wind load applied as an uniform load.

    The pressure due to wind is applied as a uniform force. It is a function
    of the pipe elevation. A pressure profile versus elevation table is used
    to determine the pressure at node i and j of a piping element. Then the
    computed average of the two pressures is applied over the entire element
    projected length as an uniform load.

    Parameters
    ----------
    name : str
        Unique name for load.

    opercase : int
        Operating case the load belongs to.

    profile : list
        Wind pressure profile. List of (elevation, pressure) tuples.

    shape : float
        The element wind shape factor. Set to typical default value of 0.65.

    direction : tuple
        The wind unit vector direction given in (x, y, z) tuple.
    """

    label = "W"

    def __init__(self, name, opercase, profile, direction=(1, 0, 0),
                 shape=0.65):
        super(Wind, self).__init__(name, opercase)
        self.profile = profile
        self.direction = direction
        self.shape = shape

    @classmethod
    def from_ASCE7(cls, name, **kwargs):
        """Generate a profile and define shape factor based on ASCE criteria.
        """
        raise NotImplementedError

    def flocal(self, element):
        f = np.zeros((12, 1), dtype=np.float64)

        L = element.length

        # the vertical direction
        vert = self.app.models.active_object.settings.vertical

        from_point = np.array(element.to_point.xyz, dtype=np.float64)
        to_point = np.array(element.from_point.xyz, dtype=np.float64)

        # angle element makes with wind direction
        a = to_point - from_point
        b = np.array(self.direction, dtype=np.float64)
        theta = arccos(a.dot(b) / (np.linalg.norm(a)*np.linalg.norm(b)))
        sint = sin(theta)
        cost = cos(theta)
        L = L*sint  # projected length

        if vert == "y":
            elev1 = from_point[1]
            elev2 = to_point[1]
            pressure = (self.profile[elev1] + self.profile[elev2]) / 2
        elif vert == "z":
            elev1 = from_point[2]
            elev2 = to_point[2]
            pressure = (self.profile[elev1] + self.profile[elev2]) / 2
        w = pressure * L    # wind pressure

        f[:, 0] = [-w*L*cost/2, -w*L*sint/2, 0, 0, 0, -w*L**2*sint/12,
                   -w*L*cost/2, -w*L*sint/2, 0, 0, 0, w*L**2*sint/12]

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
        self.Uniform = Uniform
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
