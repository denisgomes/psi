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

"""Applying loads to model elements and nodes.

Example
-------
Suppose you have a list of elements 'element_list' that you want to apply a
single or multiple loads to.

To apply a single load do the following:

>>> t = Thermal('T1', 500)
>>> t.apply(element_list)

For multiple loads use the container apply method:

>>> loads.apply([L1, L2, ..., LN], element_list)
"""

from __future__ import division

import csv
from math import pi

import numpy as np

import psi
from psi.entity import Entity, EntityContainer
from psi import units

from psi.settings import options

import sys
from contextlib import redirect_stdout


class Load(Entity):

    def __init__(self, name):
        super(Load, self).__init__(name)

    @property
    def parent(self):
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
    """The weight cases for each element.

    Parameters
    ----------
    name : str
        Unique name for load.

    gfac : float
        The gfac is set to 1 by default for earth gravity.
    """

    def __init__(self, name, gfac=1.0):
        super(Weight, self).__init__(name)
        self.gfac = gfac

    def pipe(self, element):
        """Pipe component weight

        Parameters
        ----------
        element : Element
            A piping element.

        Returns
        -------
        The weight of the element.
        """
        return element.mass() * self.gfac

    def insulation(self, element):
        """Insulation weight"""
        try:
            return element.insulation_mass() * self.gfac
        except:
            return 0.0

    def total(self, element):
        """Total weight of piping element"""
        return self.pipe(element) + self.insulation(element)

    def flocal(self, element):
        L = element.length

        # the vertical direction
        vert = options["core.vertical"]

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        if vert == "y":
            wy = -self.total(element) / L
            wz = 0.0

        elif vert == "z":
            wz = -self.total(element) / L
            wy = 0.0

        f[:, 0] = [0, wy*L/2, wz*L/2, 0, -wz*L**2/12, wy*L**2/12,
                   0, wy*L/2, wz*L/2, 0, wz*L**2/12, -wy*L**2/12]

        return f


@units.define(pres="pressure")
class Pressure(Load):
    """Internal or external pressure.

    The pressure can have a stiffening affect on the piping system called the
    Bourdon effect. For large D/t ratio pipes with large system pressures, the
    effects of ovalization can be made worse due to this effect. When a
    straight pipe is pressurized it wants to shrink axially. A pipe bend on the
    other hand wants to open up.

    The pressure on the capped ends also has the effect of pulling on the pipe
    axially.

    The sustained stress in most piping codes also use the internal pressure to
    calculate a longitudinal stress which is added to the bending stress.

    A hoop stress can also be calculated to determine if the pipe wall is sized
    properly based on code requirements.
    """

    def __init__(self, name, pres=0):
        super(Pressure, self).__init__(name)
        self.pres = pres

    def thrust(self, element):
        """Pressure thrust force.

        Force causing axial elongation due to the system being closed.
        """
        pipe = element.section

        od = pipe.od
        id_ = pipe.od - 2*pipe.thk  # inner diameter

        return self.pres * np.pi * (od - id_)**2 / 4

    def bourdon(self, element):
        """Bourdon effect.

        The axial shortening of piping due to the internal pressure. This is
        a byproduct of the poisson's ratio effect.
        """
        pipe = element.section
        mat = element.material

        od = pipe.od
        id_ = pipe.od - 2*pipe.thk
        nu = mat.nu.value
        area = pipe.area

        return 2 * nu * self.pres * area * (id_**2 / (od**2-id_**2))

    def flocal(self, element):
        """If bourdon effect is activated, a force is applied in the axial
        direction of the piping similar to the thermal loading.
        """
        f = np.zeros((12, 1), dtype=np.float64)

        # pressure thrust force - elongation
        ft = self.thrust(element)

        # bourdon effect - shortening
        fb = self.bourdon(element)

        f[:, 0] = [-ft+fb, 0, 0, 0, 0, 0,
                   ft-fb, 0, 0, 0, 0, 0]

        return f


@units.define(pres="pressure")
class Hydro(Load):
    """Hydro test pressure

    Test pressure is typically 1.5 times the design pressure. Hydro testing is
    performed to ensure that there are no leaks. Pneumatic testing can also be
    used along with RT (x-ray).
    """

    def __init__(self, name, pres=0):
        super(Pressure, self).__init__(name)
        self.pres = pres

    def flocal(self, element):
        raise NotImplementedError("implement")


@units.define(temp="temperature", tref="temperature")
class Thermal(Load):
    """Thermal expansion load"""

    def __init__(self, name, temp, tref):
        super(Thermal, self).__init__(name)
        self.temp = temp
        self.tref = tref

    def flocal(self, element):
        # NOTE: the value of E determines the elongation,
        # verify if this should be expected behavior
        delT = self.temp - self.tref

        alp = element.material.alp[self.temp]
        E = element.material.ymod[self.tref]
        A = element.section.area
        L = element.length

        # thermal axial force causing displacement
        # larger E produces larger force and displacement
        fa = (E*A) * (alp*delT)

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

    def __init__(self, name, rho, gfac=1.0):
        super(Fluid, self).__init__(name)
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

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        # note that for riser piping the fluid weigth in the column can be
        # larger because the water does not stick to the pipe
        wx = wy = wz = self.weight(element) / L

        f[:, 0] = [wx*L/2, wy*L/2, wz*L/2, 0, -wz*L**2/12, wy*L**2/12,
                   -wx*L/2, wy*L/2, wz*L/2, 0, wz*L**2/12, -wy*L**2/12]

        return f


@units.define(fx="force", fy="force", fz="force", mx="moment_input",
              my="moment_input", mz="moment_input")
class Force(Load):
    """A generic global force vector"""

    def __init__(self, name, point, fx=0, fy=0, fz=0, mx=0, my=0, mz=0):
        super(Force, self).__init__(name)
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
    """Generic uniform load"""

    def __init__(self, name, ux=0, uy=0, uz=0):
        super(Uniform, self).__init__(name)
        self.ux = ux
        self.uy = uy
        self.uz = uz

    def flocal(self, element):
        L = element.length

        # apply the weight as an uniform load
        f = np.zeros((12, 1), dtype=np.float64)

        wx = self.ux
        wy = self.uy
        wz = self.uz

        f[:, 0] = [wx*L/2, wy*L/2, wz*L/2, 0, -wz*L**2/12, wy*L**2/12,
                   -wx*L/2, wy*L/2, wz*L/2, 0, wz*L**2/12, -wy*L**2/12]

        return f


@units.define(gx="g_load", gy="g_load", gz="g_load")
class Seismic(Weight):
    """Three directional seismic load applied as an uniform g load"""

    def __init__(self, name, gx=0.0, gy=0.0, gz=0.0, gfac=1.0):
        super(Seismic, self).__init__(name, gfac)
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

        # the gfac is 1 by default for earth gravity
        wx = (self.total(element)*self.gx) / L
        wy = (self.total(element)*self.gy) / L
        wz = (self.total(element)*self.gz) / L

        f[:, 0] = [wx*L/2, wy*L/2, wz*L/2, 0, -wz*L**2/12, wy*L**2/12,
                   -wx*L/2, wy*L/2, wz*L/2, 0, wz*L**2/12, -wy*L**2/12]

        return f


class Wind(Uniform):
    """Wind load applied as an uniform load.

    The pressure due to wind is applied as a uniform force. It is a function
    of the pipe elevation.
    """

    @classmethod
    def from_ASCE7(cls, name, **kwargs):
        raise NotImplementedError


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

        A copy of the load is attached to each element, where each copy shares
        the same name. Multiple pipe runs can have a thermal load with a name
        'T1' but have varying temperatures. Each load corresponds to an in
        individual operating case.

        Parameters
        ----------
        loads : list
            A list of loads
        elements : list
            A list of elements. If elements is None, loads are applied to all
            elements.
        """
        if elements is None:
            elements = []

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.loads.update(loads)
