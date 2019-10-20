"""Applying loads to model elements.

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

import csv
import copy

import openpipe
from openpipe.core.entity import Entity, EntityContainer
from openpipe.core.units import units
from openpipe.utils.orderedset import OrderedSet


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

    def __repr__(self):
        return "%s(name='%s')" % (self.type, self.name)


@units.define(gfac="g_load")
class Weight(Load):

    def __init__(self, name, gfac=1.0):
        super(Weight, self).__init__(name)
        self.gfac = gfac


@units.define(pres="pressure")
class Pressure(Load):

    def __init__(self, name, pres):
        super(Pressure, self).__init__(name)
        self.pres = pres


@units.define(pres="pressure")
class Hydro(Load):
    """Hydro test pressure"""

    def __init__(self, name, pres):
        super(Pressure, self).__init__(name)
        self.pres = pres


@units.define(temp="temperature", tref="temperature")
class Thermal(Load):
    """Thermal expansion load"""

    def __init__(self, name, temp, tref=70):
        super(Thermal, self).__init__(name)
        self.temp = temp
        self.tref = tref


@units.define(rho="density")
class Fluid(Load):
    """Contents load"""

    def __init__(self, name, rho):
        super(Fluid, self).__init__(rho)
        self.rho = rho

    @classmethod
    def from_file(cls, name, fluid, fname=None):
        if fname is None:
            fname = openpipe.FLUID_DATA_FILE

        with open(fname, "r") as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if row["fluid"] == fluid:
                    rho = float(row["rho"])
                    return cls(name, rho)
            else:
                return None


@units.define(fx="force", fy="force", fz="force",
              mx="moment_in", my="moment_in", mz="moment_in")
class Force(Load):
    """A generic force load"""

    def __init__(self, name, point, fx=None, fy=None, fz=None,
                 mx=None, my=None, mz=None):
        super(Force, self).__init__(name)
        self.point = point
        self.fx = fx
        self.fy = fy
        self.fz = fz
        self.mx = mx
        self.my = my
        self.mz = mz


@units.define(dx="length", dy="length", dz="length",
              mx="rotation", my="rotation", mz="rotation")
class Displacement(Load):
    """A generic displacement load"""

    def __init__(self, name, point, dx=None, dy=None, dz=None,
                 rx=None, ry=None, rz=None):
        super(Displacement, self).__init__(name)
        self.point = point
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rx = rx
        self.ry = ry
        self.rz = rz


@units.define(ux="g_load", uy="g_load", uz="g_load")
class Uniform(Load):
    """Generic uniform load"""

    def __init__(self, name, ux, uy, uz):
        super(Uniform, self).__init__(name)
        self.ux = ux
        self.uy = uy
        self.uz = uz


class Seismic(Uniform):
    """One directional seismic load applied as an uniform load"""

    @classmethod
    def from_ASCE7(cls, name, **kwargs):
        raise NotImplementedError


class Wind(Uniform):
    """One directional wind load applied as an uniform load"""

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
        self.Displacement = Displacement
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
            elements.append(self.app.elements.active_objects)

        for element in elements:
            element.loads.add(loads)
