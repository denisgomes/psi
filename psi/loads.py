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

import csv

import psi
from psi.entity import Entity, EntityContainer
from psi.units import units


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

    def flocal(self):
        """Return the element local force vector for the load."""
        raise NotImplementedError("implement")

    def fglobal(self, T):
        """Return the element global force vector for the load."""
        raise NotImplementedError("implement")


@units.define(gfac="g_load")
class Weight(Load):

    def __init__(self, name, gfac=1.0):
        super(Weight, self).__init__(name)
        self.gfac = gfac


@units.define(pres="pressure")
class Pressure(Load):

    def __init__(self, name, pres=0):
        super(Pressure, self).__init__(name)
        self.pres = pres


@units.define(pres="pressure")
class Hydro(Load):
    """Hydro test pressure"""

    def __init__(self, name, pres=0):
        super(Pressure, self).__init__(name)
        self.pres = pres


@units.define(temp="temperature", tref="temperature")
class Thermal(Load):
    """Thermal expansion load"""

    def __init__(self, name, temp, tref=0):
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
            fname = psi.FLUID_DATA_FILE

        with open(fname, "r") as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if row["fluid"] == fluid:
                    rho = float(row["rho"])
                    return cls(name, rho)
            else:
                return None


@units.define(fx="force", fy="force", fz="force",
              mx="moment_input", my="moment_input", mz="moment_input")
class Force(Load):
    """A generic force load"""

    def __init__(self, name, point, fx=0, fy=0, fz=0, mx=0, my=0, mz=0):
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

    def __init__(self, name, point, dx=0, dy=0, dz=0, rx=0, ry=0, rz=0):
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

    def __init__(self, name, ux=0, uy=0, uz=0):
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

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.loads.update(loads)
