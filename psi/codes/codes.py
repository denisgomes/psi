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

"""Implementation of different piping codes"""

import math

from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.elements import (Run, Bend, Reducer, Rigid, Valve, Flange)
from psi.sifs import (Welding, Unreinforced, Reinforced, Weldolet, Sockolet,
                      Sweepolet, Weldolet, ButtWeld)
from psi.loads import Pressure, Thermal
from psi import units


class Code(Entity, ActiveEntityMixin):

    def __init__(self, name):
        super(Code, self).__init__(name)

        self.activate()

    @property
    def parent(self):
        return self.app.codes

    def activate(self):
        self.parent.activate(self)

    def apply(self, elements=None):
        self.parent.apply(self, elements)

    def sifi(self, element):
        """Element stress intensification factor.

        The code stress SIF is a fatigue factor developed by Markl and team in
        the 1950s. These factors were empiracally determined by subjecting
        piping to alternating loads and determining the effect it had over a
        number of cycles.

        By definition a SIF is the peak stress over the average stress of a
        piping component. It is a multiplier on nominal stress for typical bend
        and intersection components so that the effect of geometry and welding
        can be considered in beam analysis.

        The inplane and out of plane SIF applies to the moments and associated
        stresses in the respective directions. An inplane moment/force will
        tend to keep the component in the plane of the page. An out of plane
        moment or force will cause the component to come out of the page.
        """
        raise NotImplementedError("implement")

    def sifo(self, element):
        raise NotImplementedError("implement")

    def kfac(self, element):
        """Element flexibility factor.

        Note this factor is applied to the element for bending only. It has the
        effect of reducing the stiffness of the element in the transverse
        bending direction.
        """
        raise NotImplementedError("implement")

    def shoop(self, element, loadcase):
        """Element hoop stress due to pressure"""
        raise NotImplementedError("implement")

    def slp(self, element, loadcase):
        """Element longitudinal stress due to pressure"""
        raise NotImplementedError("implement")

    def slb(self, element, loadcase):
        """Element longitudinal stress due to bending"""
        raise NotImplementedError("implement")

    def sts(self, element, forces):
        """Element transverse shear"""
        raise NotImplementedError("implement")

    def stor(self, element, loadcase):
        """Element torsional stress"""
        raise NotImplementedError("implement")

    def sax(self, element, forces):
        """Element axial stress due to structural loading.

        F/A type stress.
        """
        raise NotImplementedError("implement")

    def sl(self, element, loadcase, point, forces):
        """Total longitudinal stress due to pressure and bending. This
        combination of stresses is also known as the code stress for most
        codes.
        """
        raise NotImplementError("implement")

    def s1(self, sl, shoop, stor, sts):
        """Maximum principal stress"""
        return (sl+shoop)*0.5 + math.sqrt(((sl-shoop)*0.5)**2 + (stor+sts)**2)

    def s2(self, sl, shoop, stor, sts):
        """Minimum principal stress"""
        return (sl+shoop)*0.5 - math.sqrt(((sl-shoop)*0.5)**2 + (stor+sts)**2)

    def max_shear(self, s1, s2):
        """Maximum shear stress"""
        return (s1-s2) * 0.5

    def sint(self, max_shear):
        """Stress intensity"""
        return 2 * max_shear

    def svon(self, s1, s2):
        """Von mises stress"""
        return math.sqrt(s1**2 + s2**2 - s1*s2)

    def sallow(self, element, loadcase):
        """Element stress allowable based on the stress type of loadcase."""
        raise NotImplementedError("implement")

    def toper(self, element, loadcase):
        """Convenience function to get the largest element temperature load for
        a particular operating case.
        """
        with units.Units(user_units="code_english"):
            thermals = []
            for loadtype, opercase in zip(loadcase.loadtypes,
                                          loadcase.opercases):
                for load in element.loads:
                    if (isinstance(load, Thermal) and
                            load.opercase == opercase):
                        thermals.append(load)

            # if length of temperature is greater than one, print warning
            # message saying multiple thermal loads exist for the same
            # operating case for the element, then proceed to take the max
            # worst pressure

            try:
                return max(thermals, key=lambda t: t.temp)
            except ValueError: # empty list
                return None

    def tmax(self, element):
        """Convenience function to get the largest element temperature from
        all operating cases.
        """
        with units.Units(user_units="code_english"):
            thermals = []
            for loadcase in element.loadcases:

                for loadtype, opercase in zip(loadcase.loadtypes,
                                              loadcase.opercases):
                    for load in element.loads:
                        if (isinstance(load, Thermal) and
                                load.opercase == opercase):
                            thermals.append(load)

            try:
                return max(thermals, key=lambda t: t.temp)
            except ValueError:
                return None

    def poper(self, element, loadcase):
        """Convenience function to get the largest element pressure load for a
        particular operating case.
        """
        with units.Units(user_units="code_english"):
            # pressure load specified for loadcase and opercase sorted by
            # maximum
            pressures = []
            for loadtype, opercase in zip(loadcase.loadtypes,
                                          loadcase.opercases):
                for load in element.loads:
                    if (isinstance(load, Pressure) and
                            load.opercase == opercase):
                        pressures.append(load)

            # if length of pressures is greater than one, print warning message
            # saying multiple pressure loads exist for the same operating case
            # for the element, then proceed to take the max worst pressure

            try:
                return max(pressures, key=lambda p: p.pres)
            except ValueError:
                return None

    def pmax(self, element):
        """Convenience function to get the largest element pressure from
        all operating cases.
        """
        with units.Units(user_units="code_english"):
            pressures = []
            for loadcase in element.loadcases:
                for loadtype, opercase in zip(loadcase.loadtypes,
                                            loadcase.opercases):
                    for load in element.loads:
                        if (isinstance(load, Pressure) and
                                load.opercase == opercase):
                            pressures.append(load)

            try:
                return max(pressures, key=lambda p: p.pres)
            except ValueError:
                return None


class CodeContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(CodeContainer, self).__init__()

        # to avoid circular import
        from .b311 import B31167

        self.B31167 = B31167

    def apply(self, code, elements=None):
        """Apply code to elements.

        Parameters
        ----------
        code : object
            A code instance.

        elements : list
            A list of elements. If elements is None, loads are applied to all
            elements.
        """
        if elements is None:
            elements = []

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.code = code
