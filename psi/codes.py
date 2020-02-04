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

"""Implementation of different piping codes"""

import math

from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.elements import (Run, Bend, Reducer, Rigid, Valve, Flange)
from psi.loads import (Weight, Pressure, Thermal)


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
        the 1950s. These factors were empiracally determined by subject piping
        to alternating loads and determining the effect it had over a number of
        cycles.

        By defination a SIF is the peak stress over the average stress of a
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

    def Shoop(self, element, loadcase):
        """Element hoop stress due to pressure"""
        raise NotImplementedError("implement")

    def Slp(self, element, loadcase):
        """Element longitudinal stress due to pressure"""
        raise NotImplementedError("implement")

    def Slb(self, element, loadcase):
        """Element longitudinal stress due to bending"""
        raise NotImplementedError("implement")

    def Stor(self, element, loadcase):
        """Element torsional stress"""
        raise NotImplementedError("implement")

    def Sax(self, element, loadcase):
        """Element axial stress due to structural loading"""
        raise NotImplementedError("implement")

    def Sallow(self, element, loadcase):
        """Element thermal stress allowable"""
        raise NotImplementedError("implement")

    def get_max_temp(self, element, loadcase):
        """Convenience function to get the largest element temperature for
        a particular operating case.
        """
        for loadtype, opercase in loadcase.loads:
            thermals = [load for load in element.loads if isinstance(load,
                        Thermal) and load.opercase==opercase]
            thermals.sort(key=lambda x: x.temp, reverse=True)

        try:
            thermal = thermals[0]
            temp = thermal.temp
        except IndexError:
            # pressure load not defined
            temp = 0

        return temp

    def get_max_pres(self, element, loadcase):
        """Convenience function to get the largest element pressure for
        a particular operating case.
        """
        # pressure load specified for loadcase and opercase sorted by maximum
        for loadtype, opercase in loadcase.loads:
            pressures = [load for load in element.loads if isinstance(load,
                         Pressure) and load.opercase==opercase]
            pressures.sort(key=lambda x: x.pres, reverse=True)

        # if length of pressures is greater than one, print warning message
        # saying multiple pressure loads exist for the same operating case for
        # the element, then proceed to take the max worst pressure

        try:
            pressure = pressures[0]
            pres = pressure.pres
        except IndexError:
            # pressure load not defined
            pres = 0

        return pres


class B311(Code):
    """B31.1 Power Piping Code Implementation.

    For each element, based on the loadcase and stress type, calculate code
    stresses.

    Each element can have one Temperature and/or Pressure load defined per one
    opercase. A load is defined to be unique using a key which checks the name,
    type and opercase.
    """

    def __init__(self, name, year='1967'):
        super(B311, self).__init__(name)
        self.year = year
        self.k = 1.15    # usage factor
        self.f = 0.80    # fatigue cycle

    def sifi(self, element):
        """In plane stress intensification factor"""
        if isinstance(element, Run):
            return 1.0

    def sifo(self, element):
        """Out of plane stress intensification factor"""
        if isinstance(element, Run):
            return 1.0

    def kfac(self, element):
        """Code flexibility factor provided for components"""
        if isinstance(element, Run):
            return 1.0

    def Shoop(self, element, loadcase):
        """Hoop stress due to pressure.

        Equal to two times the longitudinal stress due to pressure.
        """
        return 2 * self.Slp(element, loadcase)

    def Slp(self, element, loadcase):
        """Longitudinal stress due to pressure.

        Exact formulation for the pressure stress per code section 102.3.2.
        The pressure stress is a primary stress and by definition can result in
        large gross deformations and overall failure, i.e., pipe rupture.
        """
        section = element.section

        do = section.od
        di = section.id_

        p = self.get_max_pres(element, loadcase)

        # use the exact formulation
        return p*di**2 / (do**2-di**2)

    def Slb(self, element, moments):
        """Longitudinal stress due to bending moments."""
        section = element.section
        Z = section.z   # section modulus

        _, my, mz = moments
        M = math.sqrt(my**2 + mz**2)

        # note for B31.1 sifi and sifo are the same
        i = self.sifi(element)

        return M*i / Z

    def Stor(self, element, moments):
        """Shear stress due to torsion"""
        section = element.section
        do = section.od
        ip = section.ixx    # polar moment

        mx, _, _ = moments

        return mx*do / (2*ip)

    def Sax(self, element, loadcase):
        """Axial stress due to mechanical loading.

        Axial stress can be included by user defined option otherwise zero by
        default.
        """
        return 0

    def Sl(self, element, loadcase, moments):
        """Total longitudinal stress due to pressure and bending+torsion"""

        slp = self.Slp(element, loadcase)
        slb = self.Slb(element, moments)
        stor = self.Stor(element, moments)

        if loadcase.stype == "sus" or loadcase.stype == "occ":
            return slp + math.sqrt(slb**2 + 4*stor**2)

        elif loadcase.stype == "exp":
            return math.sqrt(slb**2 + 4*stor**2)

    def Sallow(self, element, loadcase):
        """Allowable stress for sustained, occasional and expansion loadcases.

        Liberal stress can be excluded for the expansion case by user defined
        option otherwise enabled by default per code.
        """
        material = element.material
        temp = self.get_max_temp(element, loadcase)

        if loadcase.stype == "sus":
            return material.sh[temp]

        elif loadcase.stype == "occ":
            # k factor is 1.15 for events lasting 8hrs or less and
            # 1.20 for 1hr or less per code para. 102.3.3, use 1.15 to
            # be conservative
            return self.k * material.sh[temp]

        elif loadcase.stype == "exp":
            pass


class CodeContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(CodeContainer, self).__init__()
        self.B311 = B311

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
