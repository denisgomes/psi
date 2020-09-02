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

"""Implementation of different piping codes"""

import math

from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.elements import (Run, Bend, Reducer, Rigid, Valve, Flange)
from psi.loads import (Weight, Pressure, Thermal)
from psi import units
from contextlib import redirect_stdout


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

    def stor(self, element, loadcase):
        """Element torsional stress"""
        raise NotImplementedError("implement")

    def sax(self, element, forces):
        """Element axial stress due to structural loading.

        F/A type stress.
        """
        raise NotImplementedError("implement")

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
                    if isinstance(load, Thermal) and load.opercase == opercase:
                        thermals.append(load)

            thermals.sort(key=lambda x: x.temp, reverse=True)

            try:
                tload = thermals[0]
            except IndexError:
                # thermal load not defined - defaults to ambient
                tload = None

            return tload

    def tmax(self, element, loadcases):
        """Convenience function to get the largest element temperature across
        all operating cases.
        """
        raise NotImplementedError("implement")

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
                    if isinstance(load, Pressure) and load.opercase == opercase:
                        pressures.append(load)

            pressures.sort(key=lambda x: x.pres, reverse=True)

            # if length of pressures is greater than one, print warning message
            # saying multiple pressure loads exist for the same operating case
            # for the element, then proceed to take the max worst pressure

            try:
                pload = pressures[0]
            except IndexError:
                # pressure load not defined - defaults to 0
                pload = None

            return pload

    def pmax(self, element, loadcases):
        """Convenience function to get the largest element pressure across
        all operating cases.
        """
        raise NotImplementedError("implement")


class B311(Code):
    """B31.1 Power Piping Code Implementation.

    For each element, based on the loadcase and stress type, calculate code
    stresses.

    Each element can have one Temperature and/or Pressure load defined per one
    opercase. A load is defined to be unique using a key which checks the name,
    type and opercase.
    """

    def __init__(self, name, year="1967"):
        super(B311, self).__init__(name)
        self.year = year
        self.k = 1.15    # usage factor
        self.f = 0.90    # fatigue reduction factor

    def h(self, element):
        """Flexibility characterisitic for fittings per the code"""
        if isinstance(element, Bend):
            # Per Appendix D of 1967 code
            R1 = element.radius                 # bend radius
            T = element.section.thk             # nominal thk
            r2 = (element.section.od - T) / 2   # mean radius
            h = (T*R1) / r2**2  # flexibility characteristic

            # stiffening effect due to flanged ends, note 3
            if element.flange == 0:
                c = 1
            elif element.flange == 1:
                c = h**(1/6)
            elif element.flange == 2:
                c = h**(1/3)

            h = c*h  # corrected h

            return h

    def sifi(self, element):
        """In plane stress intensification factor for fittings"""
        if isinstance(element, Run):
            return 1.0

        elif isinstance(element, Bend):
            h = self.h(element)
            sif = 0.9 / h**(2/3)

            return 1.0 if sif < 1 else sif

        elif isinstance(element, Reducer):
            return 1.0

        else:
            # must be 1 at a minimum
            return 1.0

    sifo = sifi     # out of plane - same value

    def kfac(self, element):
        """Code flexibility factor provided for fittings"""
        if isinstance(element, Run):
            return 1.0

        elif isinstance(element, Bend):
            h = self.h(element)
            k = 1.65 / h

            return 1.0 if k < 1 else k

        elif isinstance(element, Reducer):
            return 1.0

        else:
            return 1.0

    def shoop(self, element, loadcase):
        """Hoop stress due to pressure.

        Equal to two times the longitudinal stress due to pressure.
        """
        return 2 * self.slp(element, loadcase)

    def slp(self, element, loadcase):
        """Longitudinal stress due to pressure.

        Exact formulation for the pressure stress per code section 102.3.2.
        The pressure stress is a primary stress and by definition can result
        in large gross deformations and overall failure, i.e., pipe rupture.
        """
        with units.Units(user_units="code_english"):
            section = element.section
            do = section.od
            di = section.id_

            pload = self.poper(element, loadcase)

            try:
                # exact formulation
                return pload.pres*di**2 / (do**2-di**2)
            except AttributeError:
                return 0

    def slb(self, element, forces):
        """Longitudinal stress due to bending moments."""
        with units.Units(user_units="code_english"):
            my, mz = forces[-2:]
            M = math.sqrt(my**2 + mz**2)

            # note for B31.1 sifi equals sifo
            i = self.sifi(element)

            section = element.section
            Z = section.z   # section modulus

            return M*i / Z

    def stor(self, element, forces):
        """Shear stress due to torsion"""
        with units.Units(user_units="code_english"):
            mx = forces[3]

            section = element.section
            do = section.od
            Ip = section.ixx    # polar moment

            return mx*do / (2*Ip)

    def sax(self, element, forces):
        """Axial stress due to mechanical loading.

        Axial stress can be included by user defined option. Force is zero by
        default.
        """
        sx = 0  # code default

        if self.app.models.active_object.settings.axial_force:
            with units.Units(user_units="code_english"):
                fx = forces[0]

                section = element.section
                area = section.area

                sx = fx/area

        return sx

    def sl(self, element, loadcase, forces):
        """Total longitudinal stress due to pressure and bending+torsion. This
        combination of stresses is also known as the code stress for most
        codes.
        """
        with units.Units(user_units="code_english"):
            slp = self.slp(element, loadcase)
            slb = self.slb(element, forces)
            sax = self.sax(element, forces)
            stor = self.stor(element, forces)

            if loadcase.stype == "sus" or loadcase.stype == "occ":
                return sax + slp + math.sqrt(slb**2 + 4*stor**2)
            elif loadcase.stype == "exp":
                return math.sqrt(slb**2 + 4*stor**2)
            else:
                return 0

    def sallow(self, element, loadcase, forces):
        """Allowable stress for sustained, occasional and expansion loadcases.

        Liberal stress can be excluded for the expansion case by user defined
        option otherwise enabled by default depending on the code.
        """
        with units.Units(user_units="code_english"):
            material = element.material
            tload = self.toper(element, loadcase)

            try:
                temp = tload.temp
                tref = tload.tref
            except AttributeError:
                temp = 70   # fahrenheit
                tref = 70

            sh = material.sh[temp]
            sc = material.sh[tref]

            if loadcase.stype == "sus":
                return sh
            elif loadcase.stype == "occ":
                # k factor is 1.15 for events lasting 8hrs or less and 1.20
                # for 1hr or less per code para. 102.3.3, use 1.15 to be
                # conservative
                return self.k * sh
            elif loadcase.stype == "exp":
                liberal_stress = 0  # default per code
                if self.app.models.active_object.settings.liberal_stress:
                    liberal_stress = sh - self.sl(element, loadcase, forces)

                return self.f * (1.25*sc + 0.25*sh + liberal_stress)
            else:
                return 0


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
