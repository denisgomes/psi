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

"""Implementation of B31.1 Power Piping codes"""

import math
import sys
from contextlib import redirect_stdout

from psi.elements import (Run, Bend, Reducer, Rigid, Valve, Flange)
from psi.sifs import (Welding, Unreinforced, Reinforced, Weldolet, Sockolet,
                      Sweepolet, Weldolet, ButtWeld)
from psi.loadcase import LoadCase
from psi import units
from .codes import Code


class B31167(Code):
    """B31.1 1967 Power Piping code implementation.

    For each element, based on the loadcase and stress type, calculate code
    stresses.

    Each element can have one Temperature and/or Pressure load defined per one
    opercase. A load is defined to be unique using a key which checks the name,
    type and opercase.

    Stress intensification and flexibility factors are based on Mandatory
    Appendix D. Per code, the calculated SIFs must always be greater than or
    eqaul to 1.0.

    SIFs provided are valid for a D/t ratio of less than 100. Beyond this
    limit, the pipe behaves like large duct piping and must be modeling using
    shell elements.

    Multiple SIFs defined for a particular element for example a reducing tee
    is not properly defined by any code. Should the reducer, tee or both SIFs
    be used to determine the final stress results?
    """

    def __init__(self, name, year="1967", cycles=7000, occ_hours=8):
        super(B31167, self).__init__(name)
        self.year = year
        self.k = self.k_from_hours(occ_hours)
        self.f = self.f_from_cycles(cycles)
        self.Y = 0.0    # hoop stress factor

    @property
    def label(self):
        """Title used in report output"""
        return "B31167"

    def f_from_cycles(self, cycles):
        """Thermal stress range reduction factor due to fatigue based on the
        number of cycles per Table 102.3.2(c).
        """
        if cycles <= 7000:
            return 1
        elif 7000 < cycles <= 14000:
            return 0.9
        elif 14000 < cycles <= 22000:
            return 0.8
        elif 22000 < cycles <= 45000:
            return 0.7
        elif 45000 < cycles <= 100000:
            return 0.6
        else:
            return 0.5

    def k_from_hours(self, hours):
        """The occasional stress usage factor based on the total duration of
        the event given in hours.
        """
        if hours <= 1:
            return 1.2
        elif 1 < hours <= 8:
            return 1.15
        else:
            return 1.15

    def h(self, entity):
        """Flexibility characterisitic for fittings per the code."""
        with units.Units(user_units="code_english"):
            if isinstance(entity, Bend):
                # Per Appendix D of 1967 code
                R = entity.radius                  # bend radius
                tn = entity.section.thk            # nominal thk
                r = (entity.section.od - T) / 2    # mean radius
                h = (tn*R) / r**2                  # flexibility characteristic

                # stiffening effect due to flanged ends, note 3
                if entity.flange == 0:
                    c = 1
                elif entity.flange == 1:
                    c = h**(1/6)
                elif entity.flange == 2:
                    c = h**(1/3)

                h = c*h  # corrected h

                return h

            elif isinstance(entity, Reducer):
                return 1.0

            elif isinstance(entity, Welding):
                do = entity.do      # header pipe dia
                dob = entity.dob    # branch pipe dia
                tn = entity.tn      # nominal header thk
                rx = entity.rx      # crotch radius
                tc = entity.tc      # crotch thk
                r = (do-tn) / 2

                if rx >= dob/8 and tc >= 1.5*tn:
                    h = 4.4*tn / r
                else:
                    h = 3.1*tn / r

                return h

            elif isinstance(entity, Unreinforced):
                do = entity.do      # header pipe dia
                tn = entity.tn      # nominal header thk
                r = (do-tn) / 2
                h = tn / r

                return h

            elif isinstance(entity, Reinforced):
                do = entity.do      # header pipe dia
                tn = entity.tn      # nominal header thk
                tr = entity.tr      # pad thk
                r = (do-tn) / 2

                if tr > 1.5*tn:
                    h = 4.05*tn / r
                else:
                    h = (tn+tr/2)**(5/2) / (r*tn**(3/2))

                return h

            elif isinstance(entity, Weldolet):
                do = entity.do      # header pipe dia
                tn = entity.tn      # nominal header thk
                r = (do-tn) / 2
                h = 3.3*tn / r

                return h

            elif isinstance(entity, Sockolet):
                do = entity.do      # header pipe dia
                tn = entity.tn      # nominal header thk
                r = (do-tn) / 2
                h = 3.3*tn / r

                return h

            elif isinstance(entity, Sweepolet):
                do = entity.do      # header pipe dia
                dob = entity.dob    # branch pipe dia
                tn = entity.tn      # nominal header thk
                rx = entity.rx      # crotch radius
                tc = entity.tc      # crotch thk
                r = (do-tn) / 2

                if rx >= dob/8 and tc >= 1.5*tn:
                    h = 4.4*tn / r
                else:
                    h = 3.1*tn / r

                return h

            else:
                return 1.0

    def sifi(self, element, point, combfunc="max"):
        """In plane stress intensification factor for fittings. The sif must
        be 1 or greater.

        .. todo::
            The maximum sif at a point is taken. It should eventually be a user
            defined option to take the sum, max or average.
        """
        model = self.app.models.active_object
        section = element.section
        material = element.material
        sifs = []

        with units.Units(user_units="code_english"):
            if isinstance(element, (Run, Rigid, Valve, Flange)):
                sifs.append(1.0)

            elif isinstance(element, Bend):
                R = section.radius          # bend radius
                tn = section.thk            # nominal thickness
                r = (section.od - tn) / 2   # mean radius
                Ec = material.ymod[model.settings.tref]
                pmax = self.pmax(element)
                h = self.h(element)

                ic = 1  # corrected for pressure - note 5
                if section.is_large_bore and section.is_thin_wall:
                    ic = 1 + 3.25 * (pmax/Ec) * (r/tn)**(5/2) * (R/r)**(2/3)

                sif = (0.9 / h**(2/3)) / ic

                sifs.append(1.0 if sif < 1 else sif)

            elif isinstance(element, Reducer):
                D1 = section.od
                t1 = section.thk
                D2 = section2.od
                t2 = section2.thk
                L = element.length
                alp = element.alpha

                alpv = True if alp <= 60 else False
                conc = True if element.is_concentric else False
                d2t1 = True if section.d2t <= 100 else False
                d2t2 = True if section2.d2t <= 100 else False

                if all([alpv, conc, d2t1, d2t2]):
                    sif = 0.5 + 0.01*alp*(D2/t2)**(1/2)

                    sifs.append(1.0 if sif < 1 else sif)

                else:
                    sifs.append(2.0)

            for entity in element.sifs:

                if self.app.points(entity.point) is point:

                    if isinstance(entity, Welding):
                        # per mandatory appendix D
                        h = self.h(entity)

                        # in-plane and out-of-plane sifs are the same
                        # for B31.1 the higher is used
                        sif = 0.9 / h**(2/3)

                        # must be larger than or equal to 1.0
                        sifs.append(1.0 if sif < 1 else sif)

                    elif isinstance(entity, Unreinforced):
                        h = self.h(entity)

                        # in-plane and out-of-plane sifs are the same
                        # for B31.1 the higher is used
                        sif = 0.9 / h**(2/3)

                        sifs.append(1.0 if sif < 1 else sif)

                    elif isinstance(entity, Reinforced):
                        h = self.h(entity)

                        # in-plane and out-of-plane sifs are the same
                        # for B31.1 the higher is used
                        sif = 0.9 / h**(2/3)

                        sifs.append(1.0 if sif < 1 else sif)

                    elif isinstance(entity, Weldolet):
                        h = self.h(entity)

                        # in-plane and out-of-plane sifs are the same
                        # for B31.1 the higher is used
                        sif = 0.9 / h**(2/3)

                        sifs.append(1.0 if sif < 1 else sif)

                    elif isinstance(entity, Sockolet):
                        h = self.h(entity)

                        # in-plane and out-of-plane sifs are the same
                        # for B31.1 the higher is used
                        sif = 0.9 / h**(2/3)

                        sifs.append(1.0 if sif < 1 else sif)

                    elif isinstance(entity, Sweepolet):
                        h = self.h(entity)

                        # in-plane and out-of-plane sifs are the same
                        # for B31.1 the higher is used
                        sif = 0.9 / h**(2/3)

                        sifs.append(1.0 if sif < 1 else sif)

                    elif isinstance(entity, ButtWeld):
                        sifs.append(1.0)

                    else:
                        # must be 1 at a minimum
                        sifs.append(1.0)

            if combfunc == "max":
                return max(sifs)
            elif combfunc == "sum":
                return sum(sifs)
            elif combfunc == "avg":
                return sum(sifs) / len(sifs)

    sifo = sifi     # out of plane - same value

    def kfac(self, element):
        """Code flexibility factor for fittings."""
        model = self.app.models.active_object
        section = element.section
        material = element.material

        with units.Units(user_units="code_english"):
            if isinstance(element, Bend):
                R = section.radius          # bend radius
                tn = section.thk            # nominal thickness
                r = (section.od - tn) / 2   # mean radius
                Ec = material.ymod[model.settings.tref]
                pmax = self.pmax(element)
                h = self.h(element)

                kc = 1  # corrected for pressure - note 5
                if section.is_large_bore and section.is_thin_wall:
                    kc = 1 + 6 * (pmax/Ec) * (r/tn)**(7/3) * (R/r)**(1/3)

                k = (1.65 / h) / kc

                return 1.0 if k < 1 else k

            elif isinstance(element, Reducer):
                return 1.0

            else:
                return 1.0

    def shoop(self, element, loadcase):
        """Hoop stress due to pressure.

        Conservatively equal to P*D/2t.
        """
        with units.Units(user_units="code_english"):
            sec = element.section
            do = sec.od
            thke = sec.thke
            pload = self.poper(element, loadcase)

            try:
                return (pload.pres*do) / (2*thke)
            except AttributeError:
                return 0

    def slp(self, element, loadcase):
        """Longitudinal stress due to pressure.

        Exact formulation for the pressure stress per code section 102.3.2.
        The pressure stress is a primary stress and by definition can result
        in large gross deformations and overall failure, i.e., pipe rupture.
        """
        with units.Units(user_units="code_english"):
            section = element.section
            do = section.od
            di = section.id

            pload = self.poper(element, loadcase)

            try:
                # exact formulation
                return pload.pres*di**2 / (do**2-di**2)
            except AttributeError:
                return 0

    def slb(self, element, point, forces):
        """Longitudinal stress due to bending moments."""
        with units.Units(user_units="code_english"):
            my, mz = forces[-2:]
            M = math.sqrt(my**2 + mz**2)

            # note for B31.1 sifi equals sifo
            i = self.sifi(element, point)

            section = element.section
            Z = section.z   # section modulus

            return M*i / Z

    def sts(self, element, forces):
        """Transverse shear stress"""
        with units.Units(user_units="code_english"):
            fy, fz = forces[1:3]
            F = math.sqrt(fy**2 + fz**2)

            section = element.section
            area = section.area

            return F / area

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

    def sl(self, element, loadcase, point, forces):
        """Total longitudinal stress due to pressure and bending+torsion. This
        combination of stresses is also known as the code stress for most
        codes.
        """
        with units.Units(user_units="code_english"):
            slp = self.slp(element, loadcase)
            slb = self.slb(element, point, forces)
            sax = self.sax(element, forces)
            stor = self.stor(element, forces)

            if loadcase.stype == "sus" or loadcase.stype == "occ":
                return sax + slp + math.sqrt(slb**2 + 4*stor**2)
            elif loadcase.stype == "exp":
                return math.sqrt(slb**2 + 4*stor**2)
            else:
                return 0

    def sallow(self, element, loadcase, point, forces):
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
                    cases = self.app.models.active_object.loadcases
                    suscases = [case for case in cases if
                            isinstance(case, LoadCase) and case.stype == "sus"]
                    shsls = [sh-suscase.stresses[point].sl for suscase in suscases]
                    liberal_stresses = [shsl for shsl in shsls if shsl > 0]

                    if liberal_stresses:
                        liberal_stress = min(liberal_stresses)

                return self.f*(1.25*sc + 0.25*sh + liberal_stress)
            else:
                return 0
