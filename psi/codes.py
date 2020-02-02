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


from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.elements import (Run, Bend, Reducer, Rigid, Valve, Flange)


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

        Note this factor is applied to the element for bending only. It has
        the effect of reducing the stiffness of the element in the transverse
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

    def Sa(self, element, loadcase):
        """Element thermal stress allowable"""
        raise NotImplementedError("implement")

    def Sh(self, element, loadcase):
        """Element hot allowable"""
        raise NotImplementedError("implement")


class B311(Code):
    """B31.1 Power Piping Code Implementation.

    For each element, based on the loadcase and stress type, calculate code
    stresses.

    Each element can have one Temperature and/or Pressure load defined per one
    opercase. A load is defined to be unique using a key which checks the name,
    type and opercase.
    """

    def __init__(self, name, year='1965'):
        super(B311, self).__init__(name)
        self.year = year

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
        return 2 * self.Slp

    def Slp(self, element, loadcase):
        """Longitudinal stress due to pressure.

        Exact formulation for the pressure stress per code section 102.3.2.
        The pressure stress is a primary stress and by definition can result
        in large gross deformations and overall failure, i.e., pipe rupture.
        """
        section = element.section

        do = section.od
        di = section.id
        p = pressure.pres   # the pressure load

        # use the exact formulation
        return p*di**2 / (do**2-di**2)

    def Slb(self, element, loadcase):
        """Longitudinal stress due to bending moments."""
        pass

    def Stor(self, element, loadcase):
        pass

    def Sax(self, element, loadcase):
        pass

    def Sh(self, element, loadcase):
        """Material hot allowable"""

        # extract temperature of loadcase
        exp = loadcase.loads
        temp = exp.temp

        if loadcase.stype == "sus":
            return material.sh[temp]
        elif loadcase.stype == "occ":
            # k factor is 1.15 for events lasting 8hrs or less and
            # 1.20 for 1hr or less per code para. 102.3.3, use 1.15 to
            # be conservative
            return 1.15 * material.sh[temp]
        elif loadcase.stype == "exp":
            pass

    def Sa(self, element, loadcase):
        """Thermal stress allowable using left over sustained stress.

        Allowed only if the user selects to use liberal stress allowables.
        """
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
