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

    def sif(self, element):
        """Element stress intensification factor.

        The code stress SIF is a fatigue factor developed by Markl and his team
        in the 1950s. These factors were empiracally determined by subject
        piping to alternating loads and determining the effect it had over a
        number of cycles.

        By defination a SIF is the peak stress over the average stress of a
        piping component. It is a multiplier on nominal stress for typical
        bend and intersection components so that the effect of geometry and
        welding can be considered in beam analysis.

        The inplane and out of plane SIF applies to the moments and associated
        stresses in the respective directions. An inplane moment/force will
        tend to keep the component in the plane of the page. An out of plane
        moment or force will make the component come out of the page.
        """
        raise NotImplementedError("implement")

    def kfac(self, element):
        """Element flexibility factor.

        Note this factor is applied to the element for bending only. It has
        the effect of reducing the stiffness of the element in the transverse
        bending direction.
        """
        raise NotImplementedError("implement")


class B311(Code):

    def __init__(self, name):
        super(B311, self).__init__(name)

    def thke(self):
        """Effective thickness used for stress calculations"""
        pass

    def Shoop(self, section, pressure):
        """Hoop stress due to pressure"""
        pass

    def Slp(self, section, pressure):
        """Longitudinal stress due to pressure.

        Exact formulation for the pressure stress per code section 102.3.2.
        The pressure stress is a primary stress and by definition can result
        in large gross deformations and overall failure, i.e., pipe rupture.
        """
        do = section.od
        di = section.id
        p = pressure.pres   # the pressure load

        return p*di**2 / (do**2-di**2)

    def Slb(self, section):
        pass

    def Sh(self, material, loadcase):
        """Material hot allowable"""

        # extract temperature of loadcase
        exp = loadcase.loads
        temp = exp.temp

        if loadcase.stype == "SUS":
            return material.sh[temp]
        elif loadcase.stype == "OCC":
            # k factor is 1.15 for events lasting 8hrs or less and
            # 1.20 for 1hr or less per code para. 102.3.3, use 1.15 to
            # be conservative
            return 1.15 * material.sh[temp]
        elif loadcase.stype == "EXP":
            pass

    def Sa(self):
        """Thermal stress allowable using left over sustained stress.

        Allowed only if the user selects to use liberal stress allowables.
        """
        pass

    def sif(self, element):
        if isinstance(element, Run):
            return 1.0

    def kfac(self, element):
        if isinstance(element, Run):
            return 1.0


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
