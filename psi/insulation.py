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

"""Pipe insulation"""

import csv

import psi
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi import units


@units.define(rho="insulation_density", thk="length")
class Insulation(Entity, ActiveEntityMixin):

    def __init__(self, name, rho, thk):
        super(Insulation, self).__init__(name)
        self.rho = rho
        self.thk = thk

        self.activate()

    def activate(self):
        self.parent.activate(self)

    def apply(self, elements=None):
        self.parent.apply(elements)

    @classmethod
    def from_file(cls, name, material, thickness, fname=None):
        if fname is None:
            fname = psi.INSULATION_DATA_FILE

        with open(fname) as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if row["insul"] == material:
                    rho = float(row["rho"])
                    thk = float(thickness)
                    return cls(name, rho, thk)
            else:
                return None

    @property
    def parent(self):
        return self.app.insulation

    def __repr__(self):
        return "%s(rho=%s, thk=%s)" % ("Insulation", self.rho, self.thk)


class InsulationContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(InsulationContainer, self).__init__()
        self.Insulation = Insulation

    def apply(self, inst, elements=None):
        if elements is None:
            elements = []
            elements.append(self.app.models.active_model.active_element)

        for element in elements:
            element.insulation = inst
