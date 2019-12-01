"""Pipe insulation"""

import csv

import psi
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi import units


@units.define(rho="density", thk="length")
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
