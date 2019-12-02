# Copyright (c) 2019 Denis Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
#    * Neither the name of Pipe Stress Infinity (PSI) nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

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
