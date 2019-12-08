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

"""Main Application Class"""

import sys
import code
import psi
from psi.settings import options
from psi.entity import Entity, EntityContainer
from psi.model import ModelContainer
from psi.point import PointManager
from psi.elements import ElementContainer
from psi.sections import SectionContainer
from psi.material import MaterialContainer
from psi.sifs import SIFContainer
from psi.supports import SupportContainer
from psi.codes import CodeContainer
from psi.insulation import InsulationContainer
from psi.loads import LoadContainer
from psi.loadcase import LoadCaseContainer
from psi.results import ResultContainer
from psi.units import Units


class PSIInterpreter(code.InteractiveConsole):

    def __init__(self, locals=None, filename="<console>"):
        super(PSIInterpreter, self).__init__(locals, filename)

    def runcode(self, code):
        """Raise all expections"""
        try:
            exec(code, self.locals)
        except SystemExit:
            raise
        except:
            self.showtraceback()


class App(object):

    _app = None

    def __new__(cls, *args, **kwargs):
        if cls._app is not None:  # singleton
            return cls._app
        else:
            inst = super(App, cls).__new__(cls, *args, **kwargs)
            cls._app = inst
            return inst

    def __init__(self):
        """Initialize all managers and subsystems"""
        self.options = options

        Units._app = self
        self.units = Units()

        # pass app to objects/containers
        Entity._app = self
        EntityContainer._app = self

        self.models = ModelContainer()
        self.points = PointManager()
        self.elements = ElementContainer()
        self.sections = SectionContainer()
        self.materials = MaterialContainer()
        self.insulation = InsulationContainer()
        self.sifs = SIFContainer()
        self.supports = SupportContainer()
        self.loads = LoadContainer()
        self.loadcases = LoadCaseContainer()
        self.codes = CodeContainer()
        self.results = ResultContainer()

        self.interp = PSIInterpreter()
        self.interp.locals = self._interp_locals

    @property
    def _interp_locals(self):
        return {
                "app": self,
                "options": self.options,

                # models
                "models": self.models,
                "Model": self.models.Model,

                # branch
                "points": self.points,
                "Point": self.points.Point,

                # elements
                "elements": self.elements,
                "Run": self.elements.Run,
                "Bend": self.elements.Bend,
                "Valve": self.elements.Valve,
                "Flange": self.elements.Flange,
                "Rigid": self.elements.Rigid,
                "Reducer": self.elements.Reducer,

                # sections
                "sections": self.sections,
                "Pipe": self.sections.Pipe,
                "WideFlange": self.sections.WideFlange,

                # materials
                "materials": self.materials,
                "Material": self.materials.Material,

                # sifs
                "sifs": self.sifs,
                # tee intersections
                "Welding": self.sifs.Welding,
                "Unreinforced": self.sifs.Unreinforced,
                "Reinforced": self.sifs.Reinforced,
                "Weldolet": self.sifs.Weldolet,
                "Sockolet": self.sifs.Sockolet,
                "Sweepolet": self.sifs.Sweepolet,
                # connections
                "ButtWeld": self.sifs.ButtWeld,

                # supports
                "supports": self.supports,
                "Anchor": self.supports.Anchor,
                "GlobalX": self.supports.GlobalX,
                "GlobalY": self.supports.GlobalY,
                "GlobalZ": self.supports.GlobalZ,
                "Spring": self.supports.Spring,
                "Displacement": self.supports.Displacement,

                # codes
                "codes": self.codes,
                "B311": self.codes.B311,

                # insulation
                "insulation": self.insulation,
                "Insulation": self.insulation.Insulation,

                # loads
                "loads": self.loads,
                "Weight": self.loads.Weight,
                "Pressure": self.loads.Pressure,
                "Hydro": self.loads.Hydro,
                "Thermal": self.loads.Thermal,
                "Fluid": self.loads.Fluid,
                "Uniform": self.loads.Uniform,
                "Wind": self.loads.Wind,
                "Seismic": self.loads.Seismic,
                "Force": self.loads.Force,

                # loadcases
                "loadcases": self.loadcases,
                "LoadCase": self.loadcases.LoadCase,
                "LoadComb": self.loadcases.LoadComb,

                # results
                "results": self.results,
                "Movements": self.results.Movements,
                "Reactions": self.results.Reactions,
                "Forces": self.results.Forces,
                "Stresses": self.results.Stresses,
                }

    @property
    def banner(self):
        msg = ('Python %s\n'
               'Type "copyright", "credits" or "license" for more '
               'information.\n'
               'PSI %s -- The pipe stress analysis and design program.\n' %
               (sys.version.split('\n')[0], psi.VERSION)
               )

        return msg

    def run(self):
        """Run the PSI interpreter"""
        self.interp.interact(self.banner)


if __name__ == "__main__":
    App().run()
