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

"""Main Application Class"""

import sys
import code
import psi
from psi.settings import Configuration
from psi.entity import Entity, EntityContainer
from psi.model import ModelContainer
from psi.point import PointManager
from psi.elements import ElementContainer
from psi.sections import SectionContainer
from psi.material import MaterialContainer
from psi.sifs import SIFContainer
from psi.supports import SupportContainer
from psi.codes.codes import CodeContainer
from psi.insulation import InsulationContainer
from psi.loads import LoadContainer
from psi.loadcase import LoadCaseContainer
from psi.reports import ReportContainer
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
        Units._app = self
        self.units = Units()

        # pass app to class objects/containers
        Configuration._app = self
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
        self.reports = ReportContainer()

        self.interp = PSIInterpreter()
        self.interp.locals = self._interp_locals

    @property
    def _interp_locals(self):
        return {
                "app": self,

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
                "Bellow": self.elements.Bellow,

                # sections
                "sections": self.sections,
                "Pipe": self.sections.Pipe,
                "Beam": self.sections.Beam,

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
                "Inclined": self.supports.Inclined,
                "X": self.supports.X,
                "Y": self.supports.Y,
                "Z": self.supports.Z,
                "Spring": self.supports.Spring,
                "LimitStop": self.supports.LimitStop,
                "Lateral": self.supports.Lateral,

                # codes
                "codes": self.codes,
                "B31167": self.codes.B31167,

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
                "Wind": self.loads.Wind,
                "Seismic": self.loads.Seismic,
                "Force": self.loads.Force,
                "Displacement": self.loads.Displacement,
                "Hanger": self.loads.Hanger,

                # loadcases
                "loadcases": self.loadcases,
                "LoadCase": self.loadcases.LoadCase,
                "LoadComb": self.loadcases.LoadComb,

                # reports
                "reports": self.reports,
                "Movements": self.reports.Movements,
                "Reactions": self.reports.Reactions,
                "Forces": self.reports.Forces,
                "MForces": self.reports.MForces,
                "Stresses": self.reports.Stresses,
                "StressesExtended": self.reports.StressesExtended,
                "Codes": self.reports.Codes,
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
