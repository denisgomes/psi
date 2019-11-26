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
from psi.units import units


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
        self.units = units

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
        # self.results = ResultContainer()

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
                "Displacement": self.loads.Displacement,

                # loadcases
                "loadcases": self.loadcases,
                "LoadCase": self.loadcases.LoadCase,
                "LoadComb": self.loadcases.LoadComb,

                # # results
                # "results": self.results,
                # "Stresses": self.results.Stresses,
                # "Movements": self.results.Movements,
                }

    @property
    def banner(self):
        msg = ('Python %s\n'
               'Type "copyright", "credits" or "license" for more '
               'information.\n'
               'PSI %s -- A pipe stress analysis and design program.\n' %
               (sys.version.split('\n')[0], psi.VERSION)
               )

        return msg

    def run(self):
        """Run the PSI interpreter"""
        self.interp.interact(self.banner)


if __name__ == "__main__":
    App().run()
