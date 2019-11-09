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

        Entity._app = self  # pass app to objects
        EntityContainer._app = self   # pass app to containers

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

        self.interp = code.InteractiveConsole()
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
                "Vertical": self.supports.Vertical,
                "Lateral": self.supports.Lateral,
                "Axial": self.supports.Axial,
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
        msg = ('Pipe Stress Infinity %s - The pipe stress design and analysis program.\n'
               'Python %s\n'
               'Type "copyright", "credits" or "license" '
               'for more information.' %
               (psi.VERSION, sys.version.split('\n')[0]))

        return msg

    def run(self):
        """Run the PSI interpreter"""
        self.interp.interact(self.banner)


if __name__ == "__main__":
    App().run()
