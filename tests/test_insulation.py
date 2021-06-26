"""Pipe insulation test"""


from psi.app import App
from psi.codes import B311
from psi.elements import Run
from psi.insulation import Insulation
from psi.loads import Weight
from psi.loadcase import LoadCase
from psi.material import Material
from psi.model import Model
from psi.point import Point
from psi.sections import Pipe
from psi.supports import Anchor
from .utils import compare


def test_insulation():
    """Test combined weight of pipe and insulation with no contents"""
    # parameter
    L = 10*12

    app = App()
    mdl = Model('simple')

    # properties
    Pipe.from_file('pipe1', '10', '40')
    Insulation.from_file('insul1', 'minwool', 3)
    Material.from_file('mat1', 'A53A', 'B31.1')
    B311('B31.1')

    pt10 = Point(10)
    run20 = Run(20, L)

    anc10 = Anchor('anc10', 10)
    anc10.apply([run20])

    app.elements.select()
    dw = Weight('dw', 1)
    dw.apply()

    lc1 = LoadCase('lc1', 'sus', [Weight], [1])

    mdl.analyze()

    assert compare(lc1.reactions[pt10].fy, -521.3947)
