"""Example 1: Skewed Cantilever"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes import B311
from psi.supports import Anchor
from psi.loads import Weight
from psi.loadcase import LoadCase
from .utils import compare


@pytest.fixture()
def app():
    app = App()

    mdl = Model('loads')
    mdl.settings.vertical = "z"

    Pipe.from_file('pipe1', '10', '40')
    Material.from_file('mat1', 'A53A', 'B31.1')
    B311('code1')

    Point(10)
    run20 = Run(20, 0, 0, 10*12)
    Run(30, 0, 15*12)
    Run(40, 0, 0, -50*12)
    run50 = Run(50, 0, 25*12)

    anc1 = Anchor('anc1', 10)
    anc1.apply([run20])
    anc2 = Anchor('anc2', 50)
    anc2.apply([run50])

    app.elements.select()   # select all

    return app


def test_deadweight(app):
    """Deadweight load vector test"""
    # loads
    W1 = Weight('W1', 1)
    W1.apply()  # to selected elements

    # loadcase
    L1 = LoadCase('L1', 'sus', [Weight], [1])

    app.models('loads').analyze()

    # check displacements due to fz
    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], 0.33)
    assert compare(L1.movements[pt30][dz], -1.31)
    assert compare(L1.movements[pt30][rx], -0.4)
