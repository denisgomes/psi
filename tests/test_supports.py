"""Test all support types"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes import B311
from psi.supports import Anchor
from psi.loads import Force
from psi.loadcase import LoadCase


@pytest.fixture()
def app():
    # parameter
    L = 10*12

    app = App()
    Model('simple')

    # properties
    Pipe.from_file('PIPE1', '10', '40')
    Material.from_file('MAT1', 'A53A', 'B31.1')

    # geometry
    Point(10)
    run20 = Run(20, L)

    # code
    b311 = B311('B311')
    b311.apply([run20])

    return app


def test_anchor(app):
    # get pipe objects
    pt10 = app.points(10)
    run20 = app.elements(10, 20)

    # support
    anc1 = Anchor('A1', 10)
    anc1.apply([run20])

    # loads
    F1 = Force('F1', 1, 20, fy=-10000)
    F1.apply([run20])

    F2 = Force('F2', 2, 20, fx=-10000)
    F2.apply([run20])

    F3 = Force('F3', 3, 20, fz=-10000)
    F3.apply([run20])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])
    L2 = LoadCase('L2', 'ope', [Force], [2])
    L3 = LoadCase('L3', 'ope', [Force], [3])

    app.models('simple').analyze()

    # check reactions due to fy
    assert round(L1.reactions[pt10][1]) == 10000
    assert round(L1.reactions[pt10][5]) == 100000

    # check reaction due to fx
    assert round(L2.reactions[pt10][0]) == 10000

    # check reaction due to fz
    assert round(L3.reactions[pt10][2]) == 10000
    assert round(L3.reactions[pt10][4]) == -100000
