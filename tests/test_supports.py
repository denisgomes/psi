"""Test support types"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes import B311
from psi.supports import Anchor, GlobalX, GlobalY, GlobalZ
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
    run15 = Run(15, L/2)
    run20 = Run(20, L/2)

    # code
    b311 = B311('B311')
    b311.apply([run15, run20])

    return app


def test_anchor(app):
    """Cantilever beam with concentrated force"""

    # get pipe objects
    pt10 = app.points(10)
    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # support
    anc1 = Anchor('A1', 10)
    anc1.apply([run10])

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


def test_global_y(app):
    """Simply supported beam with central concentrated force"""

    # get pipe objects
    pt10 = app.points(10)
    pt15 = app.points(15)
    pt20 = app.points(20)

    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    gblx10 = GlobalX('GblX10', 10)
    gblx10.apply([run10])

    gbly10 = GlobalY('GblY10', 10)
    gbly10.apply([run10])

    gblz10 = GlobalZ('GblZ10', 10)
    gblz10.apply([run10])

    # constrain torsion
    gblrotx10 = GlobalX('GblRotX10', 10, is_rotational=True)
    gblrotx10.apply([run10])

    gbly20 = GlobalY('GblY20', 20)
    gbly20.apply([run20])

    gblz20 = GlobalZ('GblZ20', 20)
    gblz20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=-10000)
    F1.apply([run10])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert round(L1.reactions[pt10][1]) == 5000
    assert round(L1.reactions[pt20][1]) == 5000

    # check max moment at the middle
    assert round(L1.forces[pt15][5]) == -25000
