"""Test support types"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes.b311 import B31167
from psi.supports import Anchor, X, Y, Z, LimitStop, Inclined, Lateral, Spring
from psi.loads import Force, Weight, Thermal
from psi.loadcase import LoadCase

from .utils import compare


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
    b311 = B31167('B311')
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

    # check reactions due to fy, mz
    assert compare(L1.reactions[pt10][1], -10000)
    assert compare(L1.reactions[pt10][5], -100000)

    # check reaction due to fx
    assert compare(L2.reactions[pt10][0], -10000)

    # check reaction due to fz, fy
    assert compare(L3.reactions[pt10][2], -10000)
    assert compare(L3.reactions[pt10][4], 100000)


def test_global_y(app):
    """Simply supported beam with central concentrated force"""

    # get pipe objects
    pt10 = app.points(10)
    pt15 = app.points(15)
    pt20 = app.points(20)

    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    gblx10 = X('GblX10', 10)
    gblx10.apply([run10])

    gbly10 = Y('GblY10', 10)
    gbly10.apply([run10])

    gblz10 = Z('GblZ10', 10)
    gblz10.apply([run10])

    # constrain torsion
    gblrotx10 = X('GblRotX10', 10, is_rotational=True)
    gblrotx10.apply([run10])

    gbly20 = Y('GblY20', 20)
    gbly20.apply([run20])

    gblz20 = Z('GblZ20', 20)
    gblz20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=-10000)
    F1.apply([run10])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert compare(L1.reactions[pt10][1], -5000)
    assert compare(L1.reactions[pt20][1], -5000)

    # check max moment at the middle, mz
    assert compare(L1.forces[pt15][5], -25000)


def test_incline(app):
    """Skewed support test"""

    # get pipe objects
    pt10 = app.points(10)
    pt20 = app.points(20)

    run15 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    anc10 = Anchor('anc10', 10)
    anc10.apply([run15])

    y20 = Inclined('y20', 20, (-0.7071, 0.7071, 0))
    y20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=-10000)
    F1.apply([run15])

    W1 = Weight('W1', 1)
    W1.apply([run15, run20])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Weight, Force], [1, 1])

    app.models('simple').analyze()

    # check reactions due to fy, mx
    assert compare(L1.reactions[pt10].fy, -7108.9)
    assert compare(L1.reactions[pt10].mz, -19067.105)

    # check reactions due to fx, fy
    assert compare(L1.reactions[pt20].fx, 3295.48)
    assert compare(L1.reactions[pt20].fy, -3295.48)


def test_limitstop(app):
    """Limit stop test"""

    # get pipe objects
    pt20 = app.points(20)

    run15 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    anc10 = Anchor('anc10', 10)
    anc10.apply([run15])

    lim20 = LimitStop('lim20', 20)
    lim20.apply([run20])

    # loads
    T1 = Thermal('T1', 1, 500, 70)
    T1.apply([run15, run20])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Thermal], [1])

    app.models('simple').analyze()

    # check reactions due to fx
    assert compare(L1.reactions[pt20].fx, 1002348.1)


def test_guide(app):
    # get pipe objects
    pt10 = app.points(10)
    pt20 = app.points(20)
    run15 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    gblx10 = X('GblX10', 10)
    gblx10.apply([run15])

    gbly10 = Y('GblY10', 10)
    gbly10.apply([run15])

    gblz10 = Z('GblZ10', 10)
    gblz10.apply([run15])

    gblrotx10 = X('GblRotX10', 10, is_rotational=True)
    gblrotx10.apply([run15])

    guide20 = Lateral('guide20', 20)
    guide20.apply([run20])
    guide20.flip()  # not required for horizontal elements

    # loads
    T1 = Thermal('T1', 1, 500, 70)
    T1.apply([run15, run20])

    F1 = Force('F1', 1, 15, fz=-10000)
    F1.apply([run15])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Thermal, Force], [1, 1])

    app.models('simple').analyze()

    assert compare(L1.reactions[pt20].fx, 0)
    assert compare(L1.reactions[pt20].fy, 0)
    assert compare(L1.reactions[pt20].fz, -5000)


def test_plus_y_negative_load(app):
    """Simply supported beam with central concentrated force and a non-linear
    +Y support.
    """

    # get pipe objects
    pt10 = app.points(10)
    pt15 = app.points(15)
    pt20 = app.points(20)

    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    gblx10 = X('GblX10', 10)
    gblx10.apply([run10])

    gbly10 = Y('GblY10', 10)
    gbly10.apply([run10])

    gblz10 = Z('GblZ10', 10)
    gblz10.apply([run10])

    # constrain torsion
    gblrotx10 = X('GblRotX10', 10, is_rotational=True)
    gblrotx10.apply([run10])

    # nonlinear +Y
    gbly20 = Y('GblY20', 20, direction="+")
    gbly20.apply([run20])

    gblz20 = Z('GblZ20', 20)
    gblz20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=-10000)
    F1.apply([run10])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert compare(L1.reactions[pt10][1], -5000)
    assert compare(L1.reactions[pt20][1], -5000)
    assert compare(L1.movements[pt20][1], 0)


def test_plus_y_positive_load(app):
    """Simply supported beam with central concentrated force and a non-linear
    +Y support.
    """

    # get pipe objects
    pt10 = app.points(10)
    pt15 = app.points(15)
    pt20 = app.points(20)

    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    anc10 = Anchor('Anc10', 10)
    anc10.apply([run10])

    # nonlinear +Y
    gbly20 = Y('GblY20', 20, direction="+")
    gbly20.apply([run20])

    gblz20 = Z('GblZ20', 20)
    gblz20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=10000)
    F1.apply([run10])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert compare(L1.reactions[pt10][1], 10000)
    assert compare(L1.reactions[pt20][1], 0)


def test_minus_y_positive_load(app):
    """Simply supported beam with central concentrated force and a non-linear
    -Y support.
    """

    # get pipe objects
    pt10 = app.points(10)
    pt15 = app.points(15)
    pt20 = app.points(20)

    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    gblx10 = X('GblX10', 10)
    gblx10.apply([run10])

    gbly10 = Y('GblY10', 10)
    gbly10.apply([run10])

    gblz10 = Z('GblZ10', 10)
    gblz10.apply([run10])

    # constrain torsion
    gblrotx10 = X('GblRotX10', 10, is_rotational=True)
    gblrotx10.apply([run10])

    # nonlinear -Y
    gbly20 = Y('GblY20', 20, direction="-")
    gbly20.apply([run20])

    gblz20 = Z('GblZ20', 20)
    gblz20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=10000)
    F1.apply([run10])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert compare(L1.reactions[pt10][1], 5000)
    assert compare(L1.reactions[pt20][1], 5000)
    assert compare(L1.movements[pt20][1], 0)


def test_minus_y_negative_load(app):
    """Simply supported beam with central concentrated force and a non-linear
    -Y support.
    """

    # get pipe objects
    pt10 = app.points(10)
    pt15 = app.points(15)
    pt20 = app.points(20)

    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    anc10 = Anchor('Anc10', 10)
    anc10.apply([run10])

    # nonlinear -Y
    gbly20 = Y('GblY20', 20, direction="-")
    gbly20.apply([run20])

    gblz20 = Z('GblZ20', 20)
    gblz20.apply([run20])

    # loads
    F1 = Force('F1', 1, 15, fy=-10000)
    F1.apply([run10])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert compare(L1.reactions[pt10][1], -10000)
    assert compare(L1.reactions[pt20][1], 0)


def test_variable_spring(app):

    # get pipe objects
    pt10 = app.points(10)
    run10 = app.elements(10, 15)
    run20 = app.elements(15, 20)

    # supports
    anc1 = Anchor('A1', 10)
    anc1.apply([run10])

    spr20 = Spring('Spr20', 20, 1000, 0)
    spr20.apply([run20])

    W1 = Weight('W1', 1)
    W1.apply([run10, run20])

    # loadcase
    L1 = LoadCase('L1', 'sus', [Weight], [1])

    app.models('simple').analyze()

    # check reactions due to fy
    assert compare(L1.reactions[pt10].fy, -100)
