"""Testing various loading types"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes import B311
from psi.supports import Anchor
from psi.loads import Weight, Thermal, Pressure, Seismic
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

    anc10 = Anchor('anc10', 10)
    anc10.apply([run20])
    anc50 = Anchor('anc50', 50)
    anc50.apply([run50])

    app.elements.select()   # select all

    return app


def test_deadweight(app):
    """Check deadweight displacements due to gravity"""
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


def test_pressure(app):
    """Check longitudinal stress due to pressure loading"""
    P1 = Pressure('P1', 1, 250)
    P1.apply()

    # loadcase
    L1 = LoadCase('L1', 'sus', [Pressure], [1])

    app.models('loads').analyze()

    # check displacements due to fz
    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    # assert compare(L1.stresses[pt30][dy], 0.33)


def test_thermal(app):
    """Check displacements due to thermal loading"""
    # loads
    T1 = Thermal('T1', 1, 900, 70)
    T1.apply()

    L1 = LoadCase('L1', 'exp', [Thermal], [1])

    app.models('loads').analyze()

    # check displacements due to fz
    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], 0.758)
    assert compare(L1.movements[pt30][dz], 1.891)
    assert compare(L1.movements[pt30][rx], 0.1521)


def test_seismic_x(app):
    """Check displacements due to seismic loading in x direction"""
    # loads
    S1 = Seismic('S1', 1, gx=1.0)
    S1.apply()

    L1 = LoadCase('L1', 'occ', [Seismic], [1])

    app.models('loads').analyze()

    # check displacements due to fz
    pt30 = app.points(30)
    dx, ry, rz = 0, 4, 5

    assert compare(L1.movements[pt30][dx], 1.436)
    assert compare(L1.movements[pt30][ry], -0.1830)
    assert compare(L1.movements[pt30][rz], -0.4091)


def test_seismic_y(app):
    """Check displacements due to seismic loading in y direction"""
    # loads
    S1 = Seismic('S1', 1, gy=1.0)
    S1.apply()

    L1 = LoadCase('L1', 'occ', [Seismic], [1])

    app.models('loads').analyze()

    # check displacements due to fz
    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], 0.075)
    assert compare(L1.movements[pt30][dz], 0.251)
    assert compare(L1.movements[pt30][rx], 0.1554)


def test_seismic_z(app):
    """Check displacements due to seismic loading in z direction"""
    # loads
    S1 = Seismic('S1', 1, gz=1.0)
    S1.apply()

    L1 = LoadCase('L1', 'occ', [Seismic], [1])

    app.models('loads').analyze()

    # check displacements due to fz
    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], -0.330)
    assert compare(L1.movements[pt30][dz], 1.305)
    assert compare(L1.movements[pt30][rx], 0.3966)
