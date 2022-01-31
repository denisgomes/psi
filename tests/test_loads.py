"""Testing various loading types"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes.b311 import B31167
from psi.supports import Anchor
from psi.loads import (Weight, Thermal, Pressure, Fluid, Hydro, Seismic, Wind,
                       Force, Displacement)
from psi.loadcase import LoadCase
from .utils import compare


@pytest.fixture()
def app():
    app = App()

    mdl = Model('loads')
    mdl.settings.vertical = "z"

    Pipe.from_file('pipe1', '10', '40')
    Material.from_file('mat1', 'A53A', 'B31.1')
    B31167('code1')

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


@pytest.fixture()
def app2():
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



def test_deadweight(app):
    """Check deadweight displacements due to gravity."""
    # loads
    W1 = Weight('W1', 1)
    W1.apply()  # to selected elements

    # loadcase
    L1 = LoadCase('L1', 'sus', [Weight], [1])

    app.models('loads').analyze()

    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], 0.33)
    assert compare(L1.movements[pt30][dz], -1.31)
    assert compare(L1.movements[pt30][rx], -0.4)


def test_pressure(app):
    """Check longitudinal stress due to pressure loading, pressure thrust and
    bourdon pressure effects.


    .. note::

        The internal hoop stress formula is (P*D)/2*t. A different code
        specific formula may be used by other softwares, however, the
        common formulation used is conservative.
    """
    P1 = Pressure('P1', 1, 250)
    P1.apply()

    # loadcase
    L1 = LoadCase('L1', 'sus', [Pressure], [1])

    app.models('loads').analyze()

    # check stress due to pressure
    pt30 = app.points(30)

    assert compare(L1.stresses[pt30].shoop, 3681.51)
    assert compare(L1.stresses[pt30].slp, 1655.45)


def test_pressure_thrust(app):
    """Check axial displacement due to pressure thrust.


    .. note::

        The pressure thrust force is based on the inner diameter and the pipe
        internal pressure, which may be the different from what some other
        piping softwares use internally.
    """
    app.models('loads').settings.pressure_thrust = True

    P1 = Pressure('P1', 1, 1000)
    P1.apply()

    # loadcase
    L1 = LoadCase('L1', 'sus', [Pressure], [1])

    app.models('loads').analyze()

    # check displacement due to pressure thrust
    pt40 = app.points(40)

    assert compare(L1.movements[pt40].dz, -0.0366)


def test_pressure_bourdon(app):
    """Check axial displacement due to bourdon pressure effects."""
    pass


def test_hydro(app):
    """Check longitudinal stress due to hydrotest pressure."""
    pass


def test_thermal(app):
    """Check displacements due to thermal loading."""
    # loads
    T1 = Thermal('T1', 1, 900, 70)
    T1.apply()

    L1 = LoadCase('L1', 'exp', [Thermal], [1])

    app.models('loads').analyze()

    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], 0.758)
    assert compare(L1.movements[pt30][dz], 1.891)
    assert compare(L1.movements[pt30][rx], 0.1521)


def test_fluid(app):
    """Check deadweight displacement pipe fluid contents."""
    # loads
    W1 = Weight('W1', 1)
    W1.apply()  # to selected elements

    FL1 = Fluid.from_file('FL1', 1, 'water')
    FL1.apply()

    # loadcase
    L1 = LoadCase('L1', 'sus', [Weight, Fluid], [1, 1])

    app.models('loads').analyze()

    pt30 = app.points(30)

    assert compare(L1.movements[pt30].dz, -2.4078)


def test_force(app):
    """Check displacement due to applied force."""
    pass


def test_displacement(app2):
    """Cantilever beam with displacement support"""

    # get pipe objects
    pt10 = app2.points(10)
    pt20 = app2.points(20)
    run10 = app2.elements(10, 15)
    run20 = app2.elements(15, 20)

    # support
    anc1 = Anchor('A1', 10)
    anc1.apply([run10])

    disp1 = Displacement('D1', 1, 20, dx=1)
    disp1.apply([run20])

    # loads
    F1 = Force('F1', 1, 20, fy=-10000)
    F1.apply([run20])

    # loadcase
    L1 = LoadCase('L1', 'ope', [Force, Displacement], [1, 1])

    app2.models('simple').analyze()

    # check dx at 10
    assert compare(L1.movements[pt20].dx, 1)


def test_seismic_x(app):
    """Check displacements due to seismic loading in x direction."""
    # loads
    S1 = Seismic('S1', 1, gx=1.0)
    S1.apply()

    L1 = LoadCase('L1', 'occ', [Seismic], [1])

    app.models('loads').analyze()

    pt30 = app.points(30)
    dx, ry, rz = 0, 4, 5

    assert compare(L1.movements[pt30][dx], 1.436)
    assert compare(L1.movements[pt30][ry], -0.1830)
    assert compare(L1.movements[pt30][rz], -0.4091)


def test_seismic_y(app):
    """Check displacements due to seismic loading in y direction."""
    # loads
    S1 = Seismic('S1', 1, gy=1.0)
    S1.apply()

    L1 = LoadCase('L1', 'occ', [Seismic], [1])

    app.models('loads').analyze()

    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], 0.075)
    assert compare(L1.movements[pt30][dz], 0.251)
    assert compare(L1.movements[pt30][rx], 0.1554)


def test_seismic_z(app):
    """Check displacements due to seismic loading in z direction."""
    # loads
    S1 = Seismic('S1', 1, gz=1.0)
    S1.apply()

    L1 = LoadCase('L1', 'occ', [Seismic], [1])

    app.models('loads').analyze()

    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dy], -0.330)
    assert compare(L1.movements[pt30][dz], 1.305)
    assert compare(L1.movements[pt30][rx], 0.3966)


def test_wind_x(app):
    """Check displacements due to wind loading in x direction."""
    pass


def test_wind_y(app):
    """Check displacements due to wind loading in y direction."""
    pass


def test_wind_z(app):
    """Check displacements due to wind loading in z direction."""
    pass
