"""Testing various element types"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run, Rigid
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

    mdl = Model('elements')
    mdl.settings.vertical = "y"

    Pipe.from_file('pipe1', '10', '40')
    Material.from_file('mat1', 'A53A', 'B31.1')
    B31167('code1')

    return app


def test_rigid(app):
    """Check rigid element."""

    Point(10)
    run20 = Rigid(20, 10*12, weight=1000)

    anc10 = Anchor('anc10', 10)
    anc10.apply([run20])

    app.elements.select()   # select all

    # loads
    W1 = Weight('W1', 1)
    W1.apply()  # to selected elements

    # loadcase
    L1 = LoadCase('L1', 'sus', [Weight], [1])

    app.models('elements').analyze()

    pt10 = app.points(10)
    pt20 = app.points(20)

    rx, ry, rz = 0, 1, 2
    assert compare(L1.reactions[pt10][ry], -1000)
    # assert compare(L1.movements[pt20].dy, -0.0120)
