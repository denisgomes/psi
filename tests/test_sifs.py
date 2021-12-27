"""Testing various stress intensifications factors"""

import pytest

from psi.app import App
from psi.model import Model
from psi.point import Point
from psi.elements import Run
from psi.sections import Pipe
from psi.material import Material
from psi.codes.B311 import B31167
from psi.sifs import (Welding, Unreinforced, Reinforced, Weldolet, Sockolet,
                      Sweepolet, Weldolet, ButtWeld)
from psi.supports import Anchor
from psi.loads import Weight
from psi.loadcase import LoadCase
from .utils import compare


@pytest.fixture()
def app():
    app = App()
    points = app.points

    mdl = Model('sifs')
    mdl.settings.vertical = "z"

    Pipe.from_file('pipe1', '3.5', '40')
    Material.from_file('mat1', 'A53A', 'B31.1')
    B31167('code1')

    Point(10)
    run20 = Run(20, 5*12)
    Run(30, 5*12)

    points(20).activate()
    Run(40, 0, 0, 5*12)

    anc10 = Anchor('anc10', 10)
    anc10.apply([run20])

    app.elements.select()   # select all

    # loads
    W1 = Weight('W1', 1)
    W1.apply()  # to selected elements

    # loadcase
    L1 = LoadCase('L1', 'sus', [Weight], [1])

    app.models('sifs').analyze()

    pt30 = app.points(30)
    dy, dz, rx = 1, 2, 3

    assert compare(L1.movements[pt30][dz], -0.20901)

    return app


def test_welding_tee(app):
    """Check welding tee sif."""
    pt20 = app.points(20)
    run20 = app.elements(10, 20)

    # define welding tee
    tee20 = Welding('tee20', 20, 3.5, 0.216, 3.5, 0, 0)
    tee20.apply([run20])

    code = app.codes('code1')
    sifi20 = code.sifi(run20, pt20)

    assert compare(sifi20, 1.6366)


def test_reinforced_tee(app):
    """Check reinforced fabricated tee sif."""
    pt20 = app.points(20)
    run20 = app.elements(10, 20)

    # define reinforced tee
    tee20 = Reinforced('tee20', 20, 3.5, 0.216, 0.25)
    tee20.apply([run20])

    code = app.codes('code1')
    sifi20 = code.sifi(run20, pt20)

    assert compare(sifi20, 1.6366)


def test_unreinforced_tee(app):
    """Check unreinforced fabricated tee sif."""
    pt20 = app.points(20)
    run20 = app.elements(10, 20)

    # define reinforced tee
    tee20 = Unreinforced('tee20', 20, 3.5, 0.216)
    tee20.apply([run20])

    code = app.codes('code1')
    sifi20 = code.sifi(run20, pt20)

    assert compare(sifi20, 3.4795)


def test_weldolet(app):
    """Check weldolet sif."""
    pt20 = app.points(20)
    run20 = app.elements(10, 20)

    # define weldolet
    tee20 = Weldolet('tee20', 20, 3.5, 0.216)
    tee20.apply([run20])

    code = app.codes('code1')
    sifi20 = code.sifi(run20, pt20)

    assert compare(sifi20, 1.5698)


def test_sockolet(app):
    """Check sockolet sif."""
    pass


def test_sweepolet(app):
    """Check sweepolet sif."""
    pt20 = app.points(20)
    run20 = app.elements(10, 20)

    # define sweepolet
    tee20 = Sweepolet('tee20', 20, 3.5, 0.216, 3.5, 0, 0)
    tee20.apply([run20])

    code = app.codes('code1')
    sifi20 = code.sifi(run20, pt20)

    assert compare(sifi20, 1.6366)
