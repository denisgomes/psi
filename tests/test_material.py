"""Pipe material property tests"""

# from .test_appsetup import app
import pytest

from psi.app import App
from psi.model import Model
from psi.material import Material


@pytest.fixture()
def app():
    app = App()
    mdl = Model('simple')

    # properties
    Material.from_file('MAT1', 'A53A', 'B31.1')

    return app


def test_rho(app):
    mat = app.materials("MAT1")

    assert round(mat.rho.value, 4) == 0.2830


def test_nu(app):
    mat = app.materials("MAT1")

    assert round(mat.nu.value, 1) == 0.3


def test_alp(app):
    mat = app.materials("MAT1")

    assert round(mat.alp[-325], 6) == 0.000005
    assert round(mat.alp[70], 8) == 0.00000607
    assert round(mat.alp[1100], 8) == 0.00000812


def test_ymod(app):
    mat = app.materials("MAT1")

    assert round(mat.ymod[-325], 0) == 30000000
    assert round(mat.ymod[70], 0) == 27900000
    assert round(mat.ymod[1100], 0) == 13000000


def test_sh(app):
    mat = app.materials("MAT1")

    assert round(mat.sh[-325], 0) == 12000
    assert round(mat.sh[70], 0) == 12000
    assert round(mat.sh[1000], 0) == 9000
