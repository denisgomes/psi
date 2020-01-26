"""Pipe section property tests"""

import pytest

# from .test_appsetup import app

from psi.app import App
from psi.model import Model
from psi.sections import Pipe


@pytest.fixture()
def app():
    """A cantilevered beam model"""

    # parameter
    L = 10*12

    app = App()
    mdl = Model('simple')

    # properties
    Pipe.from_file('PIPE1', '10', '40')

    return app


def test_od(app):
    pipe = app.sections("PIPE1")

    assert round(pipe.od, 2) == 10.75


def test_thk(app):
    pipe = app.sections("PIPE1")

    assert round(pipe.thk, 3) == 0.365


def test_area(app):
    pipe = app.sections("PIPE1")

    assert round(pipe.area, 3) == 11.908


def test_ixx(app):
    pipe = app.sections("PIPE1")

    assert round(pipe.ixx, 3) == 321.468


def test_iyy(app):
    pipe = app.sections("PIPE1")

    assert round(pipe.iyy, 3) == 160.734


def test_izz(app):
    pipe = app.sections("PIPE1")

    assert round(pipe.izz, 3) == 160.734
