"""Pipe section property tests"""

import pytest

from psi.app import App


@pytest.fixture(scope="module")
def pipe():
    app = App()

    model = app.models.Model("test")
    model.units = "english"

    pipe = app.sections.Pipe.from_file("pipe", "10", "40")

    yield pipe


def test_od(pipe):
    assert round(pipe.od, 2) == 10.75


def test_thk(pipe):
    assert round(pipe.thk, 3) == 0.365


def test_area(pipe):
    assert round(pipe.area, 3) == 11.908


def test_ixx(pipe):
    assert round(pipe.ixx, 3) == 321.468


def test_iyy(pipe):
    assert round(pipe.iyy, 3) == 160.734


def test_izz(pipe):
    assert round(pipe.izz, 3) == 160.734
