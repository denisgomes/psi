"""Test different loads against a simply supported cantilevered pipe."""

import pytest

from psi.app import App


@pytest.fixture(scope="module")
def run():
    app = App()

    app.models.Model("test")

    app.sections.Pipe.from_file('pipe1', '10', '40')
    app.materials.Material.from_file('mat1', 'A53A', 'B31.1')

    L = 10*12

    app.points.Point(10)
    run = app.elements.Run(20, L)

    yield run


def test_weight_disp(self):
    """Test pipe deadweight load"""
    pass


def test_thermal_disp(self):
    """Test pipe thermal displacement"""
    pass


def test_force_disp(self):
    """Test force displacement"""
    pass
