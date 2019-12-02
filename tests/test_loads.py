"""Test different loads against a simply supported cantilevered pipe."""

import unittest

from psi.app import App


class TestLoads(unittest.TestCase):

    def setUp(self):
        app = App()

        L = 10*12

        app.models.Model("test")

        app.sections.Pipe.from_file('pipe1', '10', '40')
        app.materials.Material.from_file('mat1', 'A53A', 'B31.1')

        app.points.Point(10)
        run20 = app.elements.Run(20, L)

        anc1 = app.supports.Anchor('A1', 10)
        anc1.apply([run20])

    def test_weight_disp(self):
        """Test pipe deadweight load"""
        pass

    def test_thermal_disp(self):
        """Test pipe thermal displacement"""
        pass

    def test_force_disp(self):
        """Test force displacement"""
        pass
