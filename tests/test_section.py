"""Pipe section property tests"""

import unittest

from psi.app import App


class TestSection(unittest.TestCase):

    def setUp(self):
        app = App()
        self.model = app.models.Model("test")
        self.p10sch40 = app.sections.Pipe.from_file("P10", "10", "40")

    def test_od(self):
        """Test diameter for 10 inch schedule 40 pipe"""
        self.assertAlmostEqual(self.p10sch40.od, 10.75, places=2)

    def test_thk(self):
        """Test pipe thickness for 10 inch schedule 40 pipe"""
        self.assertAlmostEqual(self.p10sch40.thk, 0.365, places=3)

    def test_area(self):
        """Test cross-sectional area of pipe"""
        self.assertAlmostEqual(self.p10sch40.area, 11.908, places=3)

    def test_ixx(self):
        """Test polar moment of inertia of pipe"""
        self.assertAlmostEqual(self.p10sch40.ixx, 321.468, places=3)

    def test_iyy(self):
        """Test area moment of inertia about y-y axis of pipe"""
        self.assertAlmostEqual(self.p10sch40.iyy, 160.734, places=3)

    def test_izz(self):
        """Test area moment of inertia about z-z axis of pipe"""
        self.assertAlmostEqual(self.p10sch40.izz, 160.734, places=3)

