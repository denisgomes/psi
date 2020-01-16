"""Simple generic model used for testing.

This model is intended to be an easy way to setup the different pieces of PSI
for testing due to interdependance of classes. As such, it will evolve in time
as additional features and capabilities are added.

.. note::
    Additional example problems may be used to fully test the set of features.

Testing of the solver module should be used to determine the accuracy of the
results. It is not intended however to be a through check, which the V&V of the
software will achieve.
"""

import pytest

from psi.app import App
from psi.model import Model
from psi.sections import Pipe
from psi.material import Material
from psi.point import Point
from psi.elements import Run
from psi.supports import Anchor
from psi.loads import Weight, Pressure, Thermal, Force
from psi.loadcase import LoadCase, LoadComb
from psi.codes import B311


@pytest.fixture(scope="session")
def app():
    # parameter
    L = 10*12

    app = App()
    mdl = Model('simple')

    # properties
    Pipe.from_file('PIPE1', '10', '40')
    Material.from_file('MAT1', 'A53A', 'B31.1')

    # geometry
    Point(10)
    run20 = Run(20, L)

    # supports
    anc1 = Anchor('A1', 10)
    anc1.apply([run20])

    # d1 = Displacement('D1', 1, 10, dy=1.0)
    # d1.apply([run20])

    # v1 = GlobalY('V1', 20)
    # v1.apply([run20])

    # loads
    W1 = Weight('W1', 1)
    W1.apply([run20])

    P1 = Pressure('P1', 1, 250)
    P1.apply([run20])

    T1 = Thermal('T1', 1, 500, 70)
    T1.apply([run20])

    F1 = Force('F1', 1, 20, fy=-10000)
    F1.apply([run20])

    # loadcases
    L1 = LoadCase('L1', 'ope', [Weight, Pressure, Thermal, Force], [1, 1, 1, 1])
    L2 = LoadCase('L2', 'sus', [Weight, Pressure, Force], [1, 1, 1])
    LoadComb('L1-L2', 'exp', 'algebraic', [L1, L2], [1, -1])

    # code
    b311 = B311('B311')
    b311.apply([run20])

    # solve
    mdl.analyze()

    return app
