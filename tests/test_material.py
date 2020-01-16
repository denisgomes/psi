"""Pipe material property tests"""


from .test_appsetup import app


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

    assert round(mat.ymod[-325], 0) == 31400000
    assert round(mat.ymod[70], 0) == 29500000
    assert round(mat.ymod[1100], 0) == 18000000


def test_sh(app):
    mat = app.materials("MAT1")

    assert round(mat.sh[-325], 0) == 30000
    assert round(mat.sh[70], 0) == 30000
    assert round(mat.sh[1000], 0) == 17800
