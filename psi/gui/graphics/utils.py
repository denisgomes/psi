"""Various utility routines for PSI gui."""

from __future__ import division
from itertools import izip_longest
from math import sin, cos, acos, degrees, radians, sqrt


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks - from python.org
    website
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def rotz2vec(v):
    """Return the two rotation angles namely, x_ang and y_ang needed to
    rotate the object z axis to an arbitary vector.
    """
    xzlen = sqrt(v.x*v.x + v.z*v.z)
    vlen = v.magnitude()
    if xzlen == 0:
        y_ang = 0
        x_ang = radians(90)
    else:
        y_ang = acos(v.z/xzlen)
        x_ang = acos(xzlen/vlen)

    if v.y > 0:
        x_ang = -x_ang
    else:
        x_ang = x_ang

    if v.x > 0:
        y_ang = y_ang
    else:
        y_ang = -y_ang

    y_ang = degrees(y_ang)
    x_ang = degrees(x_ang)

    return x_ang, y_ang


def axis_angle_to_rota(ax, ay, az, ang):
    n = sqrt(ax**2 + ay**2 + az**2)
    uax = ax / n
    uay = ay / n
    uaz = az / n
    cos_ang = cos(radians(ang))
    sin_ang = sin(radians(ang))
    _1_cos_ang = 1 - cos_ang

    # rotation matrix elements
    r00 = _1_cos_ang * uax**2 + cos_ang
    r01 = _1_cos_ang * uax * uay - uaz * sin_ang
    r02 = _1_cos_ang * uax * uaz + uay * sin_ang
    r10 = _1_cos_ang * uax * uay + uaz * sin_ang
    r11 = _1_cos_ang * uay**2 + cos_ang
    r12 = _1_cos_ang * uay * uaz - uax * sin_ang
    r20 = _1_cos_ang * uax * uaz - uay * sin_ang
    r21 = _1_cos_ang * uay * uaz + uax * sin_ang
    r22 = _1_cos_ang * uaz**2 + cos_ang

    return (r00, r01, r02,
            r10, r11, r12,
            r20, r21, r22)


def translate(vertices, dx, dy, dz):
    """Translate a list of vertices in place"""
    new = []
    for x, y, z in grouper(vertices, 3):
        new.extend((x+dx, y+dy, z+dz))
    vertices[:] = new
