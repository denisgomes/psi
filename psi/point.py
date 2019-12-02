# Copyright (c) 2019 Denis Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
#    * Neither the name of Pipe Stress Infinity (PSI) nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Creating the global starting coordinates of a piping model.

There are two ways in which a model can be created. The first method is
commonly used to define a branch point and the elements follow the point.
The elements all derive thier attributes, such as section and material
definitions from what is currently active.

Example
-------
>>> Point(10)
>>> Run(20, 10)
>>> Bend(30, 0, 20)
>>> ...

In the sequence above, the point object has the effect of creating a vertex
in global coordinates and defining the active point for the model.  This
allows subsequent calls to Run, Bend, etc., to use the active point to define
the next vertex using incremental offsets in the x, y, and z directions. Note
that a global point is created implicitly for the 'to' point and is accessed
via the container object.

The second way to create a model is to use the Point object to create points
in global coordinates first and then to call the alternate element constructor
to define the model.

Example
-------
>>> Point(10)   # point at (0, 0, 0)
>>> Point(20, 10)
>>> Point.by_point(20, 30, 0, 20)
>>> Run.by_point(10, 20)
>>> Bend.by_point(20, 30)

This is the traditional method used to describe the connectivity of a finite
element mesh consisting of beam elements. Here, the Point object is called
multiple times to create vertices in global coordinates explicitly.
"""

from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi import units


# TODO: Check app.points to see if point exists, __new__


@units.define(_x="length", _y="length", _z="length")
class Point(Entity, ActiveEntityMixin):
    """Define the global starting point for a branch."""

    def __new__(cls, name, *args, **kwargs):
        # call base class __new__
        inst = super(Point, cls).__new__(cls, name, *args, **kwargs)

        app = cls._app
        model = app.models.active_object
        if name in model.geometry.vertices:
            raise NameError("name '%s' in use" % name)

        return inst

    def __init__(self, point, x=0, y=0, z=0):
        super(Point, self).__init__(point)
        self._x = x
        self._y = y
        self._z = z
        self._sif = set()  # sifs at intersections / connections

        app = self.app
        model = app.models.active_object

        vert = model.geometry.MV(point, self._x, self._y, self._z)
        self._vertex = vert
        vert.data["elem"] = self

        self.activate()

    @classmethod
    def by_point(cls, refpnt, point, dx=0, dy=0, dz=0):
        """Create a point by a reference point and local coordinates"""
        app = cls._app
        rp = app.points.objects.get(refpnt, None)
        if rp:
            x, y, z = rp.xyz
            return cls(point, x+dx, y+dy, z+dz)
        return

    @property
    def vertex(self):
        return self._vertex

    @property
    def xyz(self):
        """Get the global position of the branch point"""
        return self.x, self.y, self.z

    @xyz.setter
    def xyz(self, value):
        """Set the global position of the branch point"""
        self.x, self.y, self.z = value

    @property
    def x(self):
        """Get the x component of the global position of the branch point"""
        return self._x

    @x.setter
    def x(self, value):
        """Set the x component of the global position of the branch point"""
        self._x = value

        # value is in user units
        self._vertex.co[0] = value

    @property
    def y(self):
        """Get the y component of the global position of the branch point"""
        return self._y

    @y.setter
    def y(self, value):
        """Set the y component of the global position of the branch point"""
        self._y = value
        self._vertex.co[1] = value

    @property
    def z(self):
        """Get the z component of the global position of the branch point"""
        return self._z

    @z.setter
    def z(self, value):
        """Set the z component of the global position of the branch point"""
        self._z = value
        self._vertex.co[2] = value

    @property
    def parent(self):
        return self.app.points


class PointManager(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(PointManager, self).__init__()
        self.Point = Point
