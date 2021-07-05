from __future__ import division
from math import pi, sin, cos, sqrt, radians, degrees
# from itertools import izip
from functools import partial

from pyglet.gl import *
from pyglet import graphics

from psi.settings import options
from psi.gui.utils import grouper, rotz2vec, axis_angle_to_rota
from psi.gui.bounds import AABB

from psi.utils.euclid import Vector3, Point3, Matrix4


class FaceAlphaGroup(graphics.Group):

    def __init__(self, parent=None):
        super(FaceAlphaGroup, self).__init__(parent)

    def set_state(self):
        glPushAttrib(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT | GL_ENABLE_BIT)
        glColor4f(1.0, 1.0, 1.0, 0.5)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_BLEND)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)

    def unset_state(self):
        glPopAttrib()


class LineAlphaGroup(FaceAlphaGroup):

    def __init__(self, parent=None):
        super(LineAlphaGroup, self).__init__(parent)

    def set_state(self):
        glPushAttrib(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT | GL_ENABLE_BIT)
        glDisable(GL_LIGHTING)  # no lighting for lines and points
        glColor4f(1.0, 1.0, 1.0, 0.5)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_BLEND)
        glDisable(GL_DEPTH_TEST)


class PointAlphaGroup(LineAlphaGroup):

    def __init__(self, parent=None):
        super(PointAlphaGroup, self).__init__(parent)

# shared groups between all transparent primitives with similar modes
face_alpha_group = FaceAlphaGroup()
line_alpha_group = LineAlphaGroup()
point_alpha_group = PointAlphaGroup()


class PointGroup(graphics.Group):
    """Disables lighting and specifies a size for points"""

    def __init__(self, size=5, parent=None):
        super(PointGroup, self).__init__(parent)
        self._size = size

    def set_state(self):
        glPushAttrib(GL_LIGHTING_BIT | GL_ENABLE_BIT)
        glDisable(GL_LIGHTING)
        glDisable(GL_DEPTH_TEST)    # show points on top
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glPointSize(self._size)

    def unset_state(self):
        glPopAttrib()

    def _get_size(self):
        return self._size

    def _set_size(self, size):
        self._size = size

    size = property(_get_size, _set_size)


class LineGroup(graphics.Group):
    """Disables lighting and specifies a width for lines"""

    def __init__(self, width=1, parent=None):
        super(LineGroup, self).__init__(parent)
        self._width = width

    def set_state(self):
        glPushAttrib(GL_LIGHTING_BIT | GL_ENABLE_BIT)
        glDisable(GL_LIGHTING)
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glLineWidth(self._width)

    def unset_state(self):
        glPopAttrib()

    def _get_width(self):
        return self._width

    def _set_width(self, width):
        self._width = width

    width = property(_get_width, _set_width)

_point_group = PointGroup()
_line_group = LineGroup()


class Primitive(object):
    """All primitives are drawn with respect to their local object coordinate
    system, with the x axis pointing to the right, the y axis pointing up and
    the z axis pointing out of the screen.

    The anchor point is at the orgin of the object coordinate when a primitive
    is first created.  The anchor point can be moved to a different location
    with respect to the object coordinate system (ie. the local system).  All
    translations and rotations use the anchor point as the reference.  The
    position of the primitive is really the global location of the anchor
    point.  When setting the anchor point, it is important to remember that one
    must use the local object coordinate frame of reference.
    """

    __slots__ = ("_anc_x", "_anc_y", "_anc_z", "_x", "_y", "_z", "_ax", "_ay",
                 "_az", "_ang", "_mode", "_batch", "_group", "_vertex_list",
                 "_is_visible", "_is_picked", "_rgb", "_alpha", "_is_translated",
                 "_is_rotated", "_bound")

    def __init__(self, position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                 batch=None, group=None):
        self._anc_x, self._anc_y, self._anc_z = anchor
        self._x, self._y, self._z = position
        self._ax, self._ay, self._az, self._ang = (0, 0, 0, 0)
        self._mode = mode
        self._vertex_list = None
        self._is_visible = True
        self._is_picked = False
        self._is_translated = False if position == (0, 0, 0) else True
        self._is_rotated = False
        self._alpha = 255
        self._bound = AABB()

        self._batch = batch if batch is not None else graphics.Batch()

        if self._mode == GL_TRIANGLES:
            self._rgb = options["graphics.face_color"]
        elif self._mode == GL_LINES:
            self._rgb = options["graphics.line_color"]
        elif self._mode == GL_POINTS:
            self._rgb = options["graphics.point_color"]

        self._group = group

    def add_indexed(self, count, indices, *data):
        self._vertex_list = self._batch.add_indexed(count, self._mode,
                                                    self._group, indices,
                                                    *data)
        self._bound.update(self._vertex_list.vertices)
        self._update_position()
        self._update_color()

    def add(self, count, *data):
        self._vertex_list = self._batch.add(count, self._mode, self._group,
                                            *data)
        self._bound.update(self._vertex_list.vertices)
        self._update_position()
        self._update_color()

    def __del__(self):
        try:
            if self._vertex_list is not None:
                self._vertex_list.delete()
        except:
            pass

    def delete(self):
        self._vertex_list.delete()
        self._vertex_list = None
        self._group = None

    def _get_mode(self):
        return self._mode

    def _set_mode(self, mode):
        raise NotImplementedError

    mode = property(_get_mode, _set_mode)

    def _set_batch(self, batch):
        if self._batch == batch:
            return

        self._batch.migrate(self._vertex_list, self._mode, self._group, batch)
        self._batch = batch

    def _get_batch(self):
        return self._batch

    batch = property(_get_batch, _set_batch)

    def _get_vertex_list(self):
        return self._vertex_list

    def _set_vertex_list(self):
        pass

    vertex_list = property(_get_vertex_list, _set_vertex_list)

    def _set_group(self, group):
        self._group = group
        self._batch.migrate(self._vertex_list, self._mode, self._group,
                                self._batch)

    def _get_group(self):
        #return self._group.parent
        return self._group

    group = property(_get_group, _set_group)

    def _set_alpha(self, alpha):
        self._alpha = alpha
        self._update_color()

        if not isinstance(self._group, FaceAlphaGroup):
            if alpha < 255:
                if self._mode == GL_TRIANGLES:
                    self.group = face_alpha_group
                elif self._mode == GL_LINES:
                    self.group = line_alpha_group
                elif self._mode == GL_POINTS:
                    self.group = point_alpha_group
            else:
                self.group = None

    alpha = property(lambda self: self._alpha, _set_alpha)

    def _set_color(self, rgb):
        self._rgb = map(int, rgb)
        self._update_color()

    color = property(lambda self: self._rgb, _set_color)

    def _update_position(self):
        if self._is_translated:
            # object coordinate system translation given by dx, dy, dz
            vertices = []
            dx, dy, dz = (self._x - self._anc_x, self._y - self._anc_y,
                            self._z - self._anc_z)
            for x, y, z in grouper(self._vertex_list.vertices, 3):
                vertices.extend((x+dx, y+dy, z+dz))
            self._vertex_list.vertices[:] = vertices

            self._is_translated = False
            self._bound.update(self._vertex_list.vertices)
        elif self._is_rotated:
            # move to origin -> rotate -> move back
            # ([T][R][-T]) * [x, y, z, 1] - for vertices
            # ([T][R][-T]) * [nx, ny, nz, 0] - for normals
            # normal rotation does not support non-uniform scaling

            # rotate about anchor point position
            dx, dy, dz = (self._x, self._y, self._z)

            (r00, r01, r02, # rotation about anchor point
             r10, r11, r12,
             r20, r21, r22) = axis_angle_to_rota(self._ax, self._ay, self._az,
                                                    self._ang)
            vertices = []
            normals = []
            attribs = zip(self._vertex_list.vertices, self._vertex_list.normals)
            for (x, nx), (y, ny), (z, nz) in grouper(attribs, 3):
                vertices.extend((r00*x + r01*y + r02*z - dx*r00 - dy*r01 - dz*r02 + dx,
                                 r10*x + r11*y + r12*z - dx*r10 - dy*r11 - dz*r12 + dy,
                                 r20*x + r21*y + r22*z - dx*r20 - dy*r21 - dz*r22 + dz))
                normals.extend((r00*nx + r01*ny + r02*nz,
                                r10*nx + r11*ny + r12*nz,
                                r20*nx + r21*ny + r22*nz))
            self._vertex_list.vertices[:] = vertices
            self._vertex_list.normals[:] = normals

            self._is_rotated = False

            self._bound.update(self._vertex_list.vertices)

    def _update_color(self):
        if self._is_picked:
            r, g, b = options["graphics.pick_color"]
        else:
            r, g, b = self._rgb

        if self._is_visible:
            alpha = int(self._alpha)
        else:
            alpha = 0

        self._vertex_list.colors[:] = [r, g, b, alpha] * \
                                        self._vertex_list.get_size()

    def _set_anchor(self, x, y, z):
        """The anchor point is the reference point which is used as the focul
        point for translation and rotation.  It is defined with respect to the
        object coordinate system.
        """
        # object coordinate position
        obj_x, obj_y, obj_z = (self._x - self._anc_x, self._y - self._anc_y,
                                self._z - self._anc_z)
        self._anc_x, self._anc_y, self._anc_z = x, y, z

        # TODO: test!
        # update the position vector of the anchor point
        self._x, self._y, self._z = obj_x + x, obj_y + y, obj_z + z

    anchor = property(lambda self: (self._anc_x, self._anc_y, self._anc_z),
                            lambda self, t: self._set_anchor(*t))

    def _set_position(self, x, y, z):
        """Set the position of the anchor point with respect to world
        coordinates.
        """
        self._x = x
        self._y = y
        self._z = z
        self._is_translated = True
        self._update_position()

    position = property(lambda self: (self._x, self._y, self._z),
                            lambda self, t: self._set_position(*t))

    def _set_x(self, x):
        self._x = x
        self._is_translated = True
        self._update_position()

    x = property(lambda self: self._x, _set_x)

    def _set_y(self, y):
        self._y = y
        self._is_translated = True
        self._update_position()

    y = property(lambda self: self._y, _set_y)

    def _set_z(self, z):
        self._z = z
        self._is_translated = True
        self._update_position()

    z = property(lambda self: self._z, _set_z)

    def _set_rotation(self, ax, ay, az, ang):
        self._ax = ax
        self._ay = ay
        self._az = az
        self._ang = ang
        self._is_rotated = True
        self._update_position()

    rotation = property(lambda self: (self._ax, self._ay, self._az, self._ang),
                            lambda self, t: self._set_rotation(*t))

    def _set_is_visible(self, visible):
        self._is_visible = visible
        self._update_color()

    is_visible = property(lambda self: self._is_visible, _set_is_visible)

    def _set_is_picked(self, picked):
        self._is_picked = picked
        self._update_color()

    is_picked = property(lambda self: self._is_picked, _set_is_picked)

    def _set_bound(self, bound):
        self._bound = bound

    bound = property(lambda self: self._bound, _set_bound)

    def draw(self):
        if self._is_visible:
            self._group.set_state_recursive()
            self._vertex_list.draw(self._mode)
            self._group.unset_state_recursive()


class ComplexPrimitive(object):
    """A collection of primitives that share a common anchor and bounding
    volume.
    """

    __slots__ = ("_anc_x", "_anc_y", "_anc_z", "_x", "_y", "_z", "_ax", "_ay",
                 "_az", "_rgb", "_alpha", "_ang", "_mode", "is_visible",
                 "_is_picked", "_bound", "_primitives", "_batch", "_group")

    def __init__(self, position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                 batch=None, group=None):
        self._anc_x, self._anc_y, self._anc_z = anchor
        self._x, self._y, self._z = position
        self._mode = mode
        self._rgb = (255, 255, 255)
        self._alpha = 255
        self._primitives = set()
        self._is_visible = True
        self._is_picked = False
        self._bound = AABB()
        self._batch = batch
        self._group = group

    def add(self, prim):
        prim.position = self.position
        prim.anchor = self.anchor
        self._primitives.add(prim)
        self._update_bound()

    def _get_mode(self):
        return self._mode

    def _set_mode(self, mode):
        raise NotImplementedError

    mode = property(_get_mode, _set_mode)

    def _set_batch(self, batch):
        if self._batch == batch:
            return

        self._batch = batch
        for prim in self._primitives:
            prim.batch = batch

    def _get_batch(self):
        return self._batch

    batch = property(_get_batch, _set_batch)

    def _set_group(self, group):
        self._group = group

        for prim in self._primtives:
            prim.group = group

    def _get_group(self):
        return self._group

    group = property(_get_group, _set_group)

    def _get_primitives(self):
        return self._primitives

    def _set_primitives(self):
        return

    primitives = property(_get_primitives, _set_primitives)

    def _set_anchor(self, x, y, z):
        self._anc_x = x
        self._anc_y = y
        self._anc_z = z

        for prim in self._primitives:
            prim.anchor = (x, y, z)

    anchor = property(lambda self: (self._anc_x, self._anc_y, self._anc_z),
                            lambda self, t: self._set_anchor(*t))

    def _set_rotation(self, ax, ay, az, ang):
        self._ax = ax
        self._ay = ay
        self._az = az
        self._ang = ang

        for prim in self._primitives:
            prim.rotation = (ax, ay, az, ang)

        self._update_bound()

    rotation = property(lambda self: (self._ax, self._ay, self._az, self._ang),
                            lambda self, t: self._set_rotation(*t))

    def _set_position(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z

        for prim in self._primitives:
            prim.position = (x, y, z)

        self._update_bound()

    position = property(lambda self: (self._x, self._y, self._z),
                            lambda self, t: self._set_position(*t))

    def _set_x(self, x):
        self._x = x

        for prim in self._primitives:
            prim.x = x

        self._update_bound()

    x = property(lambda self: self._x, _set_x)

    def _set_y(self, y):
        self._y = y

        for prim in self._primitives:
            prim.y = y

        self._update_bound()

    y = property(lambda self: self._y, _set_y)

    def _set_z(self, z):
        self._z = z

        for prim in self._primitives:
            prim.z = z

        self._update_bound()

    z = property(lambda self: self._z, _set_z)

    def _set_is_visible(self, visible):
        self._is_visible = visible

        for prim in self._primitives:
            prim.is_visible = visible

    is_visible = property(lambda self: self._is_visible, _set_is_visible)

    def _set_is_picked(self, picked):
        self._is_picked = picked

        for prim in self._primitives:
            prim.is_picked = picked

    is_picked = property(lambda self: self._is_picked, _set_is_picked)

    def _set_bound(self, bound):
        self._bound = bound

    bound = property(lambda self: self._bound, _set_bound)

    def _set_alpha(self, alpha):
        self._alpha = alpha

        for prim in self._primitives:
            prim.alpha = alpha

    alpha = property(lambda self: self._alpha, _set_alpha)

    def _set_color(self, rgb):
        self._rgb = map(int, rgb)

        for prim in self._primitives:
            prim.color = rgb

    color = property(lambda self: self._rgb, _set_color)

    def delete(self):
        for prim in self._primitives:
            prim.delete()

    def _update_bound(self):
        self._bound.clear()
        for prim in self._primitives:
            self._bound.merge(prim.bound)


class Cylinder(Primitive):

    @classmethod
    def Cylinder1(cls, v1, v2, diameter, slices=21, stacks=1, batch=None,
                  group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        height = v.magnitude()
        radius = diameter / 2
        cyl1 = Cylinder(radius, height, slices, stacks, batch=batch,
                        group=group)
        x_ang, y_ang = rotz2vec(v)

        cyl1.position = (v1.x, v1.y, v1.z)
        cyl1.rotation = (1, 0, 0, x_ang)
        cyl1.rotation = (0, 1, 0, y_ang)

        return cyl1

    def __init__(self, radius, height, slices=21, stacks=1, position=(0, 0, 0),
                 anchor=(0, 0, 0), mode=GL_TRIANGLES, batch=None, group=None):
        super(Cylinder, self).__init__(position, anchor, mode, batch, group)

        vertices = []
        normals = []
        colors = []
        textures = []

        slices += 1
        u_step = 2*pi / (slices-1)
        v_step = height / stacks
        v = 0.
        for i in xrange(stacks+1):
            u = 0.
            for j in xrange(slices):
                cos_u = cos(u)
                sin_u = sin(u)
                x = radius*cos_u
                y = radius*sin_u
                z = v
                nx = cos_u
                ny = sin_u
                nz = 0.
                vertices.extend([x, y, z])
                normals.extend([nx, ny, nz])
                colors.extend([255]*4)
                textures.extend([z/height, z/height])
                u += u_step
            v += v_step

        indices = []
        for i in xrange(stacks):
            for j in xrange(slices-1):
                p = i*slices+j
                indices.extend([p, p+slices+1, p+slices, p, p+1, p+slices+1])

        self.add_indexed(len(vertices)//3, indices,
                         ('v3f/dynamic', vertices),
                         ('n3f/dynamic', normals),
                         ('c4B/dynamic', colors),
                         ('t2f/dynamic', textures),
                         )


class Line(Primitive):

    def __init__(self, pt1, pt2, mode=GL_LINES, batch=None, group=_line_group):
        position = (0, 0, 0)
        anchor = (0, 0, 0)
        super(Line, self).__init__(position, anchor, mode, batch, group)

        vertices = []
        for pt in (pt1, pt2):
            vertices.extend(pt)

        color = [255, 255, 255, 255,
                 255, 255, 255, 255]

        self.add(2, ('v3f/dynamic', vertices),
                    ('c4B/dynamic', color))

    @property
    def pt1(self):
        return self.vertex_list[:3]

    @pt1.setter
    def pt1(self, x, y, z):
        self.vertex_list[:3] = x, y, z
        self._is_translated = True
        self._update_position()

    @property
    def pt2(self):
        return self.vertex_list[3:6]

    @pt2.setter
    def pt2(self, x, y, z):
        self.vertex_list[3:6] = x, y, z
        self._is_translated = True
        self._update_position()


class Point(Primitive):
    def __init__(self, position, mode=GL_POINTS, batch=None, group=_point_group):
        # position and anchor must be the same
        super(Point, self).__init__(position, position, mode, batch, group)

        color = [255, 255, 255, 255]

        self.add(1, ('v3f/dynamic', position),
                    ('c4B/dynamic', color))

    def _set_position(self, x, y, z):
        """The anchor point location and the position are the same for a
        point"""
        self._x = x
        self._y = y
        self._z = z
        self._anc_x = x
        self._anc_y = y
        self._anc_z = z
        self._is_translated = True
        self._update_position()

    position = property(lambda self: (self._x, self._y, self._z),
                        lambda self, t: self._set_position(*t))


class Cone(Primitive):

    @classmethod
    def Cone1(cls, v1, v2, radius, slices=21, stacks=5, batch=None,
              group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        height = v.norm()
        cone1 = Cone(radius, height, slices, stacks, batch=batch, group=group)
        x_ang, y_ang = rotz2vec(v)
        cone1.position = (v1.x, v1.y, v1.z)
        cone1.rotation = (1, 0, 0, x_ang)
        cone1.rotation = (0, 1, 0, y_ang)
        return cone1

    def __init__(self, radius, height, slices=21, stacks=5, position=(0, 0, 0),
                    anchor=(0, 0, 0), mode=GL_TRIANGLES, batch=None,
                    group=None):
        super(Cone, self).__init__(position, anchor, mode, batch, group)

        vertices = []
        normals = []
        colors = []
        textures = []

        slices += 1
        u_step = 2*pi / (slices - 1)
        v_step = height / stacks
        r_step = -(radius / height)
        v_step_sq = v_step**2
        v = 0.0
        for i in xrange(stacks+1):
            cos_v = cos(v)
            sin_v = sin(v)
            u = 0.0
            for j in xrange(slices):
                cos_u = cos(u)
                sin_u = sin(u)
                d = r_step*v + radius
                x = d*cos_u
                y = d*sin_u
                z = v
                nx = x/sqrt(d**2 + v_step_sq)
                ny = y/sqrt(d**2 + v_step_sq)
                nz = v_step/sqrt(d**2 + v_step_sq)
                vertices.extend([x, y, z])
                normals.extend([nx, ny, nz])
                colors.extend([255]*4)
                textures.extend([z/height, z/height])
                u += u_step
            v += v_step

        indices = []
        for i in xrange(stacks):
            for j in xrange(slices - 1):
                p = i*slices + j
                indices.extend([p, p+slices+1, p+slices, p, p+1, p+slices+1])

        self.add_indexed(len(vertices)//3, indices, ('v3f/dynamic', vertices),
                                                    ('n3f/dynamic', normals),
                                                    ('c4B/dynamic', colors),
                                                    ('t2f/dynamic', textures),
                                                    )


class Disk(Primitive):

    def __init__(self, radius, slices=12, position=(0, 0, 0), anchor=(0, 0, 0),
                 mode=GL_TRIANGLES, batch=None, group=None):

        super(Disk, self).__init__(position, anchor, mode, batch, group)

        vertices = [0, 0, 0]
        normals = [0, 0, 1]
        colors = [255, 255, 255, 255]
        textures = [0, 0]
        indices = []

        u_step = radians(360/slices)
        for i in xrange(slices+1):
            x = radius * cos(u_step*i)
            y = radius * sin(u_step*i)
            z = 0
            vertices.extend([x, y, z])
            normals.extend([0, 0, 1])
            colors.extend([255, 255, 255, 255])
            textures.extend([0, 0])

        for i in xrange(slices):
            indices.extend([0, i+1, i+2])

        self.add_indexed(len(vertices)//3, indices, ('v3f/dynamic', vertices),
                                                    ('n3f/dynamic', normals),
                                                    ('c4B/dynamic', colors),
                                                    ('t2f/dynamic', textures),
                                                    )

class Reducer(Primitive):

    @classmethod
    def Reducer1(cls, v1, v2, od, od1, slices=21, stacks=1, batch=None,
                 group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        height = v.magnitude()
        red1 = Reducer(od/2, od1/2, height, slices, stacks, batch=batch,
                       group=group)
        x_ang, y_ang = rotz2vec(v)
        red1.position = (v1.x, v1.y, v1.z)
        red1.rotation = (1, 0, 0, x_ang)
        red1.rotation = (0, 1, 0, y_ang)
        return red1

    def __init__(self, radius, radius1, height, slices=21, stacks=1,
                 position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                 batch=None, group=None):
        super(Reducer, self).__init__(position, anchor, mode, batch, group)

        vertices = []
        normals = []
        colors = []
        textures = []

        slices += 1
        u_step = 2*pi / (slices - 1)
        v_step = height / stacks
        r_step = -((radius-radius1) / height)
        v_step_sq = v_step**2
        v = 0.0
        for i in xrange(stacks+1):
            cos_v = cos(v)
            sin_v = sin(v)
            u = 0.0
            for j in xrange(slices):
                cos_u = cos(u)
                sin_u = sin(u)
                d = r_step*v + radius
                x = d*cos_u
                y = d*sin_u
                z = v
                if v == 0 or v == height:
                    nx = cos_u
                    ny = sin_u
                    nz = 0.
                else:
                    nx = -x/sqrt(d**2 + v_step_sq)
                    ny = -y/sqrt(d**2 + v_step_sq)
                    nz = -v_step/sqrt(d**2 + v_step_sq)
                vertices.extend([x, y, z])
                normals.extend([nx, ny, nz])
                colors.extend([255]*4)
                textures.extend([z/height, z/height])
                u += u_step
            v += v_step

        indices = []
        for i in xrange(stacks):
            for j in xrange(slices - 1):
                p = i*slices + j
                indices.extend([p, p+slices+1, p+slices, p, p+1, p+slices+1])

        self.add_indexed(len(vertices)//3, indices, ('v3f/dynamic', vertices),
                                                    ('n3f/dynamic', normals),
                                                    ('c4B/dynamic', colors),
                                                    ('t2f/dynamic', textures),
                                                    )

class Torus(Primitive):

    def __init__(self, radius, inner_radius, slices=24, inner_slices=12,
                    position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                    batch=None, group=None):
        super(Torus, self).__init__(position, anchor, mode, batch, group)

        vertices = []
        normals = []
        colors = []
        textures = []

        slices += 1
        inner_slices += 1
        u_step = 2*pi / (slices-1)
        v_step = 2*pi / (inner_slices-1)
        u = 0.
        for i in xrange(slices):
            cos_u = cos(u)
            sin_u = sin(u)
            v = 0.
            for j in xrange(inner_slices):
                cos_v = cos(v)
                sin_v = sin(v)

                d = (radius + inner_radius * cos_v)
                x = d * cos_u
                y = d * sin_u
                z = inner_radius * sin_v

                nx = cos_u * cos_v
                ny = sin_u * cos_v
                nz = sin_v

                vertices.extend([x, y, z])
                normals.extend([nx, ny, nz])
                colors.extend([255]*4)
                textures.extend([0, 0])
                v += v_step
            u += u_step

        indices = []
        for i in xrange(slices - 1):
            for j in xrange(inner_slices - 1):
                p = i * inner_slices + j
                indices.extend([p, p+inner_slices, p+inner_slices+1,
                                    p, p+inner_slices+1, p+1])

        self.add_indexed(len(vertices)//3, indices, ('v3f/dynamic', vertices),
                                                    ('n3f/dynamic', normals),
                                                    ('c4B/dynamic', colors),
                                                    ('t2f/dynamic', textures),
                                                    )

class Valve(ComplexPrimitive):

    @classmethod
    def Valve1(cls, v1, v2, diameter, slices=21, stacks=6, batch=None,
               group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        height = v.magnitude()
        radius = diameter / 2
        valv1 = Valve(radius, height, slices=slices, stacks=stacks,
                      batch=batch, group=group)
        x_ang, y_ang = rotz2vec(v)
        valv1.position = (v1.x, v1.y, v1.z)
        valv1.rotation = (1, 0, 0, x_ang)
        valv1.rotation = (0, 1, 0, y_ang)
        return valv1

    def __init__(self, radius, length, slices=21, stacks=6, position=(0, 0, 0),
                 anchor=(0, 0, 0), mode=GL_TRIANGLES, batch=None, group=None):
        super(Valve, self).__init__(position, anchor, mode, batch, group)

        left_valve_body_cap = Disk(radius, slices=slices, mode=mode,
                                   batch=batch, group=group)

        left_valve_body_cap.rotation = (1, 0, 0, 180)
        left_valve_body = Cone(radius, length/2, slices=slices, stacks=stacks,
                               mode=mode, batch=batch, group=group)

        right_valve_body = Cone(radius, length/2, slices=slices, stacks=stacks,
                                mode=mode, batch=batch, group=group)
        right_valve_body.rotation = (1, 0, 0, 180)
        right_valve_body.z = length

        right_valve_body_cap = Disk(radius, slices=slices, mode=mode,
                                    batch=batch, group=group)
        right_valve_body_cap.z = length

        self.add(left_valve_body_cap)
        self.add(left_valve_body)
        self.add(right_valve_body)
        self.add(right_valve_body_cap)


class Box(Primitive):

    @classmethod
    def Box1(cls, v1, v2, diameter, batch=None, group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)

        v = v2 - v1
        height = v.magnitude()

        box = Box(1.25*diameter, 1.25*diameter, 0.1*height, batch=batch,
                  group=group)

        x_ang, y_ang = rotz2vec(v)

        box.position = (v1.x, v1.y, v1.z)
        box.rotation = (1, 0, 0, x_ang)
        box.rotation = (0, 1, 0, y_ang)

        return box

    def __init__(self, base, height, depth, position=(0, 0, 0), anchor=(0, 0, 0),
                    mode=GL_TRIANGLES, batch=None, group=None):
        super(Box, self).__init__(position, anchor, mode, batch, group)

        vertices = []
        normals = []
        colors = []
        textures = []
        indices = []
        x = base / 2
        y = height / 2
        z = depth / 2

        vertices.extend([ x, y, z,
                         -x, y, z,
                         -x, -y, z,
                         x, -y, z,
                         x, y, z,
                         x, -y, z,
                         x, -y, -z,
                         x, y, -z,
                         x, y, z,
                         x, y, -z,
                         -x, y, -z,
                         -x, y, z,
                         -x, y, z,
                         -x, y, -z,
                         -x, -y, -z,
                         -x, -y, z,
                         -x, -y, -z,
                         x, -y, -z,
                         x, -y, z,
                         -x, -y, z,
                         x, -y, -z,
                         -x, -y, -z,
                         -x, y, -z,
                         x, y, -z])

        normals.extend([ 0, 0, 1,
                         0, 0, 1,
                         0, 0, 1,
                         0, 0, 1,
                         1, 0, 0,
                         1, 0, 0,
                         1, 0, 0,
                         1, 0, 0,
                         0, 1, 0,
                         0, 1, 0,
                         0, 1, 0,
                         0, 1, 0,
                        -1, 0, 0,
                        -1, 0, 0,
                        -1, 0, 0,
                        -1, 0, 0,
                         0, -1, 0,
                         0, -1, 0,
                         0, -1, 0,
                         0, -1, 0,
                         0, 0, -1,
                         0, 0, -1,
                         0, 0, -1,
                         0, 0, -1])

        indices.extend([0, 1, 2,    # triangle indices
                        0, 2, 3,
                        4, 5, 6,
                        4, 6, 7,
                        8, 9, 10,
                        8, 10, 11,
                        12, 13, 14,
                        12, 14, 15,
                        16, 17, 18,
                        16, 18, 19,
                        20, 21, 22,
                        20, 22, 23])

        colors = [255, 255, 255, 255] * 24

        textures = [0, 0] * 24

        self.add_indexed(24, indices, ('v3f/dynamic', vertices),
                                      ('n3f/dynamic', normals),
                                      ('c4B/dynamic', colors),
                                      ('t2f/dynamic', textures),
                                      )

class Arrow(ComplexPrimitive):

    @classmethod
    def Arrow1(cls, v1, v2, diameter, slices=12, stacks=5, batch=None,
               group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        length = v.magnitude()
        arrow1 = Arrow(diameter, length, slices=slices, stacks=stacks,
                       batch=batch, group=group)
        x_ang, y_ang = rotz2vec(v)
        arrow1.position = (v1.x, v1.y, v1.z)
        arrow1.rotation = (1, 0, 0, x_ang)
        arrow1.rotation = (0, 1, 0, y_ang)
        return arrow1

    def __init__(self, diameter, length, slices=12, stacks=5,
                 position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                 batch=None, group=None):
        super(Arrow, self).__init__(position, anchor, mode, batch, group)

        x_arrow_body = Cylinder(diameter/3, length/2, slices=slices,
                                stacks=stacks, mode=mode, batch=batch,
                                group=group)

        x_arrow_body_cap = Disk(diameter/3, slices=slices, mode=mode,
                                batch=batch, group=group)
        x_arrow_body_cap.rotation = (1, 0, 0, 180)

        x_arrow_head = Cone(1.25*diameter, length/2, slices=slices, stacks=stacks,
                            mode=mode, batch=batch, group=group)
        x_arrow_head.z = length/2

        x_arrow_head_cap = Disk(1.25*diameter, slices=slices, mode=mode,
                                batch=batch, group=group)
        x_arrow_head_cap.rotation = (0, 1, 0, 180)
        x_arrow_head_cap.z = length/2

        self.add(x_arrow_body)
        self.add(x_arrow_body_cap)
        self.add(x_arrow_head)
        self.add(x_arrow_head_cap)


class Flange(ComplexPrimitive):

    @classmethod
    def Flange1(cls, v1, v2, diameter, slices=21, stacks=1, batch=None,
                group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        length = v.magnitude()
        flange = Flange(diameter, length, slices=slices, stacks=stacks,
                        batch=batch, group=group)
        x_ang, y_ang = rotz2vec(v)
        flange.position = (v1.x, v1.y, v1.z)
        flange.rotation = (1, 0, 0, x_ang)
        flange.rotation = (0, 1, 0, y_ang)

        return flange

    def __init__(self, diameter, length, slices=21, stacks=1,
                 position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                 batch=None, group=None):
        super(Flange, self).__init__(position, anchor, mode, batch, group)

        dsk1 = Disk(diameter/2*1.25, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk1.rotation = (1, 0, 0, 180)

        flng1 = Cylinder(diameter/2*1.25, 2*length/5, slices=slices,
                         stacks=stacks, mode=mode, batch=batch,
                         group=group)

        dsk2 = Disk(diameter/2*1.25, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk2.z = 2*length/5

        flng2 = Cylinder(diameter/2, length/5, slices=slices,
                         stacks=stacks, mode=mode, batch=batch,
                         group=group)
        flng2.z = 2*length/5

        dsk3 = Disk(diameter/2*1.25, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk3.z = 3*length/5
        dsk3.rotation = (1, 0, 0, 180)

        flng3 = Cylinder(diameter/2*1.25, 2*length/5, slices=slices,
                         stacks=stacks, mode=mode, batch=batch,
                         group=group)
        flng3.z = 3*length/5

        dsk4 = Disk(diameter/2*1.25, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk4.z = length

        self.add(dsk1)
        self.add(flng1)
        self.add(dsk2)
        self.add(flng2)
        self.add(dsk3)
        self.add(flng3)
        self.add(dsk4)


class Spring(ComplexPrimitive):

    @classmethod
    def Spring1(cls, v1, v2, diameter, slices=12, stacks=5, batch=None,
               group=None):
        v1 = Vector3(*v1)
        v2 = Vector3(*v2)
        v = v2 - v1
        length = v.magnitude()
        spring = Spring(diameter, length, slices=slices, stacks=stacks,
                        batch=batch, group=group)
        x_ang, y_ang = rotz2vec(v)
        spring.position = (v1.x, v1.y, v1.z)
        spring.rotation = (1, 0, 0, x_ang)
        spring.rotation = (0, 1, 0, y_ang)

        return spring

    def __init__(self, diameter, length, slices=21, stacks=1,
                 position=(0, 0, 0), anchor=(0, 0, 0), mode=GL_TRIANGLES,
                 batch=None, group=None):
        super(Spring, self).__init__(position, anchor, mode, batch, group)

        dsk1 = Disk(diameter/4, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk1.rotation = (1, 0, 0, 180)

        flng1 = Cylinder(diameter/4, length/5, slices=slices,
                         stacks=stacks, mode=mode, batch=batch,
                         group=group)

        dsk2 = Disk(1.5*diameter, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk2.z = length/5

        flng2 = Cylinder(1.5*diameter, length, slices=slices,
                         stacks=stacks, mode=mode, batch=batch,
                         group=group)
        flng2.z = length/5

        dsk3 = Disk(1.5*diameter, slices=slices, mode=mode, batch=batch,
                    group=group)
        dsk3.z = 6*length/5


        flng3 = Cylinder(diameter/4, length/5, slices=slices,
                         stacks=stacks, mode=mode, batch=batch,
                         group=group)
        flng3.z = 6*length/5

        box = Box(2*diameter, 2*diameter, 0.1*diameter, batch=batch,
                  group=group)
        box.z = 7*length/5

        self.add(dsk1)
        self.add(flng1)
        self.add(dsk2)
        self.add(flng2)
        self.add(dsk3)
        self.add(flng3)
        self.add(box)


class Snubber(ComplexPrimitive):

    pass
