'''The main scene to be rendered'''

from __future__ import division
from math import sin, cos, tan, radians, fabs, pi
from itertools import chain
import ctypes
from math import degrees

import numpy as np
from pyglet.gl import *
from pyglet.gl.glu import *
from pyglet.graphics import Group
from pyglet.graphics import Batch
from pyglet.text import Label
# from OpenGL.GLUT import *

# from psi.settings import options

from psi.gui.pyglettk import PygletTkFrame
from psi.gui.primitives import (Cylinder, Line, Point, Valve, Arrow,
                                Reducer, Flange, Box, Spring)
from psi.gui.primitives import ComplexPrimitive
from psi.gui.bounds import AABB
from psi.gui.bvh import BVHTree
from psi.gui.text import Label as PointLabel
from psi.gui.utils import grouper

from psi.utils.euclid import Vector3


# ctypes arrays of floats
def vec(*args):
    return (GLfloat * len(args))(*args)


class ProgLabelGroup(Group):

    def __init__(self, viewport, parent=None):
        super(ProgLabelGroup, self).__init__(parent)
        self.viewport = viewport

    def set_state(self):
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()

        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, self.viewport.width, 0, self.viewport.height, -1, 1)

    def unset_state(self):
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()


class ProgLabel(object):

    def __init__(self, viewport, font_name="Times New Roman", font_size=28,
                 bold=True, color=(255, 255, 255, 255), is_visible=True):
        self.viewport = viewport
        self._is_visible = is_visible
        w, h = viewport.width, viewport.height
        self.prog_name = Label("OpenPIPE", x=w-20, y=h-10, font_size=font_size,
                               bold=bold, color=color, font_name=font_name,
                               group=ProgLabelGroup(viewport),
                               anchor_x="right", anchor_y="top",
                               )

    @property
    def is_visible(self):
        return self._is_visible

    @is_visible.setter
    def is_visible(self, value):
        self._is_visible = value

    def draw(self):
        if self._is_visible:
            self.prog_name.x = self.viewport.width - 20
            self.prog_name.y = self.viewport.height - 10
            self.prog_name.draw()


class BackGroundGroup(Group):

    def __init__(self, parent=None):
        super(BackGroundGroup, self).__init__(parent)

    def set_state(self):
        glPushAttrib(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)
        glDepthMask(GL_FALSE)

        glClear(GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(-1, 1, -1, 1, -1, 1)

    def unset_state(self):
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopAttrib()


class BackGround(object):

    def __init__(self, top=(135, 206, 235), bottom=(209, 239, 255),
                    group=BackGroundGroup()):
        self._id = None
        self._top = top
        self._bottom = bottom
        self._group = group
        self._update()

    def _update(self):
        if self._id:
            glDeleteLists(self._id, 1)

        self._id = glGenLists(1)
        glNewList(self._id, GL_COMPILE)
        glMatrixMode(GL_MODELVIEW)
        glBegin(GL_QUADS)
        glColor3ub(*self._top)
        glVertex2f(1.0, 1.0)
        glVertex2f(-1.0, 1.0)
        glColor3ub(*self._bottom)
        glVertex2f(-1.0, -1.0)
        glVertex2f(1.0, -1.0)
        glEnd()
        glEndList()

    def _set_top(self, rgb):
        self._top = rgb
        self._update()

    top = property(lambda self: self._top, _set_top)

    def _set_bottom(self, rgb):
        self._bottom = rgb
        self._update()

    bottom = property(lambda self: self._bottom, _set_bottom)

    def draw(self):
        self._group.set_state()
        glCallList(self._id)
        self._group.unset_state()


class TriadGroup(Group):

    def __init__(self, camera, position, parent=None):
        super(TriadGroup, self).__init__(parent)
        self.camera = camera
        self.position = position
        self.transform = np.identity(4, dtype=np.float32)
        self.viewport = (GLint * 4)()

    def set_state(self):
        glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT)
        #glDisable(GL_LIGHTING)

        glGetIntegerv(GL_VIEWPORT, self.viewport)
        _, _, w, h = self.viewport
        fac = 0.15
        m = int(fac * min(w, h) * 1.25)
        if self.position == 1:
            x, y = 0, h-m
        elif self.position == 2:
            x, y = w-m, h-m
        elif self.position == 3:
            x, y = (w//2)-(m//2), (h//2)-(m//2)
        elif self.position == 4:
            x, y = 0, 0
        elif self.position == 5:
            x, y = w-m, 0
        glViewport(x, y, m, m)
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(-fac, fac, -fac, fac, -fac, fac)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        self.transform[:3,:3] = self.camera.transform[:3,:3]
        glMultMatrixf(vec(*np.ravel(self.transform)))

    def unset_state(self):
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        _, _, w, h = self.viewport
        glViewport(0, 0, w, h)
        glPopAttrib()


class Triad(object):

    def __init__(self, camera, position=4, visible=True):
        self._visible = visible

        self.id = glGenLists(1)
        glNewList(self.id, GL_COMPILE)
        glColor3ub(0, 0, 255)                           # Z axis in blue
        cyl = gluNewQuadric()
        gluCylinder (cyl, 0.005, 0.005, 0.1, 12, 1)       # Body of axis
        glColor3ub(0, 0, 255)                       # Make arrow head white
        glPushMatrix()
        glTranslatef(0.0, 0.0, 0.1)
        gluDisk(cyl, 0.01, 0.02, 12, 1)
        gluCylinder(cyl, 0.02, 0.001, 0.05, 12, 1)      # Cone at end of axis
        glPopMatrix()

        glColor3ub(0, 255, 0)                           # Y axis is green
        glPushMatrix()
        glRotatef(-90, 1, 0, 0)
        cyl = gluNewQuadric()
        gluCylinder(cyl, 0.005, 0.005, 0.1, 12, 1)        # Body of axis
        glColor3ub(0, 255, 0)                       # Make arrow head white
        glPushMatrix()
        glTranslatef(0.0, 0.0, 0.1)
        gluDisk(cyl, 0.01, 0.02, 12, 1)
        gluCylinder(cyl, 0.02, 0.001, 0.05, 12, 1)      # Cone at end of axis
        glPopMatrix()
        glPopMatrix()

        glColor3ub(255, 0, 0)                           # X axis is red.
        glPushMatrix()
        glRotatef(90, 0, 1, 0)
        cyl = gluNewQuadric()
        gluCylinder(cyl, 0.005, 0.005, 0.1, 12, 1)        # Body of axis
        glColor3ub(255, 0, 0)                       # Make arrow head white
        glPushMatrix()
        glTranslatef(0.0, 0.0, 0.1)
        gluDisk(cyl, 0.01, 0.02, 12, 1)
        gluCylinder(cyl, 0.02, 0.001, 0.05, 12, 1)      # Cone at end of axis
        glPopMatrix()
        glPopMatrix()
        glEndList()

        self._group = TriadGroup(camera, position, None)

    def draw(self):
        if self.is_visible:
            self._group.set_state()
            glCallList(self.id)
            self._group.unset_state()

    def delete(self):
        glDeleteLists(self.id, 1)

    def _get_position(self):
        return self._group.position

    def _set_position(self, value):
        self._group.position = value

    position = property(_get_position, _set_position)

    def _is_visible(self, value):
        self._visible = value

    is_visible = property(lambda self: self._visible, _is_visible)


class Light(object):

    def __init__(self, id=GL_LIGHT0, ambient=(0, 0, 0), diffuse=(1, 1, 1),
                    specular=(1, 1, 1), position=(10, 10, 10)):
        self._id = id
        self._ambient = ambient
        self._diffuse = diffuse
        self._specular = specular
        self._position = position

    @property
    def id(self):
        return self._id

    def _set_ambient(self, ax, ay, az):
        self._ambient = (ax, ay, az, 1)

    ambient = property(lambda self: self._ambient,
                       lambda self, t: self._set_ambient(*t))

    def _set_diffuse(self, dx, dy, dz):
        self._diffuse = (dx, dy, dz, 1)

    diffuse = property(lambda self: self._diffuse,
                       lambda self, t: self._set_diffuse(*t))

    def _set_specular(self, sx, sy, sz):
        self._specular = (sx, sy, sz, 1)

    specular = property(lambda self: self._specular,
                        lambda self, t: self._set_specular(*t))

    def _set_position(self, x, y, z):
        self._position = (x, y, z, 1)

    position = property(lambda self: self._position,
                        lambda self, t: self._set_position(*t))

    def draw(self):
        glEnable(self.id)
        glLightfv(self.id, GL_AMBIENT, vec(*self.ambient))
        glLightfv(self.id, GL_DIFFUSE, vec(*self.diffuse))
        glLightfv(self.id, GL_SPECULAR, vec(*self.specular))
        glLightfv(self.id, GL_POSITION, vec(*self.position))


class Camera(object):
    '''Implementation of a camera class capable of both perspective and
    orthographic viewing.
    '''

    def __init__(self, scene, znear=1, zfar=1000, fovy=45, is_persp=True):
        self._scene = scene
        self._znear = znear
        self._zfar = zfar
        self._fovy = fovy
        self._is_persp = is_persp
        self._transform = np.eye(4, dtype=np.float32)   # row major modeling transform
        self._projection = np.eye(4, dtype=np.float32)
        self.orbit_point = np.array([0, 0, 0], dtype=np.float32)   # in world space

    @property
    def scene(self):
        return self._scene

    @property
    def znear(self):
        return self._znear

    @znear.setter
    def znear(self, value):
        self._znear = value

    @property
    def zfar(self):
        return self._zfar

    @zfar.setter
    def zfar(self, value):
        self._zfar = value

    @property
    def fovy(self):
        return self._fovy

    @fovy.setter
    def fovy(self, value):
        self._fovy = value

    @property
    def is_persp(self):
        return self._is_persp

    @is_persp.setter
    def is_persp(self, value):
        self._is_persp = value

    @property
    def projection(self):
        return self._projection[:]

    @property
    def gl_projection(self):
        return vec(*self._projection.reshape(-1, order='F'))

    @projection.setter
    def projection(self, znear, zfar, fovy, is_persp):
        self._znear = znear
        self._zfar = zfar
        self._fovy = fovy
        self._is_persp = is_persp

    def general_projection(self):
        """Used for perfect panning"""
        viewport = self._scene.viewport

        half_size = 1   # pick a value
        self.half_width = self.half_height = half_size
        self.tsx = self.tsy = 0   # no stereo
        self.iez = 0    # assume ortho
        if self.is_persp:
            self.iez = tan(0.5*self.fovy) / half_size

        w, h = viewport.size
        if w > h:
            self.half_width *= w/h
        else:
            self.half_height *= h/w

        hh = self.half_height
        hw = self.half_width
        zn = self.znear
        zf = self.zfar
        iez = self.iez
        tsx = self.tsx
        tsy = self.tsy

        if self.is_persp:
            # iez = 1 / zn => ez = 1/iez
            ez = 1 / iez
            self._projection[:, :] = [[ez/hw, 0, 0, 0],
                                      [0, ez/hh, 0, 0],
                                      [-ez*tsx/hw, -ez*tsy/hh, -(2*ez-(zn+zf))/(zn-zf), -1],
                                      [0, 0, -2*(ez*(ez-(zn+zf))+zn*zf)/(zn-zf), 0],
                                     ]
            # self.translate(ez*tsx, ez*tsy, ez)
        else:
            # iez = 0
            self._projection[:, :] = [[1/hw, 0, 0, 0],
                                      [0, 1/hh, 0, 0],
                                      [-tsx/hw, -tsy/hh, -2/(zn-zf), 0],
                                      [0, 0, (zn+zf)/(zn-zf), 1],
                                     ]

    def update_projection(self, aspect):
        fovy = self._fovy
        n = self._znear
        f = self._zfar
        ar = aspect

        if self._is_persp:
            t1 = 1/tan(radians(fovy*0.5))
            self._projection[:] = [[t1/ar, 0, 0, 0],
                                   [0, t1, 0, 0],
                                   [0, 0, (f+n)/(n-f), (2*f*n)/(n-f)],
                                   [0, 0, -1, 0]]
        else:
            eye_to_orbit_point = self.orbit_point - Vector3(*self.position)
            tan_fovy_2 = tan(radians(fovy*0.5))

            dz = eye_to_orbit_point.magnitude()

            # dz = abs(eye_to_orbit_point.z)

            if dz <= 0:
                dz = 1

            if ar < 1:
                # w = tan_fovy_2 * n * dz
                # h = w * (1/ar)
                h = tan_fovy_2 * n * dz
                w = h * ar
            else:
                w = tan_fovy_2 * n * dz
                h = w * (1/ar)
                # h = tan_fovy_2 * n * dz
                # w = h * ar

            l, r, b, t = -w, w, -h, h

            self._projection[:] = [[2/(r-l), 0, 0, -(r+l)/(r-l)],
                                   [0, 2/(t-b), 0, -(t+b)/(t-b)],
                                   [0, 0, -2/(f-n), -(f+n)/(f-n)],
                                   [0, 0, 0, 1]]

    @property
    def transform(self):
        '''Returns the camera viewing transform'''
        return np.linalg.inv(self._transform)

    def get_inv_proj_tf(self):
        """Get the inverse projection matrix for selection"""
        return np.linalg.inv(self._projection)

    def get_inv_tf(self):
        """Same as tranform"""
        return np.linalg.inv(self._transform)

    @property
    def gl_transform(self):
        '''Returns the modeling transform that OpenGL requires in linear column
        major format.
        '''
        return vec(*self._transform.reshape(-1, order='F'))

    def translate(self, x, y, z):
        '''Translating the camera effectively moves the world in the
        inverse direction.  Note that the homogeneous transform matrix
        is in row major form as in many reference texts. Also, the camera
        transform is pre-multiplied as its a modeling transform.
        '''
        transform = np.array([[1, 0, 0, -x],
                              [0, 1, 0, -y],
                              [0, 0, 1, -z],
                              [0, 0, 0, 1]], dtype=np.float32)
        self._transform = np.dot(transform, self._transform)

    def rotatez(self, ang):
        '''Rotate the camera about it's z-axis by a user defined angle.'''
        c = cos(radians(-ang))
        s = sin(radians(-ang))
        transform = np.array([[c, -s, 0, 0],
                              [s, c, 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]], dtype=np.float32)
        self._transform = np.dot(transform, self._transform)

    def rotatey(self, ang):
        '''Rotate the camera about it's y-axis by a user defined angle.'''
        c = cos(radians(-ang))
        s = sin(radians(-ang))
        transform = np.array([[c, 0, s, 0],
                              [0, 1, 0, 0],
                              [-s, 0, c, 0],
                              [0, 0, 0, 1]], dtype=np.float32)
        self._transform = np.dot(transform, self._transform)

    def rotatex(self, ang):
        '''Rotate the camera about it's x-axis by a user defined angle.'''
        c = cos(radians(-ang))
        s = sin(radians(-ang))
        transform = np.array([[1, 0, 0, 0],
                              [0, c, -s, 0],
                              [0, s, c, 0],
                              [0, 0, 0, 1]], dtype=np.float32)
        self._transform = np.dot(transform, self._transform)

    def axisangle(self, x, y, z, ang):
        '''Rotate the camera about an arbitrary axis by a user defined
        angle.
        '''
        c = cos(radians(-ang))
        s = sin(radians(-ang))
        t = 1 - c
        transform = np.array([[t*x*x+c, t*x*y-z*s, t*x*z+y*s, 0],
                              [t*x*y+z*s, t*y*y+c, t*y*z-x*s, 0],
                              [t*x*z-y*s, t*y*z+x*s, t*z*z+c, 0],
                              [0, 0, 0, 1]], dtype=np.float32)
        self._transform = np.dot(transform, self._transform)

    def _get_position(self):
        '''Return the global position of the camera.'''
        inv_tf = np.linalg.inv(self._transform)
        return inv_tf[:3, 3]

    def _set_position(self, x, y, z):
        '''Set the global position of the camera.  Note that the negatives
        sign defines the viewing transform of the camera.
        '''
        new_pos = np.array([x, y, z])
        # get world position of camera
        inv_tf = np.linalg.inv(self._transform)
        cam_pos = inv_tf[:3, 3]
        new_cam_pos = cam_pos - new_pos
        delta_pos = np.dot(inv_tf, np.append(new_cam_pos, 1.0))
        dx, dy, dz = delta_pos[:3]
        self.translate(-dx, -dy, -dz)

    position = property(_get_position,
                        lambda self, args: self._set_position(*args))

    def lookat(self, eyex, eyey, eyez, ctrx, ctry, ctrz):
        '''Position the camera such that it is positioned towards a
        target point in world space.  The camera origin is defined as
        being the eye point in world space.
        '''
        ctr = np.array([ctrx, ctry, ctrz], dtype=np.float32)
        eye = np.array([eyex, eyey, eyez], dtype=np.float32)

        z = ctr - eye
        zunit = z / np.linalg.norm(z)
        # in case 'lookat' vector is also 'up' vector
        if (fabs(zunit[0]) < self.model.settings.epsilon
                and fabs(zunit[2]) < self.model.settings.epsilon):
            if zunit[1] > 0:
                up = np.array([0, 0, -1], dtype=np.float32)
            else:
                up = np.array([0, 0, 1], dtype=np.float32)
        else:
            up = np.array([0, 1, 0], dtype=np.float32)
        x = np.cross(up, zunit)
        xunit = x / np.linalg.norm(x)
        yunit = np.cross(xunit, zunit)  # make sure 'up' is orthogonal

        xunit = xunit.reshape((3, 1))   # make into column arrays
        yunit = yunit.reshape((3, 1))
        zunit = zunit.reshape((3, 1))
        eye = -eye.reshape((3, 1))      # convert to viewing transform
        self._transform[:3, :4] = np.hstack((xunit, yunit, zunit, eye))

    def clear_transform(self):
        self._transform[:] = np.eye(4, dtype=np.float32)

    def distance_from_selection(self):
        pass

    def front_view(self):
        pass

    def top_view(self):
        pass

    def bottom_view(self):
        pass

    def left_view(self):
        pass

    def right_view(self):
        pass

    def back_view(self):
        pass

    def fit(self):
        pass


class Ray(object):

    def __init__(self):
        self.ori = Vector3(0, 0, 0)
        self.dir = Vector3(0, 0, 0)

    def hits_prim(self, prim):
        """Since a primitive may be a ComplexPrimitive ie. consist of other
        primitives or even other ComplexPrimitives, a recursive test must be
        performed to determine whether the ray will interect it.

        For a simple primitive, the primitive bound is first checked.  If the
        bound test passes, then each individual triangle, line or point is
        checked to find a hit.

        If the primitive is complex in nature, a similar test is performed as
        in the case of a simple primitive, however, this test is performed for
        all the constituent primitives.  As soon as a hit is found, the testing
        stops and the complex primtive is deemed to be intersected in its
        entirity.
        """
        if isinstance(prim, ComplexPrimitive):  # bottom-top
            for prim in prim.primitives:
                if self.hits_bound(prim):
                    return self.hits_prim(prim)

        if prim.mode == GL_TRIANGLES:
            if self.hits_bound(prim):
                return self._hits_triangles(prim)
        elif prim.mode == GL_LINES:
            return self._hits_lines(prim)
        elif prim.mode == GL_POINTS:
            return self._hits_vertices(prim)

        return False    # sentinal

    def _hits_triangles(self, prim):
        """Iterate all the triangles of a primitive and return true if there
        is a hit.  When the first instance of ray intersection a triangle is
        found, a boolean True value is returned (ie. it is not necessary to
        iterate over all the triangles.

        http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
        """
        # TODO: Refactor without making unnecessary function call
        # TODO: See if slicing operations can be reduced to one slice

        # get the full list of vertices from the domain by typecasting bytes array
        start = prim.vertex_list.start
        count = prim.vertex_list.count
        domain = prim.vertex_list.get_domain()
        attribute = domain.attribute_names['vertices']
        buffer = attribute.buffer
        ptr_type = ctypes.POINTER(attribute.c_type * buffer.size)
        vertices = ctypes.cast(buffer.data_ptr, ptr_type).contents

        try:
            indices = prim.vertex_list.indices

            for i, j, k in grouper(indices, 3):
                v0 = Vector3(*vertices[i*3:i*3+3])
                v1 = Vector3(*vertices[j*3:j*3+3])
                v2 = Vector3(*vertices[k*3:k*3+3])

                e1 = v1 - v0
                e2 = v2 - v0
                h = self.dir.cross(e2)
                a = e1.dot(h)
                if a > -0.00001 and a < 0.00001:
                    continue    # False so move to next iteration

                f = 1 / a
                s = self.ori - v0
                u = f * s.dot(h)
                if u < 0.0 or u > 1.0:
                    continue

                q = s.cross(e1)
                v = f * self.dir.dot(q)
                if v < 0.0 or u + v > 1.0:
                    continue

                t = f * e2.dot(q)
                if t > 0.00001:
                    return True # exit search - hit found
                else:
                    continue

            return False
        except AttributeError as e:   # not indexed
            for v0, v1, v2 in grouper(vertices, 3):
                e1 = v1 - v0
                e2 = v2 - v0
                h = self.dir.cross(e2)
                a = e1.dot(h)
                if a > -0.00001 and a < 0.00001:
                    continue    # False so move to next iteration

                f = 1 / a
                s = self.ori - v0
                u = f * s.dot(h)
                if u < 0.0 or u > 1.0:
                    continue

                q = s.cross(e1)
                v = f * self.dir.dot(q)
                if v < 0.0 or u + v > 1.0:
                    continue

                t = f * e2.dot(q)
                if t > 0.00001:
                    return True
                else:
                    continue

            return False

    def _hits_lines(self, prim):
        """Ray line intersection testing"""
        tol = 0.05  # tolerance

        ro = self.ori
        rd = self.dir
        u = rd.normalized()
        ev1 = Vector3(*prim.vertex_list.vertices[:3])
        ev2 = Vector3(*prim.vertex_list.vertices[3:6])
        v = ev2 - ev1
        w = ro - ev1

        a = u.dot(u)
        b = u.dot(v)
        c = v.dot(v)
        d = u.dot(w)
        e = v.dot(w)

        D = a*c - b*b
        sD = D
        tD = D

        if D < tol:
            sN = 0.0
            sD = 1.0
            tN = e
            tD = c
        else:
            sN = (b*e - c*d)
            tN = (a*e - b*d)
            if sN < 0.0:
                sN = 0.0
                tN = e
                tD = c

        if tN < 0.0:
            tN = 0.0
            if -d < 0.0:
                sN = 0.0
            else:
                sN = -d
                sD = a
        elif (tN > tD):
            tN = tD
            if (-d + b) < 0.0:
                sN = 0
            else:
                sN = (-d + b)
                sD = a

        if abs(sN) < tol:
            sc = 0.0
        else:
            sc = sN / sD

        if abs(tN) < tol:
            tc = 0.0
        else:
            tc = tN / tD

        dP = w + (sc * u) - (tc * v)

        if dP.magnitude() < tol:
            return True

        return False

    def _hits_vertices(self, prim):
        """Ray-vertex intersection"""
        tol = 0.35  # tolerance

        ro = self.ori
        rd = self.dir
        rp = Vector3(*prim.vertex_list.vertices)
        rdu = rd.normalized()
        rpo = rp - ro
        if rpo.dot(rdu) < 0.0:  # point behind ray
            return False
        else:
            if degrees(rdu.angle(rpo)) < tol:
                return True
            else:
                return False
        return False

    def hits_bound(self, prim):
        minb = Vector3(*prim.bound.min)
        maxb = Vector3(*prim.bound.max)

        dir_x = 1.0 / self.dir.x
        dir_y = 1.0 / self.dir.y
        dir_z = 1.0 / self.dir.z

        t1 = (minb.x - self.ori.x) * dir_x
        t2 = (maxb.x - self.ori.x) * dir_x
        t3 = (minb.y - self.ori.y) * dir_y
        t4 = (maxb.y - self.ori.y) * dir_y
        t5 = (minb.z - self.ori.z) * dir_z
        t6 = (maxb.z - self.ori.z) * dir_z

        tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6))
        tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6))

        if tmax < 0:
            t = tmax
            return False

        if tmin > tmax:
            t = tmax
            return False

        t = tmin

        return True


class Scene(object):

    def __init__(self, viewport, model):
        self._viewport = viewport
        self._model = model

        # accelerator data structure
        self.bvhtree = BVHTree()

        # scene objects
        self.light = Light()
        self.camera = Camera(self)
        self.background = BackGround(
                            top=self.model.settings.background_top_color,
                            bottom=self.model.settings.background_bottom_color
                            )
        self.triad = Triad(self.camera)
        self.prog_label = ProgLabel(viewport)

        # selection ray and bound
        self.ray = Ray()
        self.hitlist = set()
        self.picked_bound = AABB()
        self.model_bound = AABB()

        # primitive mapping
        self.pipes = {}         # pipe -> (line, pt1, pt2)
        self.bends = {}         # bends -> (curve, pt1, pt20
        self.reducers = {}      # reducers -> (line, pt1, pt2)
        self.valves = {}        # valves -> (line, pt1, pt2)
        self.flanges = {}       # flange -> (line, pt1, pt2)
        self.anchors = {}       # anchor -> pt
        self.supports = {}      # support -> pt
        self.springs = {}       # spring -> pt
        self.beams = {}         # beam -> (line, pt1, pt2)

        # centerline model
        self.lines = {}
        self.points = {}

        # batches
        self.face_batch = Batch()
        self.line_batch = Batch()
        self.point_batch = Batch()
        self.point_label_batch = Batch()

        # drawing flags
        self.draw_mode = 1      # 1=default, 2=line, 3=transparent
        self.draw_logo = True
        self.draw_point_labels = True
        self.draw_triad = True

    @property
    def viewport(self):
        return self._viewport

    @property
    def model(self):
        return self._model

    @property
    def mouse(self):
        return self._viewport.mouse

    @property
    def keyboard(self):
        return self._viewport.keyboard

    def init_gl(self):
        '''One-time GL setup parameters'''
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LEQUAL)

        glShadeModel(GL_SMOOTH)
        glEnable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, vec(0.2, 0.2, 0.2, 1.0))

        # material parameters
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, vec(1, 1, 1, 1))
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, vec(0, 0, 0, 1))
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 75)

    def create_objects(self):
        '''create opengl objects when opengl is initialized'''
        pass

    def update_object_resize(self, width, height):
        '''Called when the client window is resized'''
        self.width, self.height = width, height
        self.aspect = width / height
        self.camera.update_projection(self.aspect)

    @property
    def primitives(self):
        return chain(self.pipes, self.bends, self.reducers, self.valves,
                     self.flanges, self.anchors, self.supports, self.springs,
                     self.beams, self.lines, self.points)

    def select(self, norm_x, norm_y):
        """Select primitives based on normalized mouse xpos and ypos"""
        inv_proj = self.camera.get_inv_proj_tf()
        dv = np.dot(inv_proj, np.array([norm_x, norm_y, -1, 1]))
        nx, ny, nz = dv[:3]

        if self.camera.is_persp:
            ray_origin = np.array([0, 0, 0, 1])
            ray_direction = np.array([nx, ny, nz, 0])
        else:
            # for ortho selection the origin of the ray is varient
            ray_origin = np.array([nx, ny, nz, 1])
            ray_direction = np.array([0, 0, -1, 0])

        # invert camera modelview matrix - viewing to camera transform
        mv_inv = self.camera.get_inv_tf()

        # multiply to get ray in world coordinates
        ray_origin = np.dot(mv_inv, ray_origin)
        ray_direction = np.dot(mv_inv, ray_direction)

        self.ray.ori.x, self.ray.ori.y, self.ray.ori.z = ray_origin[:3]
        self.ray.dir.x, self.ray.dir.y, self.ray.dir.z = ray_direction[:3]

        for prim in self.bvhtree.get_picked_set(self.ray):
            if prim not in self.hitlist and prim.is_visible:
                self.hitlist.add(prim)
                self.picked_bound.merge(prim.bound)
                prim.is_picked = True

        if self.hitlist:
            self.camera.orbit_point = Vector3(*self.picked_bound.center())
        else:
            # self.bvhtree.get_model_bound(self.model_bound)
            self.model_bound = self.bvhtree.root.bound
            self.camera.orbit_point = Vector3(*self.model_bound.center())

        self._viewport.draw_objects()   # highlight or unhighlight

    def clear_hits(self):
        for prim in self.hitlist:
            prim.is_picked = False
        self.hitlist.clear()
        self.picked_bound.clear()

    def build(self):
        """Build the bounding volume heirarchy tree"""
        proxies = [prim for prim in self.primitives]
        self.bvhtree = self.bvhtree.build(proxies, mode=1)
        self.camera.orbit_point = Vector3(*self.bvhtree.root.bound.center())

    def add_point(self, pt, font_size=8, color=(0, 0, 0), rebuild=True):
        p1 = Point(pt, batch=self.point_batch)
        lbl1 = PointLabel("", pt=pt, font_size=8, color=(0, 0, 0),
                          batch=self.point_label_batch)

        self.points[p1] = lbl1

        if rebuild:
            self.build()

        return p1

    def add_pipe(self, pt1, pt2, od, rebuild=True):
        """Add a primitive to the scene. For each primitive add the points,
        line and faces. Draw point and lines together and point and faces
        together.
        """
        # pipe from point
        pos1 = pt1.position
        pos2 = pt2.position

        pipe = Cylinder.Cylinder1(pos1, pos2, od, batch=self.face_batch)
        pipe.color = self.model.settings.pipe_color

        line = Line(pos1, pos2, batch=self.line_batch)

        # create relations
        self.pipes[pipe] = (line, pt1, pt2)
        self.lines[line] = pipe

        # add to tree
        if rebuild:
            self.build()

        return pipe

    def add_valve(self, pt1, pt2, od, rebuild=True):
        pos1 = pt1.position
        pos2 = pt2.position

        valve = Valve.Valve1(pos1, pos2, od, batch=self.face_batch)
        valve.color = self.model.settings.valve_color

        line = Line(pos1, pos2, batch=self.line_batch)

        # create relations for later querying
        self.valves[valve] = (line, pt1, pt2)
        self.lines[line] = valve

        if rebuild:
            self.build()

        return valve

    def add_reducer(self, pt1, pt2, od, od1, rebuild=True):
        pos1 = pt1.position
        pos2 = pt2.position

        red = Reducer.Reducer1(pos1, pos2, od, od1, batch=self.face_batch)
        red.color = self.model.settings.reducer_color

        line = Line(pos1, pos2, batch=self.line_batch)

        # create relations for later querying
        self.reducers[red] = (line, pt1, pt2)
        self.lines[line] = red

        if rebuild:
            self.build()

        return red

    def add_flange(self, pt1, pt2, od, rebuild=True):
        pos1 = pt1.position
        pos2 = pt2.position

        flng = Flange.Flange1(pos1, pos2, od, batch=self.face_batch)
        flng.color = self.model.settings.flange_color

        line = Line(pos1, pos2, batch=self.line_batch)

        # create relations for later querying
        self.flanges[flng] = (line, pt1, pt2)
        self.lines[line] = flng

        if rebuild:
            self.build()

        return flng

    def add_anchor(self, pt, od, dircos=(1, 0, 0), size=1.0, rebuild=True):
        """Add an anchor support to the scene"""
        pt1 = pt.position

        pt2 = Vector3(*dircos).normalized()[:]

        anchor = Box.Box1(pt1, pt2, od*1.25, batch=self.face_batch)
        anchor.color = self.model.settings.anchor_color

        # create relations
        self.supports[anchor] = (pt,)

        if rebuild:
            self.build()

        return anchor

    def add_support(self, pt, od, dircos=(1, 0, 0), axial=False, size=1.0,
                    rebuild=True):
        """Add a rigid support given the centerpoint, pipe od and the direction
        cosine.
        """
        unit = Vector3(*dircos).normalized()
        if axial:
            pt2 = Vector3(*pt.position)
            pt1 = pt2 + 1.25 * od * unit
        else:
            pt2 = Vector3(*pt.position) + (od/2) * unit
            pt1 = pt2 + 1.25 * od * unit

        top_supp = Arrow.Arrow1((pt1)[:], (pt2)[:], 0.15*od, batch=self.face_batch)
        top_supp.color = self.model.settings.support_color

        unit = -Vector3(*dircos).normalized()
        if axial:
            pt2 = Vector3(*pt.position)
            pt1 = pt2 + 1.25 * od * unit
        else:
            pt2 = Vector3(*pt.position) + (od/2) * unit
            pt1 = pt2 + 1.25 * od * unit

        bot_supp = Arrow.Arrow1(pt1[:], pt2[:], 0.15*od, batch=self.face_batch)
        bot_supp.color = self.model.settings.support_color


        self.supports[top_supp] = (pt,)
        self.supports[bot_supp] = (pt,)

        if rebuild:
            self.build()

        return top_supp, bot_supp


    def add_spring(self, pt, od, dircos=(0, 1, 0), base_type=False, size=1.0,
                   rebuild=True):
        unit = Vector3(*dircos).normalized()

        pt2 = Vector3(*pt.position) + (od/2) * unit
        pt1 = pt2 + 1.25 * od * unit

        top_supp = Spring.Spring1((pt2)[:], (pt1)[:], 0.25*od,
                                  batch=self.face_batch)
        top_supp.color = self.model.settings.spring_color

        self.supports[top_supp] = (pt,)

        if rebuild:
            self.build()

        return top_supp

    def add_snubber(self, pt, od, dircos=(0, 1, 0), size=1.0, rebuild=True):
        pass

    def draw_objects(self):
        """Render all objects in the scene"""
        # gather mouse input first - this resolves the ortho related panning
        # issue
        if self.mouse.tx or self.mouse.ty or self.mouse.tz:
            # perfect pan and zoom to be added
            self.camera.translate(-self.mouse.tx, self.mouse.ty, self.mouse.tz)
        elif self.mouse.rx or self.mouse.ry:
            camtf = self.camera.transform
            opwrtcp = self.camera.orbit_point - Vector3(*camtf[:3, 3])
            new_pos = np.dot(np.append(opwrtcp, 1), camtf)
            x, y, z = new_pos[:3]
            self.camera.translate(x, y, z)  # move cam to orbit points
            self.camera.rotatex(-self.mouse.rx)
            self.camera.rotatey(-self.mouse.ry)
            # move cam back with side effect from rotation above
            self.camera.translate(-x, -y, -z)

        # call before draw - avoids buggy draws
        self.mouse.clear()

        # start OGL drawing stuff
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # viewport transform
        glViewport(0, 0, self.width, self.height)

        # projection transform - gl accepts column major format
        glMatrixMode(GL_PROJECTION)
        # update required every frame for ortho zooming
        glLoadIdentity()
        self.camera.update_projection(self.aspect)
        # self.camera.general_projection()

        glLoadMatrixf(self.camera.gl_projection)

        # start rendering
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.light.draw()
        self.background.draw()

        if self.draw_triad:
            self.triad.draw()

        glLoadMatrixf(self.camera.gl_transform)

        # draw stuff
        if self.draw_mode == 1:
            self.face_batch.draw()
        elif self.draw_mode == 2:
            self.line_batch.draw()

        self.point_batch.draw()

        if self.draw_logo:
            self.prog_label.draw()

        if self.draw_point_labels:
            self.point_label_batch.draw()

        if self.model.settings.debug:
            # draw bvhtree nodes
            glPushAttrib(GL_CURRENT_BIT | GL_LIGHTING_BIT)
            glDisable(GL_LIGHTING)
            glColor3ub(0, 0, 0)
            for node in self.bvhtree.iter_nodes():
                sx, sy, sz = node.bound.size()
                x, y, z = sx*0.5, sy*0.5, sz*0.5
                vertices = [x, y, z, -x, y, z, -x, -y, z, x, -y, z, x, y, -z,
                            -x, y, -z, -x, -y, -z, x, -y, -z]
                indices = [0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0,
                            4, 4, 5, 5, 1, 1, 0, 3, 7, 7, 6, 6, 2, 2, 3, 0, 3,
                            3, 7, 7, 4, 4, 0, 1, 2, 2, 6, 6, 5, 5, 1]
                glPushMatrix()
                cx, cy, cz = node.bound.center()
                glTranslatef(cx, cy, cz)
                pyglet.graphics.draw_indexed(len(vertices)//3, GL_LINES,
                                                indices, ("v3f", vertices))
                glPopMatrix()
            glPopAttrib()
