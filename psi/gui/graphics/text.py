"""Implementation of 2d text in a 3d environment"""

from __future__ import division
from pyglet.gl import *
from pyglet import text
from pyglet import graphics


class _MVPGroup(graphics.Group):
    """Queries the modelview, projection and viewport transforms once every
    draw.
    """

    def __init__(self, parent=None):
        super(_MVPGroup, self).__init__(parent)
        self.modelview = (GLdouble * 16)()
        self.projection = (GLdouble * 16)()
        self.viewport = (GLint * 4)()
        self.depthrange = (GLdouble * 2)()
        self.mvp = (GLdouble * 16)()

    def set_state(self):
        glGetDoublev(GL_MODELVIEW_MATRIX, self.modelview)
        glGetDoublev(GL_PROJECTION_MATRIX, self.projection)
        glGetIntegerv(GL_VIEWPORT, self.viewport)
        glGetDoublev(GL_DEPTH_RANGE, self.depthrange)

        # multiply to generate MVP matrix
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadMatrixd(self.projection)
        glMultMatrixd(self.modelview)
        glGetDoublev(GL_MODELVIEW_MATRIX, self.mvp)
        glPopMatrix()

_mvp_group = _MVPGroup()    # query once

class LabelGroup(graphics.Group):
    """Each label is transformed from 3d space to 2d screen space"""

    def __init__(self, parent=_mvp_group):
        super(LabelGroup, self).__init__(parent)
        self._x, self._y, self._z = 0, 0, 0

    def set_state(self):
        # TODO: Set far clipping plane for text

        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT)
        glDisable(GL_LIGHTING)
        #glDisable(GL_BLEND)
        glDisable(GL_DEPTH_TEST)

        (a00, a01, a02, a03,
         a10, a11, a12, a13,
         a20, a21, a22, a23,
         a30, a31, a32, a33) = self.parent.mvp[:16] # OpenGL column major
        a, b, width, height = self.parent.viewport
        near, far = self.parent.depthrange

        # perspective divide
        try:
            w = a03*self.x + a13*self.y + a23*self.z + a33
            xndc = (a00*self.x + a10*self.y + a20*self.z + a30)/w
            yndc = (a01*self.x + a11*self.y + a21*self.z + a31)/w
            zndc = (a02*self.x + a12*self.y + a22*self.z + a32)/w
        except ZeroDivisionError:   # TODO: ok? revisit to optimize
            xndc = (a00*self.x + a10*self.y + a20*self.z + a30)
            yndc = (a01*self.x + a11*self.y + a21*self.z + a31)
            zndc = (a02*self.x + a12*self.y + a22*self.z + a32)

        # screen coordinates
        xs = 0.5*width*xndc + (a + 0.5*width)
        ys = 0.5*height*yndc + (b + 0.5*height)
        zs = 0.5*zndc*(far-near) + 0.5*(far+near)

        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()

        glOrtho(a, width, b, height, -1, 1)

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glTranslated(xs, ys, zs)

    def unset_state(self):
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopAttrib()

    def _get_x(self):
        return self._x

    def _set_x(self, x):
        self._x = x

    x = property(_get_x, _set_x)

    def _get_y(self):
        return self._y

    def _set_y(self, y):
        self._y = y

    y = property(_get_y, _set_y)

    def _get_z(self):
        return self._z

    def _set_z(self, z):
        self._z = z

    z = property(_get_z, _set_z)


class Label(text.Label):

    def __init__(self, text='',
                 font_name=None, font_size=None, bold=False, italic=False,
                 color=(255, 255, 255),
                 pt=(0, 0, 0), width=None, height=None,
                 anchor_x='left', anchor_y='baseline',
                 align='left',
                 multiline=False, dpi=None, batch=None):
        self._group = LabelGroup()
        self.x, self.y, self.z = pt    # follows self._group assignment
        super(Label, self).__init__(text=text,
                font_name=font_name, font_size=font_size, bold=bold,
                italic=italic, color=color+(255,),
                x=0, y=0,  # note x and y of base class set to 0
                width=width, height=height,
                anchor_x=anchor_x, anchor_y=anchor_y,
                align=align,
                multiline=multiline, dpi=dpi, batch=batch,
                group=self._group)

    def _get_x(self):
        return self._x

    def _set_x(self, x):
        self._x = x
        self._group.x = x

    x = property(_get_x, _set_x)

    def _get_y(self):
        return self._y

    def _set_y(self, y):
        self._y = y
        self._group.y = y

    y = property(_get_y, _set_y)

    def _get_z(self):
        return self._z

    def _set_z(self, z):
        self._z = z
        self._group.z = z

    z = property(_get_z, _set_z)
