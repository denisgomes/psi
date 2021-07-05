import sys


if sys.version_info[0] < 3:
    import Tkinter as tk
    import Dialog as dialog
else:
    import tkinter as tk
    from tkinter import dialog as dialog

from pyopengltk import OpenGLFrame
import pyglet

from pyglet.text import Label

# this line is very important, we're tricking pyglet into thinking there is a
# context avalible but we can't make it work with the shadow window that allows
# sharing of object between contexts
pyglet.options['shadow_window'] = False

# now that that is set we can import gl and get on our way
from pyglet import gl


class PygletTkFrame(OpenGLFrame):

    def __init__(self, *args, **kw):
        # Set background to empty string to avoid
        # flickering overdraw by Tk
        kw['bg'] = ""
        tk.Frame.__init__(self, *args, **kw)
        self.bind('<Map>', self.tkMap)
        self.bind('<Configure>', self.tkResize)
        self.bind('<Expose>', self.tkExpose)
        self.gl_initialized = False

    def tkMap(self, evt):
        """" Called when frame goes onto the screen """
        self._wid = self.winfo_id()
        self.tkCreateContext()
        self.prepare_gl()

    def tkExpose(self, evt):
        self.prepare_gl()
        self._display()

    def prepare_gl(self):
        super(PygletTkFrame, self).tkMakeCurrent()

        if not self.gl_initialized:
            self.initgl()
            self.gl_initialized = True
            self.on_reshape()

        self.pygletcontext.set_current()

    def on_reshape(self):
        if self.width <= 0:
            self.width = 1
        if self.height <= 0:
            self.height = 1

        if self.gl_initialized:
            self.pygletcontext.set_current()

        # gl.glViewport(0, 0, self.width, self.height)
        # gl.glMatrixMode(gl.GL_PROJECTION)
        # gl.glLoadIdentity()
        # gl.glOrtho(0, self.width, 0, self.height, 1, -1)
        # gl.glMatrixMode(gl.GL_MODELVIEW)

        if self.gl_initialized:
            self.update_object_resize(self.width, self.height)

    def tkResize(self, evt):
        """
        Things to do on window resize:
        Adjust viewport:
            glViewPort(0,0, width, height)
        Adjust projection matrix:
            glFrustum(left * ratio, right * ratio, bottom, top, nearClip,farClip)
        or
            glOrtho(left * ratio, right * ratio, bottom, top, nearClip,farClip)
        or
            gluOrtho2D(left * ratio, right * ratio, bottom, top)
        (assuming that left, right, bottom and top are all equal and
         ratio=width/height)
        """
        self.width, self.height = evt.width, evt.height
        # tkResize -> tkMap -> tkExpose  - order of method call
        if self.winfo_ismapped():
            self.prepare_gl()
            self.tkMakeCurrent()

    @property
    def size(self):
        return self.width, self.height

    def _display(self):
        """on draw"""
        self.update_idletasks()
        self.tkMakeCurrent()
        self.draw_objects()
        self.tkSwapBuffers()

    def initgl(self):
        # normal gl init
        if pyglet.version > "1.1.4":
            self.pygletcontext = PygletTkContext()
        else:
            self.pygletcontext = gl.Context()
        self.pygletcontext.set_current()

        # gl.glEnable(gl.GL_BLEND)
        # gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        # gl.glEnable(gl.GL_TEXTURE_2D)
        # gl.glClearColor(0, 0, 0, 1)

        # create objects to draw
        self.create_objects()

    def redraw(self):
        # For the user code
        self.draw_objects()

    # to implement by sub class
    def create_objects(self):
        '''create opengl objects when opengl is initialized'''
        self.label = pyglet.text.Label('Hello, World Pyglet-Tkinter!',
                                       font_name='Times New Roman',
                                       font_size=12,
                                       color=(255, 225, 0, 255),
                                       italic=True,
                                       x=self.width//2, y=self.height//2,
                                       anchor_x='center', anchor_y='center')

    def update_object_resize(self, width, height):
        '''called when the window receives only if opengl is initialized'''
        pass

    def draw_objects(self):
        '''called in the middle of ondraw after the buffer has been cleared'''
        self.label.draw()


class PygletTkContext(gl.Context):

    def __init__(self, config=None, context_share=None):
        self.config = config
        self.context_share = context_share
        self.canvas = None

        if context_share:
            self.object_space = context_share.object_space
        else:
            self.object_space = gl.ObjectSpace()

    def attach(self, canvas=None):
        pass

    def detach(self):
        pass

    def set_current(self):
        # XXX not per-thread
        gl.current_context = self

        # XXX
        gl.gl_info.set_active_context()
        gl.glu_info.set_active_context()

        # Implement workarounds
        if not self._info:
            self._info = gl.gl_info.GLInfo()
            self._info.set_active_context()
            for attr, check in self._workaround_checks:
                setattr(self, attr, check(self._info))

        # Release textures and buffers on this context scheduled for deletion.
        # Note that the garbage collector may introduce a race condition,
        # so operate on a copy of the textures/buffers and remove the deleted
        # items using list slicing (which is an atomic operation)
        if self.object_space._doomed_textures:
            textures = self.object_space._doomed_textures[:]
            textures = (gl.GLuint * len(textures))(*textures)
            gl.glDeleteTextures(len(textures), textures)
            self.object_space._doomed_textures[0:len(textures)] = []
        if self.object_space._doomed_buffers:
            buffers = self.object_space._doomed_buffers[:]
            buffers = (gl.GLuint * len(buffers))(*buffers)
            gl.glDeleteBuffers(len(buffers), buffers)
            self.object_space._doomed_buffers[0:len(buffers)] = []


if __name__ == "__main__":
    from tkinter import *

    root = Tk()
    app = PygletTkFrame(root, width=320, height=200)
    app.pack(fill=BOTH, expand=YES)
    app.mainloop()
