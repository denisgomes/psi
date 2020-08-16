from __future__ import division

from pyglet.graphics import Batch
from psi.gui.scene import Scene
from psi.gui.pygletk import PygletTkFrame


class Viewport(PygletTkFrame):

    def __init__(self, master, model, *args, mouse=TKMouse,
                 keyboard=TKKeyboard, **kw):
        super(Viewport, self).__init__(*args, **kwargs)
        self._master = master
        self._scene = Scene(self, model)
        self._mouse = Mouse(self, model)
        self._keyboard = Keyboard(self, model)

    @property
    def mouse(self):
        return self._mouse

    @mouse.setter
    def mouse(self, value):
        self._mouse = value
        self._mouse.bind()

    @property
    def keyboard(self):
        return self._keyboard

    @keyboard.setter
    def keyboard(self, value):
        self._keyboard = value
        self._keyboard.bind()

    @property
    def scene(self):
        return self._scene

    def init_gl(self):
        Viewport.init_gl(self)
        self.scene.init_gl()

    def create_objects(self):
        """Create opengl objects when opengl is initialized"""
        self.scene.create_objects()

    def update_object_resize(self, width, height):
        """Called when the window receives only if opengl is initialized"""
        self.scene.update_object_resize(width, height)

    def draw_objects(self):
        """Called in the middle of on_draw after the buffer has been cleared"""
        # self.update_idletasks()
        self.scene.draw_objects()
        # self.swap_buffers()


class Mouse(object):
    """A three buttom mouse.

    Left Button - Rotate
    Middle Button - Zoom (Fine/Coarse)
    Right Button - Pan
    """

    def __init__(self, master, model):
        self.master = master    # parent window
        self.model = model
        self.xpos, self.ypos = 0, 0
        self.old_xpos, self.old_ypos = 0, 0
        self.viewx, self.viewy, self.viewz = 0, 0, 0
        self.old_viewx, self.old_viewy, self.old_viewz = 0, 0, 0

        self.tx, self.ty, self.tz = 0, 0, 0
        self.rx, self.ry = 0, 0

        # bind events on master
        self.bind()

    def get_norm_coords(self):
        """Get the normalized screen coordinates wrt to a coordinate system
        centered at the middle of the screen.  This system has a range from
        -1 < x < 1 and -1 < y < 1.
        """
        w, h = self.master.size
        win_y = (h - self.ypos) - h / 2
        norm_y = win_y / (h / 2)
        win_x = self.xpos - (w / 2)
        norm_x = win_x / (w / 2)
        return (norm_x, norm_y)

    def clear(self):
        self.tx, self.ty, self.tz, self.rx, self.ry = 0, 0, 0, 0, 0

    def bind(self):
        # left mouse button
        self.master.bind("<Button-1>", self.on_left_click)
        self.master.bind("<ButtonRelease-1>", self.on_left_release)
        self.master.bind("<B1-Motion>", self.on_left_drag)
        self.master.bind("<Double-1>", self.on_left_double_click)

        # middle mouse button
        self.master.bind("<Button-2>", self.on_middle_click)
        self.master.bind("<ButtonRelease-2>", self.on_middle_release)
        self.master.bind("<B2-Motion>", self.on_middle_drag)
        self.master.bind("<Double-2>", self.on_middle_double_click)

        # right mouse button
        self.master.bind("<Button-3>", self.on_right_click)
        self.master.bind("<ButtonRelease-3>", self.on_right_release)
        self.master.bind("<B3-Motion>", self.on_right_drag)
        self.master.bind("<Double-3>", self.on_right_double_click)

        # motion
        self.master.bind("<Motion>", self.on_motion)

    def on_left_click(self, event):
        # set keyboard focus
        self.master.focus_set()
        # do something else
        # print "clicked at", event.x, event.y
        self.xpos, self.ypos = self.old_xpos, self.old_ypos = event.x, event.y

    def on_left_release(self, event):
        self.xpos, self.ypos = event.x, event.y
        self.master.config(cursor="arrow")

        # selection
        self.master.scene.clear_hits()

        norm_x, norm_y = self.get_norm_coords()
        self.master.scene.select(norm_x, norm_y)

    def on_left_drag(self, event):
        self.master.config(cursor="exchange")

        x, y = event.x, event.y
        self.old_xpos, self.old_ypos = self.xpos, self.ypos
        self.xpos, self.ypos = x, y
        dx, dy = self.xpos-self.old_xpos, self.ypos-self.old_ypos

        self.tx, self.ty, self.tz = 0, 0, 0
        self.ry, self.rx = (dx*self.model.settings.mouse_rotate_factor,
                            dy*self.model.settings.mouse_rotate_factor)

        self.master.draw_objects()

    def on_left_double_click(self, event):
        # initiate scene selection
        pass

    def on_middle_click(self, event):
        # set keyboard focus
        self.master.focus_set()

    def on_middle_release(self, event):

        self.xpos, self.ypos = event.x, event.y
        self.master.config(cursor="arrow")

    def on_middle_drag(self, event):
        self.master.config(cursor="double_arrow")
        raise NotImplementedError("implement middle motion event")

    def on_middle_double_click(self, event):
        pass

    def on_right_click(self, event):
        # set keyboard focus
        self.master.focus_set()

        self.xpos, self.ypos = self.old_xpos, self.old_ypos = event.x, event.y

    def on_right_release(self, event):

        self.xpos, self.ypos = event.x, event.y
        self.master.config(cursor="arrow")

    def on_right_drag(self, event):
        self.master.config(cursor="fleur")
        x, y = event.x, event.y
        self.old_xpos, self.old_ypos = self.xpos, self.ypos
        self.xpos, self.ypos = x, y
        dx, dy = self.xpos-self.old_xpos, self.ypos-self.old_ypos

        self.tx, self.ty = (dx*self.model.settings.mouse_translate_factor,
                            dy*self.model.settings.mouse_translate_factor)

        self.tz, self.ry, self.rx = 0, 0, 0

        self.master.draw_objects()

    def on_right_double_click(self, event):
        pass

    def on_motion(self, event):
        """Shoot rays down for tentative selection"""
        pass


class Keyboard(object):

    def __init__(self, master, model):
        self.master = master
        self.model = model
        self.bind()

    def bind(self):
        # keypress
        self.master.bind("<Key>", self.on_keypress)

    def on_keypress(self, event):
        print("pressed", repr(event.char))


if __name__ == "__main__":
    import tkinter as tk
    from psi.model import Model

    win = tk.Tk()
    mdl = Model('test')

    viewport = Viewport(win, mdl, width=800, height=500)

    # pipe = viewport.scene.add_pipe((0, 0, 0), (0, 2, 0), 1)
    # viewport.scene.camera.translate(0, 0, 15)

    # print cyl.alpha

    win.title("Tkinter Pyglet")
    win.mainloop()
