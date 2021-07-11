"""The PSI Graphical User Interface (GUI)"""

import tkinter as tk
import tkinter.ttk as ttk

from psi.gui.ribbon import Ribbon
from psi.gui.graphics.pyglettk import PygletTkFrame


class GUI(ttk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self._create_widgets()
        self._do_layout()

    def _create_widgets(self):
        # make ribbon
        self.ribbon = Ribbon(self)

        # pyglettk area
        self.graphics = PygletTkFrame(self)

        # icon frame
        self.icon_bar = IconBar(self)

        # make tables, interp, etc...
        self.footer = ttk.Frame(self, height=144)

        # make status
        self.status_bar = StatusBar(self)

    def _do_layout(self):
        self.ribbon.pack(fill=tk.X)
        self.icon_bar.pack(fill=tk.X)
        self.graphics.pack(expand=True, fill=tk.BOTH)
        self.footer.pack(fill=tk.X)
        self.status_bar.pack(fill=tk.X)

        self.pack(expand=True, fill=tk.BOTH)


class Footer(ttk.Frame):

    def __init__(self, parent):
        super(Footer, self).__init__(parent)
        self.parent = parent

        self._create_widgets()
        self._do_layout()

    def _create_widgets(self):
        self.data_notebook_frame = ttk.Frame(self, height=144)
        self.data_notebook = ttk.Notebook(self.data_notebook_frame)

        self.data_notebook_element_frame = ttk.Frame(self.data_notebook)
        self.data_notebook.add(self.data_notebook_element_frame, text="Element")

    def _do_layout(self):
        self.data_notebook_element_frame.pack(fill=tk.X)

        self.data_notebook.pack(fill=tk.X)
        self.data_notebook_frame.pack(fill=tk.X)

        self.pack(expand=True, fill=tk.X)


class StatusBar(ttk.Frame):

    def __init__(self, parent):
        super(StatusBar, self).__init__(parent)
        self.parent = parent

        self._create_widgets()
        self._do_layout()

    def _create_widgets(self):
        self.notify_label = ttk.Label(self)

    def _do_layout(self):
        self.notify_label.pack(expand=True, fill=tk.X)


class IconBar(ttk.Frame):

    def __init__(self, parent):
        super(IconBar, self).__init__(parent)
        self.parent = parent

        self._create_widgets()
        self._do_layout()

    def _create_widgets(self):
        self.single_select_button = ttk.Button(self)
        self.sep1 = ttk.Separator(self, orient="vertical")
        self.range_select_button = ttk.Button(self)

    def _do_layout(self):
        self.single_select_button.pack(side=tk.LEFT, padx=1, pady=1)
        self.sep1.pack(side=tk.LEFT, fill=tk.Y, padx=3, pady=3)
        self.range_select_button.pack(side=tk.LEFT, padx=1, pady=1)


if __name__ == "__main__":
    root = tk.Tk()
    app = GUI(master=root)
    root.title("PS8")
    root.mainloop()
