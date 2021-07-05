"""The PSI Graphical User Interface (GUI)"""

import tkinter as tk
import tkinter.ttk as ttk


class Application(ttk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.pack()

        self.create_widgets()
        self.do_layout()

    def create_widgets(self):
        # make ribbon
        # pyglettk area
        # make tables
        # make status
        pass

    def do_layout(self):
        pass


if __name__ == "__main__":
    root = tk.Tk()
    app = Application(master=root)
    root.title("PS8")
    root.mainloop()
