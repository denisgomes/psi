import tkinter as tk
import tkinter.ttk as ttk


PAGE_HEIGHT = 144
RELIEF = "raised"


style = ttk.Style()


class Ribbon(ttk.Notebook):

    def __init__(self, parent):
        super(Ribbon, self).__init__(parent)
        self.parent = parent

        self._create_widgets()
        self._bind_widgets()

    def _create_widgets(self):
        pages = [FilePage, GeometryPage, ViewPage, SpecPage, LoadingPage,
                 AnalysisDesignPage, HelpPage]
        page_handles = ["file", "geometry", "view", "spec", "loading",
                        "analysisdesign", "help"]
        pages_text = ["FILE", "GEOMETRY", "VIEW", "SPECIFICATION", "LOADING",
                      "ANALYSIS & DESIGN", "HELP"]

        for handle, page, text in zip(page_handles, pages, pages_text):
            setattr(self, handle, page(self))
            self.add(getattr(self, handle), text=text)

    def _bind_widgets(self):
        pass


class FilePage(ttk.Frame):

    def __init__(self, parent):
        super(FilePage, self).__init__(parent, height=PAGE_HEIGHT,
                                       relief=RELIEF)
        self.parent = parent


class GeometryPage(ttk.Frame):

    def __init__(self, parent):
        super(GeometryPage, self).__init__(parent, height=PAGE_HEIGHT,
                                           relief=RELIEF)
        self.parent = parent


class ViewPage(ttk.Frame):

    def __init__(self, parent):
        super(ViewPage, self).__init__(parent, height=PAGE_HEIGHT,
                                       relief=RELIEF)
        self.parent = parent


class SpecPage(ttk.Frame):

    def __init__(self, parent):
        super(SpecPage, self).__init__(parent, height=PAGE_HEIGHT,
                                       relief=RELIEF)
        self.parent = parent


class LoadingPage(ttk.Frame):

    def __init__(self, parent):
        super(LoadingPage, self).__init__(parent, height=PAGE_HEIGHT,
                                          relief=RELIEF)
        self.parent = parent


class AnalysisDesignPage(ttk.Frame):

    def __init__(self, parent):
        super(AnalysisDesignPage, self).__init__(parent, height=PAGE_HEIGHT,
                                                 relief=RELIEF)
        self.parent = parent


class HelpPage(ttk.Frame):

    def __init__(self, parent):
        super(HelpPage, self).__init__(parent, height=PAGE_HEIGHT,
                                       relief=RELIEF)
        self.parent = parent


if __name__ == "__main__":
    win = tk.Tk()

    ribbon = Ribbon(win)
    ribbon.pack(side=tk.TOP, expand=True, fill=tk.X)

    win.mainloop()
