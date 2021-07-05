import tkinter.ttk as ttk


class Ribbon(ttk.Notebook):

    def __init__(self, parent):
        super(Ribbon, self).__init__(parent)


class FilePage(ttk.Frame):
    pass


class GeometryPage(ttk.Frame):
    pass


class ViewPage(ttk.Frame):
    pass


class SpecPage(ttk.Frame):
    pass


class Loading(ttk.Frame):
    pass


class AnalysisDesignPage(ttk.Frame):
    pass


class HelpPage(ttk.Frame):
    pass
