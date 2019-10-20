from __future__ import division

from collections import OrderedDict
import gzip
import cPickle
import copy

from openpipe.settings import options
from openpipe.core.entity import (Entity, EntityContainer, ActiveEntityMixin,
                                  ActiveEntityContainerMixin)
from openpipe.core.topology import Geometry
from openpipe.core.units import units


# TODO: Raise exception if model is not active on attribute access
# TODO: Store settings locally, validate and merge on model open


class Model(Entity, ActiveEntityMixin):
    """The model object contains all other objects"""

    def __init__(self, name):
        self._settings = copy.deepcopy(options)
        self._geometry = Geometry()

        # internal containers
        self._points = OrderedDict()
        self._elements = OrderedDict()
        self._sections = OrderedDict()
        self._materials = OrderedDict()
        self._insulation = OrderedDict()
        self._codes = OrderedDict()
        self._sifs = OrderedDict()
        self._supports = OrderedDict()
        self._loads = OrderedDict()
        self._loadcases = OrderedDict()

        # active model objects
        self._active_point = None
        self._active_section = None
        self._active_material = None
        self._active_insulation = None
        self._active_code = None

        super(Model, self).__init__(name)   # call last
        self.activate()     # activate on init

    @property
    def settings(self):
        return self._settings

    # def __getattribute__(self, name):
    #     # to avoid recursion
    #     app = super(Model, self).__getattribute__("_app")

    #     if self is app.models.active_object:
    #         return super(Model, self).__getattribute__(name)

    #     raise AttributeError("model must be active")

    @property
    def geometry(self):
        return self._geometry

    @property
    def active_point(self):
        return self.app.points.active_object

    @active_point.setter
    def active_point(self, point):
        self.app.points.active_object = point

    @property
    def active_section(self):
        return self.app.sections.active_object

    @active_section.setter
    def active_section(self, section):
        self.app.sections.active_object = section

    @property
    def active_material(self):
        return self.app.materials.active_object

    @active_material.setter
    def active_material(self, material):
        self.app.materials.active_object = material

    @property
    def active_insulation(self):
        return self.app.insulation.active_object

    @active_insulation.setter
    def active_insulation(self, insulation):
        self.app.insulation.active_object = insulation

    @property
    def active_code(self):
        return self.app.codes.active_object

    @active_code.setter
    def active_code(self, code):
        self.app.codes.active_object = code

    @property
    def parent(self):
        return self.app.models

    @property
    def units(self):
        """Assign model units derived from application units"""
        # return options["core.units"]
        return self._settings["core.units"]

    @units.setter
    def units(self, name):
        units.set_user_units(name)
        self._settings["core.units"] = name
        # options["core.units"] = name

    @property
    def vertical(self):
        return self._settings["core.vertical"]
        # return options["core.vertical"]

    @vertical.setter
    def vertical(self, value):
        self._settings["core.vertical"] = value
        # options["core.vertical"] = value

    @property
    def tref(self):
        return self._settings["core.tref"]
        # return options["core.tref"]

    @tref.setter
    def tref(self, value):
        self._settings["core.tref"] = value
        # options["core.tref"] = value

    def plot(self, size=(800, 600)):
        """In interactive mode, the user is allowed to visually inspect and
        interact with the model.
        """
        import wx
        from openpipe.ui.viewer import OpenpipeModelViewer
        app = wx.App(redirect=False)
        frm = OpenpipeModelViewer(self, size=size)
        frm.Show()
        app.MainLoop()

    def close(self):
        self.parent.close(self)

    def save(self):
        self.parent.save(self)

    def save_as(self, fname):
        self.parent.save_as(self, fname)


class ModelContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(ModelContainer, self).__init__()
        self.Model = Model

    def open(self, fname, merge=False):
        """Open a model file. Load the gunzipped pickled file.

        If merge is set to True, try to merge model settings with application
        settings even if versions are not the same.
        """
        with gzip.open(fname, 'rb') as fp:
            name, inst = cPickle.load(fp)   # model instance
            self.new(name, inst)
            inst.activate()

            return inst

    def activate(self, inst):
        app = self.app

        app.models._active_object = inst

        app.points._objects = inst._points
        app.points._active_object = inst._active_point

        app.sections._objects = inst._sections
        app.sections._active_object = inst._active_section

        app.materials._objects = inst._materials
        app.materials._active_object = inst._active_material

        app.insulation._objects = inst._insulation
        app.insulation._active_object = inst._active_insulation

        app.codes._objects = inst._codes
        app.codes._active_object = inst._active_code

        app.sifs._objects = inst._sifs

        app.supports._objects = inst._supports

        app.loads._objects = inst._loads

        app.loadcases._objects = inst._loadcases

    def close(self, inst):
        """Closes a model"""
        self.delete(inst)

    def save(self, inst=None):
        """Saves a gunzipped pickled model to disk.  This will overwrite any
        existing file with the same name inside the directory.  If an instance
        is not provided, the active model is saved.
        """
        if inst is None:
            inst = self._active_object
        fname = inst.name
        self.save_as(inst, fname)

    def save_as(self, inst, fname):
        """Save a model with a different filename"""
        with gzip.GzipFile(fname, 'wb') as fp:
            fp.write(cPickle.dumps((inst.name, inst), 1))
        return inst
