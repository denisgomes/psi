# Pipe Stress Infinity (PSI) - The pipe stress design and analysis software.
# Copyright (c) 2019 Denis Gomes

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import division

from collections import OrderedDict
import gzip
import pickle

# from psi.settings import options
from psi.settings import Configuration
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.topology import Geometry
from psi.solvers import static, modal


# TODO: Raise exception if model is not active on attribute access
# TODO: Store settings locally, validate and merge on model open


class Model(Entity, ActiveEntityMixin):
    """The model object contains all internal objects."""

    def __init__(self, name):
        """Create a model instance.

        Example
        -------

        .. code-block:: python

            >>> mdl = Model("mdl1")
        """
        self._config = Configuration()
        self._jobname = None

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
        self._reports = OrderedDict()

        # active model objects
        self._active_point = None
        self._active_elements = None
        self._active_section = None
        self._active_material = None
        self._active_insulation = None
        self._active_code = None
        self._active_report = None

        super(Model, self).__init__(name)   # call last
        self.activate()     # activate on init


    @property
    def jobname(self):
        """Get the jobname"""
        return self._jobname

    @jobname.setter
    def jobname(self, name):
        """Set the jobname"""
        self._jobname = name

    @property
    def settings(self):
        """Get the settings dictionary"""
        return self._config

    @property
    def geometry(self):
        """Get the geometry instance"""
        return self._geometry

    @property
    def points(self):
        """Get the list of point"""
        return self._points.values()

    @property
    def elements(self):
        """Get the list of elements"""
        return self._elements.values()

    @property
    def sections(self):
        """Get the list of sections"""
        return self._sections.values()

    @property
    def materials(self):
        """Get the list of materials"""
        return self._materials.values()

    @property
    def insulation(self):
        """Get the list of insulation"""
        return self._insulation.values()

    @property
    def codes(self):
        """Get the list of codes"""
        return self._codes.values()

    @property
    def sifs(self):
        """Get the list of sifs"""
        return self._sifs.values()

    @property
    def supports(self):
        """Get the list of supports"""
        return self._supports.values()

    @property
    def loads(self):
        """Get the list of loads"""
        return self._loads.values()

    @property
    def loadcases(self):
        """Get the list of loadcases"""
        return self._loadcases.values()

    @property
    def active_point(self):
        """Get the active point instance"""
        return self.app.points.active_object

    @active_point.setter
    def active_point(self, point):
        """Set the active point instance"""
        self.app.points.active_object = point

    @property
    def active_elements(self):
        """Get the list of active elements"""
        return self.app.elements.active_objects

    @active_elements.setter
    def active_elements(self, elements):
        """Set the list of active elements"""
        self.app.elements.active_objects = elements

    @property
    def active_section(self):
        """Get the active section"""
        return self.app.sections.active_object

    @active_section.setter
    def active_section(self, section):
        """Set the active section"""
        self.app.sections.active_object = section

    @property
    def active_material(self):
        """Get the active material"""
        return self.app.materials.active_object

    @active_material.setter
    def active_material(self, material):
        """Set the active material"""
        self.app.materials.active_object = material

    @property
    def active_insulation(self):
        """Get the active insulation"""
        return self.app.insulation.active_object

    @active_insulation.setter
    def active_insulation(self, insulation):
        """Set the active insulation"""
        self.app.insulation.active_object = insulation

    @property
    def active_code(self):
        """Get the active code"""
        return self.app.codes.active_object

    @active_code.setter
    def active_code(self, code):
        """Set the active code"""
        self.app.codes.active_object = code

    @property
    def active_report(self, code):
        """Get the active report"""
        return self.app.reports.active_object

    @active_report.setter
    def active_report(self, report):
        """Set the active report"""
        self.app.reports.active_object = report

    @property
    def parent(self):
        """Return the parent ModelContainer instance"""
        return self.app.models

    def close(self):
        """Close the model instance"""
        self.parent.close(self)

    def save(self):
        """Save the model instance persistantly"""
        self.parent.save(self)

    def save_as(self, fname):
        """Save the model using the given filename"""
        self.parent.save_as(self, fname)

    def analyze(self, mode="static"):
        """Run the analysis"""
        self.parent.analyze(self, mode)


class ModelContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(ModelContainer, self).__init__()
        self.Model = Model

    def open(self, fname, merge=False):
        """Open a model file.

        Load the gunzipped pickled file. If merge is set to True, try to merge
        model settings with application settings even if versions are not the
        same.
        """
        with gzip.open(fname, 'rb') as fp:
            name, inst = pickle.load(fp)   # model instance

            # the user units should be set here
            ##
            self.new(name, inst)
            inst.activate()

            return inst

    def activate(self, inst):
        app = self.app

        app.models._active_object = inst

        app.points._objects = inst._points
        app.points._active_object = inst._active_point

        app.elements._objects = inst._elements
        app.elements._active_objects = inst._active_elements

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

        app.reports._objects = inst._reports
        app.reports._active_object = inst._active_report

        # add others here

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
            fp.write(pickle.dumps((inst.name, inst), 1))
        return inst

    def analyze(self, inst, mode="static"):
        """Run the analysis"""
        if mode == "static":
            static(inst)
        elif mode == "modal":
            modal(inst)

    def check(self, inst):
        """Check the model input parameters before analyzing"""
        pass
