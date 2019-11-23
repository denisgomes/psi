from __future__ import division

from collections import OrderedDict
from contextlib import redirect_stdout
import gzip
import logging
import pickle
import copy
import sys

import numpy as np

from psi.settings import options
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi.topology import Geometry
from psi.units import units


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
        self._active_elements = None
        self._active_section = None
        self._active_material = None
        self._active_insulation = None
        self._active_code = None

        super(Model, self).__init__(name)   # call last
        self.activate()     # activate on init

    @property
    def settings(self):
        return self._settings

    @property
    def geometry(self):
        return self._geometry

    @property
    def points(self):
        return self._points.values()

    @property
    def elements(self):
        return self._elements.values()

    @property
    def sections(self):
        return self._sections.value()

    @property
    def materials(self):
        return self._materials.values()

    @property
    def insulation(self):
        return self._insulation.values()

    @property
    def codes(self):
        return self._codes.values()

    @property
    def sifs(self):
        return self._sifs.values()

    @property
    def supports(self):
        return self._supports.values()

    @property
    def loads(self):
        return self._loads.values()

    @property
    def loadcases(self):
        return self._loadcases.values()

    @property
    def active_point(self):
        return self.app.points.active_object

    @active_point.setter
    def active_point(self, point):
        self.app.points.active_object = point

    @property
    def active_elements(self):
        return self.app.elements.active_objects

    @active_elements.setter
    def active_elements(self, elements):
        self.app.elements.active_objects = elements

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

    def close(self):
        self.parent.close(self)

    def save(self):
        self.parent.save(self)

    def save_as(self, fname):
        self.parent.save_as(self, fname)

    def analyze(self):
        self.parent.analyze(self)

    def render(self):
        self.parent.render(self)


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
            name, inst = pickle.load(fp)   # model instance
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

    def analyze(self, inst):
        """Go through all loadcases, solve each one and then combine them to
        generate the final results.

        For each element calculate the element stiffness matrix and assemble
        the system stiffness matrix using the nodal degree of freedom matrix.

        The loads defined for each element loadcase is combined into a single
        load vector. The load vector for each element loadcase is then
        assembled into a global load vector. Multiple global load vectors are
        solved for at once using gaussian elimination.

        Note that to simplify the FEA solution, all elements are essentially
        beams, including bends and reducers which are approximations made by
        chaining multiple beams together.

        1. Create NDOF matrix.

        2. For each element
            a. Construct the local stiffness matrix.
            b. Construct local force matrix, one for each load case.
            c. Transform local stiffness and force matrix.
            d. Add global element stiffness matrix and the global element
               force vector to the global system stiffness and force vectors,
               respectively.

        3. Apply the boundary conditions using the penalty method.

        4. Solve the global system using guassian elimination techniques.

            AX=B

            Where A is the global system matrix.
            B is the matrix of force vectors, one for each loadcase, and
            X is the matrix of solutions vectors for each loadcase.

        5. Use the calculated displacements to calculate the element force and
           moments.

        6. Use the results from step 5 to calculate code stresses.
        """
        tqdm = logging.getLogger("tqdm")
        tqdm.info("*** Starting analysis...")

        # Note: All internal object data is stored in SI once loaded from
        # external files, disabling unit consersion allows for working with
        # only SI
        units.disable()
        tqdm.info("*** Switching to base units.")

        # debuging code
        # with redirect_stdout(sys.__stdout__):
            # print(inst.points.values())
            # print(inst.elements.values())
            # print(inst.supports.values())
            # print(inst.loads.values())
            # print(element.klocal(200))

            # for element in inst.elements.values():
            #     print(element.klocal(273))

        # do stuff here
        tqdm.info("*** Assembling system stiffness and force matrix.")

        ndof = 6    # nodal degrees of freedom
        en = 2      # number of nodes per element
        nn = len(inst.points)
        lc = len(inst.loadcases)

        # similar to nodal dof matrix
        points = list(inst.points)

        # global system stiffness matrix
        Ks = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)

        # system force matrix, one loadcase per column where a loadcase
        # consists of one or more loads
        Fs = np.zeros((nn*ndof, lc), dtype=np.float64)

        for element in inst.elements:
            idxi = points.index(element.from_point)
            idxj = points.index(element.to_point)

            # node and corresponding dof (start, finish)
            niqi, niqj = idxi*ndof, idxi*ndof + ndof
            njqi, njqj = idxj*ndof, idxj*ndof + ndof

            # element stiffness at room temp, conservative stress
            ke = element.kglobal(273)

            # assemble global stiffness matrix
            Ks[niqi:niqj, niqi:niqj] += ke[niqi:niqj, niqi:niqj]
            Ks[niqi:niqj, njqi:njqj] += ke[niqi:niqj, njqi:njqj]
            Ks[njqi:njqj, niqi:niqj] += ke[njqi:njqj, niqi:niqj]
            Ks[njqi:njqj, njqi:njqj] += ke[njqi:njqj, njqi:njqj]

            # iterate each loadcase adding loads
            for i, loadcase in enumerate(inst.loadcases):
                with redirect_stdout(sys.__stdout__):
                    print(loadcase)

                # sum of all loads in a loadcase
                fe = np.zeros((en*ndof, 1), dtype=np.float64)

                for load in loadcase.loads:

                    if load in element.loads:
                        fe += load.fglobal(element)

                # assembly global system force matrix
                Fs[niqi:niqj, i] += fe[niqi:niqj]
                Fs[njqi:njqj, i] += fe[njqi:njqj]

        with redirect_stdout(sys.__stdout__):
            print(Fs)

        # supports and hangers using penalty method

        # solve using gaussian elimination - solve for displacements

        # extract element global forces and moments

        # calculate code stresses

        units.enable()

    def check(self, inst):
        """Check the model input parameters before analyzing"""
        pass
