# Pipe Stress Infinity (PSI) - The pipe stress analysis and design software.
# Copyright (c) 2021 Denis Gomes <denisgomes@consultant.com>

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

"""Base entity and entity container for all psi objects.

Entity manages objects created by user defined unique name. After an object is
created all functions work with objects only. A specific object can be called
its user defined name.
"""

import logging
from collections import OrderedDict


class Entity(object):
    """Base class for psi objects"""

    _app = None

    def __new__(cls, name, *args, **kwargs):
        # before creation and initialization
        app = cls._app

        # catch exceptions that raise cryptic messages
        if isinstance(name, (int, str)) is False:
            raise NameError("name must be integer or string")

        # avoid recursive import
        from psi.model import Model
        if cls is not Model:
            model = app.models.active_object
            assert model is not None, "create or activate a model"

        inst = super(Entity, cls).__new__(cls)

        # objects should not be replaced
        if name in inst.parent.objects:
            del(inst)
            raise NameError("name '%s' in use" % name)

        return inst

    def __init__(self, name, register=True):
        if register:
            self.parent.new(name, self)

    def outlog(self, message, level=logging.INFO):
        """Log message to output file"""
        logger = logging.getLogger("stdout")
        logger.log(level, message)

    def errlog(self, message, level=logging.INFO):
        """Log message to error file"""
        logger = logging.getLogger("stderr")
        logger.log(level, message)

    def statlog(self, message, level=logging.INFO):
        """Log message to the tqdm progress bar"""
        logger = logging.getLogger("tqdm")
        logger.log(level, message)

    def _key(self):
        return (self.type, self.name)

    def __eq__(self, other):
        return self._key() == other._key()

    def __hash__(self):
        return hash(self._key())

    @property
    def app(self):
        return self.__class__._app

    def parent(self):
        raise NotImplementedError("abstract method")

    @property
    def type(self):
        return self.__class__.__name__

    @property
    def name(self):
        return self.parent.name(self)

    @name.setter
    def name(self, value):
        self.parent.name(self, value)

    def delete(self):
        self.parent.delete(self)

    def __repr__(self):
        return "%s %s" % (self.type, self.name)


class EntityContainer(object):
    """Base container object for an entity"""

    _app = None

    def __init__(self):
        self._objects = OrderedDict()

    @property
    def objects(self):
        return self._objects

    @property
    def app(self):
        return self.__class__._app

    def __len__(self):
        return self.count()

    def outlog(self, message, level=logging.INFO):
        """Log message to output file"""
        logger = logging.getLogger("stdout")
        logger.log(level, message)

    def errlog(self, message, level=logging.INFO):
        """Log message to error file"""
        logger = logging.getLogger("stderr")
        logger.log(level, message)

    def statlog(self, message, level=logging.INFO):
        """Log message to the tqdm progress bar"""
        logger = logging.getLogger("tqdm")
        logger.log(level, message)

    def list(self):
        return list(self._objects.values())

    def count(self):
        return len(self._objects)

    def new(self, name, inst):
        self._objects[name] = inst

    def name(self, inst, new_name=None):
        """Given an instance return the name or assign a new name if the
        name does not already exist.
        """
        if new_name:
            assert new_name not in self._objects.keys(), "name is taken"

        for key, val in self._objects.items():
            if val is inst:
                if new_name is None:
                    return key
                else:
                    del(self._objects[key])
                    self.new(new_name, val)

    def delete(self, inst):
        for key, val in self._objects.items():
            if val is inst:
                del self._objects[key]
                return val

    def update(self, inst, **kwargs):
        for key, val in self._objects.items():
            if val is inst:
                for k, v in kwargs.item():
                    setattr(inst, k, v)
                return val

    def __iter__(self):
        for obj in self._objects.values():
            yield obj

    def check(self):
        raise NotImplementedError("implement")

    def __call__(self, name):
        """Return an object given its name"""
        return self._objects[name]

    def __contains__(self, obj):
        """Return true if the object exists"""
        return obj in self._objects.values()

    def __getitem__(self, key):
        # return self.__call__(name)
        return list(self._objects.values())[key]

    def __repr__(self):
        return str(self.list())


class ActiveEntityMixin(object):

    @property
    def is_active(self):
        return self.parent.is_active(self)

    def activate(self):
        self.parent.activate(self)


class ActiveEntityContainerMixin(object):
    """The parent container keeps track of the active object"""

    def __init__(self):
        self._active_object = None

    def is_active(self, inst):
        if inst is self._active_object:
            return True
        else:
            return False

    @property
    def active_object(self):
        return self._active_object

    @active_object.setter
    def active_object(self, inst):
        self._active_object = inst

    def activate(self, inst):
        self._active_object = inst

    def delete(self, inst):
        if self._active_object is inst:
            self._active_object = None

        self.delete(inst)
