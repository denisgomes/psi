# Copyright (c) 2019 Denis Gomes
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
#    * Neither the name of Pipe Stress Infinity (PSI) nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Defining units for object attributes that make up the data model.

By default the core data model does not use units explicitly. It has implicit
'internal' base SI units associated with each of the managed quantaties. When
the user changes or defines a new unit system, the internal units and values
remain changed. The various data model quantaties are converted from user
defined units to 'internal' units and back again on the fly at runtime. As a
result, the core data model quantaties can remain consistent when the model is
finally analyzed in base units.
"""

import os
import csv

from pint import UnitRegistry

from psi import UNITS_DIRECTORY
from psi.settings import options

UREG = UnitRegistry()
Q_ = UREG.Quantity


class Quantity(object):
    """Descriptor for attributes with units.

    Parameters
    ----------
    name : str
        Attribute name

    utype : str
        Default internal type of unit, i.e. length, pressure, etc.

    is_active : bool
        If classmethod is set to True, unit conversion is set active.
    """

    base_units = {}
    user_units = {}
    is_active = True

    def __init__(self, name, utype):
        self.name = name
        self.utype = utype

    def __get__(self, inst, cls):
        if inst is None:
            # if attribute called from the class
            return self
        else:
            value = inst.__dict__[self.name]

            # conversion is off
            if self.is_active is False:
                return value

            try:
                # do conversion on value
                qty = Q_(value, self.base_units[self.utype])
                usr_qty = qty.to(self.user_units[self.utype])
                return usr_qty.magnitude
            except TypeError as e:
                if "NoneType" in str(e):
                    # ie. for a temp independent property
                    # temp is None, unit conversion skipped
                    return value
                else:
                    raise e

    def __set__(self, inst, value):
        # conversion is off
        if self.is_active is False:
            inst.__dict__[self.name] = value

        else:

            try:
                # do conversion on value
                qty = Q_(value, self.user_units[self.utype])
                def_qty = qty.to(self.base_units[self.utype])
                inst.__dict__[self.name] = def_qty.magnitude
            except TypeError as e:
                if "NoneType" in str(e):
                    # ie. for a temp independent property
                    # temp is None, unit conversion skipped
                    inst.__dict__[self.name] = value
                else:
                    raise e


class UnitsContextManagerMixin(object):
    """Allows for any units to be used within a specific code section"""

    def __enter__(self):
        self.set_user_units(self.user_units)

    def __exit__(self, type, value, traceback):
        name = self.app.models.active_object.units
        self.set_user_units(name)


class Units(UnitsContextManagerMixin):
    """Allows for units to be defined for managed attributes"""

    _app = None

    def __init__(self, base_units="base", user_units="english"):
        self.set_base_units(base_units)
        self.set_user_units(user_units)

    @property
    def app(self):
        return self._app

    def load_units_file(self, name):
        with open(os.path.join(UNITS_DIRECTORY, name + ".csv")) as csvfile:
            reader = csv.DictReader(csvfile)
            units = next(reader)

        return units

    def set_base_units(self, name=None):
        """The user may set the base units. Use at your own risk!"""
        if name is None:
            name = "base"

        self.base_units = name
        Quantity.base_units.update(self.load_units_file(name))

    def set_user_units(self, name=None):
        """Load user defined units"""
        if name is None:
            name = options["core.units"]

        self.user_units = name
        Quantity.user_units.update(self.load_units_file(name))

    def disable(self):
        """Turn off unit conversion so that only values in base units are
        returned. Automatic unit conversion is turned off before the model is
        analyzed so that base SI units are used.
        """
        Quantity.is_active = False

    def enable(self):
        """Turn on automatic unit conversion which is the default."""
        Quantity.is_active = True


def define(**kwargs):
    """Class decorator that helps define quantities for managed attributes.
    These attributes have a Quantity descriptor that stores a value and
    dimensionality. For example, a 10in sch 40 pipe has a diameter value of
    10.75in and a dimentionality of 'length'. When the attribute value is set
    using the current user units, the value is converted to base units and
    stored.  When the attribute is accessed, the value is converted from base
    units to user units on the fly at runtime.
    """
    def decorate(cls):
        for name, unit in kwargs.items():
            setattr(cls, name, Quantity(name, unit))
        return cls
    return decorate
