"""Defining units for object attributes that make up the data model.

By default the core data model does not use explicit units, though it has
implicit 'internal' base SI units associated with each of the managed
quantaties. When the user changes or defines a new unit system, the internal
units are unchanged. The various data model quantaties are converted from user
defined units to 'internal' units and back again on the fly. As a result, the
core data model quantaties can remain consistent. All calculations are also
therefore consistent.
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


class Units(object):
    """Allows for units to be defined for managed attributes"""

    def __init__(self, base_units="base", user_units="english"):
        self.set_base_units(base_units)
        self.set_user_units(user_units)

    def load_units_file(self, name):
        with open(os.path.join(UNITS_DIRECTORY, name + ".csv")) as csvfile:
            reader = csv.DictReader(csvfile)
            units = next(reader)

        return units

    def set_base_units(self, name=None):
        """The user may set the base units. Use at own risk!"""
        if name is None:
            name = "base"

        Quantity.base_units.update(self.load_units_file(name))

    def set_user_units(self, name=None):
        """Load user defined units"""
        if name is None:
            name = options["core.units"]

        Quantity.user_units.update(self.load_units_file(name))

    def disable(self):
        """Turn off unit conversion so that only values in default units are
        returned. Automatic unit conversion is turned off before the model is
        analyzed so that base SI units are used.
        """
        Quantity.is_active = False

    def enable(self):
        """Turn on automatic unit conversion which is the default."""
        Quantity.is_active = True

    @staticmethod
    def define(**kwargs):
        """Class decorator that helps define quantities for managed attributes.
        These attributes have a Quantity descriptor that stores a value and
        dimensionality. For example, a 10in sch 40 pipe has a diameter value of
        10.75in and a dimentionality of 'length'. When the attribute value is
        set using the current user units, the value is converted to base units
        and stored.  When the attribute is accessed, the value is converted
        from base units to user units on the fly.
        """
        def decorate(cls):
            for name, unit in kwargs.items():
                setattr(cls, name, Quantity(name, unit))
            return cls
        return decorate


units = Units()
