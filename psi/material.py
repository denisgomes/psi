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

from __future__ import division
import csv

import psi
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi import units


@units.define(_temps="temperature")
class Property(object):

    def __init__(self, name):
        self._name = name
        self._temps = []
        self._values = []

    @property
    def name(self):
        return self._name

    @property
    def _data(self):
        return list(zip(self._temps, self._values))

    @_data.setter
    def _data(self, data):
        self._temps, self._values = list(zip(*data))

    @property
    def value(self):
        """Return a single temperature independent property value"""
        if self._data and self._data[0][0] is None:
            return self._data[0][1]

    @value.setter
    def value(self, val):
        """Used to set a single value for a temperature independent
        property.
        """
        self._data = [(None, val)]

    @property
    def table(self):
        """Return tabular list of (temp, value) pairs"""
        # temp independent property
        if len(self._data) == 1 and self._data[0][0] is None:
            return self._data[0][1]

        return self._data

    @table.setter
    def table(self, data):
        """Set tabular list of (temp, value) pairs"""
        self._data = sorted(data, key=lambda val: val[0])

    def __getitem__(self, item):
        """Temperature dependent interpolation"""
        tvals = self._data  # temp vs. vals

        # not a func of temp
        if len(tvals) == 1 and tvals[0][0] is None:
            return tvals[0][1]
        else:
            # if table is not empty
            if tvals:
                sorted(tvals, key=lambda p: p[0])   # by temp
                temps, vals = zip(*tvals)

                # input temp in range?
                mintemp = min(temps)
                maxtemp = max(temps)
                if item < mintemp or item > maxtemp:
                    return None  # temp out of range

                # find bounding indices and interpolate
                for idx, temp in enumerate(temps):
                    if temp >= item:
                        uidx = idx
                        lidx = uidx - 1
                        break
                if temp == item:    # exact match
                    return vals[uidx]
                else:
                    x1, y1 = temps[lidx], vals[lidx]
                    x2, y2 = temps[uidx], vals[uidx]

                    try:
                        return y1 - (y2-y1)/(x2-x1) * (x1-item)
                    except TypeError:   # None value
                        return None

    def clear(self):
        self._temps = []
        self._values = []

    def __repr__(self):
        if len(self._data) == 1:
            return str(self._data[0][1])

        return str(self.table)


@units.define(_values="expansion_coefficient")
class Alp(Property):
    """Coefficient of thermal expansion"""

    def __init__(self):
        super(Alp, self).__init__("Thermal Expansion Coefficient")


@units.define(_values="elastic_modulus")
class YMod(Property):
    """Hot Young's modulus"""

    def __init__(self):
        super(YMod, self).__init__("Young's Modulus")


@units.define(_values="pipe_density")
class Rho(Property):
    """Density"""

    def __init__(self):
        super(Rho, self).__init__("Density")


class Nu(Property):
    """Poisson's ratio"""

    def __init__(self):
        super(Nu, self).__init__("Poisson's Ratio")


@units.define(_values="pressure")
class Sh(Property):
    """Hot material allowable"""

    def __init__(self):
        super(Sh, self).__init__("Hot Allowable")


class Material(Entity, ActiveEntityMixin):

    def __init__(self, name):
        super(Material, self).__init__(name)
        self._rho = Rho()
        self._nu = Nu()
        self._alp = Alp()
        self._sh = Sh()
        self._ymod = YMod()

        self.activate()

    def apply(self, elements=None):
        self.parent.apply(self, elements)

    @property
    def parent(self):
        return self.app.materials

    def activate(self):
        self.parent.activate(self)

    @classmethod
    def from_file(cls, name, material, code, fname=None,
                  default_units="english"):
        """Load a material from the database"""
        if fname is None:
            fname = psi.MATERIAL_DATA_FILE

        def num(st):
            """Convert a string to an int or float"""
            try:
                return int(st)
            except ValueError:
                try:
                    return float(st)
                except ValueError:
                    return None     # value missing

        def convert(prop):
            vals = prop.split(",")
            if len(vals) == 1:
                return num(*vals)
            else:
                return [num(val) for val in vals]

        with open(fname) as csvfile:
            reader = csv.DictReader(csvfile)

            (NAME, CODE, RHO, NU, TEMP, ALP, YMOD, SH) = \
                ("name", "code", "rho", "nu", "temp", "alp", "ymod", "sh")
            for row in reader:
                if row[NAME] == material and row[CODE] == code:
                    rho = row[RHO]
                    nu = row[NU]
                    temp = row[TEMP]
                    alp = row[ALP]
                    ymod = row[YMOD]
                    sh = row[SH]

                    rho = convert(rho)
                    nu = convert(nu)
                    tempc = convert(temp)
                    alpc = convert(alp)
                    ymodc = convert(ymod)
                    shc = convert(sh)

                    with units.Units(user_units=default_units):
                        mat = cls(name)

                        # rho and nu not func of temp
                        mat.rho.value = rho
                        mat.nu.value = nu

                        # remove temps with no values given
                        mat.alp.table = [(t, v) for t, v in zip(tempc, alpc)
                                         if v is not None]
                        mat.sh.table = [(t, v) for t, v in zip(tempc, shc)
                                        if v is not None]
                        mat.ymod.table = [(t, v) for t, v in zip(tempc, ymodc)
                                          if v is not None]

                        return mat
            else:
                return None

    @property
    def sh(self):
        return self._sh

    @property
    def alp(self):
        return self._alp

    @property
    def ymod(self):
        return self._ymod

    @property
    def rho(self):
        return self._rho

    @property
    def nu(self):
        return self._nu

    def __repr__(self):
        return "%s %s" % (self.type, self.name)


class MaterialContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(MaterialContainer, self).__init__()
        self.Material = Material

    def apply(self, inst, elements=None):
        """Apply a material to elements"""
        if elements is None:
            elements = []

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.material = inst