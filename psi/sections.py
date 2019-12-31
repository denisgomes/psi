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

from abc import abstractproperty
import csv
from math import pi

import psi
from psi.entity import (Entity, EntityContainer, ActiveEntityMixin,
                        ActiveEntityContainerMixin)
from psi import units


class Section(Entity, ActiveEntityMixin):

    def __init__(self, name):
        super(Section, self).__init__(name)
        self.activate()

    @property
    def parent(self):
        return self.app.sections

    def apply(self, elements=None):
        self.parent.apply(self, elements)

    @abstractproperty
    def area(self):
        pass

    @abstractproperty
    def ixx(self):
        pass

    @abstractproperty
    def iyy(self):
        pass

    @abstractproperty
    def izz(self):
        pass


@units.define(od="length", thk="length")
class Pipe(Section):

    @classmethod
    def from_file(cls, name, nps, sch, corra=None, milltol=None, fname=None,
                  default_units='english'):
        """Create a pipe object from a csv data file.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        nps : str
            Nominal pipe size.

        sch : str
            Pipe schedule designation.

        corro : float
            Corrosion allowance.

        milltol : float
            Mill tolerance.

        fname : str
            Full path to the csv data file used to do the lookup.

        default_units : str
            The units used for the data. Must be one of the units defined in
            the psi.UNITS_DIRECTORY path.

        Example
        -------
        Create a 10" schedule 40 pipe and activate the section.

        .. code-block:: python

            >>> p1 = Pipe.from_file("p1", "10", "40")
        """
        if fname is None:
            fname = psi.PIPE_DATA_FILE

        with open(fname, "r") as csvfile:
            reader = csv.DictReader(csvfile)

            # allows a schedule to have multiple names
            sch_map = {}
            schedules = reader.fieldnames[2:]
            for schedule in schedules:
                schs = schedule.split("/")
                for s in schs:
                    sch_map[s] = schedule

            for row in reader:
                if row["nps"] == nps:
                    try:
                        od = float(row["od"])
                        thk = float(row[sch_map[sch]])

                        # using data file units to initialize
                        with units.Units(user_units=default_units):
                            return cls(name, od, thk, corra=corra,
                                       milltol=milltol)

                    except ValueError:  # calling float on empty string
                        return None
            else:
                return None

    def __init__(self, name, od, thk, corra=None, milltol=None):
        """Create a pipe object.

        Parameters
        ----------
        name : str
            Unique name for pipe instance.

        od : float
            Actual outer diameter of pipe, not to be confused with the nomimal
            od. Note that for 14 inch pipe and above, the actual and nominal
            diameters are the same.

        thk : float
            Pipe wall thickness.

        corro :  float
            Corrosion allowance (CA).

            For systems that transport corrosive substances, the pipe wall may
            errode over the course of time. Typically, this type of corrosion
            is added on when the minimum design thickness is determined for
            pressure retention. Systems that experience significant corrosion
            are far more prone to fatigue related failures. A typical value for
            CA is 0.062 inches for piping.

        milltol : float
            Mill tolerance.

            When a mandrel is pushed through a hot billet to create the "hole"
            in a seamless pipe, a variation in thickness may occur if the
            mandrel is off-course. This tolerance, which is (+/-) 12.5% is
            typically accounted for in the min wall calculation per the
            applicable piping code so the user need not put in a value unless
            otherwise required. It does not apply to welded pipes which need to
            consider the weld joint efficiency factor since steel plates are
            easier to manufacture with high accuracy.

            The specified mill tolerance is typically used to calculate the
            hoop stress of the pipe depending on the code and has not bearing
            on the stiffeness calculation of the element

        Example
        -------
        Create a 10" schedule 40 pipe and activate the section.

        .. code-block:: python

            >>> p1 = Pipe("p1", 10.75, 0.365)
        """
        super(Pipe, self).__init__(name)
        self.od = od
        self.thk = thk
        self.corra = corra
        self.milltol = milltol

    @property
    def thke(self):
        """The effective wall thickness used for stress calcultions when
        corrosion and mill tolerance are taken into consideration. Both the
        pressure and bending stresses are affected, if the code calls for it.
        For example B31.1 does not require that the effective thickness be used
        for stress calculation.

        Only used for code stress calculations and as a result, the stress will
        be higher. This wall is not, however, used for the stiffness
        calculations, since a lighter wall will result in lower loads and
        ultimately give less conservative stress values. In other words, the
        actual thickness is used for the moment of inertia and area of pipe
        formulations.

        Depending on the code, the mill tolerance may be included for hoop
        stress calculation only to ensure the proper min wall was specified by
        the designer.
        """
        nomthk = effthk = self.thk

        if self.corra:
            effthk -= self.corra

        if self.milltol:   # pos or neg
            effthk += (self.milltol/100.0) * nomthk

        return effthk

    @property
    def nps(self):
        """The nominal pipe diameter is used to determine the bend radius for
        a pipe bend.

        Note that the actual and nominal diameter are the same number for 14
        inch pipe and large.
        """
        pass

    @property
    def id_(self):
        """The pipe inner diameter, not to be confused with the object id."""
        return self.od - 2*self.thk

    @property
    def area(self):
        """Cross sectional area of the pipe."""
        do = self.od
        di = self.od - 2*self.thk
        return (pi/4) * (do**2 - di**2)

    @property
    def izz(self):
        """Area moment of inertia about the local z-z axis. The horizontal
        bending axis.
        """
        do = self.od
        di = self.od - 2*self.thk
        return (pi/64) * (do**4 - di**4)

    @property
    def iyy(self):
        """Area moment of inertia about the local y-y axis. The vertical
        bending axis.
        """
        return self.izz     # same due to symmetry

    @property
    def ixx(self):
        """Polar moment of inertia about the longitudinal local x-axis"""
        do = self.od
        di = self.od - 2*self.thk
        return (pi/32) * (do**4 - di**4)

    @property
    def is_thin_wall(self):
        """Check to see if the pipe is thin wall.

        Per the applicable code, the pipe pressure will influence the
        stiffness of a thin wall pipe. This is usually taken into account
        by modifying the flexibility factor of the element matrix by a
        pressure factor.

        If the od to thk ratio of a pipe is greater than 100 then the code
        computed sifs and flexibility factors no longer apply. This is usually
        true of large diameter duct piping, for which local effects have a
        significant influence on the overall behavior of the pipe and results.
        """
        od_over_thk = self.od / self.thk
        return od_over_thk > 10 and od_over_thk <= 100

    @property
    def is_heavy_wall(self):
        """Check to determine if the pipe is heavy wall.

        Pipes with a diameter to thickness ratio less than 10 are considered
        having heavy walls.
        """
        return (self.od / self.thk) <= 10

    @property
    def is_small_bore(self):
        """Pipe sizes 2" and below are considered small bore.

        Imperial units are used here for the relational test.
        """
        with units.Units(user_units="english"):
            return self.od <= 2

    @property
    def is_large_bore(self):
        """Pipe sizes larger than 2" are considered large bore."""
        return not self.is_small_bore

    @property
    def d2t(self):
        """The D/t ratio of the section.

        If a pipe has a D/t ratio of larger than 100 it behaves more like a
        shell than a beam.

        Code based SIF and flexibility factors do not apply for D/t ratios
        above 100 due to limitations in the testing performed by Markl and
        company.
        """
        return self.od / self.thk


@units.define(d="length", tw="length", bf="length", tf="length")
class WideFlange(Section):

    @classmethod
    def from_file(cls, name, size="W4x13", fname=None,
                  default_units="english"):
        if fname is None:
            fname = psi.BEAM_DATA_FILE

        with open(fname, "r") as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                if row["size"] == size:
                    d = float(row["d"])
                    tw = float(row["tw"])
                    bf = float(row["bf"])
                    tf = float(row["tf"])

                    with units.Units(user_units=default_units):
                        return cls(name, d, tw, bf, tf)

    def __init__(self, name, d, tw, bf, tf):
        super(WideFlange, self).__init__(name)
        self.d = d      # depth
        self.tw = tw    # web thickness
        self.bf = bf    # flange width
        self.tf = tf    # flange thickness

    @property
    def area(self):
        d = self.d
        tw = self.tw
        bf = self.bf
        tf = self.tf

        return 2*tf*bf + (d-2*tf)*tw

    @property
    def izz(self):
        """About strong axis"""
        d = self.d
        tw = self.tw
        bf = self.bf
        tf = self.tf

        return (1/12)*tw*(d-2*tf)**3 + (1/2)*(tf*bf)*(d-tf)**2

    @property
    def iyy(self):
        """About weak axis"""
        d = self.d
        tw = self.tw
        bf = self.bf
        tf = self.tf

        return (1/12)*(d-2*tf)*tw**3 + (1/6)*tf*bf**3

    @property
    def ixx(self):
        """About longitudinal axis"""
        return self.izz + self.iyy


class SectionContainer(EntityContainer, ActiveEntityContainerMixin):

    def __init__(self):
        super(SectionContainer, self).__init__()
        self.Pipe = Pipe
        self.WideFlange = WideFlange

    def apply(self, inst, elements=None):
        """Apply a section to elements.

        Example
        ------
        Create a 10" schedule 40 pipe and assign it to all active elements.

        .. code-block:: python

            >>> p1 = Pipe.from_file("p1", "10", "40")
            >>> sections.apply(p1)
        """
        if elements is None:
            elements = []

            for element in self.app.elements.active_objects:
                elements.append(element)

        for element in elements:
            element.section = inst
