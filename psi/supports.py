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

"""Implementation of different types of theoretical pipe supports.

Supports are implemented using the penalty approach where the global system
stiffness and force matrices are modified by the support stiffness value and
displacement respectively.

Axis aligned supports directly scale the diagonal elements of the global
stiffness matrix. Inclined or skewed supports are implemented using constraint
equations and modify more than the diagonal terms as the displacements are
coupled via the direction cosines.

For support displacements, the support displacement vector is multiplied by the
support stiffness and added to the corresponding force vector. Refer to
"Introduction to Finite Elements in Engineering" by Chandrupatla and Belegundu
for additional details. The user specified global support displacement is
projected in the support direction to determine the displacement components. If
a component is 0, the support will not be active along that axis, and the pipe
will be able to move in that direction.

A stiffness value of 1000*K, where K is the largest stiffness in the global
stiffness matrix has shown to produce good results based on textbook examples.
Reasonable default values for translation and rotation stiffness are specified
for each support. The user can change the default values via model settings.

X, Y and Z supports are inherited from AbstractSupport and can take several
different forms. They can be snubbers (only active in occasional load cases),
single or bi-directional, translational or rotational and/or define friction
and gaps.
"""

import csv
import os
import warnings
import sys
from contextlib import redirect_stdout

import numpy as np

import psi
from psi import units
from psi.entity import Entity, EntityContainer
from psi.loads import Displacement
from psi.units import DEFAULT_UNITS


@units.define(translation_stiffness="translation_stiffness",
              rotation_stiffness="rotation_stiffness")
class Support(Entity):

    def __init__(self, name, point):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        translation_stiffnesss : float
            Stiffness in the translational directions.

            .. note::
               The default value is based on english units.

        rotation_stiffness : float
            Stiffness in the rotational directions.

            .. note::
               The default value is based on english units.
        """
        super(Support, self).__init__(name)
        self.point = point
        self.elememt = None

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness

    def apply(self, element):
        """Apply the support to the element.

        Parameters
        ----------
        element : Element
            An element object.

        Supports are applied to an element at a node on the element.

        Example
        -------
        Create a support at node 20 of run element 20.

        .. code-block:: python

            >>> run20 = elements(10, 20)
            >>> anc = Anchor("anc20", 20)
            >>> anc.apply([run20])
        """
        self.parent.apply([self], element)

    @property
    def parent(self):
        return self.app.supports

    def kglobal(self, element):
        """The support stiffness matrix used to modify the system stiffness
        matrix.

        .. note::

            kglobal is a matrix. For axis aligned supports only the diagonal
            terms consist of the penalty terms. This is not the case for
            skewed supports which modify more than just the diagonal elements.
        """
        k = np.zeros((12, 12), dtype=np.float64)

        c = self.cglobal(element)   # penalty vector
        di = np.diag_indices(6)

        if self.point == element.from_point.name:
            k[:6, :6][di] = c[:6, 0]

        elif self.point == element.to_point.name:
            k[6:12, 6:12][di] = c[6:12, 0]

        return k

    def cglobal(self, element):
        """The support penalty vector consisting of translation and rotation
        terms. The first 6 DOFs apply to node i and last 6 to node j. The first
        3 DOFs are translational DOFs for node i. DOFs 3 to 6 are the
        rotational DOFs for node i.

        .. note::

            The cglobal vector appears in the kglobal diagonal terms when the
            support is axis aligned.
        """
        raise NotImplementedError("implement")

    def dglobal(self, element, loadcase):
        """Nodal displacement vector used for penalty method.

        By default all supports except for a support with displacement load is
        assumed to have 0 imposed displacement.

        The displacement is given as a vector given with dx, dy, dz, rx, ry and
        rz at the support node.
        """
        a = np.zeros((12, 1), dtype=np.float64)

        disps = [disp for disp in element.loads if isinstance(disp,
                 Displacement)]
        for disp in disps:
            if disp in loadcase and self.element.from_point == disp.point:
                a[:6, 0] += disp.to_list()

            elif disp in loadcase and self.element.to_point == disp.point:
                a[6:12, 0] += disp.to_list()

        return a


class Anchor(Support):
    """Support with all 6 degrees of freedom at a node fixed."""

    def __init__(self, name, point):
        """Create an anchor support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        translation_stiffnesss : float
            Stiffness in the translational direction.

        rotation_stiffness : float
            Stiffness in the rotational direction.
        """
        super(Anchor, self).__init__(name, point)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.point == element.from_point.name:
            c[:3, 0] = [self.translation_stiffness] * 3
            c[3:6, 0] = [self.rotation_stiffness] * 3

        elif self.point == element.to_point.name:
            c[6:9, 0] = [self.translation_stiffness] * 3
            c[9:12, 0] = [self.rotation_stiffness] * 3

        return c


@units.define(_gap="length")
class AbstractSupport(Support):    # Abstract Base Class
    """A generic rigid support with different configurations.

    1. Translational or rotational (X, Y, Z, RX, RY, RZ)
    2. Single or double-acting, (+/-) (X, Y, Z, RX, RY, RZ)
    3. Gap and friction,
    4. Snubber (+/-) (XSNB, YSNB, ZSNB).

    The default support is a rigid double-acting translational support with no
    friction or gap, in other words a linear support.

    .. note::
        Snubbers are translational linear supports, so gaps and friction are
        not valid. They are active only for occasional cases.

    .. note::
        Rotational supports do not have friction.

    X, Y, Z - support for friction

    +X, +Y, +Z, +RX, +RY, +RZ - support for gaps
    +X, +Y, +Z - support for friction

    -X, -Y, -Z, -RX, -RY, -RZ - support for gaps
    -x, -Y, -Z - support for friction

    SNB - translational linear support
    """

    def __init__(self, name, point, direction=None, is_rotational=False):
        """Create a support instance at a node point.

        Parameters
        ----------
        name : str
            Unique name for pipe object.

        point : Point
            Point instance where support is located.

        direction : str
            Support direction. Default is None. "+" and "-" is used to
            specify a directional support (non-linear).

        mu : float
            Support friction (non-linear).

        is_rotational : bool
            True if the support is a rotational restraint.

        is_snubber : bool
            True if the support is snubber.

        translation_stiffnesss : float
            Stiffness in the translational direction.

        rotation_stiffness : float
            Stiffness in the rotational direction.

        gap : float
            Support gap (non-linear).
        """
        super(AbstractSupport, self).__init__(name, point)
        self.direction = direction
        self.is_rotational = is_rotational
        self._mu = 0
        self._gap = 0
        self._is_snubber = False

        model = self.app.models.active_object
        with units.Units(user_units=DEFAULT_UNITS):
            self.translation_stiffness = model.settings.translation_stiffness
            self.rotation_stiffness = model.settings.rotation_stiffness

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, value):
        if self._is_rotational or self._is_snubber:
            raise ValueError("cannot set friction")

        self._mu = value

    @property
    def gap(self):
        return self._gap

    @gap.setter
    def gap(self, value):
        if self._is_snubber:
            raise ValueError("cannot set gap")

        self._gap = value

    @property
    def is_rotational(self):
        return self._is_rotational

    @is_rotational.setter
    def is_rotational(self, value):
        assert isinstance(value, bool), "value must be boolean"

        self._is_rotational = value
        if value is True:
            self._is_snubber = False
            self._mu = 0.0

    @property
    def is_snubber(self):
        return self._is_snubber

    @is_snubber.setter
    def is_snubber(self, value):
        assert isinstance(value, bool), "value must be boolean"

        self._is_snubber = value
        if value is True:
            self._is_rotational = False
            self._mu = 0.0
            self._gap = 0.0

    @property
    def is_nonlinear(self):
        return self.has_friction or self.has_gap

    @property
    def has_friction(self):
        return self.mu > 0.0

    @property
    def has_gap(self):
        return self.gap > 0.0 or str(self.direction) in "+-"

    @property
    def type(self):
        supdir = self.direction if self.direction else ""

        if self.is_snubber:
            _type = supdir + super(AbstractSupport, self).type + "SNB"
        elif self.is_rotational:
            _type = supdir + "R" + super(AbstractSupport, self).type
        else:
            _type = super(AbstractSupport, self).type

        return _type


class Inclined(AbstractSupport):
    """A skewed roller support not aligned with a global axis.

    Parameters
    ----------
    dircos : tuple
        Support direction cosine. Also pointing in the positive support
        direction.

        .. note::
            The support direction gives the sense of the support, "+" or "-".
    """

    def __init__(self, name, point, dircos, direction=None,
                 is_rotational=False):
        super(Inclined, self).__init__(name, point, direction, is_rotational)
        self.dircos = dircos

    def kglobal(self, element):
        k = np.zeros((12, 12), dtype=np.float64)
        B1, B2, B3 = self.dircos
        B = np.array([[B1], [B2], [B3]], dtype=np.float64) @  \
            np.array([[B1, B2, B3]], dtype=np.float64)

        if self.is_rotational is False:
            k_trans = self.translation_stiffness * B

            if self.point == element.from_point.name:
                k[:3, :3] = k_trans[:, :]
            elif self.point == element.to_point.name:
                k[6:9, 6:9] = k_trans[:, :]

        else:
            k_rot = self.rotation_stiffness * B

            if self.point == element.from_point.name:
                k[3:6, 3:6] = k_rot[:, :]
            elif self.point == element.to_point.name:
                k[9:12, 9:12] = k_rot[:, :]

        return k

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        # does c need to be multiplied by dircos?
        # does d need to be multiplied by dircos?

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[:3, 0] = [self.translation_stiffness if dc != 0 else 0
                            for dc in self.dircos]
            elif self.point == element.to_point.name:
                c[6:9, 0] = [self.translation_stiffness if dc != 0 else 0
                             for dc in self.dircos]
        else:
            if self.point == element.from_point.name:
                c[3:6, 0] = [self.rotation_stiffness if dc != 0 else 0
                             for dc in self.dircos]
            elif self.point == element.to_point.name:
                c[9:12, 0] = [self.rotation_stiffness if dc != 0 else 0
                              for dc in self.dircos]

        return c

    @property
    def type(self):
        supdir = self.direction if self.direction else ""

        if self.is_snubber:
            _type = supdir + "INC" + "SNB"
        elif self.is_rotational:
            _type = supdir + "R" + "INC"
        else:
            _type = "INC"

        # add on the dircos
        _type += repr(self.dircos)

        return _type


class X(AbstractSupport):
    """Support aligned with the global x direction."""

    @property
    def dircos(self):
        return (1, 0, 0)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[0, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[6, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[3, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[9, 0] = self.rotation_stiffness

        return c


class Y(AbstractSupport):
    """Support aligned with the global y direction."""

    @property
    def dircos(self):
        return (0, 1, 0)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[1, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[7, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[4, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[10, 0] = self.rotation_stiffness

        return c


class Z(AbstractSupport):
    """Support aligned with the global z direction."""

    @property
    def dircos(self):
        return (0, 0, 1)

    def cglobal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[11, 0] = self.rotation_stiffness

        return c


class LimitStop(AbstractSupport):
    """Support aligned with the axial direction of the pipe.

    LimitStop supports are used to redirect thermal movement. They are commonly
    used for rack piping with expansion loops.
    """

    @property
    def dircos(self):
        # dircos depends on element direction
        r = (np.array(self.element.to_point.xyz, dtype=np.float64) -
             np.array(self.element.from_point.xyz, dtype=np.float64))
        dc = tuple(r / np.linalg.norm(r))

        return dc

    def clocal(self, element):
        """The local element x direction is the inline/axial direction"""
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                c[0, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                c[6, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                c[3, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                c[9, 0] = self.rotation_stiffness

        return c

    def cglobal(self, element):
        """Transform such that support is always inline with respect to the
        element coordinate system.
        """
        T = element.T()

        return T.transpose() @ self.clocal(element)

    @property
    def type(self):
        supdir = self.direction if self.direction else ""

        if self.is_snubber:
            _type = supdir + "LIM" + "SNB"
        elif self.is_rotational:
            _type = supdir + "R" + "LIM"
        else:
            _type = "LIM"

        return _type


class Lateral(AbstractSupport):
    """Support perpendicular to the pipe run direction.

    The support direction can be toggled between the two local element lateral
    directions of the pipe.

    .. note::
        If a pipe is horizontal, there is only one lateral direction. Any other
        configuration will result in two lateral directions (perpendicular to
        the pipe run direction), including for skewed piping.
    """

    _dir1 = True

    def flip(self):
        """Select alternate support lateral direction.

        .. note::
            This method does not affect the support direction for a support on
            horizontal piping.
        """
        if self._dir1:
            self._dir1 = False
        else:
            self._dir1 = True

    @property
    def dircos(self):
        vert = self.app.models.active_object.settings.vertical
        dc = self.element.dircos()

        dc1 = tuple(dc[1, :])
        dc2 = tuple(dc[2, :])

        local_x = dc[0, :]
        if vert == "y":
            up = np.array([[0], [1], [0]], dtype=np.float32)

            if np.dot(local_x, up) == 0:
                dc1 = dc2 = tuple(dc[2, :])

        elif vert == "z":
            up = np.array([[0], [0], [1]], dtype=np.float32)

            if np.dot(local_x, up) == 0:
                dc1 = dc2 = tuple(dc[1, :])

        if self._dir1:
            return dc1
        else:
            return dc2

    def clocal(self, element):
        c = np.zeros((12, 1), dtype=np.float64)

        if self.is_rotational is False:
            if self.point == element.from_point.name:
                if element.is_horizontal:
                    # element local y is always aligned with global y
                    # for horizontal elements, thus lateral in local z
                    c[2, 0] = self.translation_stiffness
                elif self._dir1:
                    c[1, 0] = self.translation_stiffness
                else:
                    c[2, 0] = self.translation_stiffness
            elif self.point == element.to_point.name:
                if element.is_horizontal:
                    c[8, 0] = self.translation_stiffness
                elif self._dir1:
                    c[7, 0] = self.translation_stiffness
                else:
                    c[8, 0] = self.translation_stiffness

        else:
            if self.point == element.from_point.name:
                if element.is_horizontal:
                    c[5, 0] = self.rotation_stiffness
                elif self._dir1:
                    c[4, 0] = self.rotation_stiffness
                else:
                    c[5, 0] = self.rotation_stiffness
            elif self.point == element.to_point.name:
                if element.is_horizontal:
                    c[11, 0] = self.rotation_stiffness
                elif self._dir1:
                    c[10, 0] = self.rotation_stiffness
                else:
                    c[11, 0] = self.rotation_stiffness

        return c

    def cglobal(self, element):
        T = element.T()

        return T.transpose() @ self.clocal(element)

    @property
    def type(self):
        supdir = self.direction if self.direction else ""

        if self.is_snubber:
            _type = supdir + "LAT" + "SNB"
        elif self.is_rotational:
            _type = supdir + "R" + "LAT"
        else:
            _type = "LAT"

        return _type


@units.define(_spring_rate="translation_stiffness", _cold_load="force",
              component_weight="force")
class Spring(Support):
    """Define a spring support.

    The hanger selection algorithm requires a hot load and movement. If the
    variability is less than 3% or greater than the input variability or 25% by
    default a constant spring is selected.

    Parameters
    ----------
    name : str
        Unique name for pipe object.

    point : Point
        Point instance where support is located.

    spring_rate : float
        Spring stiffness

    cold_load : float
        The spring load in the cold installation position.

    variability : float
        The allowable difference in cold to hot load given in percentage.
        Set to default value at 25 percent per MSS-SP-58.

    is_constant : bool
        If set to True, spring is set to constant effort.

        .. note::
            Variability is set to zero and the hot and cold loads are the same
            value. The spring_rate is also set to None since a constant force
            is applied.

    is_ftype : bool
        Is set to bottom supporting F type if set to True.

    component_weight : float
        The weight of support component such as pipe clamp and rod that should
        be accounted for in spring sizing.

    is_program_designed : bool
        Set to True by default.
    """

    def __init__(self, name, point, spring_rate=None, cold_load=None,
                 variability=25, is_constant=False, is_ftype=False,
                 component_weight=0, is_program_designed=True):

        super(Spring, self).__init__(name, point)
        self._spring_rate = spring_rate
        self._cold_load = cold_load
        self._variability = variability     # 25% per MSS-SP-58
        self._is_constant = is_constant     # true if constant effort
        self.is_ftype = is_ftype            # bottom supported
        self.component_weight = component_weight
        self._size = None    # size
        self._stype = None   # type

        self._is_program_designed = is_program_designed
        if spring_rate and cold_load:
            self._is_program_designed = False

    @classmethod
    def variable(cls, name, point, spring_rate, cold_load):
        inst = cls(name, point, spring_rate, cold_load)

        return inst

    @classmethod
    def constant(cls, name, point, load):
        inst = cls(name, point, None, load, is_constant=True)

        return inst

    @property
    def is_program_designed(self):
        return self._is_program_designed

    @is_program_designed.setter
    def is_program_designed(self, value):
        self._is_program_designed = value

    @property
    def size(self):
        return self._size

    @property
    def stype(self):
        return self._stype

    @property
    def spring_rate(self):
        if self.is_constant:
            return None
        else:
            return self._spring_rate

    @spring_rate.setter
    def spring_rate(self, value):
        if self.is_program_designed:
            raise ValueError("cannont set value for program designed spring")
        self._spring_rate = value

    @property
    def cold_load(self):
        return self._cold_load

    @cold_load.setter
    def cold_load(self, value):
        if self.is_program_designed:
            raise ValueError("cannont set value for program designed spring")
        self._cold_load = value

    @staticmethod
    def _nan_helper(y):
        """Helper for missing data used to interpolate values. Code courtesy of
        stackoverflow question 6518811.
        """
        return np.isnan(y), lambda z: z.nonzero()[0]

    @staticmethod
    def _pick_variable(udata, ldata, top_out_load_row, bottom_out_load_row,
                       kdata, hot_load, movement, sizes=None, stypes=None):
        """Helper function which returns a spring stiffness and other data if
        available.

        The algorithm finds the hot load at the middle of the table such that
        the spring does not top or bottom out. It then uses the movement to
        go across the corresponding row to determine the first bounding value.
        It uses the load row and the displacment column to determine the spring
        stiffness.

        Parameters
        ----------
        udata : numpy array
            Spring displacement data.

        ldata : numpy array
            Spring load data.

        top_out_load_row : int
            Spring load above which the spring tops out.

        bottom_out_load_row : int
            Spring load row below which the spring bottoms out.

        kdata : numpy array
            Spring stiffness data.

        hot_load : float
            The hot load to design for.

        movement : float
            The movement to design for.

        sizes : numpy array
            Sizes if available.

        stypes : numpy array
            Hanger types if available.
        """
        assert (ldata[0, 0] <= hot_load <= ldata[-1, -1]), ("hot load out of "
                "range, use rigid support instead")
        # find all applicable columns from variable table
        idxs = []
        for i, col in enumerate(ldata.transpose()): # down column
            break_outer = False
            for j, load in enumerate(col):
                if j == 0 and load >= hot_load:
                    break_outer = True
                    break

                if load >= hot_load:
                    if top_out_load_row <= j <= bottom_out_load_row:
                        idxs.append((j, i)) # store row and column
                    break

            if break_outer:     # break nested for
                break

        # print(idxs)

        # get the load closest to middle of table
        if not idxs:
            raise IndexError("spring load not found")
        else:
            dist = []
            midrow = 0.5 * (top_out_load_row + bottom_out_load_row)
            for row, _ in idxs:
                dist.append(abs(midrow - row))
            # get index of load closest to middle
            mididx = dist.index(min(dist))

        # print(mididx)

        # use load from above to find displacement
        lrow, lcol = idxs[mididx]   # load row and col
        assert (udata[lrow, -1] <= movement <= udata[lrow, 0]), ("movement "
                "out of range, use constant spring or rigid support")
        ucol = 0
        for i, mvt in enumerate(reversed(udata[lrow])): # right to left
            if movement <= mvt:
                ucol = i    # spring displacement column
                break

        # use load and displacement to find spring stiffness
        spring_rate = kdata[ucol, lcol]
        stype = stypes[ucol] if stypes.tolist() else None
        size = sizes[lcol] if sizes.tolist() else None

        return spring_rate, stype, size

    @staticmethod
    def _pick_constant(cudata, cldata, hot_load, movement, csizes):
        """Design a constant spring based on hot load and movement.

        1. Determine total travel.

            TT = AT + greater of (1.20*AT or 1")

            Round TT to nearest 1/2".

        2. Find the load that most closely bounds the hot load.

        3. Determine the corresponding constant spring size.
        """
        total_travel = max(1.2*movement, 1+movement)
        assert (cudata[0] <= total_travel <= cudata[-1]), ("movement out of",
                "range.")

        ucol = 0    # mvt column
        for i, udat in enumerate(cudata):
            if udat >= total_travel:
                ucol = i
                break

        lrow = 0    # load row
        for i, cldat in enumerate(cldata[:, ucol]):
            if cldat >= hot_load:
                lrow = i
                break

        load = cldata[lrow, ucol]
        travel = cudata[ucol]
        size = csizes[lrow]

        return load, travel, size

    @classmethod
    def from_file(cls, name, point, hot_load, movement, variability=25,
                  is_ftype=False, component_weight=0, vendor="anvil"):
        """Pick a spring from a vendor hanger table based on the hot load and
        movement.

        If the movement is small for the corresponding hot load a rigid support
        in the vertical direction is returned.

        If however, the movement is beyond the maximum range prescribed in the
        vendor table for the hot load, a constant spring load with the required
        movement is selected.
        """
        vsfname = vendor + "_variable.csv"
        csfname = vendor + "_constant.csv"
        with open(os.path.join(psi.SPRING_DIRECTORY, vsfname), "r") as vsfile:
            vsdata = np.array(list(csv.reader(vsfile)))

        with open(os.path.join(psi.SPRING_DIRECTORY, csfname), "r") as csfile:
            csdata = np.array(list(csv.reader(csfile)))

        # reset per vendor as applicable
        default_catalog_units = "english"

        if vendor == "anvil":
            stypes = np.asarray(vsdata[0, :5])
            sizes = np.asarray(vsdata[0, 5:])

            # interpolate missing spring deflection data
            udata = np.array(vsdata[1:30, :5])
            udata[udata==""] = np.nan
            udata = udata.astype(dtype=np.float)
            # fill in nan values with interpolated data col by col
            for i, col in enumerate(udata.transpose()):
                nans, x = Spring._nan_helper(col)
                udata[:, i][nans] = np.interp(x(nans), x(~nans),
                                              udata[:, i][~nans])

            ldata = np.asarray(vsdata[1:30, 5:], dtype=np.float)
            top_out_load_row = 4        # going above tops out support
            bottom_out_load_row = 24    # going below bottoms out support

            kdata = np.array(vsdata[31:36, 5:])
            kdata[kdata=='-'] = np.nan
            kdata = kdata.astype(dtype=np.float)

            # constant spring data
            csizes = np.asarray(csdata[1:, 0])
            cudata = np.asarray(csdata[0, 1:], dtype=np.float)

            cldata = np.array(csdata[1:, 1:])
            cldata[cldata==""] = np.nan
            cldata = cldata.astype(dtype=np.float)

        hot_load += component_weight

        # use load and displacement to find spring stiffness
        spring_rate, stype, size = Spring._pick_variable(udata, ldata,
                top_out_load_row, bottom_out_load_row, kdata, hot_load,
                movement, sizes=sizes, stypes=stypes)

        # calculate cold load based on the hot load
        if movement < 0:
            cold_load = hot_load - (abs(movement)*spring_rate)
        else:
            cold_load = hot_load + (abs(movement)*spring_rate)

        # TODO: calculate cold load and move to next or previous load column if
        # spring tops or bottoms out - recalculate corresponding hot load

        # actual variation between cold and hot load
        var = abs(((hot_load-cold_load) / hot_load) * 100)

        spring_selection_failed = False
        try:
            assert (var <= variability and var < 25), "variability over limit"

            # create spring instance
            with units.Units(user_units=default_catalog_units):
                inst = cls(name, point, spring_rate, cold_load,
                        variability=variability, is_constant=False,
                        is_ftype=is_ftype, component_weight=component_weight)
                inst._stype = stype
                inst._size = size
        except AssertionError:
            # design a constant spring
            load, travel, size = Spring._pick_constant(cudata, cldata,
                    hot_load, movement, csizes)

            with units.Units(user_units=default_catalog_units):
                inst = cls(name, point, None, load, variability=variability,
                        is_constant=True, is_ftype=is_ftype,
                        component_weight=component_weight)
                inst._stype = None
                inst._size = size

        return inst

    @property
    def variability(self):
        if self.is_constant:
            return 0
        else:
            return self._variability

    @variability.setter
    def variability(self, variability):
        """The allowable change in load from the hot load to cold load"""
        self._variability = variability
        if variability == 0:
            self.is_constant = True

    @property
    def is_constant(self):
        return self._is_constant

    @is_constant.setter
    def is_constant(self, value):
        self._is_constant = value

    def hot_load(self, movement):
        """Return the spring hot load based on the movement."""
        if self.is_constant:
            return self.cold_load
        else:
            if movement < 0:
                hl = (self.cold_load + self.component_weight +
                        abs(movement)*self.spring_rate)
            else:
                hl = (self.cold_load + self.component_weight -
                        abs(movement)*self.spring_rate)

        return hl

    def cglobal(self, element):
        """The spring rate of a variable spring will alter the diagonal of the
        stiffness matrix in the vertical direction (y or z). A constant spring
        will not have any stiffness effect.
        """
        c = np.zeros((12, 1), dtype=np.float64)

        vert = self.app.models.active_object.settings.vertical
        rate = self.spring_rate
        if self.is_constant:
            rate = 0

        if vert == "y":
            if self.point == element.from_point.name:
                c[1, 0] = rate
            elif self.point == element.to_point.name:
                c[7, 0] = rate

        elif vert == "z":
            if self.point == element.from_point.name:
                c[2, 0] = rate
            elif self.point == element.to_point.name:
                c[8, 0] = rate

        return c

    def fglobal(self, element):
        """For both a variable and a constant spring a load is applied to the
        point to simulate the installation load the support provides.  Note the
        sign of the load is positive (upward).
        """
        f = np.zeros((12, 1), dtype=np.float64)

        vert = self.app.models.active_object.settings.vertical
        cl = self.cold_load

        if vert == "y":
            if self.point == element.from_point.name:
                f[:6, 0] = [0, cl, 0, 0, 0, 0]
            elif self.point == element.to_point.name:
                f[6:12, 0] = [0, cl, 0, 0, 0, 0]

        elif vert == "z":
            if self.point == element.from_point.name:
                f[:6, 0] = [0, 0, cl, 0, 0, 0]
            elif self.point == element.to_point.name:
                f[6:12, 0] = [0, 0, cl, 0, 0, 0]

        return f


class SupportContainer(EntityContainer):

    def __init__(self):
        super(SupportContainer, self).__init__()
        self.Anchor = Anchor
        self.Inclined = Inclined
        self.LimitStop = LimitStop
        self.Lateral = Lateral
        self.Spring = Spring
        self.X = X; self.Y = Y; self.Z = Z

    def apply(self, supports=[], elements=[]):
        """Apply supports to elements.

        A reference of the support is attached to each element, a one to one
        assignment.

        Parameters
        ----------
        supports : list
            A list of supports

        elements : list
            A list of elements. If elements is None, loads are applied to all
            elements.
        """
        for support, element in zip(supports, elements):
            support.element = element
            element.supports.add(support)
