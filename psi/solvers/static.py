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

"""Static solver"""


import math
import logging

import numpy as np
from tqdm import trange

from psi.supports import Inclined, AbstractSupport, Spring
from psi.loadcase import LoadCase, LoadComb
from psi.loads import Displacement
from psi import units
from .codecheck import element_codecheck


def order_of_mag(n):
    """Return the order of magnitude of a number n"""
    return math.floor(math.log10(n))


def static(model):
    """Run a static analysis of the system."""

    tqdm = logging.getLogger("tqdm")
    tqdm.info("*** Preprocessing, initializing analysis...")

    # Note: All internal object data is stored in SI once loaded from an
    # external files, disabling unit conversion allows for working with
    # consistent SI
    model.app.units.disable()
    tqdm.info("*** Switching to base units.")

    # do stuff here
    tqdm.info("*** Assembling system stiffness and force matrices.")

    ndof = 6    # nodal degrees of freedom
    en = 2      # number of nodes per element

    # similar to nodal dof matrix
    points = list(model.points)
    nn = len(model.points)

    # for unsolvable loadcases
    nan = np.empty((nn*ndof, 1))
    nan[:, 0] = np.nan

    # iterate each primary loadcase
    loadcases = [lc for lc in model.loadcases if isinstance(lc, LoadCase)]

    for i, loadcase in enumerate(loadcases):
        tqdm.info("*** Solving loadcase %s" % loadcase.name)

        # global system stiffness matrix
        Ks = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)

        # global system force matrix consisting of summation of loads
        Fs = np.zeros((nn*ndof, 1), dtype=np.float64)

        # solution matrices for primary load cases
        R = np.zeros((nn*ndof, 1), dtype=np.float64)    # reactions
        Fi = np.zeros((nn*ndof, 1), dtype=np.float64)   # internal forces

        # pre-processing elements
        for element in model.elements:
            idxi = points.index(element.from_point)
            idxj = points.index(element.to_point)

            # node and corresponding dof (start, finish), used to define the
            # elements of the system stiffness and force matrices
            niqi, niqj = idxi*ndof, idxi*ndof + ndof
            njqi, njqj = idxj*ndof, idxj*ndof + ndof

            # element stiffness at room temp, conservative stresses
            keg = element.kglobal(model.settings.tref)

            # assemble global stiffness matrix, quadrant 1 to 4
            Ks[niqi:niqj, niqi:niqj] += keg[:6, :6]         # 2nd
            Ks[niqi:niqj, njqi:njqj] += keg[:6, 6:12]       # 1st
            Ks[njqi:njqj, niqi:niqj] += keg[6:12, :6]       # 3rd
            Ks[njqi:njqj, njqi:njqj] += keg[6:12, 6:12]     # 4th

            # stiffness matrix modified at nodes with supports
            # supports with displacments also modify global force vector
            for support in element.supports:
                # modify diagonal elements, penalty method, by adding large
                # stiffnesses to the diagonals where a support is located
                # if support is skewed off diagonal elements are modified too
                ksup = support.kglobal(element)

                Ks[niqi:niqj, niqi:niqj] += ksup[:6, :6]        # 2nd
                Ks[njqi:njqj, njqi:njqj] += ksup[6:12, 6:12]    # 4th

                # large stiffness added to each force matrix component with
                # non-zero support displacements based on loadcase
                # NOTE: this only applies to supports with displacement, for
                # all others dsup is 0 and so added equivalent force is also 0
                csup = support.cglobal(element)
                dsup = support.dglobal(element, loadcase)   # support disp

                Fs[niqi:niqj, 0] += (csup[:6, 0] * dsup[:6, 0])
                Fs[njqi:njqj, 0] += (csup[6:12, 0] * dsup[6:12, 0])

                # add cold load if spring support or constant
                # note spring rate for variables is added above
                if isinstance(support, Spring):
                    fspring = support.fglobal(element)

                    Fs[niqi:niqj, 0] += fspring[:6, 0]
                    Fs[njqi:njqj, 0] += fspring[6:12, 0]

            # load vector summed at nodes with forces
            feg = np.zeros((en*ndof, 1), dtype=np.float64)
            for load in element.loads:
                if load in loadcase:
                    feg += load.fglobal(element)

                    # modify stiffness matrix for nodal displacement
                    if isinstance(load, Displacement):
                        kforc = load.kglobal(element)

                        Ks[niqi:niqj, niqi:niqj] += kforc[:6, :6]
                        Ks[njqi:njqj, njqi:njqj] += kforc[6:12, 6:12]
                else:
                    # otherwise print a warning
                    pass
            # assemble global system force matrix per loadcase
            Fs[niqi:niqj, 0] += feg[:6, 0]
            Fs[njqi:njqj, 0] += feg[6:12, 0]

        tqdm.info("    --> Solving system equations for displacements.")
        if model.settings.weak_springs:
            tqdm.info("    --> Turning on weak springs.")

            di = np.diag_indices_from(Ks)   # Ks is square

            # support stiffness is reduced by 75% order of magnitude
            trans_order = order_of_mag(model.settings.translation_stiffness) * 0.75
            rota_order = order_of_mag(model.settings.rotation_stiffness) * 0.75

            weak_spring_arr = np.zeros(Ks.shape[0], dtype=np.float64)
            weak_spring_arr[::3] = model.settings.translation_stiffness * 10**-trans_order
            weak_spring_arr[3::3] = model.settings.rotation_stiffness * 10**-rota_order

            # apply small stiffness to diagonals for numerical stability
            Ks[di] += weak_spring_arr

        try:
            # all supports are initially assumed linear
            Xs = np.linalg.solve(Ks, Fs)
        except np.linalg.LinAlgError:
            loadcase.movements.results = nan
            loadcase.reactions.results = nan
            loadcase.forces.results = nan

            if not model.settings.weak_springs:
                raise np.linalg.LinAlgError("a solver error occured, ",
                                            "try turning on weak springs.")
            else:
                raise np.linagl.LinAlgError("a solver error occured, ",
                                            "check model for errors.")

        # one directional (i.e. gaps) and/or friction is nonlinear
        nonlinear = {support: "active" for support in model.supports
                     if isinstance(support, AbstractSupport)
                     and support.is_nonlinear}  # only handles gaps
        if nonlinear:
            converge_status = [False] * len(nonlinear)

            # stiffness matrix and force vector modified at nodes with supports
            t = trange(model.settings.nonlinear_iteration)
            for _ in t:
                t.set_description("Iterating loadcase %s..." % loadcase.name)
                t.refresh()

                # iterate non-linear supports to determine final state
                for i, support in enumerate(nonlinear.keys()):
                    element = support.element

                    idxi = points.index(element.from_point)
                    idxj = points.index(element.to_point)

                    # node and corresponding dof (start, finish), used to
                    # define the elements of the system stiffness and force
                    # matrices
                    niqi, niqj = idxi*ndof, idxi*ndof + ndof
                    njqi, njqj = idxj*ndof, idxj*ndof + ndof

                    # calculate support reactions for non-linear supports
                    if isinstance(support, Inclined):
                        csup = support.cglobal(element)             # stiffness
                        dsup = support.dglobal(element, loadcase)   # displacement
                        B1, B2, B3 = support.dircos
                        B = np.array([[B1, B2, B3]], dtype=np.float64)
                        B2 = np.append(B, B)

                        R[niqi:niqj, 0] += (-(csup[:6, 0]*B2) *
                                            ((Xs[niqi:niqj, 0].dot(B2)) -
                                                dsup[:6, 0]))
                        R[njqi:njqj, 0] += (-(csup[6:12, 0]*B2) *
                                            ((Xs[njqi:njqj, 0].dot(B2)) -
                                                dsup[6:12, 0]))
                    else:
                        csup = support.cglobal(element)             # stiffness
                        dsup = support.dglobal(element, loadcase)   # displacement

                        R[niqi:niqj, 0] += (-csup[:6, 0] *
                                            (Xs[niqi:niqj, 0] - dsup[:6, 0]))
                        R[njqi:njqj, 0] += (-csup[6:12, 0] *
                                            (Xs[njqi:njqj, 0] - dsup[6:12, 0]))

                    # support can be "+/-", "+", or "-"
                    # note that a dircos always point in the positive direction
                    # the solver is linear, assumes bi-directional supports,
                    # i.e. restricts movement in both directions

                    suppsign = support.direction    # may be +, - or None
                    suppdir = np.array(support.dircos, dtype=np.float64)

                    # extract reaction, displacement and set gap
                    gsup = np.zeros((12, 1), dtype=np.float64)
                    gapvec = support.gap * np.array(support.dircos,
                                                    dtype=np.float64)
                    if support.point == element.from_point.name:
                        if support.is_rotational is False:
                            reac = -R[niqi:niqj, 0][:3]
                            disp = Xs[niqi:niqj, 0][:3]
                            gsup[:3, 0] = gapvec
                        else:
                            reac = -R[niqi:niqj, 0][3:6]
                            disp = Xs[niqi:niqj, 0][3:6]
                            gsup[3:6, 0] = gapvec
                    elif support.point == element.to_point.name:
                        if support.is_rotational is False:
                            reac = -R[njqi:njqj, 0][:3]
                            disp = Xs[njqi:njqj, 0][:3]
                            gsup[6:9, 0] = gapvec
                        else:
                            reac = -R[njqi:njqj, 0][3:6]
                            disp = Xs[njqi:njqj, 0][3:6]
                            gsup[9:12, 0] = gapvec

                    # component in support direction
                    reac = np.dot(reac, suppdir)
                    disp = np.dot(disp, suppdir)
                    gap = np.dot(gapvec, suppdir)

                    if suppsign is None:
                        # bi-directional supports components always against
                        # support direction
                        reac = -abs(reac)
                        disp = -abs(disp)
                        gap = abs(gap)
                    elif suppsign == "-":
                        reac *= -1
                        disp *= -1
                        gap *= -1

                    # check support reaction in support direction
                    if nonlinear[support] == "active":
                        if reac < 0 and support.gap == 0:
                            nonlinear[support] = "active"
                            converge_status[i] = True
                        else:
                            # remove inactive support
                            nonlinear[support] = "inactive"
                            converge_status[i] = False

                            ksup = support.kglobal(element)
                            Ks[niqi:niqj, niqi:niqj] -= ksup[:6, :6]      # 2nd
                            Ks[njqi:njqj, njqi:njqj] -= ksup[6:12, 6:12]  # 4th

                            csup = support.cglobal(element)
                            dsup = support.dglobal(element, loadcase) + gsup

                            Fs[niqi:niqj, 0] -= (csup[:6, 0] * dsup[:6, 0])
                            Fs[njqi:njqj, 0] -= (csup[6:12, 0] * dsup[6:12, 0])

                    # check support displacement in support direction
                    elif nonlinear[support] == "inactive":
                        if disp >= 0:
                            nonlinear[support] = "inactive"
                            converge_status[i] = True
                        else:
                            # add active support
                            nonlinear[support] = "active"
                            converge_status[i] = False

                            ksup = support.kglobal(element)
                            Ks[niqi:niqj, niqi:niqj] += ksup[:6, :6]      # 2nd
                            Ks[njqi:njqj, njqi:njqj] += ksup[6:12, 6:12]  # 4th

                            csup = support.cglobal(element)
                            dsup = support.dglobal(element, loadcase) + gsup

                            Fs[niqi:niqj, 0] += (csup[:6, 0] * dsup[:6, 0])
                            Fs[njqi:njqj, 0] += (csup[6:12, 0] * dsup[6:12, 0])

                if all(converge_status):    # check convergence
                    break

                try:
                    R[:, 0] = 0                     # R cleared
                    Xs = np.linalg.solve(Ks, Fs)    # Ks and Fs modified
                except np.linalg.LinAlgError:
                    loadcase.movements.results = nan
                    loadcase.reactions.results = nan
                    loadcase.forces.results = nan

                    raise np.linagl.LinAlgError("a solver error occured, ",
                                                "check model for errors.")
            else:
                loadcase.movements.results = nan
                loadcase.reactions.results = nan
                loadcase.forces.results = nan

                tqdm.info("    --> Loadcase %s failed to converge." % loadcase.name)

        tqdm.info("    --> Post processing elements...")
        tqdm.info("    --> Calculating support reactions and internal forces""")
        R[:, 0] = 0     # reset reaction vector
        for element in model.elements:
            idxi = points.index(element.from_point)
            idxj = points.index(element.to_point)

            # node and corresponding dof (start, finish)
            niqi, niqj = idxi*ndof, idxi*ndof + ndof
            njqi, njqj = idxj*ndof, idxj*ndof + ndof

            # element local stiffness and transformation matrix
            kel = element.klocal(model.settings.tref)
            T = element.T()

            for support in element.supports:
                if support not in nonlinear or nonlinear[support] == "active":
                    if isinstance(support, Inclined):
                        csup = support.cglobal(element)     # stiffness
                        dsup = support.dglobal(element, loadcase)   # displacement
                        B1, B2, B3 = support.dircos
                        B = np.array([[B1, B2, B3]], dtype=np.float64)
                        B2 = np.append(B, B)

                        R[niqi:niqj, 0] += (-(csup[:6, 0]*B2) *
                                            ((Xs[niqi:niqj, 0].dot(B2)) -
                                                dsup[:6, 0]))
                        R[njqi:njqj, 0] += (-(csup[6:12, 0]*B2) *
                                            ((Xs[njqi:njqj, 0].dot(B2)) -
                                                dsup[6:12, 0]))
                    else:
                        csup = support.cglobal(element)     # stiffness
                        dsup = support.dglobal(element, loadcase)   # displacement

                        R[niqi:niqj, 0] += (-csup[:6, 0] *
                                            (Xs[niqi:niqj, 0] - dsup[:6, 0]))
                        R[njqi:njqj, 0] += (-csup[6:12, 0] *
                                            (Xs[njqi:njqj, 0] - dsup[6:12, 0]))
                else:
                    # inactive nonlinear support
                    R[niqi:niqj, 0] += 0
                    R[njqi:njqj, 0] += 0

            # calculate element local forces and moments using the local
            # element stiffness matrix and fi = kel*x where x is the local
            # displacement vector given by (T * _x)
            _x = np.zeros((12, 1), dtype=np.float64)
            _x[:6, 0] = Xs[niqi:niqj, 0]
            _x[6:12, 0] = Xs[njqi:njqj, 0]

            # x is the element local displacements
            fi = kel @ (T @ _x)

            # internal force and moment matrix at node i and j
            Fi[niqi:niqj, 0] = fi[:6, 0]
            Fi[njqi:njqj, 0] = fi[6:12, 0]

        tqdm.info("    --> Writing loadcase results data.")
        # X[:, 0], R[:, 0] and Fi[:, 0] return row vectors
        loadcase.movements.results = Xs[:, 0]
        loadcase.reactions.results = -R[:, 0]   # action force on support
        loadcase.forces.results = Fi[:, 0]

    # switch back to user units - analysis is complete
    model.app.units.enable()

    tqdm.info("*** Performing element code checking.")

    S = np.zeros((nn, 10), dtype=np.float64)        # stresses
    for i, loadcase in enumerate(loadcases):
        C = []  # code used at each node
        for element in model.elements:
            element_codecheck(points, loadcase, element, S)
            C.extend(2 * [element.code.label])
        loadcase.stresses.results = (S[:, :], C)

    tqdm.info("*** Code checking complete.")
    tqdm.info("*** Analysis complete!\n")
