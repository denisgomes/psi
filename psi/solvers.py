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

"""Static, modal and dynamic solvers."""

import math
import sys
import logging
from contextlib import redirect_stdout

import numpy as np

from psi.loadcase import LoadCase, LoadComb
from psi import units


def order_of_mag(n):
    """Return the order of magnitude of a number n"""
    return math.floor(math.log10(n))


def calculate_element_code_stresses(points, loadcase, element, S, lcasenum):
    ndof = 6

    idxi = points.index(element.from_point)
    idxj = points.index(element.to_point)

    # node and corresponding dof (start, finish)
    niqi, niqj = idxi*ndof, idxi*ndof + ndof
    njqi, njqj = idxj*ndof, idxj*ndof + ndof

    with units.Units(user_units="code_english"):
        # Note: units are changed to code_english for the moments
        # to have units of inch*lbf per code requirement
        # code equations are units specific, i.e. imperial or si

        fori = loadcase.forces.results[niqi:niqj, 0]
        forj = loadcase.forces.results[njqi:njqj, 0]

        # code stresses per element node i and j for each loadcase
        shoop = element.code.shoop(element, loadcase)

        saxi = element.code.sax(element, loadcase, fori)
        saxj = element.code.sax(element, loadcase, forj)

        stori = element.code.stor(element, fori)
        storj = element.code.stor(element, forj)

        # pressure stress is same at both nodes
        slp = element.code.slp(element, loadcase)

        slbi = element.code.slb(element, fori)
        slbj = element.code.slb(element, forj)

        # total code stress
        sli = element.code.sl(element, loadcase, fori)
        slj = element.code.sl(element, loadcase, forj)

        # fitting and nodal sifs, sum together, take max or average?
        sifi = element.code.sifi(element)
        sifo = element.code.sifo(element)

        sallowi = element.code.sallow(element, loadcase, fori)
        sallowj = element.code.sallow(element, loadcase, forj)

        try:
            sratioi = sli / sallowi  # code ratio at node i
            sratioj = slj / sallowj  # code ratio at node j
        except ZeroDivisionError:
            sratioi = 0
            sratioj = 0

        # hoop, sax, stor, slp, slb, sl, sifi, sifj, sallow, ir
        # take the worst code stress at node
        if sratioi > S[idxi, -1, lcasenum]:
            S[idxi, :10, lcasenum] = (shoop, saxi, stori, slp, slbi, sli,
                                      sifi, sifo, sallowi, sratioi)

        if sratioj > S[idxj, -1, lcasenum]:
            S[idxj, :10, lcasenum] = (shoop, saxj, storj, slp, slbj, slj,
                                      sifi, sifo, sallowj, sratioj)

        # with redirect_stdout(sys.__stdout__):
        #     print(S)

        # TODO : Implement Ma, Mb and Mc calculation loads
        # for each loadcase where Ma is for sustained, Mb is
        # for occasional and Mc is for expansion type loads
        # This applies to code stress calculations only


def static(model):
    """Run a static analysis of the system.

    Go through all loadcases, solve each one and then combine them to generate
    the final results.

    For each element calculate the element stiffness matrix and assemble the
    system stiffness matrix using the nodal degree of freedom matrix.

    The loads defined for each element loadcase is combined into a single load
    vector. The load vector for each element loadcase is then assembled into a
    global load vector. Multiple global load vectors are solved for at once
    using gaussian elimination.

    Note that to simplify the FEA solution, all elements are essentially beams,
    including bends and reducers which are approximations made by chaining
    multiple beams together.

    1. Create NDOF matrix.

    2. For each element
        a. Construct the local stiffness matrix.
        b. Construct local force matrix, one for each load case.
        c. Transform local stiffness and force matrix.
        d. Add global element stiffness matrix and the global element force
           vector to the global system stiffness and force vectors,
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

    Non-linear supports

    For each non-linear support and each loadcase:

    Initially a non-linear (+Y) support is converted to a linear full Y
    support.

    If the solution shows a (-Y) load on the support, the support is modeled
    correctly. If not, the program gets rid of the support altogether and
    tracks the displacement at that point.

    If the point shows a (+Y) deflection, this support is modeled correclty for
    this loadcase. Else, the stiffness matrix is reset to a full Y.

    Loading sequence and non-linear supports.

    Non-linear supports use an iterative approach to determine the final
    support loads and pipe configuration. The sequence in which loads are
    applied matters as superposition is not valid for non-linear analysis.
    Each load is applied and the displacements extracted. These movements are
    then used for the next step, ie. the model mesh is modified to incorporate
    the new point locations. As a result the system stiffness matrix is updated
    each iteration.
    """
    tqdm = logging.getLogger("tqdm")
    tqdm.info("*** Preprocessing, initializing analysis...")

    # Note: All internal object data is stored in SI once loaded from an
    # external files, disabling unit consersion allows for working with
    # only SI
    model.app.units.disable()
    tqdm.info("*** Switching to base units.")

    # do stuff here
    tqdm.info("*** Assembling system stiffness and force matrices.")

    ndof = 6    # nodal degrees of freedom
    en = 2      # number of nodes per element
    nn = len(model.points)
    lc = len(model.loadcases)

    # similar to nodal dof matrix
    points = list(model.points)

    # global system stiffness matrix
    Ks = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)

    # global system force matrix, one loadcase per column
    # a loadcase consists of one or more loads
    Fs = np.zeros((nn*ndof, lc), dtype=np.float64)

    # pre-processing elements
    for element in model.elements:
        idxi = points.index(element.from_point)
        idxj = points.index(element.to_point)

        # node and corresponding dof (start, finish), used to define the
        # elements of the system stiffness and force matrices
        niqi, niqj = idxi*ndof, idxi*ndof + ndof
        njqi, njqj = idxj*ndof, idxj*ndof + ndof

        # element stiffness at room temp, conservative stresses
        keg = element.kglobal(294.2611)

        # with redirect_stdout(sys.__stdout__):
        #     print(keg)

        # assemble global stiffness matrix, quadrant 1 to 4
        Ks[niqi:niqj, niqi:niqj] += keg[:6, :6]         # 2nd
        Ks[niqi:niqj, njqi:njqj] += keg[:6, 6:12]       # 1st
        Ks[njqi:njqj, niqi:niqj] += keg[6:12, :6]       # 3rd
        Ks[njqi:njqj, njqi:njqj] += keg[6:12, 6:12]     # 4th

        # with redirect_stdout(sys.__stdout__):
        #     print(Ks)

        # modify diagonal elements, penalty method, by adding large
        # stiffnesses to the diagonals where a support is located
        di = np.diag_indices(6)     # diagonal indices for 6x6 matrix
        for support in element.supports:
            ksup = support.kglobal(element)
            # assign to diagonal elements of system matrix
            Ks[niqi:niqj, niqi:niqj][di] += ksup[:6, 0]     # 2nd
            Ks[njqi:njqj, njqi:njqj][di] += ksup[6:12, 0]   # 4th

        # with redirect_stdout(sys.__stdout__):
        #     print(Ks)

        # iterate each loadcase adding loads
        for i, loadcase in enumerate(model.loadcases):
            # with redirect_stdout(sys.__stdout__):
            #     print(loadcase)

            # sum of all loads in a primary/primitive loadcase
            if isinstance(loadcase, LoadCase):
                feg = np.zeros((en*ndof, 1), dtype=np.float64)

                for loadtype, opercase in zip(loadcase.loadtypes,
                                              loadcase.opercases):

                    for load in element.loads:
                        if (isinstance(load, loadtype) and
                                load.opercase == opercase):
                            feg += load.fglobal(element)
                        else:
                            pass
                            # otherwise print a warning

                # assemble global system force matrix per loadcase
                Fs[niqi:niqj, i] += feg[:6, 0]
                Fs[njqi:njqj, i] += feg[6:12, 0]

            # large stiffness added to each force matrix component with
            # non-zero support displacements, again per loadcase
            # NOTE: this only applies to Displacement supports, for all
            # others dsup is 0 and so added stiffness is also 0
            for support in element.supports:

                ksup = support.kglobal(element)
                dsup = support.dglobal(element)     # support displacement

                Fs[niqi:niqj, i] += (ksup[:6, 0] * dsup[:6, 0])
                Fs[njqi:njqj, i] += (ksup[6:12, 0] * dsup[6:12, 0])

    # with redirect_stdout(sys.__stdout__):
    #     print(Ks)

    tqdm.info("*** Solving system equations for displacements.")

    if model.settings.weak_springs:
        tqdm.info("*** Turning on weak springs.")

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
        X = np.linalg.solve(Ks, Fs)
    except np.linalg.LinAlgError:
        if not model.settings.weak_springs:
            raise np.linalg.LinAlgError("Solver error, try turning on",
                                        "weak springs")
        raise

    # with redirect_stdout(sys.__stdout__):
    #     print(X)

    tqdm.info("*** Post processing elements...")
    R = np.zeros((nn*ndof, lc), dtype=np.float64)   # reactions
    Fi = np.zeros((nn*ndof, lc), dtype=np.float64)  # internal forces

    tqdm.info("*** Calculating support reactions and internal forces.")
    for element in model.elements:
        idxi = points.index(element.from_point)
        idxj = points.index(element.to_point)

        # node and corresponding dof (start, finish)
        niqi, niqj = idxi*ndof, idxi*ndof + ndof
        njqi, njqj = idxj*ndof, idxj*ndof + ndof

        # element local stiffness and transformation matrix
        kel = element.klocal(294.2611)

        T = element.T()

        # nodal displacement vector per loadcase
        for i, loadcase in enumerate(model.loadcases):
            # reaction forces and moments
            for support in element.supports:
                ksup = support.kglobal(element)
                dsup = support.dglobal(element)     # support displacement

                R[niqi:niqj, i] += (-ksup[:6, 0] * (X[niqi:niqj, i] -
                                                    dsup[:6, 0]))
                R[njqi:njqj, i] += (-ksup[6:12, 0] * (X[njqi:njqj, i] -
                                                      dsup[6:12, 0]))

            # calculate element local forces and moments using the local
            # element stiffness matrix and fi = kel*x where x is the local
            # displacement vector given by (T * _x)
            _x = np.zeros((12, 1), dtype=np.float64)
            _x[:6, 0] = X[niqi:niqj, i]
            _x[6:12, 0] = X[njqi:njqj, i]

            # x is the element local displacements
            fi = kel @ (T @ _x)

            # internal force and moment matrix at node i and j
            Fi[niqi:niqj, i] = fi[:6, 0]
            Fi[njqi:njqj, i] = fi[6:12, 0]

    tqdm.info("*** Writing loadcase results data.")
    for i, loadcase in enumerate(model.loadcases):
        # load combs are combined later
        if isinstance(loadcase, LoadCase):
            # X[:, i], R[:, i] and Fi[:, i] return row vectors
            loadcase.movements.results = X[:, i]
            loadcase.reactions.results = R[:, i]
            loadcase.forces.results = Fi[:, i]

    # with redirect_stdout(sys.__stdout__):
    #     print(X)
    #     print(R)
    #     print(Fi)

    # switch back to user units - analysis is complete
    model.app.units.enable()

    tqdm.info("*** Element code checking.")

    S = np.zeros((nn, 10, lc), dtype=np.float64)
    C = []  # code listing per element node

    for element in model.elements:

        for i, loadcase in enumerate(model.loadcases):
            if isinstance(loadcase, LoadCase):
                calculate_element_code_stresses(points, loadcase, element, S, i)
                C.extend(2 * [element.code.name])

        for i, loadcase in enumerate(model.loadcases):
            if isinstance(loadcase, LoadComb) and loadcase.method == "algebraic":
                # load combs are combined at runtime except below for special
                # case where combination is algebraic stresses are combined
                # based on combined element forces
                pass

    tqdm.info("*** Code checking complete.")

    tqdm.info("*** Writing code checking results data.")
    for i, loadcase in enumerate(model.loadcases):

        if isinstance(loadcase, LoadCase):
            # X[:, i], R[:, i] and Fi[:, i] return row vectors
            loadcase.stresses.results = (S[:, :, i], C)

        elif isinstance(loadcase, LoadComb) and loadcase.method == "algebraic":
            # load combs are combined at runtime except below
            # for special case where combination is algebraic
            pass

    tqdm.info("*** Analysis complete!\n\n\n")


def modal(model):
    """Run a modal analysis of the system extracting frequencies and mode
    shapes. Solving the eigenvalue problem.

    All nonlinearity is neglected for this type of analysis. In other words,
    all support types are considered double acting and friction is not
    considered.
    """
    pass
