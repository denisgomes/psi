# Pipe Stress Infinity (PSI) - The pipe stress analysis and design software.
# Copyright (c) 2020 Denis Gomes

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

"""Static, modal and dynamic solvers.

Static
======
Go through all loadcases, solve each one and then combine them to generate the
final results.

For each element calculate the element stiffness matrix and assemble the system
stiffness matrix using the nodal degree of freedom. The DOFs of a node is based
on the position/index of the point in the points list

The loads defined for each element primary/primitive loadcase is combined into
a single load vector. The load vector for each element loadcase is then
assembled into a global load matrix. Multiple global load vectors are solved
for at once using gaussian elimination.

Note that to simplify the FEA solution, all elements are simple beams,
including bends and reducers which are approximated by multiple beams chained
together.

General Procedure
-----------------
1. Define NDOF indices for element nodes.

2. For each element
    a. Construct the local stiffness matrix.
    b. Construct local force vector, one for each primitive loadcase.
    c. Transform local stiffness and force matrix.
    d. Add global element stiffness matrix and the global element force vectors
       to the global system stiffness and forces matrix, respectively.

3. Apply the boundary conditions using the penalty method.

4. Solve the global system using guassian elimination techniques.

   AX=B

   Where A is the global system matrix.

   B is the matrix of force vectors, one for each primitive loadcase, and
   X is the matrix of solutions vectors for each primitive loadcase.

5. Use the calculated displacements for each primitive loadcase to calculate
   the element force and moments.

6. Use the results from step 5 to calculate nodal reactions and element forces.
   Finally use the element forces to calculate code stresses.

Reducers and Bends
------------------
Both reducers and bends are topologically curves and therefore approximations.
The midpoint of a bend is created as a side effect of defining the bend when
it is created. All the other vertices and inaccessible from a nodal standpoint.
In other words, the user does not have control of the other underlying vertices
because only the vertex that corresponds to the midnode is referenced by the
midpoint.

The static solver preprocesses all reducers and bends such that a point is made
for each underlying vertex at runtime. Once the solution is generated all the
temporary points are deleted along with the nodal data corresponding to each
point. For the case of a bend, the solution for the midpoint and end points are
kept. For the reducer, the results for the end points are kept only similar to
all the other element.

Tee Flexibility
---------------
...

Skewed Supports
---------------
Constraint equations.

Spring Algorithm
----------------
The program tries to select a rigid support, variable or constant spring where
a spring support is specified by the user. The algorithm must satisfy all the
operating cases such that it does not bottom or top out. For each operating
case a linear deadweight analysis (W+P+D+F) is performed and a +Y support is
placed at each supports. Remove the force from the weight case at each spring
support location with the equivalent upward force and run the operating case
with the temperature. Determine the upward movement of the pipe at the spring
support. Use the applied force (i.e. the hot load) and the movement to pick a
hanger from the manufacturer catalog.

Master Slave (CNode)
--------------------
Constraint equations.

Non-linear supports
-------------------
For each loadcase, for each nonlinear support, per each solve:

Initially a non-linear (+Y) support is converted to a linear full Y support. If
the solution shows a (-Y) load on the support, the support is modeled correctly
and "active" for the loadcase. If not, the program marks it as "inactive" and
gets rid of the support altogether and tracks the displacement at that point.
If the point shows a (+Y) deflection, this support is modeled correclty for
this loadcase. Else, the stiffness matrix is reset to a full Y. If any of the
nonlinear assumption proves incorrect reanalyze with the updated stiffness and
force vector (if required). Continue the process checking all the "active"
supports for loads and all the "inactive" supports for displacement. When all
"active" supports have a (-Y) and the "inactive" suports show a positive up
displacement. A status change flag variable is used to indicate a support that
that initially shows a (-Y) load but then starts to show an uplift.

Friction
--------
...

Loading sequence and non-linear supports
----------------------------------------
Non-linear supports use an iterative approach to determine the final support
loads and pipe configuration. The sequence in which loads are applied matters
as superposition is not valid for non-linear analysis. Each load is applied
and the displacements extracted. These movements are then used for the next
step, ie. the model mesh is modified to incorporate the new point locations.
As a result the system stiffness matrix is updated each iteration.
"""

from collections import defaultdict, Counter
from itertools import zip_longest
import math
import sys
import logging
from contextlib import redirect_stdout
from functools import partial

import numpy as np
from scipy.sparse import linalg as splinalg

from psi.loadcase import LoadCase, LoadComb
from psi import units


def order_of_mag(n):
    """Return the order of magnitude of a number n"""
    return math.floor(math.log10(n))


def element_codecheck(points, loadcase, element, S, lcasenum):
    """Perform element code checking based on the assigned code for primitive
    load cases and load combinations.
    """
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

        if isinstance(loadcase, LoadCase):
            fori = loadcase.forces.results[niqi:niqj, 0]
            forj = loadcase.forces.results[njqi:njqj, 0]

            # code stresses per element node i and j for each loadcase
            shoop = element.code.shoop(element, loadcase)

            # pressure stress is same at both nodes
            slp = element.code.slp(element, loadcase)

            saxi = element.code.sax(element, fori)
            saxj = element.code.sax(element, forj)

            stori = element.code.stor(element, fori)
            storj = element.code.stor(element, forj)

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
                sratioi = sli / sallowi     # code ratio at node i
                sratioj = slj / sallowj     # code ratio at node j
            except ZeroDivisionError:
                sratioi = 0
                sratioj = 0

        elif isinstance(loadcase, LoadComb):
            # fitting and nodal sifs, sum together, take max or average?
            sifi = element.code.sifi(element)
            sifo = element.code.sifo(element)

            loadcomb = loadcase
            shoop_list, slp_list = [], []
            saxi_list, stori_list, slbi_list, sli_list = [], [], [], []
            saxj_list, storj_list, slbj_list, slj_list = [], [], [], []
            for factor, loadcase in zip_longest(loadcomb.factors,
                                                loadcomb.loadcases,
                                                fillvalue=1):
                fori = loadcase.forces.results[niqi:niqj, 0]
                forj = loadcase.forces.results[njqi:njqj, 0]

                shoop_list.append(factor * element.code.shoop(element, loadcase))
                slp_list.append(factor * element.code.slp(element, loadcase))

                saxi_list.append(factor * element.code.sax(element, fori))
                saxj_list.append(factor * element.code.sax(element, forj))

                stori_list.append(factor * element.code.stor(element, fori))
                storj_list.append(factor * element.code.stor(element, forj))

                slbi_list.append(factor * element.code.slb(element, fori))
                slbj_list.append(factor * element.code.slb(element, forj))

                # total code stress
                sli_list.append(element.code.sl(element, loadcase, fori))
                slj_list.append(element.code.sl(element, loadcase, forj))

            if loadcomb.method == "scaler":
                shoop = sum(shoop_list)
                slp = sum(slp_list)

                saxi = sum(saxi_list)
                saxj = sum(saxj_list)

                stori = sum(stori_list)
                storj = sum(storj_list)

                slbi = sum(slbi_list)
                slbj = sum(slbj_list)

                sli = sum(sli_list)
                slj = sum(slj_list)

            elif loadcomb.method == "algebraic":
                # stress per algebraic combination of forces
                fori = loadcomb.forces.results[niqi:niqj, 0]
                forj = loadcomb.forces.results[njqi:njqj, 0]

                shoop = sum(shoop_list)
                slp = sum(slp_list)

                saxi = element.code.sax(element, fori)
                saxj = element.code.sax(element, forj)

                stori = element.code.stor(element, fori)
                storj = element.code.stor(element, forj)

                slbi = element.code.slb(element, fori)
                slbj = element.code.slb(element, forj)

                # total code stress
                sli = element.code.sl(element, loadcase, fori)
                slj = element.code.sl(element, loadcase, forj)

            elif loadcomb.method == "srss":
                # note: sign of factor has no effect, always positive
                shoop = sum([s**2 for s in shoop_list])
                slp = sum([s**2 for s in slp_list])

                saxi = sum([s**2 for s in saxi_list])
                saxj = sum([s**2 for s in saxj_list])

                stori = sum([s**2 for s in stori_list])
                storj = sum([s**2 for s in storj_list])

                slbi = sum([s**2 for s in slbi_list])
                slbj = sum([s**2 for s in slbj_list])

                sli = sum([s**2 for s in sli_list])
                slj = sum([s**2 for s in slj_list])

            elif loadcomb.method == "abs":
                shoop = sum([abs(s) for s in shoop_list])
                slp = sum([abs(s) for s in slp_list])

                saxi = sum([abs(s) for s in saxi_list])
                saxj = sum([abs(s) for s in saxj_list])

                stori = sum([abs(s) for s in stori_list])
                storj = sum([abs(s) for s in storj_list])

                slbi = sum([abs(s) for s in slbi_list])
                slbj = sum([abs(s) for s in slbj_list])

                sli = sum([abs(s) for s in sli_list])
                slj = sum([abs(s) for s in slj_list])

            elif loadcomb.method == "signmax":
                shoop = max(shoop_list)
                slp = max(slp_list)

                saxi = max(saxi_list)
                saxj = max(saxj_list)

                stori = max(stori_list)
                storj = max(storj_list)

                slbi = max(slbi_list)
                slbj = max(slbj_list)

                sli = max(sli_list)
                slj = max(slj_list)

            elif loadcomb.method == "signmin":
                shoop = min(shoop_list)
                slp = min(slp_list)

                saxi = min(saxi_list)
                saxj = min(saxj_list)

                stori = min(stori_list)
                storj = min(storj_list)

                slbi = min(slbi_list)
                slbj = min(slbj_list)

                sli = min(sli_list)
                slj = min(slj_list)

            # take the sqrt last
            if loadcomb.method == "srss":
                shoop = math.sqrt(shoop)
                slp = math.sqrt(slp)

                saxi = math.sqrt(saxi)
                saxj = math.sqrt(saxj)

                stori = math.sqrt(stori)
                storj = math.sqrt(storj)

                slbi = math.sqrt(slbi)
                slbj = math.sqrt(slbj)

                sli = math.sqrt(sli)
                slj = math.sqrt(slj)

            # allowable loadcomb stress
            sallowi_list = []
            sallowj_list = []
            # determine the loadcomb code stress allowable
            for loadcase in loadcomb.loadcases:
                stype = loadcase.stype              # save type
                loadcase.stype = loadcomb.stype     # change to loadcomb type

                fori = loadcase.forces.results[niqi:niqj, 0]
                forj = loadcase.forces.results[njqi:njqj, 0]

                # calculate loadcomb allowable
                sallowi = element.code.sallow(element, loadcase, fori)
                sallowi_list.append(sallowi)

                sallowj = element.code.sallow(element, loadcase, forj)
                sallowj_list.append(sallowj)

                # revert to loadcase stype
                loadcase.stype = stype
            sallowi = min(sallowi_list)
            sallowj = min(sallowj_list)

            try:
                sratioi = sli / sallowi     # code ratio at node i
                sratioj = slj / sallowj     # code ratio at node j
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

        # TODO : Implement Ma, Mb and Mc calculation loads
        # for each loadcase where Ma is for sustained, Mb is
        # for occasional and Mc is for expansion type loads
        # This applies to code stress calculations only


def static(model):
    """Run a static analysis of the system."""

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
            raise np.linalg.LinAlgError("solver error, try turning on",
                                        "weak springs.")
        else:
            raise np.linagl.LinAlgError("solver error, make sure the model",
                                        "has no error.")

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
            loadcase.reactions.results = -R[:, i]   # action force on support
            loadcase.forces.results = Fi[:, i]

    # switch back to user units - analysis is complete
    model.app.units.enable()

    tqdm.info("*** Element code checking.")
    S = np.zeros((nn, 10, lc), dtype=np.float64)
    C = []  # code per element node
    for element in model.elements:
        for i, loadcase in enumerate(model.loadcases):
            element_codecheck(points, loadcase, element, S, i)
            C.extend(2 * [element.code.label])
            loadcase.stresses.results = (S[:, :, i], C)

    tqdm.info("*** Code checking complete.")
    tqdm.info("*** Analysis complete!\n")


def modal(model, nmodes=3):
    """Run a modal analysis of the system extracting frequencies and mode
    shapes. Solving the eigenvalue problem.

    All nonlinearity is neglected for this type of analysis. In other words,
    all support types are considered double acting and friction is not
    considered.

    The global mass and stiffness matrix are first calculated. Then all the
    rows and columns corresponding to a support location and direction are
    removed (i.e. effectively fixed). This is true for supports with gaps and
    even for springs. Finally, the dynamical matrix is calculated and the
    eigenvalue solution is solved.

    For each mode the following table is stored:

    Mode #, Mode (Freq), Mode(Period), Participation Factor (X, Y, Z),
    Mass Participation (X, Y, Z)

    Total captured modal mass (X, Y, Z)
    Total system pass

    The participation factors is directly from the eigvector.
    The mass participation is calculated by multiplying by system Mass matrix.

    From 'Theory of Matrix Structural Analysis' by J.S. Przemienieski.

    Default internal units: m, kg, sec

    Dynamical matrix

    D = np.linalg.inv(Ks) @ Ms
    eigvals, eigvecs = np.linalg.eig(D)

    inveigval - array of inverse eigenvalues (1 / omega)
    eigvecmat - matrix of eigenvectors where the column [:, i] corresponds
    to eigenvalue [i]

    Where the eigenvalues equals 1 there is a support at that nodal degree of
    freedom and the eigenvector is zero when using the penalty approach.

    The eigenvalue array polulates right to left, opposite of how a list
    populates using append.

    AutoPIPE uses lumped mass, does not consider rotational inertia except
    for nodes with eccentric weights, and considers shear stiffness.
    """
    assert nmodes >= 3, "number of modes must be greater than 3"

    tqdm = logging.getLogger("tqdm")
    tqdm.info("*** Preprocessing, initializing analysis...")
    model.app.units.disable()
    tqdm.info("*** Switching to base units.")

    # do stuff here
    tqdm.info("*** Assembling system mass and stiffness matrices.")

    ndof = 6    # nodal degrees of freedom
    en = 2      # number of nodes per element
    nn = len(model.points)
    lc = len(model.loadcases)

    # similar to nodal dof matrix
    points = list(model.points)

    # global system stiffness matrix
    Ks = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)
    Ms = np.zeros((nn*ndof, nn*ndof), dtype=np.float64)

    zeros = partial(np.zeros, (12, 1))
    fixed_dof = defaultdict(zeros)

    for element in model.elements:
        idxi = points.index(element.from_point)
        idxj = points.index(element.to_point)

        # node and corresponding dof (start, finish), used to define the
        # elements of the system stiffness and force matrices
        niqi, niqj = idxi*ndof, idxi*ndof + ndof
        njqi, njqj = idxj*ndof, idxj*ndof + ndof

        meg = element.mglobal()
        keg = element.kglobal(294.2611)     # room temp

        # assemble global mass matrix, quadrant 1 to 4
        Ms[niqi:niqj, niqi:niqj] += meg[:6, :6]         # 2nd
        Ms[niqi:niqj, njqi:njqj] += meg[:6, 6:12]       # 1st
        Ms[njqi:njqj, niqi:niqj] += meg[6:12, :6]       # 3rd
        Ms[njqi:njqj, njqi:njqj] += meg[6:12, 6:12]     # 4th

        # assemble global stiffness matrix, quadrant 1 to 4
        Ks[niqi:niqj, niqi:niqj] += keg[:6, :6]         # 2nd
        Ks[niqi:niqj, njqi:njqj] += keg[:6, 6:12]       # 1st
        Ks[njqi:njqj, niqi:niqj] += keg[6:12, :6]       # 3rd
        Ks[njqi:njqj, njqi:njqj] += keg[6:12, 6:12]     # 4th

        # modify diagonal elements, penalty method, by adding large
        # stiffnesses to the diagonals where a support is located
        di = np.diag_indices(6)     # diagonal indices for 6x6 matrix
        for support in element.supports:
            ksup = support.kglobal(element)

            # bump up all support stiffness and nodal mass
            # ksup[np.where(ksup > 0)] = np.inf
            fixed_dof[support.point] += ksup

            Ks[niqi:niqj, niqi:niqj][di] += ksup[:6, 0]     # 2nd
            Ks[njqi:njqj, njqi:njqj][di] += ksup[6:12, 0]   # 4th

            Ms[niqi:niqj, niqi:niqj][di] += ksup[:6, 0]     # 2nd
            Ms[njqi:njqj, njqi:njqj][di] += ksup[6:12, 0]   # 4th

    tqdm.info("*** Reducing system mass and stiffness matrices.")
    # reduce assembled stiffness and mass matrices
    # Krows, Kcols = np.where(Ks == np.inf)
    # Mrows, Mcols = np.where(Ms == np.inf)

    # with redirect_stdout(sys.__stdout__):
    #     print(Krows, Kcols)

    # delete rows and columns
    # Kr = np.delete(np.delete(Ks, Krows, 0), Kcols, 1)
    # Mr = np.delete(np.delete(Ms, Mrows, 0), Mcols, 1)

    tqdm.info("*** Calculating dynamical matrix.")
    D = np.linalg.inv(Ks) @ Ms

    # total fixed dofs
    total_fixed_dofs = 0
    for kmat in fixed_dof.values():
        zrow, _ = np.where(kmat != 0)
        total_fixed_dofs += len(zrow)

    tqdm.info("*** Solving for eigenvalues and eigenvectors.")
    # inveigvals, eigvecmat = splinalg.eigsh(Mr, 12, M=Kr)
    inveigvals, eigvecmat = splinalg.eigs(D, nmodes+total_fixed_dofs)

    # divide by 2*pi to get the circular freq.
    eigvals = np.sqrt(1/inveigvals.real) / (2*np.pi)
    eigvecmat = np.real(eigvecmat)  # take real part

    # sort eigvals and corresponding vectors smallest -> largest
    idx = eigvals.argsort()[::1]
    eigvals = eigvals[idx]
    eigvecmat = eigvecmat[:, idx]

    # remove the trivial results in the front
    eigvals = eigvals[total_fixed_dofs:]
    eigvectmat = eigvecmat[:, total_fixed_dofs:]

    with redirect_stdout(sys.__stdout__):
        # print(total_fixed_dofs)
        print(eigvals)
        print(eigvecmat[:, 0])

    # switch back to user units - analysis is complete
    model.app.units.enable()

    # create output table

    tqdm.info("*** Analysis complete!\n\n\n")

    return eigvals, eigvecmat
