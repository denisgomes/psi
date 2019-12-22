"""Static, modal and dynamic solvers."""

import sys
import logging
from contextlib import redirect_stdout

import numpy as np

from psi.loadcase import LoadCase
from psi import units


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

    Non linear supports use an iterative approach to determine the final
    support loads and pipe configuration. The sequence in which loads are
    applied matters as superposition is not valid for non-linear analysis.
    Each load is applied and the displacements extracted. These movements are
    then used for the next step, ie. the model mesh is modified to incorporate
    the new point locations. As a result the system stiffness matrix is updated
    each iteration.
    """
    tqdm = logging.getLogger("tqdm")
    tqdm.info("*** Preprocessing, initializing analysis...")

    # Note: All internal object data is stored in SI once loaded from
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

            if isinstance(loadcase, LoadCase):
                # sum of all loads in a loadcase
                feg = np.zeros((en*ndof, 1), dtype=np.float64)

                for loadtype, opercase in loadcase.loads:

                    for load in element.loads:

                        if (isinstance(load, loadtype) and
                                load.opercase == opercase):
                            feg += load.fglobal(element)

                # assemble global system force matrix
                Fs[niqi:niqj, i] += feg[:6, 0]
                Fs[njqi:njqj, i] += feg[6:12, 0]

            # large stiffness added to each force matrix component with
            # non-zero support displacements
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
    X = np.linalg.solve(Ks, Fs)

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

            # internal force and moment matrix
            Fi[niqi:niqj, i] = fi[:6, 0]
            Fi[njqi:njqj, i] = fi[6:12, 0]

    tqdm.info("*** Writing loadcase results data.")
    for i, loadcase in enumerate(model.loadcases):
        # load combs are combined later
        if isinstance(loadcase, LoadCase):
            # X[:, i], R[:, i] and Fi[:, i] are row vectors
            loadcase.movements.results = X[:, i]
            loadcase.reactions.results = R[:, i]
            loadcase.forces.results = Fi[:, i]

    # with redirect_stdout(sys.__stdout__):
    #     print(X)
    #     print(R)
    #     print(Fi)

    tqdm.info("*** Analysis complete!\n")
    model.app.units.enable()


def modal(model):
    """Run a modal analysis of the system extracting frequencies and mode
    shapes.

    All nonlinearity is neglected for this type of analysis. In other words,
    all support types are considered double acting and friction is not
    considered.
    """
    pass
