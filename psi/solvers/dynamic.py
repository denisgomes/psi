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

"""Dynamic solvers"""

from collections import defaultdict
import sys
import logging
from contextlib import redirect_stdout
from functools import partial

import numpy as np
from scipy.sparse import linalg as splinalg


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
        keg = element.kglobal(model.settings.tref)     # room temp

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
    eigvecmat = eigvecmat[:, total_fixed_dofs:]

    with redirect_stdout(sys.__stdout__):
        # print(total_fixed_dofs)
        print(eigvals)
        print(eigvecmat[:, 0])

    # switch back to user units - analysis is complete
    model.app.units.enable()

    # create output table

    tqdm.info("*** Analysis complete!\n\n\n")

    return eigvals, eigvecmat
