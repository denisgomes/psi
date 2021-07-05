Solvers
=======

Static
------
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
~~~~~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~
Both reducers and bends are topologically curves and therefore approximations.
The midpoint of a bend is created as a side effect of defining the bend when
it is created. All the other vertices are inaccessible from a nodal standpoint.
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
~~~~~~~~~~~~~~~
Under Construction!

Skewed Supports
~~~~~~~~~~~~~~~
Implemented using local stiffness transformation.

Spring Algorithm
~~~~~~~~~~~~~~~~
The program tries to place a rigid support, variable or constant spring where
a spring support is specified by the user. The algorithm must satisfy all the
operating cases such that it does not bottom or top out. For each operating
case a linear deadweight analysis (W+P+D+F) is performed and a +Y support is
placed at each supports. Remove the force from the weight case at each spring
support location with the equivalent upward force and run the operating case
with the temperature. Determine the upward movement of the pipe at the spring
support. Use the applied force (i.e. the hot load) and the movement to pick a
hanger from the manufacturer catalog.

Master Slave (CNode)
~~~~~~~~~~~~~~~~~~~~
Constraint equations.

Non-Linear Supports
~~~~~~~~~~~~~~~~~~~
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

.. note::

    A *static* nonlinear analysis is performed by iterating until the solution
    converges.

Friction
~~~~~~~~
Under Construction!

.. note::

    A *static* nonlinear analysis is performed by iterating until the solution
    converges.

Loading Sequence and Non-Linear Supports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Non-linear supports use an iterative approach to determine the final support
loads and pipe configuration. The sequence in which loads are applied matters
as superposition is not valid for non-linear analysis. Each load is applied
and the displacements extracted. These movements are then used for the next
step, ie. the model mesh is modified to incorporate the new point locations.
As a result the system stiffness matrix is updated each iteration.


Modal
-----
Under Construction!


Harmonic
--------
Under Construction!


Spectra
-------
Under Construction!


Time History
------------
Under Construction!
