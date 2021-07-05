Modeling
========
Under the hood, PSI uses finite element beam elements to model various piping
components typically encountered in piping system.


Components
----------
The geometry of piping components are defined by the coordinates of the 'from'
and 'to' points. Each component must also have a section and material assigned
to be fully defined. Extra data such as insulation/cladding/refractory, piping
code, and sifs/connections may also be applied based on design requirements.


Run
~~~
The most generic piping component is a pipe Run; all other elements are derived
from a Run or approximated using a Run. A Run element is simply a 3D beam
element consisting of two nodes, each with 6 Degrees of Freedom (DOF),
resulting in a total of 12 DOFs per element. The user must define the global
(x,y,z) coordinate positions of the 'from' and 'to' nodal points for a Run.
Doing so effectively set the length (i.e. the extents) of the element.

The stiffness matrix derivation for a beam element is commonly documented in
various books and in literature. PSI uses the Timoshenko beam derivation from
*"Theory of Matrix Structural Analysis"* by J.S. Przemieniecki, which accounts
for shear deflection, an effect more pronounced for thick wall piping. Note
that shear effects are turned on my default but may be turned off by the user
from model settings.

The element (local) stiffness matrix is a 12x12 matrix:

.. math::

    \begin{align}
    k_{m,n} =
    \begin{pmatrix}
    k_{1,1} & k_{1,2} & \cdots & k_{1,n} \\
    k_{2,1} & k_{2,2} & \cdots & k_{2,n} \\
    \vdots  & \vdots  & \ddots & \vdots  \\
    k_{m,1} & k_{m,2} & \cdots & k_{m,n}
    \end{pmatrix}
    \end{align}

with m=n=12 and where,

.. math::

    \begin{align*}
    k_{1,1} & = \dfrac{EA}{L}\\
    k_{1,7} = k_{7,1} & = \dfrac{-EA}{L}\\
    k_{2,2} & = \dfrac{12EI_z}{L^{3}(1+\Phi_y)}\\
    k_{2,6} = k_{6,2} & = \dfrac{6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{2,8} = k_{8,2} & = \dfrac{-12EI_z}{L^{3}(1+\Phi_y)}\\
    \end{align*}

.. math::

    \begin{align*}
    k_{2,12} = k_{12,2} & = \dfrac{6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{3,3} & = \dfrac{12EI_y}{L^{3}(1+\Phi_z)}\\
    k_{3,5} = k_{5,3} & = \dfrac{-6EI_y}{L^{2}(1+\Phi_z)}\\
    k_{3,9} = k_{9,3} & = \dfrac{-12EI_y}{L^{3}(1+\Phi_z)}\\
    k_{3,11} = k_{11,3} & = \dfrac{-6EI_y}{L^{2}(1+\Phi_z)}\\
    \end{align*}

.. math::

    \begin{align*}
    k_{4,4} & = \dfrac{GJ}{L}\\
    k_{10,4} = k_{4,10} & = \dfrac{-GJ}{L}\\
    k_{5,5} & = \dfrac{(4+\Phi_z)EI_y}{L(1+\Phi_z)}\\
    k_{5,9} = k_{9,5} & = \dfrac{6EI_y}{L^{2}(1+\Phi_z)}\\
    k_{5,11} = k_{11,5} & = \dfrac{(2-\Phi_z)EI_y}{L(1+\Phi_z)}\\
    \end{align*}

.. math::

    \begin{align*}
    k_{6,6} & = \dfrac{(4+\Phi_y)EI_z}{L(1+\Phi_y)}\\
    k_{6,8} = k_{8,6} & = \dfrac{-6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{6,12} = k_{12,6} & = \dfrac{(2-\Phi_y)EI_z}{L(1+\Phi_y)}\\
    k_{7,7} & = \dfrac{EA}{L}\\
    k_{8,8} & = \dfrac{12EI_z}{L^{3}(1+\Phi_y)}\\
    \end{align*}

.. math::

    \begin{align*}
    k_{8,12} = k_{12,8} & = \dfrac{-6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{9,9} & = \dfrac{12EI_y}{L^{3}(1+\Phi_z)}\\
    k_{9,11} = k_{11,9} & = \dfrac{6EI_y}{L^{2}(1+\Phi_z)}\\
    k_{10,10} & = \dfrac{GJ}{L}\\
    k_{11,11} & = \dfrac{(4+\Phi_z)EI_y}{L(1+\Phi_z)}\\
    k_{12,12} & = \dfrac{(4+\Phi_y)EI_z}{L(1+\Phi_y)}\\
    \end{align*}

and where,

.. math::

    \begin{align*}
    \Phi_y = \dfrac{12EI_z}{GA_yL^{2}}\\
    \Phi_z = \dfrac{12EI_y}{GA_zL^{2}}
    \end{align*}

The element global stiffness matrix is calculated by pre and post multiplying
the element local stiffness matrix by the transpose of the transformation
matrix and transformation matrix respectively, as shown below:

.. math::

    \begin{align}
    k_{global} = T^{T} * k_{local} * T
    \end{align}

The local to global transformation matrix is used to convert quantaties that
are defined with respect to element coordinates to global coordinates. It
consists of the direction cosines of the element local axes given in matrix
format as shown below where m=n=12:

.. math::

    \begin{align}
    T_{m,n} =
    \begin{pmatrix}
    dc_{1,1} & dc_{1,2} & dc_{1,3} & 0 & \cdots & 0 \\
    dc_{2,1} & dc_{2,2} & dc_{2,3} & 0 & \cdots & 0 \\
    dc_{3,1} & dc_{3,2} & dc_{3,3} & 0 & \cdots & 0 \\
    \vdots  & \vdots  & \ddots & \vdots  \\
    T_{m,1} & T_{m,2} & \cdots & T_{m,n}
    \end{pmatrix}
    \end{align}

The 3x3 direction cosine matrix (dc) is given below:

.. math::

    \begin{align}
    dc_{3,3} =
    \begin{pmatrix}
    localx_{1} & localy_{1} & localz_{1}
    localx_{2} & localy_{2} & localz_{2}
    localx_{3} & localy_{3} & localz_{3}
    \end{pmatrix}
    \end{align}

For a Run element, the local x direction is given by the vector from the 'from'
point to the 'to' point. The local y is parallel to the global vertical
direction and the local z is the cross product of the local x and y axes.

When the local x is parallel to the vertical direction, the local y is aligned
with global x (arbitrarily) and the local z is the cross product of local x
with local y.


Bend
~~~~
Pipe bends are approximated using multiple Run elements strung together based
on an underlying rational bezier curve used to locate each point coordinate.
Similar to a Run, a Bend has a 'from' and 'to' point. It also has a 'near',
'far' and 'mid' point. The 'near' and 'far' points are the start and end points
for the Bend. The 'mid' point of the bend is the physical center of the Bend
located on the arc length. Physical quantaties are calculated at these three
point locations and directly exposed to the user. The results for all other
points are solved but not directly accessible.

.. note::

    The 'to' point is the corner of the bend element and does not fall on the
    Bend arc length.

Based on the geometry of the Bend and the governing piping code, the Bend
flexibility factor is calculated, and the stiffness matrices of all the
approximating Run elements are divided by the same factor.

.. note::

   Only the bending stiffnesses are affected by the increase in flexibility and
   therefore altered.


Reducer
~~~~~~~
Similar to a pipe Bend, a Reducer is approximated using multiple Run elements
with reducing cross-sections.


Rigid
~~~~~
A Rigid element is used to model relatively stiff components in a piping system
such as equipment for example to implicitly determine the nozzle movements and
nozzle loads at the interface points. These elements can be weightless or have
a mass assigned to them. They can also have a temperature assigned to them or
have the fluid contents turned off.

.. note::

   The relative stiffness is achieved by multiplying the thickness of the
   piping section by a factor of 10 times.


Valve
~~~~~
A Valve element is basically a Rigid element with a mass assigned to it. The
mass can be user defined or pulled for a data file. From a stress standpoint
a Valve locally stiffens a piping system. It may also have Flanges defined
that effectively increases the overall weight of the Valve.


Flange
~~~~~~
Similar to a Valve, a Flange is also a Rigid element. The mass of the Flange
can be user defined or pulled for a data file. Leak testing of Flanges is
possible dependent on the code.


Bellow
~~~~~~
Under Construction!


Component Data
--------------

Sections
~~~~~~~~
Under Construction!


Materials
~~~~~~~~~
Under Construction!


Insulation / Refractory / Cladding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Under Construction!


Codes
~~~~~
Under Construction!


SIFs/Connections
~~~~~~~~~~~~~~~~
Under Construction!
