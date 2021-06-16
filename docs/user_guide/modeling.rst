Modeling
========
Under the hood, PSI uses beam elements to model the various piping components
encountered in a typically piping system. The piping components are assigned a
section and material based on the industry standard and design requirements for
the project.

Components
----------
The most generic piping component is a pipe run. A run element consists of two
nodes with 6 degrees of freedom (DOF) at each node, resulting in a total of 12
DOFs each.

The stiffness matrix derivation for a beam element is commonly documented in
various books and in literature. PSI uses the Timoshenko beam derivation from
*"Theory of Matrix Structural Analysis"* by J.S. Przemieniecki, which accounts
for shear deflection, an effect more pronounced for thick wall piping. Note
that shear effects may be turned off by the user from model settings.

The element stiffness matrix is given as:

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
    k_{1,7} & = k_{7,1} = \dfrac{-EA}{L}\\
    k_{2,2} & = \dfrac{12EI_z}{L^{3}(1+\Phi_y)}\\
    k_{2,6} & = k_{6,2} = \dfrac{6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{2,8} & = k_{8,2} = \dfrac{-12EI_z}{L^{3}(1+\Phi_y)}\\
    k_{2,12} & = k_{12,2} = \dfrac{6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{3,3} & = \dfrac{12EI_y}{L^{3}(1+\Phi_z)}\\
    k_{3,5} & = k_{5,3} = \dfrac{-6EI_y}{L^{2}(1+\Phi_z)}\\
    k_{3,9} & = k_{9,3} = \dfrac{-12EI_y}{L^{3}(1+\Phi_z)}\\
    k_{3,11} & = k_{11,3} = \dfrac{-6EI_y}{L^{2}(1+\Phi_z)}\\
    k_{4,4} & = \dfrac{GJ}{L}\\
    k_{10,4} & = k_{4,10} = \dfrac{-GJ}{L}\\
    k_{5,5} & = \dfrac{(4+\Phi_z)EI_y}{L(1+\Phi_z)})\\
    k_{5,9} & = k_{9,5} = \dfrac{6EI_y}{(L^{2}(1+\Phi_z)}\\
    k_{5,11} & = k_{11,5} = \dfrac{(2-\Phi_z)EI_y}{L(1+\Phi_z)}\\
    k_{6,6} & = \dfrac{(4+\Phi_y)EI_z}{L(1+\Phi_y)}\\
    k_{6,8} & = k_{8,6} = \dfrac{-6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{6,12} & = k_{12,6} = \dfrac{(2-\Phi_y)EI_z}{L(1+\Phi_y)}\\
    k_{7,7} & = \dfrac{EA}{L}\\
    k_{8,8} & = \dfrac{12EI_z}{L^{3}(1+\Phi_y)}\\
    k_{8,12} & = k_{12,8} = \dfrac{-6EI_z}{L^{2}(1+\Phi_y)}\\
    k_{9,9} & = \dfrac{12EI_y}{L^{3}(1+\Phi_z)}\\
    k_{9,11} & = k_{11,9} = \dfrac{6EI_y}{L^{2}(1+\Phi_z)}\\
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


Sections
--------
Under Construction!


Materials
---------
Under Construction!
