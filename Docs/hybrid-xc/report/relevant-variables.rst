Relevant variables
==================

Scalars
-------

The calculation of Hartree-Fock-related quantities heavily relies on the
following scalars:

- maxnh: maximum number of non-zero H matrix elements
- na_s: number of atoms in the supercell
- na_u: number of atoms in the unit cell
- nlhmax: maximum number of basis orbitals interacting with any orbital
- no_l: number of orbitals local to a node (usually 1)
- no_s: number of orbitals in the supercell
- no_u: number of orbitals in the unit cell
- spin%Grid: number of spin components

.. note:: The nspin variable was replaced by a derived type tSpin between
   SIESTA 3 and 4, hence the spin%Grid field.

Here are the specifications of all these scalars, as well as their possible
aliases within HONPAS:

============  ====  ====================  ============
Name          Type  Module                Aliases
============  ====  ====================  ============
maxnh         int   sparse_matrices       maxnd
na_s          int   siesta_geom           na
na_u          int   siesta_geom           nua
nlhmax        int   hsparse               N/A
no_l          int   siesta_geom           N/A
no_s          int   atomlist              norb, nuotot
no_u          int   atomlist              N/A
spin%Grid     int   m_spin                nspin
============  ====  ====================  ============

.. warning:: The value of na_s and no_s may vary during some calculations.
   This may happen in the atomlist and m_iostruct modules, as well as
   in the coor subroutine.


Arrays
------

The calculation of Hartree-Fock-related quantities heavily relies on the
following arrays:

- iaorb: atomic index of each orbital
- indxua: index of equivalent atoms in the unit cell
- indxuo: index of equivalent orbitals in the unit cell
- iphorb: index of each orbital within the same atom
- isa: species index of each atom
- listh:
- listhptr:
- nsc: diagonal elements of the supercell
- scell: supercell vectors, column-wise
- ucell: unit cell vectors, column-wise
- xa: atomic positions in the supercell

Here are the specifications of all these arrays, as well as their possible
aliases within HONPAS:

============  ====  ====================  ============  ====================
Name          Type  Module                Aliases       Dimensions
============  ====  ====================  ============  ====================
iaorb         int   atomlist                            no_u
indxua        int   atomlist                            na_s
indxuo        int   atomlist                            no_u
iphorb        int   atomlist                            no_u
isa           int   siesta_geom           N/A           na_s
listh         int   sparse_matrices       listd         max(1,nlhmax)
listhptr      int   sparse_matrices       listdptr      max(1,no_l)
nsc           int   siesta_geom           N/A           3
numh          int   sparse_matrices       numd          no_l
scell         dble  siesta_geom           cell          3, 3
ucell         dble  siesta_geom           N/A           3, 3
xa            dble  siesta_geom           N/A           3, na_s
============  ====  ====================  ============  ====================

