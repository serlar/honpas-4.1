SIESTA is both a method and its computer program implementation, to
perform efficient electronic structure calculations and ab initio
molecular dynamics simulations of molecules and solids. SIESTA's
efficiency stems from the use of strictly localized basis sets and from
the implementation of linear-scaling algorithms which can be applied to
suitable systems. A very important feature of the code is that its
accuracy and cost can be tuned in a wide range, from quick exploratory
calculations to highly accurate simulations matching the quality of
other approaches, such as plane-wave and all-electron methods.

The possibility of treating large systems with some first-principles
electronic-structure methods has opened up new opportunities in many
disciplines. The SIESTA program is distributed freely to academics and
has become quite popular, being increasingly used by researchers in
geosciences, biology, and engineering (apart from those in its natural
habitat of materials physics and chemistry). Currently there are several
thousand users all over the world, and the paper describing the method
(J. Phys. Cond. Matt. 14, 2745 (2002)) has received more than 5000
citations so far.

SIESTA's main characteristics are:

-   It uses the standard Kohn-Sham self-consistent density functional
    method in the local density (LDA-LSD) or generalized gradient
    (GGA) approximations. Recent versions implement a functional capable
    of describing van der Waals interactions.
-   It employs norm-conserving pseudopotentials in their fully
    nonlocal (Kleinman-Bylander) form.
-   It uses atomic orbitals with finite support as a basis set, allowing
    unlimited multiple-zeta and angular momenta, polarization and
    off-site orbitals. Finite-support basis sets are the key for
    calculating the Hamiltonian and overlap matrices in O(N) operations.
-   Projects the electron wavefunctions and density onto a real-space
    grid in order to calculate the Hartree and exchange-correlation
    potentials and their matrix elements.

SIESTA can be compiled for serial or parallel execution (under MPI), and
can provide (the list is continuously expanding):

-   Total and partial energies.
-   Atomic forces.
-   Stress tensor.
-   Electric dipole moment.
-   Atomic, orbital and bond populations (Mulliken).
-   Electron density.
-   Geometry relaxation, fixed or variable cell.
-   Constant-temperature molecular dynamics (Nose thermostat).
-   Variable cell dynamics (Parrinello-Rahman).
-   Spin polarized calculations (collinear or not).
-   k-sampling of the Brillouin zone.
-   Local and orbital-projected density of states.
-   COOP and COHP curves for chemical bonding analysis.
-   Dielectric polarization.
-   Vibrations (phonons).
-   Band structure.

Starting from version 3.0, SIESTA includes the TranSIESTA module, which
provides the ability to model open-boundary systems where ballistic
electron transport is taking place. Using TranSIESTA one can compute
electronic transport properties, such as the zero-bias conductance and
the I-V characteristic, of a nanoscale system in contact with two
electrodes at different electrochemical potentials.
