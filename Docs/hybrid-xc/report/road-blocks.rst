Road blocks
===========

dhscf.F
-------

The *dhscf.F* file has been modified significantly between versions 3.2 and 4.0
of SIESTA, making the compilation of its HONPAS version fail with SIESTA
4.1-b3. Most changes were brought in at revision 368, aka "Incorporation of
BSC code", which had a strong impact on the SCF-related source code and
removed a few variables.

The compilation errors with the pristine version of *dhscf.F* are the
following:

.. code::

   ~siesta-honpas/Src/dhscf.F:1749:28:

      .                 listh, samexa, Hmat, Exc)
                            1
   Error: Symbol ‘listh’ at (1) has no IMPLICIT type
   ~siesta-honpas/Src/dhscf.F:1748:67:

      .                 listdptr, listd, Dscf, maxnh, numh, listhptr, 
                                                                   1
   Error: Symbol ‘listhptr’ at (1) has no IMPLICIT type
   ~siesta-honpas/Src/dhscf.F:1741:10:

       do n = 1,nXCfunc
          1
   Error: Symbol ‘n’ at (1) has no IMPLICIT type
   ~siesta-honpas/Src/dhscf.F:1747:53:

      .                 na, isa, xa, indxua, cell, nsc, maxnd, numd, 
                                                     1
   Error: Symbol ‘nsc’ at (1) has no IMPLICIT type
   ~siesta-honpas/Src/dhscf.F:1748:57:

      .                 listdptr, listd, Dscf, maxnh, numh, listhptr, 
                                                         1
   Error: Symbol ‘numh’ at (1) has no IMPLICIT type
   ~siesta-honpas/Src/dhscf.F:1741:22:

       do n = 1,nXCfunc
                      1
   Error: Symbol ‘nxcfunc’ at (1) has no IMPLICIT type
   ~siesta-honpas/Src/dhscf.F:1749:36:

      .                 listh, samexa, Hmat, Exc)
                                    1
   Error: Symbol ‘samexa’ at (1) has no IMPLICIT type
   arch.make:203: recipe for target 'dhscf.o' failed
   make: *** [dhscf.o] Error 1

A manual diff3 between the 3.2, HONPAS and 4.1-b3 versions of dhscf.F showed
that the Git rebase upon SIESTA 4.1-b3 had been performed correctly, since the
resulting manually merged file was only differing regarding whitespace on one
single line.

Analysing the changes between revision 367 and 368 of the trunk of SIESTA
revealed that the removed variables, *listh*, *listhptr*, *n*, *nsc*, *numh*,
*nxcfunc*, and *samexa*, had already been unused for some time and were
logically pruned out of the source code.

The Hartree-Fock block in the routine was the following:

.. code::

   C-----------------------------------------------------------------------
   C add Hartree-Fock exchange potential to Hamiltonian (nuo, norb)
   C xmqin, January 2014
   C----------------------------------------------------------------------
         call timer( 'DHFX', 1 )
   !      if(node.eq.0) write(6,*) 'Hybrid DFT calculations'
         do n = 1,nXCfunc
            if (leqi(XCauth(n),'HSE06') .or. leqi(XCauth(n),'hse06')
        .       .or. leqi(XCauth(n),'PBE0') .or . leqi(XCauth(n),'pbe0'))
        .   then
   !! Both cell and supercell have been diagonalized.
         call update_hfx( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
        .                 na, isa, xa, indxua, cell, nsc, maxnd, numd, 
        .                 listdptr, listd, Dscf, maxnh, numh, listhptr, 
        .                 listh, samexa, Hmat, Exc)
            endif
         enddo
   
         call timer( 'DHFX', 2 )
   C------------------------------------------------------------------------


Dealing with numh, listh, and listhptr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coming soon ...


Dealing with samexa
~~~~~~~~~~~~~~~~~~~


