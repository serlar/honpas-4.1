Discussions
===========

2017-10-18
----------

Yann Pouillon
~~~~~~~~~~~~~

I'm now in the middle of adapting the implementation of hybrid XC functionals
into SIESTA. Part of the process is going well, but I'm constantly hitting the
same wall for the rest. It can be summarized in one commit::

    SIESTA trunk revno #368: Incorporation of BSC code

To be able to make the hybrid XC code work with the newest version of SIESTA, I
need to understand exactly what happened at that time, since there were broad
changes to both the structure of the code and the flow of data.

Digging into the available SIESTA documentation, I found a reference to these
modifications in Docs/CHANGES::

    An overview of the rationale and implementation of the new
    functionality can be found in

    An efficient implementation of a QM-MM method in SIESTA,
    C. F. Sanz-Navarro, R. Grima, A. Garcia, E. A. Bea, A. Soba,
    J.M. Cela, and P. Ordejon,
    Theoretical Chemistry Accounts. Online pub Sep 2010.
    DOI: 10.1007/s00214-010-0816-5

    We are working on a set of developer-oriented notes to facilitate
    future enhancements.

I downloaded and read the article. Although it helped me understand what has
been done and why, I still miss important information about how it was done, in
particular because many variables used by the hybrid XC implementation were
removed / replaced by something else.

The core of the problem is dhscf.F, unsurprisingly mentioned in the article,
which is also where the hybrid XC implementation plugs-in into SIESTA.

I thus would like to turn to the second document, the above-mentioned
"developer-oriented notes", but I have no idea where to look for them.

Could you give me access to these notes and/or the possibly available reports
about the implementation from BSC? That would allow me to save a lot of time.


Alberto Garcia
~~~~~~~~~~~~~~

The main BSC change to dhscf is to switch from one parallel distribution to
another depending on the task at hand. I spent quite a lot of time documenting
the flow of data and the distribution changes in the routine itself, but I am
afraid that we do not have any more documents on the issue.

You seem to be basing your work on a previous implementation of the hybrid
code. Which one is that? Is the code accessible somewhere? Maybe I can help you
if we look together at more concrete places in the code.


José María Soler
~~~~~~~~~~~~~~~~

The unfortunate fact was that Rogeli Grima and I developed simultaneously
the code for redistributing the density over processors for the xc part.
Both codes were based on the same ideas of Jose M Cela, but mine was more
complex because it was prepared also for vdW functionals. On the other
hand, Georg Huhs found that Grima's code was faster sometimes, so that
Alberto included both.

However, as far as I understand, they run alternatively and they have no
interactions at all. So I guess you will have to work with one of them
only, since I think it is not worth to implement the hybrid functionals for
both versions.

By the way, the siestaXC/gridXC library includes only my version.


2018-07-12
----------

Javier Junquera
~~~~~~~~~~~~~~~

We have made a lot of progress in the merge of your subroutines within the
trunk of SIESTA. Indeed, there is only one piece of information that is
missing.

In the subroutine RLYLM that computes the product of r^l times Y_lm, where r is
a radial distance and Y_lm an spherical harmonic, you have changed some of the
signs that were in all the SIESTA versions.

Can you remember what is the origin of this change?


Honghui Shang
~~~~~~~~~~~~~

I have found the reason of changing RLYLM in my notebooks.

In a test of Si2H6 molecule (sz basis set), the four-center electron repulsion
integral (6 6 6 4) was found to be positive when using NAO2GTO method, however,
when using my previous direct-NAO method (J. Phys. Chem. A 2010, 114,
1039–1043), such integral was negative, in order to get the same result as the
benchmark (direct-NAO method), the subroutine RLYLM was changed.

Coming back to the merging, could that be possible to keep these two version of
RLYLM in SIESTA? One (original version) for normal DFT calculation, and the
other one (my changed sign version) for Hatree-Fock calculations? Or do you
have any other suggestions?

