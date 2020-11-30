#
# Single-test makefile template
#
# You can edit the SIESTA macro here, or pass it on the command line

MPI=mpirun -np 2
SIESTA=../../../siesta

# Example for BSC runs
#
#MPI= mpirun -np 2 
#SIESTA= ../../../siesta

# Make compatibility layer for old test-runs
ifeq ($(strip $(firstword $(SIESTA))),mpirun)
MPI=
endif
ifeq ($(strip $(firstword $(SIESTA))),mpiexec)
MPI=
endif

#----------------------------------------------------------------------------
REFERENCE_DIR?=../../../Tests/Reference
REFERENCE_CHECKER?=../cmp_digest.sh

label=work

.PHONY: completed
completed: completed_$(label)

completed_$(label):
	@echo ">>>> Running $(name)1 test..."
	@if [ -d $(label) ] ; then rm -rf $(label) ; fi; mkdir $(label)
	@if [ -n "$(EXTRAFILES)" ] ; then cp -f $(EXTRAFILES) $(label) ; fi
	@for i in `cat $(name).pseudos` ; do \
          echo "    ==> Copying pseudopotential file for $$i..." ;\
          ln ../Pseudos/$$i.psf $(label)/$$i.psf ;\
         done
	@echo "    ==> Running SIESTA as $(MPI) ${SIESTA}"
	@(cd $(label) ; $(MPI) $(SIESTA) 2>&1 > $(name)1.out < ../$(name)1.fdf) \
          && touch completed_$(label)
	@if [ -f completed_$(label) ] ; then cp $(label)/$(name)1.out . ;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name)1 did not complete successfully";\
         fi
	@if [ -f Ge.TDXV ] ; then cp Ge.TDXV $(label)/ ; fi
	@if [ -f completed_$(label) ] ; then rm completed_$(label) ; fi
	@echo ">>>> Running $(name)2 test..."
	@echo "    ==> Running SIESTA as $(MPI) ${SIESTA}"
	@(cd $(label) ; $(MPI) $(SIESTA) 2>&1 > $(name)2.out < ../$(name)2.fdf) \
          && touch completed_$(label)
	@if [ -f completed_$(label) ] ; then cp $(label)/$(name)2.out . ;\
           echo "    ===> SIESTA finished successfully";\
         else \
           echo " **** Test $(name)2 did not complete successfully";\
         fi
#
check: completed check-only

check-only:
	@echo "    ==> Running check for system $(name)1"
	@REFERENCE_DIR=$(REFERENCE_DIR) sh $(REFERENCE_CHECKER) $(name)1.out
	@echo "    ==> Running check for system $(name)2"
	@REFERENCE_DIR=$(REFERENCE_DIR) sh $(REFERENCE_CHECKER) $(name)2.out
#
clean:
	@echo ">>>> Cleaning $(name) test..."
	rm -rf $(label) completed_$(label) *.out *.xml
