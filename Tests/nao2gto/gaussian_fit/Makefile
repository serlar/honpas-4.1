#
# Makefile for the Gaussian fit tester
#

FC = gfortran
FCFLAGS = -g3 -ggdb -Og -march=native -I.

                    # ------------------------------------ #

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) $(FCFLAGS) -o $@ -c $<

                    # ------------------------------------ #

all_targets all: test_gaussian_fit test_gaussian_fit.dat
	rm -f test_gaussian_fit.out
	-./test_gaussian_fit
	test -s test_gaussian_fit.out && gnuplot test_gaussian_fit.gnu

test_gaussian_fit: test_gaussian_fit.f90 nao2gto_fit.o nao2gto_fit_fcn.o
	$(FC) $(FCFLAGS) -o $@ $^ -L. -lminpack -lblas

nao2gto_fit.o: nao2gto_fit_fcn.o

clean:
	rm -f nao2gto_fit.mod nao2gto_fit.o
	rm -f nao2gto_fit_fcn.mod nao2gto_fit_fcn.o
	rm -f test_gaussian_fit test_gaussian_fit.o test_gaussian_fit.out
