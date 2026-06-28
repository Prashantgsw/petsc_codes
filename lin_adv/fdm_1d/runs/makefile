!FC = gfortran

verify_convergence: verify_convergence.f90
	$(FC) -O3 verify_convergence.f90 -o verify_convergence -llapack -lblas

clean:
	rm -f verify_convergence *.o *.mod