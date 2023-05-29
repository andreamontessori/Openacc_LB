CC=mpicc 
CFLAGS=-O2 -cpp
FC=mpif90 -DSERIAL
FC=nvfortran


CUDAFLAGS = -cuda -gpu=cc70 -O0 -g -Mbounds -Mchkptr -Mchkstk
CUDAFLAGS = -cuda -gpu=cc70,keepptx -O0 -g
CUDAFLAGS = -cuda -fast -gpu=cc80 -DMYDIMESION=128 -DTILE1=128 -DTILE2=1 -DTILE3=1
CUDAFLAGS = -gpu=cc75,cuda11.0,lineinfo,ptxinfo -cpp -cuda -O3



main.x: cudavars_mod.o pbc_kernels_mod.o streamcoll_bulk_kernels_mod.o \
	streamcoll_bc_kernels_mod.o correct_press_kernels_mod.o prints_mod.o cuLLB_regularized_3d.o 
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -o $@ cuLLB_regularized_3d.o cudavars_mod.o \
	pbc_kernels_mod.o streamcoll_bulk_kernels_mod.o streamcoll_bc_kernels_mod.o \
	correct_press_kernels_mod.o prints_mod.o

prints_mod.o: Makefile defines.h prints_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c prints_mod.f90

cudavars_mod.o: Makefile cudavars_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c cudavars_mod.f90

pbc_kernels_mod.o: Makefile pbc_kernels_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c pbc_kernels_mod.f90

streamcoll_bulk_kernels_mod.o: Makefile streamcoll_bulk_kernels_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c streamcoll_bulk_kernels_mod.f90

streamcoll_bc_kernels_mod.o: Makefile streamcoll_bc_kernels_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c streamcoll_bc_kernels_mod.f90

correct_press_kernels_mod.o: Makefile correct_press_kernels_mod.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c correct_press_kernels_mod.f90

cuLLB_regularized_3d.o: Makefile defines.h cuLLB_regularized_3d.f90
	$(FC) $(CUDAFLAGS) $(F90FLAGS) -c cuLLB_regularized_3d.f90

clean:
	@rm -rf main.x *.o *.mod *.ptx
