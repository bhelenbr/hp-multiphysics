DIRS = symbolic_function utilities input_map quadtree spline++
TRI_DIRS = tri_basis tri_mesh tri_hp
TET_DIRS = tet_basis tet_mesh tet_hp

#DEFINES = -DBZ_DEBUG
#DEFINES +=-Df2cFortran
#DEFINES += -DIBMR2Fortran
export DEFINES

OPT = -O3 -fpic -Wno-deprecated-declarations -Wuninitialized
#OPT += -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment
#OPT += -fsanitize=memory
#OPT += -fsanitize=thread
#OPT += -funroll-loops -falign-loops=16
#OPT += -Wno-long-double  #(OS X)
#OPT += -g  #(Gnu Debug)
#OPT = -debug #(for intel debug)
#OPT = -Ofast=IP30  #(sauter)
#OPT = -Ofast=IP27  #(ctr-sgi1)
#OPT = -Ofast=IPrk5 #(origin)
export OPT

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

PACKAGES=${HOME}/Packages
LIBDIRS=-L${PACKAGES}/lib -L${ROOT_DIR}/lib
#LIBDIRS += -framework veclib (OS X)
export LIBDIRS

INCLUDEDIRS=-I${PACKAGES}/include -I${ROOT_DIR}/include
export INCLUDEDIRS

#LIBS += -lasan -lubsan
export LIBS

if iseq(${CXXMPI},"mpiicpc")
then
	LIBBLAS = -lmkl -lmkl_lapack -lmkl_intel_thread -qopenmp #(for intel acres)
else
	LIBBLAS  = -lopenblas #(acres gnu)
endif
#LIBBLAS = -framework veclib #(OS X)
#LIBBLAS = -llapack -lblas #(cares)
#LIBBLAS = -latlas -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
#LIBBLAS = -lmkl -lmkl_lapack #(for intel)
#LIBBLAS = -lessl /afs/cu/software/Development/lapack-3.0/rs_aix52/LAPACK/lapack_RS6K.a -lxlf -lxlf90 #(jupiter)
export LIBBLAS


#AR = libtool
#ARFLAGS = -o #(darwin)
#ARFLAGS = --mode=link ${CC} -o #(cares)
#export AR
#export ARFLAGS

all: dirs tri_hp tet_hp 

tet_hp: tet_mesh tet_basis tri_mesh $(DIRS) force_look
	cd $@; $(MAKE) $(MFLAGS)

tet_mesh: tri_mesh $(DIRS) force_look
	cd $@; $(MAKE) $(MFLAGS)

tet_basis: utilities force_look
	cd $@; $(MAKE) $(MFLAGS)

tri_hp: tri_mesh tri_basis force_look
	cd $@; $(MAKE) $(MFLAGS)

tri_mesh: $(DIRS) utilities input_map quadtree spline++ symbolic_function force_look
	cd $@; $(MAKE) $(MFLAGS)

tri_basis: utilities force_look
	cd $@; $(MAKE) $(MFLAGS)

quadtree: utilities input_map force_look
	cd $@; $(MAKE) $(MFLAGS)

spline++: utilities input_map force_look
	cd $@; $(MAKE) $(MFLAGS)

symbolic_function: utilities input_map force_look
	cd $@; $(MAKE) $(MFLAGS)

input_map: utilities force_look
	cd $@; $(MAKE) $(MFLAGS)

utilities: force_look
	cd $@; $(MAKE) $(MFLAGS)
	
dirs: 
	+@[ -d  bin ] || mkdir -p bin
	+@[ -d  lib ] || mkdir -p lib
	+@[ -d  include ] || mkdir -p include

clean:
	for d in $(TRI_DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done
	for d in $(TET_DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done
	for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

force_look:
	true
