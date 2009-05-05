ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: planar axi

planar: tri_hp tri_hp_mpi
axi: tri_hp_axi tri_hp_axi_mpi

tri_hp:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp
	cp _$@/tri_hp $(HOME)/bin/$@

.PHONY: tri_hp

tri_hp_axi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp EXTRAFLAGS=-DAXISYMMETRIC
	cp _$@/tri_hp $(HOME)/bin/$@

.PHONY: tri_hp_axi

tri_hp_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_mpi EXTRAFLAGS=-DMPISRC
	cp _$@/tri_hp_mpi $(HOME)/bin/$@

.PHONY: tri_hp_mpi

tri_hp_axi_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_mpi EXTRAFLAGS='-DMPISRC -DAXISYMMETRIC'
	cp _$@/tri_hp_mpi $(HOME)/bin/$@

.PHONY: tri_hp_axi_mpi

tri_hp_drop:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp EXTRAFLAGS='-DDROP -DAXISYMMETRIC'
	cp _$@/tri_hp $(HOME)/bin/$@

.PHONY: tri_hp_drop

tri_hp_drop_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_mpi EXTRAFLAGS='-DDROP -DAXISYMMETRIC -DMPISRC'
	cp _$@/tri_hp_mpi $(HOME)/bin/$@

.PHONY: tri_hp_drop_mpi

tri_hp_drop_stokes:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp EXTRAFLAGS='-DDROP -DAXISYMMETRIC -DINERTIALESS'
	cp _$@/tri_hp $(HOME)/bin/$@

.PHONY: tri_hp_drop_stokes

tri_hp_drop_stokes_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_mpi EXTRAFLAGS='-DDROP -DAXISYMMETRIC -DINERTIALESS -DMPISRC'
	cp _$@/tri_hp_mpi $(HOME)/bin/$@

.PHONY: tri_hp_drop_stokes_mpi

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _tri_hp*

else

VPATH = ../

BASE = adapt.o allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o minvrt.o movco.o movfin.o \
nstage.o output.o ptprobe.o tadvance.o tobasis.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o ins/surface.o
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/length.o cd/rsdl.o cd/tstep.o
OBJPS = planestrain/allocate.o planestrain/bdry.o planestrain/getnewbdry.o planestrain/length.o planestrain/rsdl.o planestrain/tstep.o
OBJSWIRL = swirl/allocate.o swirl/bdry.o swirl/getnewbdry.o swirl/getnewibc.o swirl/length.o swirl/rsdl.o swirl/tstep.o
OBJBUOY = buoyancy/allocate.o buoyancy/tstep.o buoyancy/rsdl.o
OBJSW = swe/allocate.o swe/getnewbdry.o swe/bdry.o swe/getnewibc.o swe/rsdl.o swe/tstep.o
OBJLVLSET = lvlset/allocate.o  lvlset/bdry.o lvlset/rsdl.o  lvlset/tstep.o
OBJEXPLICIT = explicit/allocate.o  explicit/minvrt.o explicit/rsdl.o  explicit/tstep.o


OBJ = $(BASE)
OBJ += $(OBJCD) $(OBJINS) $(OBJPS) $(OBJSWIRL) $(OBJBUOY) $(OBJSW) $(OBJLVLSET) $(OBJEXPLICIT)

DEFINES = -Df2cFortran $(EXTRAFLAGS)
CPPFLAGS = -g -O3 -I$(HOME)/Codes/include/ -I$(HOME)/Packages/include
LIB = -L$(HOME)/Codes/lib -L$(HOME)/Packages/lib -lquad -lspline -ltri_basis -lmyblas -lutil -linput_map -lbinio -lmuparser -lmetis -lm -lpth -lblitz
CXXMPI = ${CXX}

### Begin machine dependent

## Darwin:
## Machine specific for OS X
#CPPFLAGS += -O3 -Wno-long-double 
#LIB += -framework veclib

## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
#CPPFLAGS += -O3 -cxxlib-icc -I/opt/mpich/myrinet/intel/include #(cares)
#LIB += -L$(HOME) -L/usr/lib/gcc-lib/i386-pc-linux/3.2.3/ -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
#LIB +=    -L/opt/mpich/myrinet/intel/lib -lpmpich++ -lmpich -L/opt/gm/lib -lgm -lpthread #(cares)
#LDFLAGS += -static-libcxa #(cares)
CPPFLAGS += -I/opt/mpich/myrinet/gnu/include #(cares)
LIB += -L$(HOME) -latlas -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
CXXMPI = /opt/mpich/myrinet/gnu/bin/mpicxx

## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)
#LIB += -lcomplib.sgimath

## Aix:
## Machine specific for Jupiter
#CPPFLAGS += -O3
#LIB += -lessl \
/afs/cu/software/Development/lapack-3.0/rs_aix52/LAPACK/lapack_RS6K.a -lxlf -lxlf90
#DEFINES = -DIBMR2Fortran

### End machine dependent
CPPFLAGS += $(DEFINES)

tri_hp: dirs $(OBJ)
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh $(LIB)

tri_hp_mpi: dirs $(OBJ)
	${CXXMPI} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh_mpi $(LIB)

dirs:
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cd ] || mkdir -p cd
	+@[ -d  planestrain ] || mkdir -p planestrain
	+@[ -d  swirl ] || mkdir -p swirl
	+@[ -d  buoyancy ] || mkdir -p buoyancy
	+@[ -d  swe ] || mkdir -p swe
	+@[ -d  lvlset ] || mkdir -p lvlset
	+@[ -d  explicit ] || mkdir -p explicit

endif
