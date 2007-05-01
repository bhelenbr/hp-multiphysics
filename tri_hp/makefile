ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: planar axi

planar: hp hp_mpi
axi: hp_axi hp_axi_mpi

hp:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp
	cp _$@/hp $(HOME)/bin/$@

.PHONY: hp

hp_axi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp EXTRAFLAGS=-DAXISYMMETRIC
	cp _$@/hp $(HOME)/bin/$@

.PHONY: hp_axi

hp_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp_mpi EXTRAFLAGS=-DMPISRC
	cp _$@/hp_mpi $(HOME)/bin/$@

.PHONY: hp_mpi

hp_axi_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp_mpi EXTRAFLAGS='-DMPISRC -DAXISYMMETRIC'
	cp _$@/hp_mpi $(HOME)/bin/$@

.PHONY: hp_axi_mpi

hp_drop:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp EXTRAFLAGS='-DDROP -DAXISYMMETRIC'
	cp _$@/hp $(HOME)/bin/$@

.PHONY: hp_drop

hp_drop_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp_mpi EXTRAFLAGS='-DDROP -DAXISYMMETRIC -DMPISRC'
	cp _$@/hp_mpi $(HOME)/bin/$@

.PHONY: hp_drop_mpi

hp_drop_stokes:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp EXTRAFLAGS='-DDROP -DAXISYMMETRIC -DINERTIALESS'
	cp _$@/hp $(HOME)/bin/$@

.PHONY: hp_drop_stokes

hp_drop_stokes_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile hp_mpi EXTRAFLAGS='-DDROP -DAXISYMMETRIC -DINERTIALESS -DMPISRC'
	cp _$@/hp_mpi $(HOME)/bin/$@

.PHONY: hp_drop_stokes_mpi

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _serial _mpi _axi_serial _axi_mpi

else

VPATH = ../

BASE = adapt.o allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o minvrt.o movco.o movfin.o \
nstage.o output.o ptprobe.o tadvance.o tobasis.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o ins/surface.o
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/length.o cd/rsdl.o cd/tstep.o
OBJPS = planestrain/allocate.o planestrain/bdry.o planestrain/getnewbdry.o planestrain/length.o planestrain/rsdl.o planestrain/tstep.o
OBJSWIRL = swirl/allocate.o swirl/bdry.o swirl/getnewbdry.o swirl/getnewibc.o swirl/length.o swirl/rsdl.o swirl/tstep.o
OBJBUOY = buoyancy/allocate.o buoyancy/tstep.o buoyancy/rsdl.o
OBJSW = swe/allocate.o swe/getnewbdry.o swe/getnewibc.o swe/rsdl.o swe/tstep.o
OBJLVLSET = lvlset/allocate.o  lvlset/bdry.o lvlset/rsdl.o  lvlset/tstep.o


OBJ = $(BASE)
OBJ += $(OBJCD) $(OBJINS) $(OBJPS) $(OBJSWIRL) $(OBJBUOY) $(OBJSW) $(OBJLVLSET)

DEFINES = -Df2cFortran $(EXTRAFLAGS)
CPPFLAGS = -I$(HOME)/Codes/include/ -I$(HOME)/Packages/include
LIB = -L$(HOME)/Codes/lib -L$(HOME)/Packages/lib -lquad -lbasis -lmyblas -lutil -linput_map -lmuparser -lmetis -lm
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
#LIB +=   -L/opt/mpich/myrinet/intel/lib -lpmpich++ -lmpich -L/opt/gm/lib -lgm -lpthread #(cares)
#LDFLAGS += -static-libcxa #(cares)
CPPFLAGS += -O3 -I/opt/mpich/myrinet/gnu/include #(cares)
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

hp: dirs $(OBJ)
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -lmesh $(LIB)

hp_mpi: dirs $(OBJ)
	${CXXMPI} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -lmesh_mpi $(LIB)

dirs:
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cd ] || mkdir -p cd
	+@[ -d  planestrain ] || mkdir -p planestrain
	+@[ -d  swirl ] || mkdir -p swirl
	+@[ -d  buoyancy ] || mkdir -p buoyancy
	+@[ -d  swe ] || mkdir -p swe
	+@[ -d  lvlset ] || mkdir -p lvlset

endif
