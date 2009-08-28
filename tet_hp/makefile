ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: planar 

planar: tet_hp

tet_hp:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile tet_hp
	cp _serial/tet_hp $(HOME)/bin/$@

.PHONY: tet_hp

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _serial _mpi _axi_serial _axi_mpi

else

VPATH = ../

BASE = allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o mg_prolongate.o mg_restrict.o minvrt.o  \
tadvance.o output.o ptprobe.o particle.o tobasis.o update.o test.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o 
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/rsdl.o cd/setup_preconditioner.o

OBJ = $(BASE)
OBJ += $(OBJCD) $(OBJINS) 

DEFINES = -Df2cFortran $(EXTRAFLAGS)
CPPFLAGS = -I$(HOME)/Codes/include/ -I$(HOME)/Packages/include
LIB = -L$(HOME)/Codes/lib -L$(HOME)/Packages/lib -ltri_mesh -lquad -ltet_basis -lmyblas -lutil 
-linput_map -lmuparser -lmetis -lm
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

tet_hp: dirs $(OBJ)
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltet_mesh $(LIB)


dirs:
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cd ] || mkdir -p cd
endif
