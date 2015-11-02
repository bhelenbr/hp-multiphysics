OBJ = instances.o

CPPFLAGS = -I../include/ -I$(PACKAGES)/include
CPPFLAGS += -fno-inline
#CPPFLAGS += -debug #(for icc)
DEFINES = -Df2cFortran
LIB = -L../lib -L$(PACKAGES)/lib -lmyblas -lm



### Begin machine dependent

## Darwin:
## Machine specific for OS X
#AR = libtool
#ARFLAGS = -o
#CPPFLAGS += -O3 -funroll-loops -falign-loops=16
#LIB += -framework veclib
 
## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
#CXX=g++3 #(shelob)
CPPFLAGS += -O3
#CPPFLAGS += -O3 -cxxlib-icc
#LIB += /usr/lib/liblapack.so.3gf /usr/lib/libatlas.so.3gf  #(ubuntu)
#LIB += -L. -latlas -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
#DEFINES = -DpgiFortran
#LIB += -L. -latlas -llapack -lcblas -lf77blas -latlas -lg2c #(cares)
LIB += -llapack -lblas
#LIB += -lmkl -lmkl_lapack #(for icc)


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

all: tri_basis libtri_basis.a 

clean:
	rm -f *.o *~

tri_basis: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIB)

libtri_basis.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtri_basis.a ../lib/
	cp -f tri_basis.h ../include/

$(OBJ): tri_basis.h

