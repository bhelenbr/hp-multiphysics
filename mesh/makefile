OBJ = coarsen.o copy.o getnewblock.o insert.o out_mesh.o r_mesh.o \
tprims.o collapse.o findmatch.o getnewr_bdry.o rebay.o setup.o\
triangulate.o boundary.o connect.o getnewbdry.o in_mesh.o meshutilities.o \
refineby2.o swap.o yaber.o adapt.o blocks.o

MAINOBJ = main.o 

DEFINES = -DMETIS $(EXTRAFLAGS)

LIB = -L$(HOME)/Codes/lib -L$(HOME)/Packages/lib -lutil -lquad -lparseargs -linput_map -lm -lmuparser -lmetis
CPPFLAGS = -I$(HOME)/Codes/include/ -I$(HOME)/Packages/include

### Begin machine dependent

## Darwin:
## Machine specific for Darwin
#AR = libtool
#ARFLAGS = -o
#CPPFLAGS += -O3

## Linux:
##  Machine specific for Linux
CXX=icpc #(cares)
#CXX=g++3 #(shelob)
#CPPFLAGS += -O3
CPPFLAGS += -O3 -cxxlib-icc -I/opt/mpich/myrinet/intel/include #(cares)
LIB +=   -L/opt/mpich/myrinet/intel/lib -lpmpich++ -lmpich -L/opt/gm/lib -lgm -lpthread #(cares)

## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter
#CPPFLAGS += -O3

### End machine dependent

CPPFLAGS += $(DEFINES)

all: libmesh.a mesh
all_mpi: libmesh_mpi.a mesh_mpi

libmesh.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libmesh.a $(HOME)/Codes/lib/
	cp -f *.h $(HOME)/Codes/include/
	
mesh: $(OBJ) $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) $(OBJ) $(LIB)
	mv mesh $(HOME)/bin
	
libmesh_mpi.a: $(OBJ) 
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libmesh_mpi.a $(HOME)/Codes/lib/
	cp -f *.h $(HOME)/Codes/include/

mesh_mpi: $(OBJ) $(MAINOBJ) 
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) $(OBJ) $(LIB)
	mv mesh_mpi $(HOME)/bin
