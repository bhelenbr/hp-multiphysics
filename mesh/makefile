OBJ = coarsen.o copy.o getnewblock.o insert.o out_mesh.o r_mesh.o \
tprims.o collapse.o findmatch.o getnewr_bdry.o rebay.o setup.o\
triangulate.o boundary.o connect.o getnewbdry.o in_mesh.o meshutilities.o \
refineby2.o swap.o yaber.o

MAINOBJ = blocks.o main.o 

LIB = -L$(HOME)/Codes/lib -lutil -lquad -lparseargs -linput_map -lm
CPPFLAGS = -I$(HOME)/Codes/include/

### Begin machine dependent

## Darwin:
## Machine specific for Darwin
#AR = libtool
#ARFLAGS = -o
#CPPFLAGS += -O3

## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
CXX=g++3 #(shelob)
CPPFLAGS += -O3

## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter
#CPPFLAGS += -O3

### End machine dependent

all: libmesh.a mblock

libmesh.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libmesh.a $(HOME)/Codes/lib/
	cp -f mesh.h $(HOME)/Codes/include/
	cp -f r_mesh.h $(HOME)/Codes/include/
	
mblock: $(OBJ) $(MAINOBJ)
	${CXX} -o $@ $(MAINOBJ) $(OBJ) $(LIB)

$(OBJ): mesh.h r_mesh.h


