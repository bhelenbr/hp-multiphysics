ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: libtet_mesh.a  tet_mesh 

libtet_mesh.a:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@

.PHONY: libtet_mesh.a

tet_mesh:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@

.PHONY: tet_mesh

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _serial _mpi
else

VPATH = ../

OBJ = connect.o copy.o findmatch.o getnewbdry.o input.o \
meshutilities.o output.o refineby2.o setup.o setupbdry.o tet_boundary.o\
tet_prims.o vertexcircle.o test.o

MAINOBJ = main.o
DEFINES = -DMETIS -Df2cFortran  $(EXTRAFLAGS)
#DEFINES = -DMETIS -Df2cFortran -DUSING_MADLIB $(EXTRAFLAGS)

LIB = -L$(HOME)/Codes/lib -L$(HOME)/Packages/lib -ltri_mesh -lutil -lquad -lparseargs -linput_map -lbinio -lm -lmuparser -lmetis -lspline -lg2c -lpth -lblitz 
CPPFLAGS = -O3 -g -I$(HOME)/Codes/include/ -I$(HOME)/Packages/include
CXXMPI = ${CXX}

### Begin machine dependent

## Darwin:
## Machine specific for Darwin
#AR = libtool
#ARFLAGS = -o

## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
#CXX=g++3 #(shelob)
#CPPFLAGS += -cxxlib-icc -I/opt/mpich/myrinet/intel/include #(cares)
#LIB +=    -L/opt/mpich/myrinet/intel/lib -lpmpich++ -lmpich -L/opt/gm/lib -lgm -lpthread #(cares)
CPPFLAGS += -I/opt/mpich/myrinet/gnu/include #(cares)
CXXMPI = /opt/mpich/myrinet/gnu/bin/mpicxx

## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter

### End machine dependent

CPPFLAGS += $(DEFINES)

libtet_mesh.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libtet_mesh.a $(HOME)/Codes/lib/
	cp -f ../*.h $(HOME)/Codes/include/

tet_mesh: $(OBJ) $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) $(OBJ) $(LIB)
	mv tet_mesh $(HOME)/bin


 endif
