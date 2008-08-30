ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: libtri_mesh.a libtri_mesh_mpi.a tri_mesh tri_mesh_mpi

libtri_mesh.a:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@

.PHONY: libtri_mesh.a

tri_mesh:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@

.PHONY: tri_mesh

libtri_mesh_mpi.a:
	+@[ -d  _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS=-DMPISRC

.PHONY: libtri_mesh_mpi.a

tri_mesh_mpi:
	+@[ -d  _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS=-DMPISRC

.PHONY: tri_mesh_mpi

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _serial _mpi
else

VPATH = ../

OBJ = coarsen.o copy.o getnewblock.o insert.o output.o r_mesh.o \
tprims.o collapse.o findmatch.o getnewr_bdry.o rebay.o setup.o\
triangulate.o tri_boundary.o connect.o getnewbdry.o input.o meshutilities.o \
refineby2.o swap.o yaber.o adapt.o blocks.o

MAINOBJ = main.o 

DEFINES = -DMETIS -Df2cFortran $(EXTRAFLAGS)

LIB = -L$(HOME)/Codes/lib -L$(HOME)/Packages/lib -lutil -lquad -lparseargs -linput_map -lbinio -lm -lmuparser -lmetis -lspline -lg2c -lpth 
CPPFLAGS = -I$(HOME)/Codes/include/ -I$(HOME)/Packages/include
CXXMPI = ${CXX}

### Begin machine dependent

## Darwin:
## Machine specific for Darwin
#AR = libtool
#ARFLAGS = -o
CPPFLAGS += -O3

## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
#CXX=g++3 #(shelob)
#CPPFLAGS += -O3
#CPPFLAGS += -O3 -cxxlib-icc -I/opt/mpich/myrinet/intel/include #(cares)
#LIB +=    -L/opt/mpich/myrinet/intel/lib -lpmpich++ -lmpich -L/opt/gm/lib -lgm -lpthread #(cares)
CPPFLAGS += -O3 -I/opt/mpich/myrinet/gnu/include #(cares)
CXXMPI = /opt/mpich/myrinet/gnu/bin/mpicxx

## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter
#CPPFLAGS += -O3

### End machine dependent

CPPFLAGS += $(DEFINES)

libtri_mesh.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libtri_mesh.a $(HOME)/Codes/lib/
	cp -f ../*.h $(HOME)/Codes/include/
	
tri_mesh: $(OBJ) $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) $(OBJ) $(LIB)
	mv tri_mesh $(HOME)/bin
	
libtri_mesh_mpi.a: $(OBJ) 
	$(AR) $(ARFLAGS) $@ $(OBJ)
	mv -f libtri_mesh_mpi.a $(HOME)/Codes/lib/
	cp -f ../*.h $(HOME)/Codes/include/

tri_mesh_mpi: $(OBJ) $(MAINOBJ) 
	${CXXMPI} $(CPPFLAGS) -o $@ $(MAINOBJ) $(OBJ) $(LIB)
	mv tri_mesh_mpi $(HOME)/bin

 endif
