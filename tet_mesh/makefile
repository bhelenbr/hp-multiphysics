ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: tet_mesh  tet_mesh_mpi 

libtet_mesh.a:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DPTH'

.PHONY: libtet_mesh.a

tet_mesh:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DPTH'

.PHONY: tet_mesh

libtet_mesh_mpi.a:
	+@[ -d _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DMPISRC'

.PHONY: libtet_mesh_mpi.a

tet_mesh_mpi:
	+@[ -d _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DMPISRC'

.PHONY: tet_mesh_mpi

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _serial _mpi
else

VPATH = ../

OBJ = connect.o copy.o findmatch.o getnewbdry.o input.o \
meshutilities.o output.o refineby2.o setup.o setupbdry.o tet_boundary.o\
tet_prims.o vertexcircle.o test.o MAdLibInterface.o

MAINOBJ = main.o
DEFINES = -DMETIS -Df2cFortran  $(EXTRAFLAGS)
DEFINES = -DMETIS -Df2cFortran -DUSING_MADLIB $(EXTRAFLAGS)

LIB = -L$(HOME)/Codes/lib -L$(PACKAGES)/lib -ltri_mesh -lutil -lquad -lparseargs -linput_map -lspline -lmyblas
LIB += -lbinio -lm -lmuparser -lmetis -lpth -lblitz -lparseargs -lMAdLib

CPPFLAGS = -O3 -g -I$(HOME)/Codes/include/ -I$(PACKAGES)/include
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
#CPPFLAGS += -I/opt/mpich/myrinet/gnu/include #(cares)
CPPFLAGS += -I/opt/openmpi/include #(cares)
#CXXMPI = /opt/mpich/myrinet/gnu/bin/mpicxx
CXXMPI = /opt/openmpi/bin/mpicxx
LIB += -lblas -llapack #(cares)

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
	cp -f libtet_mesh.a $(HOME)/Codes/lib/
	cp -f ../*.h $(HOME)/Codes/include/

tet_mesh: libtet_mesh.a $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) -ltet_mesh $(LIB)
	cp tet_mesh $(HOME)/bin

libtet_mesh_mpi.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtet_mesh_mpi.a $(HOME)/Codes/lib/
	cp -f ../*.h $(HOME)/Codes/include/

tet_mesh_mpi: libtet_mesh_mpi.a $(MAINOBJ)
	${CXXMPI} $(CPPFLAGS) -o $@ $(MAINOBJ) -ltet_mesh_mpi $(LIB)
	cp tet_mesh_mpi $(HOME)/bin

%.o: %.cpp
	$(CXX) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

 endif
