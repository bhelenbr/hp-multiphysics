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
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DMPISRC' CXX=mpic++

.PHONY: libtet_mesh_mpi.a

tet_mesh_mpi:
	+@[ -d _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DMPISRC' CXX=mpic++

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

DEFINES += -DMETIS -DUSING_MADLIB $(EXTRAFLAGS)
CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -ltri_mesh -lutil -lquad -linput_map -lspline -lmyblas
LIBS += -lm -lmuparser -lmetis -lpth -lblitz -lMAdLib
LIBS += ${LIBBLAS}

libtet_mesh.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtet_mesh.a ../../lib/
	cp -f ../*.h ../../include/

tet_mesh: libtet_mesh.a $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) -ltet_mesh $(LIBS)
	cp tet_mesh ../bin

libtet_mesh_mpi.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtet_mesh_mpi.a ../../lib/
	cp -f ../*.h ../../include/

tet_mesh_mpi: libtet_mesh_mpi.a $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) -ltet_mesh_mpi $(LIBS)
	cp tet_mesh_mpi ../bin

%.o: %.cpp
	$(CXX) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

 endif
