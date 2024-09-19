ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: tri_mesh tri_mesh_mpi

libtri_mesh.a:
	+@[ -d  _serial ] || mkdir -p _serial 
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DPTH'

.PHONY: libtri_mesh.a

tri_mesh:
	+@[ -d  _serial ] || mkdir -p _serial
	+@$(MAKE) -C _serial -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DPTH'
	cp _serial/tri_mesh ../bin/$@

.PHONY: tri_mesh

libtri_mesh_mpi.a:
	+@[ -d  _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DMPISRC' CXX=${CXXMPI}

.PHONY: libtri_mesh_mpi.a

tri_mesh_mpi:
	+@[ -d  _mpi ] || mkdir -p _mpi
	+@$(MAKE) -C _mpi -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DMPISRC' CXX=${CXXMPI}
	cp _mpi/tri_mesh_mpi ../bin/$@

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
refineby2.o swap.o yaber.o adapt.o blocks.o mapped_mesh.o

MAINOBJ = main.o

DEFINES += -DMETIS -DMY_SPARSE $(EXTRAFLAGS)
CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lutil -lquad -linput_map -lspline -lmyblas
LIBS += -lm -lmuparser -lmetis -lparmetis -lpth -lblitz -lnetcdf
LIBS += ${LIBBLAS}

libtri_mesh.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtri_mesh.a ../../lib/
	cp -f ../*.h ../../include/

tri_mesh: libtri_mesh.a $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) -ltri_mesh $(LIBS)

libtri_mesh_mpi.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtri_mesh_mpi.a ../../lib/
	cp -f ../*.h ../../include/

tri_mesh_mpi: libtri_mesh_mpi.a $(MAINOBJ)
	${CXX} $(CPPFLAGS) -o $@ $(MAINOBJ) -ltri_mesh_mpi $(LIBS)

%.o: %.cpp
	$(CXX) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

 endif

