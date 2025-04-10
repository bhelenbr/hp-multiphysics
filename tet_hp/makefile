ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: tet_hp tet_hp_mpi

tet_hp:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tet_hp EXTRAFLAGS='-DPTH'
	cp _$@/tet_hp ../bin/$@

.PHONY: tet_hp

tet_hp_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tet_hp_mpi EXTRAFLAGS='-DMPISRC' CXX=mpic++
	cp _$@/tet_hp_mpi ../bin/$@

.PHONY: tet_hp_mpi

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _tet_hp*

else

VPATH = ../

BASE = allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o mg_prolongate.o mg_restrict.o \
tadvance.o output.o ptprobe.o particle.o tobasis.o update.o test.o hp_getnewbdry.o setup_preconditioner.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/setup_preconditioner.o 
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/rsdl.o cd/setup_preconditioner.o 
#cd/minvrt.o
OBJCD_MULTI = cd_multi/allocate.o cd_multi/rsdl.o cd_multi/setup_preconditioner.o
OBJCNS = cns/allocate.o cns/bdry.o cns/getnewbdry.o cns/getnewibc.o cns/update.o cns/rsdl.o cns/setup_preconditioner.o 
OBJCNS_EXPLICIT = cns_explicit/allocate.o cns_explicit/bdry.o cns_explicit/getnewbdry.o cns_explicit/getnewibc.o cns_explicit/rsdl.o cns_explicit/setup_preconditioner.o 
OBJMASS = mass/allocate.o mass/bdry.o mass/getnewbdry.o mass/getnewibc.o mass/rsdl.o mass/setup_preconditioner.o

OBJ = $(BASE)
OBJ += $(OBJCD) 
OBJ += $(OBJCD_MULTI) 
OBJ += $(OBJINS) 
OBJ += $(OBJCNS)
#OBJ += $(OBJCNS_EXPLICIT)
#OBJ += $(OBJMASS)


DEFINES += $(EXTRAFLAGS)
CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lquad -lspline -ltet_basis -lmyblas -lutil -linput_map -lmuparser -lmetis -lm -lpth -lblitz -lMAdLib
LIBS += ${LIBBLAS}

tet_hp: $(OBJ) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh -ltet_mesh $(LIB)

tet_hp_mpi: $(OBJ) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh_mpi -ltet_mesh_mpi $(LIB)
	
dirs:
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cns ] || mkdir -p cns
	+@[ -d  cns_explicit ] || mkdir -p cns_explicit
	+@[ -d  cd ] || mkdir -p cd
	+@[ -d  cd_multi ] || mkdir -p cd_multi
	+@[ -d  mass ] || mkdir -p mass
	touch dirs
	
	
%.o: %.cpp dirs
	$(CXX) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

endif
