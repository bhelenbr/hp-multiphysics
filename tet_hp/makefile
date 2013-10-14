ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: tet_hp tet_hp_mpi

tet_hp:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tet_hp EXTRAFLAGS='-DPTH'
	cp _$@/tet_hp $(HOME)/bin/$@

.PHONY: tet_hp

tet_hp_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tet_hp_mpi EXTRAFLAGS='-DMPISRC'
	cp _$@/tet_hp_mpi $(HOME)/bin/$@

.PHONY: tet_hp_mpi

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _tet_hp*

else

VPATH = ../

BASE = allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o mg_prolongate.o mg_restrict.o minvrt.o  \
tadvance.o output.o ptprobe.o particle.o tobasis.o update.o test.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o 
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/rsdl.o cd/setup_preconditioner.o 
#cd/minvrt.o
OBJCNS = cns/allocate.o cns/bdry.o cns/getnewbdry.o cns/getnewibc.o cns/update.o cns/rsdl.o cns/tstep.o 
OBJCNS_EXPLICIT = cns_explicit/allocate.o cns_explicit/bdry.o cns_explicit/getnewbdry.o cns_explicit/getnewibc.o cns_explicit/rsdl.o cns_explicit/tstep.o 
OBJMASS = mass/allocate.o mass/bdry.o mass/getnewbdry.o mass/getnewibc.o mass/rsdl.o mass/setup_preconditioner.o

OBJ = $(BASE)
OBJ += $(OBJCD) 
OBJ += $(OBJINS) 
OBJ += $(OBJCNS)
#OBJ += $(OBJCNS_EXPLICIT)
#OBJ += $(OBJMASS)

CXX=mpic++
DEFINES = -Df2cFortran $(EXTRAFLAGS)
CPPFLAGS = -O3 -g -I$(HOME)/Codes/include/ -I$(PACKAGES)/include
LIB = -L$(HOME)/Codes/lib -L$(PACKAGES)/lib -lquad -lspline -ltet_basis -lmyblas -lutil -linput_map -lbinio -lmuparser -lmetis -lm -lpth -lblitz -lMAdLib

### Begin machine dependent

## Linux:
LIB += -llapack -lblas 

## Darwin:
#LIB += -framework veclib

### End machine dependent

CPPFLAGS += $(DEFINES)

tet_hp: $(OBJ) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh -ltet_mesh $(LIB)

tet_hp_mpi: $(OBJ) 
	${CXXMPI} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh_mpi -ltet_mesh_mpi $(LIB)
	
dirs:
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cns ] || mkdir -p cns
	+@[ -d  cns_explicit ] || mkdir -p cns_explicit
	+@[ -d  cd ] || mkdir -p cd
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
