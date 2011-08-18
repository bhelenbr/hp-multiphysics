ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: planar 

planar: tet_hp tet_hp_mpi

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
	rm -rf _tri_hp*

else

VPATH = ../

BASE = allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o mg_prolongate.o mg_restrict.o minvrt.o  \
tadvance.o output.o ptprobe.o particle.o tobasis.o update.o test.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o 
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/rsdl.o cd/setup_preconditioner.o
OBJCNS = cns/allocate.o cns/bdry.o cns/getnewbdry.o cns/getnewibc.o cns/update.o cns/rsdl.o cns/tstep.o 
OBJCNS_EXPLICIT = cns_explicit/allocate.o cns_explicit/bdry.o cns_explicit/getnewbdry.o cns_explicit/getnewibc.o cns_explicit/rsdl.o cns_explicit/tstep.o 

OBJ = $(BASE)
OBJ += $(OBJCD) 
OBJ += $(OBJINS) 
OBJ += $(OBJCNS)
OBJ += $(OBJCNS_EXPLICIT)

DEFINES = -Df2cFortran $(EXTRAFLAGS)
CPPFLAGS = -O3 -g -Wno-trigraphs -Wreturn-type -Wunused-variable -Wuninitialized -I$(HOME)/Codes/include/ -I$(PACKAGES)/include
LIB = -L$(HOME)/Codes/lib -L$(PACKAGES)/lib -lquad -lspline -ltet_basis -lmyblas -lutil -linput_map -lbinio -lmuparser -lmetis -lm -lpth -lblitz

CXXMPI = ${CXX}
### Begin machine dependent

## Darwin:
## Machine specific for OS X
#CPPFLAGS += -O3 -Wno-long-double 
#LIB += -framework veclib

## Linux:
LIB += -L$(HOME)  -llapack -lblas #(cares)
CPPFLAGS += -I/opt/mpich/myrinet/gnu/include #(cares)
CXXMPI = /opt/mpich/myrinet/gnu/bin/mpicxx

## To use openmpi
#CPPFLAGS += -I/opt/openmpi/include
#CXXMPI = /opt/openmpi/bin/mpic++

# (ubuntu) should really make symbolic links to liblapack and libatlas
#LIB += /usr/lib/liblapack.so.3gf /usr/lib/libatlas.so.3gf  #(ubuntu)

### End machine dependent

CPPFLAGS += $(DEFINES)

tet_hp: dirs $(OBJ) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh -ltet_mesh $(LIB)

tet_hp_mpi: dirs $(OBJ) 
	${CXXMPI} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh_mpi -ltet_mesh_mpi $(LIB)
	
dirs:
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cns ] || mkdir -p cns
	+@[ -d  cns_explicit ] || mkdir -p cns_explicit
	+@[ -d  cd ] || mkdir -p cd

%.o: %.cpp
	$(CXX) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

endif