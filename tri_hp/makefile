ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: planar axi

planar: tri_hp tri_hp_mpi tri_hp_petsc
axi: tri_hp_axi tri_hp_axi_mpi tri_hp_axi_petsc

tri_hp:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp EXTRAFLAGS='-DPTH'
	cp _$@/tri_hp $(HOME)/bin/$@

.PHONY: tri_hp

tri_hp_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_mpi EXTRAFLAGS='-DMPISRC' CXX=mpic++
	cp _$@/tri_hp_mpi $(HOME)/bin/$@

.PHONY: tri_hp_mpi

tri_hp_petsc:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_petsc EXTRAFLAGS='-DMPISRC -Dpetsc' CXX=mpic++
	cp _$@/tri_hp_petsc $(HOME)/bin/$@

.PHONY: tri_hp_petsc

tri_hp_axi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp EXTRAFLAGS='-DAXISYMMETRIC -DPTH'
	cp _$@/tri_hp $(HOME)/bin/$@

.PHONY: tri_hp_axi

tri_hp_axi_mpi:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_mpi EXTRAFLAGS='-DMPISRC -DAXISYMMETRIC' CXX=mpic++
	cp _$@/tri_hp_mpi $(HOME)/bin/$@

.PHONY: tri_hp_axi_mpi

tri_hp_axi_petsc:
	+@[ -d  _$@ ] || mkdir -p _$@
	+@$(MAKE) -C _$@ -f $(CURDIR)/makefile tri_hp_petsc EXTRAFLAGS='-DAXISYMMETRIC -Dpetsc -DMPISRC' CXX=mpic++
	cp _$@/tri_hp_petsc $(HOME)/bin/$@

.PHONY: tri_hp_axi_petsc

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _tri_hp*

else

VPATH = ../

BASE = adapt.o allocate.o copy.o getnewblock.o gtol.o hp_boundary.o hp_getnewibc.o l2error.o main.o minvrt.o movco.o movfin.o \
nstage.o output.o ptprobe.o tadvance.o tobasis.o hp_coupled_boundary.o
OBJINS = ins/allocate.o ins/bdry.o ins/getnewbdry.o ins/getnewibc.o ins/length.o ins/rsdl.o ins/tstep.o ins/surface.o ins/surface2.o
OBJCNS = cns/allocate.o cns/bdry.o cns/getnewbdry.o cns/getnewibc.o cns/length.o cns/rsdl.o cns/tstep.o cns/update.o
OBJCNSEXPLICIT = cns_explicit/allocate.o cns_explicit/bdry.o cns_explicit/getnewbdry.o cns_explicit/getnewibc.o cns_explicit/length.o cns_explicit/rsdl.o cns_explicit/tstep.o cns/minvrt.o
OBJCD = cd/allocate.o cd/bdry.o cd/getnewbdry.o cd/getnewibc.o cd/length.o cd/rsdl.o cd/tstep.o # cd/minvrt.o
OBJPS = planestrain/allocate.o planestrain/bdry.o planestrain/getnewbdry.o planestrain/length.o planestrain/rsdl.o planestrain/tstep.o
OBJSWIRL = swirl/allocate.o swirl/bdry.o swirl/getnewbdry.o swirl/getnewibc.o swirl/length.o swirl/rsdl.o swirl/tstep.o
OBJBUOY = buoyancy/allocate.o buoyancy/tstep.o buoyancy/rsdl.o buoyancy/length.o buoyancy/bdry.o buoyancy/melt.o buoyancy/surface.o buoyancy/melt_kinetics.o
OBJSW = swe/allocate.o swe/getnewbdry.o swe/bdry.o swe/getnewibc.o swe/rsdl.o swe/tstep.o
OBJLVLSET = lvlset/allocate.o  lvlset/bdry.o lvlset/rsdl.o  lvlset/tstep.o lvlset/reinit.o
OBJEXPLICIT = explicit/allocate.o  explicit/minvrt.o explicit/rsdl.o  explicit/tstep.o
OBJNONNEWTONIAN = nonnewtonian/tri_hp_nonnewtonian.o
OBJPETSC = petsc.o

# COMMENT OUT LINES HERE FOR FASTER COMPILING AND THEN UNDEFINE IN GETNEWBLOCK.CPP
OBJ = $(BASE)
OBJ += $(OBJCD)
OBJ += $(OBJINS)
#OBJ += $(OBJCNS) 
#OBJ += $(OBJCNSEXPLICIT)
#OBJ += $(OBJLVLSET)
#OBJ += $(OBJBUOY)
#OBJ += $(OBJSW)   
#OBJ += $(OBJEXPLICIT)
#OBJ += $(OBJPS)
#OBJ += $(OBJSWIRL)
#OBJ += $(OBJNONNEWTONIAN)

DEFINES = -Df2cFortran -DMY_SPARSE $(EXTRAFLAGS)
CPPFLAGS = -O3 -g -I$(HOME)/Codes/include/ -I$(PACKAGES)/include

LIB = -L$(HOME)/Codes/lib -L$(PACKAGES)/lib -lquad -lspline -ltri_basis -lmyblas -lutil -linput_map -lbinio -lmuparser -lmetis -lm -lpth -lblitz
LIBPETSC = -lpetsc -ldmumps -lmumps_common -lpord -lHYPRE -lscalapack -lsuperlu_dist_3.3 -lparmetis -ldl

### Begin machine dependent

## Darwin:
## Machine specific for OS X
#CPPFLAGS += -O3 -Wno-long-double 
#LIB += -framework veclib

## Linux:
LIB += -llapack -lblas #(cares)

## Intel Compilers & MKL
#LIB += -lmkl -lmkl_lapack #(hawk)

### End machine dependent

CPPFLAGS += $(DEFINES)

tri_hp: $(OBJ) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh $(LIB)

tri_hp_mpi: $(OBJ) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) -ltri_mesh_mpi $(LIB)

tri_hp_petsc: $(OBJ) $(OBJPETSC) 
	${CXX} $(LDFLAGS) -o $@ $(CPPFLAGS) $(OBJ) $(OBJPETSC) -ltri_mesh_mpi $(LIB) $(LIBPETSC) 

dirs: 
	+@[ -d  ins ] || mkdir -p ins
	+@[ -d  cns ] || mkdir -p cns
	+@[ -d  cns_explicit ] || mkdir -p cns_explicit
	+@[ -d  cd ] || mkdir -p cd
	+@[ -d  planestrain ] || mkdir -p planestrain
	+@[ -d  swirl ] || mkdir -p swirl
	+@[ -d  buoyancy ] || mkdir -p buoyancy
	+@[ -d  swe ] || mkdir -p swe
	+@[ -d  lvlset ] || mkdir -p lvlset
	+@[ -d  explicit ] || mkdir -p explicit
	+@[ -d  nonnewtonian ] || mkdir -p nonnewtonian
	touch dirs


%.o: %.cpp dirs
	$(CXX) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

endif
