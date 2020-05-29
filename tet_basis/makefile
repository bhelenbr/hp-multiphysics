OBJ = initialize.o intgrt.o intgrt2d.o intgrt1d.o proj.o \
proj1d.o proj2d.o ptvalues.o hpprobe.o outputstuff.o 

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS} 
LIBS += -lmyblas -lm
LIBS += ${LIBBLAS}

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

all: tet_basis libtet_basis.a 

clean:
	rm tet_basis *.o

tet_basis: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIBS)

libtet_basis.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtet_basis.a ../lib/
	cp -f tet_basis.h ../include/

$(OBJ): tet_basis.h
