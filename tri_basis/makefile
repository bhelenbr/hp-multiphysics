OBJ = instances.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lmyblas -lm
LIBS += ${LIBBLAS}

all: tri_basis libtri_basis.a 

clean:
	rm -f *.o *~

tri_basis: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIBS)

libtri_basis.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libtri_basis.a ../lib/
	cp -f tri_basis.h ../include/

$(OBJ): tri_basis.h

