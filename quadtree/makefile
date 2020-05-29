OBJ = quad.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lquad -lm

all: libquad.a quadtest

clean:
	rm *.o quadtest

libquad.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libquad.a ../lib/
	cp -f quadtree.h ../include/

quadtest: $(OBJ) main.o
	${CXX} -o $@ main.o ${LIBS}

$(OBJ): quadtree.h quad_impl.h
