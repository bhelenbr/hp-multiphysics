OBJ = spline.o spline_functions2D.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lspline -lm
LIBS += ${LIBBLAS}

all: libspline.a splinetest

clean:
	rm *.o splinetest

libspline.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libspline.a ../lib/
	cp -f spline.h ../include/

splinetest: $(OBJ) main.o libspline.a
	${CXX} -o $@ main.o ${LIBS}

$(OBJ): spline.h spline_impl.h
