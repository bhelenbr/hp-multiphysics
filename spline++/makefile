OBJ = spline.o spline_functions2D.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lspline -lm
LIBS += ${LIBBLAS}

all: libspline.a spline

clean:
	rm *.o spline

libspline.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libspline.a ../lib/
	cp -f spline.h ../include/

spline: $(OBJ) main.o libspline.a
	${CXX} -o $@ main.o ${LIBS}
	cp $@ ../bin/$@

$(OBJ): spline.h spline_impl.h
