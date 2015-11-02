OBJ = spline.o

CPPFLAGS = -Df2cFortran -I${HOME}/Codes/include/ -I$(PACKAGES)/include/

# CPFLAGS += -debug -DBZ_DEBUG #(for icc)

### Begin machine dependent

## Darwin:
## Machine specific for Darwin
#AR = libtool
#ARFLAGS = -o
#CPPFLAGS += -O3

## Linux:
##  Machine specific for Linux
#CXX=icpc #(cares)
#CXX=g++3 #(shelob)
CPPFLAGS += -O3

## IRIX64:
#CPPFLAGS += -Ofast=IP30  #(sauter)
#CPPFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter
#CPPFLAGS += -O3

### End machine dependent

all: libspline.a splinetest

clean:
	rm *.o splinetest

libspline.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libspline.a ../lib/
	cp -f spline.h ../include/

splinetest: $(OBJ) main.o libspline.a
	${CXX} -o $@ main.o -L../lib -lspline -llapack -lblas -lm

$(OBJ): spline.h spline_impl.h
