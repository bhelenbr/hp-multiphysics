OBJ = quad.o

CPPFLAGS = -I../include/ -I$(PACKAGES)/include

# CPPFLAGS += -debug -DBZ_DEBUG #(for icc)

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

all: libquad.a quadtest

clean:
	rm *.o quadtest

libquad.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f libquad.a ../lib/
	cp -f quadtree.h ../include/

quadtest: $(OBJ) main.o
	${CXX} -o $@ main.o -L../lib -lquad -lm

$(OBJ): quadtree.h quad_impl.h
