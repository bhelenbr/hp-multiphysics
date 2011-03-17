OBJ = quad.o

CPPFLAGS = -I$(HOME)/Codes/include/

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
	cp -f libquad.a $(HOME)/Codes/lib/
	cp -f quadtree.h $(HOME)/Codes/include/

quadtest: $(OBJ) main.o
	${CXX} -o $@ main.o -L$(HOME)/Codes/lib -lquad -lm

$(OBJ): quadtree.h quad_impl.h
