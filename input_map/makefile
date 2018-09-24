OBJ = input_map.o

CPPFLAGS = -O3

#CPPFLAGS = -debug -DBZ_DEBUG #(for intel procs)

CPPFLAGS += -I$(PACKAGES)/include
LIBS = -L$(PACKAGES)/lib -L../lib -lparseargs -lmuparser

### Begin machine dependent

## Darwin:
## Machine specific for Darwin
#AR = libtool
#ARFLAGS = -o
 
## Linux:
##  Machine specific for Linux
#CXX=g++3 #(shelob)
#CXX=icpc #(cares)
#CPPFLAGS += -cxxlib-icc
#AR=libtool #(cares)
#ARFLAGS = --mode=link ${CC} -o #(cares)

## IRIX64:
#CPPFLAGS = -Ofast=IP30  #(sauter)
#CPPFLAGS = -Ofast=IP27  #(ctr-sgi1)
#CPPFLAGS = -Ofast=IPrk5 #(origin)

### End machine dependent
  
all: mod_map libinput_map.a

clean:
	rm mod_map *.o

mod_map: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIBS)
	cp mod_map.bash ../bin/mod_map

libinput_map.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f input_map.h ../include/
	cp -f libinput_map.a ../lib/
	
$(OBJ): input_map.h
