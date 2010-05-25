OBJ = input_map.o

CPPFLAGS = -O3

CPPFLAGS += -I$(HOME)/Packages/include
LIBS = -L$(HOME)/Packages/lib -lparseargs -lmuparser

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

mod_map: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIBS)
	cp mod_map ~/bin/

libinput_map.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f input_map.h $(HOME)/Codes/include/
	mv -f libinput_map.a $(HOME)/Codes/lib/
	
$(OBJ): input_map.h
