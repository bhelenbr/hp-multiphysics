OBJU = utilities.o
OBJM =  myblas.o
OBJP = parseargs.o

### Begin machine dependent

## Darwin:
## Machine specific for OS X
#AR = libtool
#ARFLAGS = -o
#CFLAGS += -O3
 
## Linux:
##  Machine specific for Linux
#CC=icc #(cares)
CC=gcc3 #(shelob)
CFLAGS += -O3

## IRIX64:
#CFLAGS += -Ofast=IP30  #(sauter)
#CFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter
#CFLAGS += -O3

### End machine dependent


all:libutil.a libmyblas.a libparseargs.a
	rm *.o

libutil.a: utilities.h $(OBJU)
	$(AR) $(ARFLAGS) $@ $(OBJU)
	cp -f utilities.h $(HOME)/Codes/include/
	mv -f libutil.a $(HOME)/Codes/lib/

libmyblas.a: myblas.h $(OBJM)
	$(AR) $(ARFLAGS) $@ $(OBJ) $(OBJM)
	cp -f myblas.h $(HOME)/Codes/include/
	mv -f libmyblas.a $(HOME)/Codes/lib/
	
libparseargs.a: parseargs.h $(OBJP)
	$(AR) $(ARFLAGS) $@ $(OBJ) $(OBJP)
	cp -f parseargs.h $(HOME)/Codes/include/
	cp -f gtypes.h $(HOME)/Codes/include/
	mv -f libparseargs.a $(HOME)/Codes/lib/
