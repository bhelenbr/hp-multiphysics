OBJU = utilities.o
OBJM = myblas.o matrix.o

CPPFLAGS += -Df2cFortran -I$(PACKAGES)/include

### Begin machine dependent

## Darwin:
## Machine specific for OS X
#AR = libtool
#ARFLAGS = -o
#CFLAGS += -O3
 
## Linux:
##  Machine specific for Linux
# CC=icc #(cares)
#CC=gcc3 #(shelob)
CFLAGS += -O3

## IRIX64:
#CFLAGS += -Ofast=IP30  #(sauter)
#CFLAGS += -Ofast=IP27  #(ctr-sgi1)
#CFLAGS += -Ofast=IPrk5 #(origin)

## Aix:
## Machine specific for Jupiter
#CFLAGS += -O3

### End machine dependent


all:libutil.a libmyblas.a

clean:
	rm *.o

libutil.a: utilities.h $(OBJU)
	$(AR) $(ARFLAGS) $@ $(OBJU)
	cp -f utilities.h ../include/
	cp -f libutil.a ../lib/

libmyblas.a: myblas.h $(OBJM)
	$(AR) $(ARFLAGS) $@ $(OBJ) $(OBJM)
	cp -f myblas.h ../include/
	cp -f libmyblas.a ../lib/
