OBJU = utilities.o
OBJM =  myblas.o
CFLAGS = -O3

all:libutil.a libmyblas.a
	rm *.o

libutil.a: utilities.h $(OBJU)
	libtool --mode=link ${CC} -o $@ $(OBJU)
#	ar -r $@ $(OBJU)
	cp -f utilities.h $(HOME)/Codes/include/
	mv -f libutil.a $(HOME)/Codes/lib/

libmyblas.a: myblas.h $(OBJM)
	libtool --mode=link ${CC} -o $@ $(OBJM)
#	ar -r $@ $(OBJM)
	cp -f myblas.h $(HOME)/Codes/include/
	mv -f libmyblas.a $(HOME)/Codes/lib/
