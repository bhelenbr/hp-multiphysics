OBJU = utilities.o
OBJM =  myblas.o
CFLAGS = -O3

all:libutil.a libmyblas.a
	rm *.o

libutil.a: utilities.h $(OBJU)
	libtool --mode=link ${CC} -o $@ $(OBJU)
	cp -f utilities.h /home/helenbrk/codes/include/
	mv -f libutil.a /home/helenbrk/codes/lib/

libmyblas.a: myblas.h $(OBJM)
	libtool --mode=link ${CC} -o $@ $(OBJM)
	cp -f myblas.h /home/helenbrk/codes/include/
	mv -f libmyblas.a /home/helenbrk/codes/lib/
