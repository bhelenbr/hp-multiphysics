BASISOBJ = initialize.o intgrt.o intgrt1d.o proj.o \
proj1d.o ptvalues.o hpprobe.o

CPPFLAGS = -O3 -I/home/helenbrk/codes/include/

libbasis.a: $(BASISOBJ)
	libtool --mode=link ${CC} -o $@ $(BASISOBJ)
	mv -f libbasis.a /home/helenbrk/codes/lib/
	cp -f basis.h /home/helenbrk/codes/include/
	rm *.o

$(BASISOBJ): hpbasis.h
