BASISOBJ = initialize.o intgrt.o intgrt1d.o proj.o \
proj1d.o ptvalues.o hpprobe.o

CPPFLAGS = -Df2cFortran -O3 -I$(HOME)/Codes/include/

libbasis.a: $(BASISOBJ)
	libtool --mode=link ${CC} -o $@ $(BASISOBJ)
	mv -f libbasis.a $(HOME)/Codes/lib/
	cp -f hpbasis.h $(HOME)/Codes/include/
	rm *.o

$(BASISOBJ): hpbasis.h
