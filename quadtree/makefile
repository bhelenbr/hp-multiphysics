OBJ = quad.o
CPPFLAGS = -O3 -I$(HOME)/Codes/include/

#CXX = g++3

libquad.a: quadtree.h $(OBJ)
	libtool --mode=link ${CC} -o $@ $(OBJ)
	mv -f libquad.a $(HOME)/Codes/lib/
	cp -f quadtree.h $(HOME)/Codes/include/
	rm *.o
