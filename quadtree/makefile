OBJ = quad.o
CPPFLAGS = -O3 -I/home/helenbrk/codes/include/

libquad.a: quadtree.h $(OBJ)
	libtool --mode=link ${CC} -o $@ $(OBJ)
	mv -f libquad.a /home/helenbrk/codes/lib/
	cp -f quadtree.h /home/helenbrk/codes/include/
	rm *.o
