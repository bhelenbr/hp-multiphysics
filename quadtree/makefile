OBJ = quad.o
CPPFLAGS = -O3 -I$(HOME)/Codes/include/
CCFLAGS = -O3 -I$(HOME)/Codes/include/

#CXX = g++3

libquad.a: quadtree.h $(OBJ)
#	libtool --mode=link ${CC} -o $@ $(OBJ)
	ar -r $@ $(OBJ)
	mv -f libquad.a $(HOME)/Codes/lib/
	cp -f quadtree.h $(HOME)/Codes/include/

quadtest: libquad.a main.o
	${CCC} -o $@ main.o -L$(HOME)/Codes/lib -lquad -lm
