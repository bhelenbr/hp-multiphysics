MESHOBJ = boundary.o coarsen.o collapse.o connect.o copy.o findmatch.o in_mesh.o insert.o \
main.o out_mesh.o rboundary.o rebay.o refineby2.o r_mesh.o setup.o swap.o tprims.o \
triangulate.o yaber.o

CPPFLAGS = -g -I/home/helenbrk/codes/include/

mesh: $(MESHOBJ)
	${CXX} ${CPPFLAGS} -o $@ $(MESHOBJ) ../lib/libquad.a ../lib/libutil.a 

libmesh.a: $(MESHOBJ)
	libtool --mode=link ${CC} -o $@ $(MESHOBJ)
	mv -f libmesh.a /home/helenbrk/codes/lib/
	cp -f mesh.h /home/helenbrk/codes/include/
	cp -f r_mesh.h /home/helenbrk/codes/include/
	rm *.o

$(MESHOBJ): mesh.h r_mesh.h
