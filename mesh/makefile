MESHOBJ = coarsen.o collapse.o copy.o findmatch.o in_mesh.o insert.o out_mesh.o rebay.o  \
	  setup.o swap.o tprims.o triangulate.o yaber.o mvpttobdry.o connect.o r_mesh.o 

CPPFLAGS = -O3 -I/home/helenbrk/codes/include/

libmesh.a: $(MESHOBJ)
	libtool --mode=link ${CC} -o $@ $(MESHOBJ)
	mv -f libmesh.a /home/helenbrk/codes/lib/
	cp -f mesh.h /home/helenbrk/codes/include/
	cp -f r_mesh.h /home/helenbrk/codes/include/
	rm *.o

$(MESHOBJ): mesh.h r_mesh.h
