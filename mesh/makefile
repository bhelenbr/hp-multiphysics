MESHOBJ = coarsen.o ldmesh.o metric.o triangulate.o insert.o main.o out_mesh.o tprims.o
QTREEOBJ = quad.o
UTOBJ = xmalloc.o

unstructured.exe: $(QTREEOBJ) $(MESHOBJ) $(UTOBJ)
	c++ -o $@ $(MESHOBJ) $(QTREEOBJ) $(UTOBJ)

$(MESHOBJ): mesh.h

$(QTREEOBJ): quadtree.h

$(UTOBJ): utilities.h
