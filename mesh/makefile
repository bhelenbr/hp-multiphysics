OBJ = blocks.o connect.o insert.o out_mesh.o r_mesh.o tprims.o \
copy.o setup.o triangulate.o \
coarsen.o findmatch.o rebay.o swap.o yaber.o \
collapse.o in_mesh.o refineby2.o tprims3d.o generate.o

MAINOBJ = main.o locals.o dummy.o 

CXX = g++3

CPPFLAGS = -DCAPRI -DUNDERSCORE $(CAPRIOPTS) -DBO_DEC -DLIN -DLININTEL -I/Users/helenbrk/Codes/include/ -I/usr/local/OCCapri/include
LIBPATH = -L/home/helenbrk/Codes/lib -L$(CASROOT)/Linux/lib -L/usr/local/OCCapri/lib
LIBS = -lquad -lutil -linput_map -lcapri -lFWOSPlugin -lPTKernel -lTKBool -lTKCAF -lTKCDF -lTKernel -lTKFeat -lTKGeomAlgo -lTKMath -lTKOffset -lTKPCAF -lTKPrim -lTKPShape -lTKService -lTKTopAlgo -lTKBRep -lTKG2d -lTKG3d -lTKGeomBase -lTKShHealing

mesh: $(OBJ) $(MAINOBJ) block.h  blocks.h  boundary.h  mboundary.h  mesh.h  rblock.h  rboundary.h  r_mesh.h  sblock.h
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) $(MAINOBJ) $(LIBPATH) $(LIBS) 

libmesh.a: $(MESHOBJ)
	libtool --mode=link ${CC} -o $@ $(MESHOBJ)
	mv -f libmesh.a /home/helenbrk/Codes/lib/
	cp -f *.h /home/helenbrk/Codes/include/

