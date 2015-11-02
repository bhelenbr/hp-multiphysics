OBJ=parseargs.o

# FOR INTEL COMPILER
# CC=icc
# CFLAGS="-cxxlib-icc"

libparseargs.a: $(OBJ)
	cp *.h ../include
	ar rv $@ $(OBJ)
	ranlib $@
	cp $@ ../lib
	
clean:
	rm *.o
	rm libparseargs.a
