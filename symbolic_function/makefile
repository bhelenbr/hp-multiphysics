OBJ=main.o

LIB = -L$(PACKAGES)/lib -lmuparser -L../lib -linput_map
CPPFLAGS += -I$(PACKAGES)/include -I../include

# FOR CLUSTER
#CXX=icpc
#CPPFLAGS += -cxxlib-icc
#CPPFLAGS += -debug (for icc)

main: $(OBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(OBJ) $(LIB)
	cp symbolic_function.h ../include

clean:
	rm main

$(OBJ): symbolic_function.h

