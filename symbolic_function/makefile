OBJ=main.o

LIB = -L$(HOME)/Packages/lib -lmuparser -L$(HOME)/Codes/lib -linput_map
CPPFLAGS += -I$(HOME)/Packages/include -I$(HOME)/Codes/include

# FOR CLUSTER
#CXX=icpc
#CPPFLAGS += -cxxlib-icc

all: header main

clean:
	rm main

header:
	cp symbolic_function.h ../include
	
main: $(OBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(OBJ) $(LIB)
	
$(OBJ): symbolic_function.h

