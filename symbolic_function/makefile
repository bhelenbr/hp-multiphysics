OBJ=main.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lmuparser -linput_map

main: $(OBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(OBJ) $(LIBS)
	cp symbolic_function.h ../include

clean:
	rm main

$(OBJ): symbolic_function.h

