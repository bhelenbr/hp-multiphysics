OBJ = input_map.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)

LIBS += ${LIBDIRS}
LIBS += -lmuparser

all: mod_map libinput_map.a

clean:
	rm mod_map *.o

mod_map: $(OBJ) main.o
	${CXX} ${CPPFLAGS} -o $@ $(OBJ) main.o $(LIBS)

libinput_map.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)
	cp -f input_map.h ../include/
	cp -f libinput_map.a ../lib/
	
$(OBJ): input_map.h
