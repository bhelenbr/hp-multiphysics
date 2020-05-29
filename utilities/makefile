OBJU = utilities.o
OBJM = myblas.o matrix.o

CPPFLAGS += ${OPT}
CPPFLAGS += ${INCLUDEDIRS}
CPPFLAGS += $(DEFINES)
CFLAGS += ${OPT} ${INCLUDIRS} ${DEFINES}

all:libutil.a libmyblas.a

clean:
	rm *.o

libutil.a: utilities.h $(OBJU)
	$(AR) $(ARFLAGS) $@ $(OBJU)
	cp -f utilities.h ../include/
	cp -f libutil.a ../lib/

libmyblas.a: myblas.h $(OBJM)
	$(AR) $(ARFLAGS) $@ $(OBJ) $(OBJM)
	cp -f myblas.h ../include/
	cp -f libmyblas.a ../lib/
