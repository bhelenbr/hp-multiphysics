ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

# ShowMesh.App is an automator application and can be edited with
# Apple's automator.app.  The script has a hardwired path ${HOME}/Codes which may
# need to be changed
# Also, to make this the default way to open meshes, double click the .grd file
# then choose application and pick ShowMesh.app.

MAKETARGET = $(MAKE) -C $(OBJDIR) -f $(CURDIR)/makefile $(MAKECMDGOALS)

all: showmesh showgrid

showmesh:
	+@[ -d  _showmesh ] || mkdir -p _showmesh
	+@$(MAKE) -C _showmesh -f $(CURDIR)/makefile $@ EXTRAFLAGS='-DEASYMESH'
	cp _showmesh/showmesh ../bin/$@

.PHONY: showmesh

showgrid:
	+@[ -d  _showgrid ] || mkdir -p _showgrid
	+@$(MAKE) -C _showgrid -f $(CURDIR)/makefile $@
	cp _showgrid/showgrid ../bin/$@

.PHONY: showgrid

makefile : ;
%.mk :: ;

.PHONY: clean
clean:
	rm -rf _showmesh _showgrid
else

VPATH = ../

OBJ = showmesh2.o

DEFINES = $(EXTRAFLAGS)

LIB = -L/opt/X11/lib -lX11

CPPFLAGS = -O3 -g -I/opt/X11/include 
CPPFLAGS += $(DEFINES)

showmesh: $(OBJ)
	${CC} $(CPPFLAGS) -o $@ $(OBJ) $(LIB)

showgrid: $(OBJ)
	${CC} $(CPPFLAGS) -o $@ $(OBJ) $(LIB)

%.o: %.cpp
	$(CC) ${CPPFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

 endif
