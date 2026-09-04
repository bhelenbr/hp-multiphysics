ifeq (,$(filter _%,$(notdir $(CURDIR))))

.SUFFIXES:

# ShowMesh.App is an automator application and can be edited with
# Apple's automator.app.  The script has a hardwired path ${HOME}/Codes which may
# need to be changed
# Also, to make this the default way to open meshes, double click the .grd file
# then choose application and pick ShowMesh.app.

# Makefile for the SDL2 port of ShowMesh
# Works on both macOS (Homebrew) and Linux, as long as `pkg-config`
# can find sdl2 and SDL2_ttf.
#
# macOS setup:   brew install sdl2 sdl2_ttf pkg-config
# Debian/Ubuntu: sudo apt install build-essential pkg-config libsdl2-dev libsdl2-ttf-dev
# Fedora:        sudo dnf install gcc pkgconf-pkg-config SDL2-devel SDL2_ttf-devel

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

# OBJ = showmesh.o
# LIB = -L/opt/X11/lib -lX11
# CFLAGS = -O3 -g -I/opt/X11/include 
# CFLAGS += $(DEFINES) ${EXTRAFLAGS}

OBJ = showmesh_sdl.o
PKGS    = sdl2 SDL2_ttf
CFLAGS  = -O2 -Wall $(shell pkg-config --cflags $(PKGS))
LIB  = $(shell pkg-config --libs $(PKGS)) -lm
CFLAGS  += $(DEFINES) ${EXTRAFLAGS}

showmesh: $(OBJ)
	${CC} $(CFLAGS) -o $@ $(OBJ) $(LIB)

showgrid: $(OBJ)
	${CC} $(CFLAGS) -o $@ $(OBJ) $(LIB)

%.o: %.cpp
	$(CC) ${CFLAGS} -MD -c -o $@ $<
	@cp $*.d $*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	rm -f $*.d

-include $(OBJ:.o=.P)

 endif
