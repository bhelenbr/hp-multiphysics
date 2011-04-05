export PACKAGES = /export/apps

DIRS = utilities input_map quadtree spline++ symbolic_function
TRI_DIRS = tri_basis tri_mesh tri_hp
TET_DIRS = tet_basis tet_mesh tet_hp

all: dirs tri_hp tet_hp 

tet_hp: tet_mesh tet_basis tri_mesh $(DIRS) force_look
	cd $@; $(MAKE) $(MFLAGS)

tet_mesh: tri_mesh $(DIRS) force_look
	cd $@; $(MAKE) $(MFLAGS)

tet_basis: utilities force_look
	cd $@; $(MAKE) $(MFLAGS)

tri_hp: tri_mesh tri_basis force_look
	cd $@; $(MAKE) $(MFLAGS)

tri_mesh: $(DIRS) utilities input_map quadtree spline++ symbolic_function force_look
	cd $@; $(MAKE) $(MFLAGS)

tri_basis: utilities force_look
	cd $@; $(MAKE) $(MFLAGS)

quadtree: utilities input_map force_look
	cd $@; $(MAKE) $(MFLAGS)

spline++: utilities input_map force_look
	cd $@; $(MAKE) $(MFLAGS)

symbolic_function: utilities input_map force_look
	cd $@; $(MAKE) $(MFLAGS)

input_map: utilities force_look
	cd $@; $(MAKE) $(MFLAGS)

utilities: force_look
	cd $@; $(MAKE) $(MFLAGS)
	
dirs: 
	+@[ -d  lib ] || mkdir -p lib
	+@[ -d  include ] || mkdir -p include

clean :
	for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done
	for d in $(TRI_DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done
	for d in $(TET_DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

force_look:
	true
