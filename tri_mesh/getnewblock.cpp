#include "blocks.h"
#include "block.h"
#include "r_tri_mesh.h"
#include "mapped_mesh.h"

class btype {
    public:
        const static int ntypes = 4;
        enum ids {r_tri_mesh,spline_mapped_mesh,polar_mapped_mesh,polar_log_mapped_mesh};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i)
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};
const char btype::names[ntypes][40] = {"r_tri_mesh","spline_mapped_mesh","polar_mapped_mesh","polar_log_mapped_mesh"};

multigrid_interface* block::getnewlevel(input_map& inmap) {
    std::string keyword,val,ibcname;
    std::istringstream data;
    int type;
    
    /* FIND BLOCK TYPE */
    if (inmap.get(idprefix+"_type",val)) {
        type = btype::getid(val.c_str());
    }
    else {
        if (inmap.get("blocktype",val)) {
            type = btype::getid(val.c_str());
        }
        else {
            type = btype::r_tri_mesh;
        }
    }
    
    switch(type) {
        case btype::r_tri_mesh: {
            r_tri_mesh *temp = new r_tri_mesh();
            return(temp);
        }
        case btype::spline_mapped_mesh: {
            mapped_mesh *temp = new mapped_mesh();
            temp->map = make_shared<spline_mapping>();
            return(temp);
        }
        case btype::polar_mapped_mesh: {
            mapped_mesh *temp = new mapped_mesh();
            temp->map = make_shared<polar_mapping>();
            return(temp);
        }
        case btype::polar_log_mapped_mesh: {
            mapped_mesh *temp = new mapped_mesh();
            temp->map = make_shared<polar_log_mapping>();
            return(temp);
        }
        default: {
            r_tri_mesh *temp = new r_tri_mesh();
            return(temp);
        }
    }
    return(0);
}

/* This routine waits for everyone to exit nicely */
void sim::finalize(int line,const char *file, std::ostream *log) {
    *log << "Exiting at line " << line << " of file " << file << std::endl;
#ifdef PTH
    pth_exit(NULL);
#endif
#ifdef BOOST
    throw boost::thread_interrupted();
#endif
#ifdef PTH
    pth_kill();
#endif
#ifdef petsc
    PetscFinalize();
#endif
#ifdef MPISRC
    MPI_Finalize();
#endif
    
    std::exit(0);
}

/* This routine forces everyone to die */
void sim::abort(int line,const char *file, std::ostream *log) {
    *log << "Exiting at line " << line << " of file " << file << std::endl;
    for (int b=0;b<blks.myblock;++b) {
        sim::blks.blk(b)->output("aborted_solution", block::display);
        sim::blks.blk(b)->output("aborted_solution", block::restart);
    }
#ifdef petsc
    PetscFinalize();
#endif
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    
    /* Terminates all threads */
    std::exit(1);
}
