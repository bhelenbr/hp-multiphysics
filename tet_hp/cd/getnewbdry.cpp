#include "tet_hp_cd.h"
#include "bdry_cd.h"
#include <input_map.h>

using namespace bdry_cd;

/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tet_hp_cd_stype {
    public:
        static const int ntypes = 3;
        enum ids {unknown=-1,plain,dirichlet,adiabatic};
        static const char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i)
                if (!strcmp(nin,names[i])) return(i);
            return(unknown);
        }
};

const char tet_hp_cd_stype::names[ntypes][40] = {"plain","dirichlet","adiabatic"};

hp_face_bdry* tet_hp_cd::getnewfaceobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_face_bdry *temp;  
    

    keyword =  fbdry(bnum)->idprefix + "_cd_type";
    if (bdrydata.get(keyword,val)) {
        type = tet_hp_cd_stype::getid(val.c_str());
        if (type == tet_hp_cd_stype::unknown)  {
            *gbl->log << "unknown side type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = tet_hp_cd_stype::unknown;
    }

    switch(type) {
        case tet_hp_cd_stype::plain: {
            temp = new hp_face_bdry(*this,*fbdry(bnum));
            break;
        }
        case tet_hp_cd_stype::dirichlet: {
            temp = new dirichlet(*this,*fbdry(bnum));
            break;
        }
        case tet_hp_cd_stype::adiabatic: {
            temp = new neumann(*this,*fbdry(bnum));
            break;
        }

        default: {
            temp = tet_hp::getnewfaceobject(bnum,bdrydata);
            break;
        }
    }
    
    return(temp);
}
