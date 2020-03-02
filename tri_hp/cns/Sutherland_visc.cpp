//
//  Sutherland_visc.cpp
//  tri_hp
//
//  Created by Luke D'Aquila on 1/7/20.
//

#include "tri_hp_cns.h"

#idef Sutherland
void tri_hp_cns::Sutherland_visc(FLT RT){
    
    gbl->mu = gbl->s1*pow(RT,1.5)/(RT+gbl->s2);
    
    return;
    
}
#endif
