#include "blocks.h"
#include "block.h"
#include "mgblock.h"
#include "r_mesh.h"

class btype {
   public:
      enum ids {plain=1};
};

block* blocks::getnewblock(int idnum, input_map blockdata) {
   std::string keyword,val;
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   char idntystring[10];
   int type;        
   block *temp;  
   

   sprintf(idntystring,"b%d",idnum);
   keyword = std::string(idntystring) + ".type";
   if (blockdata.get(keyword,val)) {
      data.str(val);
      data >> type;  
      data.clear(); 
   }
   else {
      if (!blockdata.get("blocktype",val)) {
         *sim::log << "couldn't find block type" << std::endl;
      }
   }
   
   switch(type) {
      case btype::plain: {
         temp = new mgrid<r_mesh>(idnum);
         break;
      }
      default: {
         std::cout << "unrecognizable block type: " <<  type << std::endl;
         temp = new mgrid<r_mesh>(idnum);
         break;
      }
   } 
      
   return(temp);
}
