#include "blocks.h"
#include "block.h"
#include "mgblock.h"
#include "r_mesh.h"

block* blocks::getnewblock(int idnum, std::map<std::string,std::string> *blockdata) {
   std::string keyword;
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   char idntystring[10];
   int type;        
   block *temp;  
   
   type = idnum&0xffff;

   if (blockdata) {
      sprintf(idntystring,"b%d",idnum);
      keyword = std::string(idntystring) + ".type";
      mi = (*blockdata).find(keyword);
      if (mi != (*blockdata).end()) {
         data.str((*blockdata)[keyword]);
         data >> type;  
         data.clear(); 
      }
      else {
         keyword = "blocktype";
         data.str((*blockdata)[keyword]);
         if (!(data >> type)) {
            *log << "couldn't find block type" << std::endl;
         }
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
