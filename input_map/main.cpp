#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include "input_map.h"

#define NO_TESTING

#ifndef TESTING
#include <parseargs.h>

static GBool Delete = gFalse;
static char NewFile[100] = "";
GBool printHelp = gFalse;


static ArgDesc argDesc[] = {
  {"-h",        argFlag,      &printHelp,      0,
    "print usage information"},
  {"-o", argString,    NewFile,     100,
    "output file name"},
  {NULL}
};
#endif

int main (int argc, char *argv[]) {
    input_map mymap;
    std::string key,value;
    std::istringstream data;
    std::string idprefix("prefix");
    
	
	
#ifndef TESTING
    GBool ok;

    // parse args
    ok = parseArgs(argDesc, &argc, argv);
    if (!ok || printHelp  || argc != 4) {
        fprintf(stderr, "mod_map utility ");
        printUsage("mod_map", "<inputfile> <keyword> <value>]", argDesc);
        exit(1);
    }
    
    mymap.input(argv[1]);
    mymap[argv[2]] = argv[3];
    
    std::ofstream fout;
    if (NewFile[0] != '\0') {
        fout.open(NewFile);
    }
    else {
        fout.open(argv[1]);
    }
    fout << mymap;
    fout.close();
    return 0;

#else
    /* INPUT MAP FROM FILE */
    mymap.input(argv[1]);
    

    /* SET MAP VALUE FROM FLOAT */
    std::ostringstream myout;
    myout.str("");
    myout << 3.72569 << std::flush;
    mymap["a"]  = myout.str();
    
    mymap["formula"] = "a + 3.0";
    
    std::cout << std::endl << std::endl << "OUTPUTING MAP AT BEGINNING" << std::endl;
    std::cout << mymap;
    mymap.echo = true;
    mymap.echoprefix = "#";
    std::cout << std::endl << std::endl;

        
    /* LOAD DOUBLE USING FORMULA */
    double b;
    mymap.get("formula",b);    
    
    /* ERASE ENTRY */
    std::map<std::string,std::string>::iterator mi;
    mi = mymap.find("formula");
    mymap.erase(mi);
        
    /* LOAD INTEGER */
    int itercrsn;
    mymap.getwdefault("itercrsn",itercrsn,1);
        
    /* LOAD INTEGER WITH DEFAULT VALUE */
    int log2p;
    if (!mymap.get(idprefix + ".log2p",log2p)) mymap.getwdefault("log2p",log2p,0);
    
    /* LOAD STRING */
    std::string filename;
    if (!mymap.get("mesh",filename)) {std::cout << "no mesh name\n"; exit(1);}
    
    /* LOAD ARRAY OF FLOATS */
    double dx[3], dflt[3] = {0.0, 1.0, 4.0};
    mymap.getwdefault("dx",dx,3,dflt);
    
    /* USE ONE KEYWORD TO FIND ANOTHER */
    std::string val;
    if (mymap.get(idprefix+".matching_block",val)) {
        if (!mymap.get(val,log2p)) {
            std::cout << "couldn't find matching blocks log2p" << std::endl;
            exit(1);
        }
    }
      
    /* OUTPUT MAP */
    std::cout << std::endl << std::endl << "OUTPUTING MAP AT END" << std::endl;
    std::cout << mymap;
#endif

}
