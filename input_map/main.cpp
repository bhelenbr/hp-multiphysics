#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include "input_map.h"
#include <muParser.h>
#include <unistd.h>

//#define TESTING

int main (int argc, char *argv[]) {
	input_map mymap;

#ifndef TESTING
    std::string deleteKeyword;
    std::string file,newFile;
    bool useNewFile = false, deleteKey = false;
    int opt;
    while ((opt = getopt (argc, argv, "d:o:h")) != -1) {
        switch (opt) {
            case 'd': {
                deleteKey = true;
                deleteKeyword = optarg;
                break;
            }
            case 'o': {
                useNewFile = true;
                newFile = optarg;
                break;
            }
            case 'h': {
                std::cout << "mod_map utility" << std::endl;
                std::cout << "mod_map [-odh] inputfile [keyword] [value]" << std::endl;
                std::cout << " -h prints usage information" << std::endl;
                std::cout << "\"mod_map -d keyword inputfile\" deletes keyword form inputfile" << std::endl;
                std::cout << "\"-o filename\" outputs to file instead of replacing current file" << std::endl;
                return 1;
            }
                
            case '?': {
                std::cerr << "Unknown option character " << optopt << std::endl;
                std::cerr << "Use mod_map -h for usage" << std::endl;
                return 1;
            }
            default: {
                std::cerr << "Use mod_map -h for usage" << std::endl;
                return 1;
            }
        }
    }

    int index = optind;
    if (argc -index < 1) {
        std::cerr << "Missing input filename?" << std::endl;
        return 1;
    }
    else {
        file = argv[index++];
        mymap.input(file);
    }
    
	std::ofstream fout;
	if (useNewFile) {
		fout.open(newFile);
	}
	else {
		fout.open(file);
	}
	
    if (deleteKey) {
		std::map<std::string,std::string>::iterator mi;
        mi = mymap.find(deleteKeyword);
		if (mi != mymap.end()) {
			mymap.erase(mi);
		}
	}
	else {
		// add or update
        if (argc -index <2) {
            std::cerr << "Missing keyword and values" << std::endl;
            return 1;
        }
        else {
            mymap[argv[index]] = argv[index+1];
        }
	}

	fout << mymap;
	fout.close();
	return 0;

#else
	mu::Parser P;
	try {
		P.SetExpr("23685*3");
		std::cout <<  P.Eval() << std::endl;
//		mu::varmap_type variables = P.GetUsedVar();
//		mu::varmap_type::const_iterator item = variables.begin();
//		for (; item!=variables.end(); ++item) {
//			std::cout << item->first << std::endl;
//		}
	}
	catch (mu::Parser::exception_type &e) {
		std::cout << "Message:  " << e.GetMsg() << std::endl;
		std::cout << "Formula:  " << e.GetExpr() << std::endl;
		std::cout << "Token:     " << e.GetToken() << std::endl;
		std::cout << "Position: " << e.GetPos() << std::endl;
		std::cout << "Errc:      " << e.GetCode() << std::endl;
		return 0;
	}
	catch(std::exception &e) {
		std::cout << e.what() << std::endl;
	}
	catch(...) {
		std::cout << "I am confused" << std::endl;
	}
	return 0;

	/* INPUT MAP FROM FILE */
	mymap.input(argv[1]);

    std::string key,value;
    std::istringstream data;
    std::string idprefix("prefix");
    
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

	/* LOAD INT USING FORMULA */
	int c;
	std::cout << mymap.get("formula",c) << ' ' << c << std::endl;

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
