/*
 *  input_map.cpp
 *  input map
 *
 *  Created by Brian Helenbrook on Fri Sep 06 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "input_map.h"
#include <fstream>
#include <muParser.h>
#include <tgmath.h>
#include <regex>

using namespace std;

string trim(string str) {
    char const* delims = " \t\r\n";
    string::size_type pos, pos2;
    
    // trim leading whitespace
    string::size_type  notwhite = str.find_first_not_of(delims);
    str.erase(0,notwhite);
    
    // trim trailing whitespace
    notwhite = str.find_last_not_of(delims);
    str.erase(notwhite+1);
    
    // compress middle white space
    pos = 0;
    while( (pos = str.find_first_of(delims,pos)) < str.length()) {
        pos2 = str.find_first_not_of(delims,pos);
        str.replace(pos,pos2-pos," ");
        ++pos;
    }
    
    return(str);
}

bool input_map::input(const std::string &filename) {
    ifstream infile(filename.c_str());
    if (!infile) {
        return(false);
    }
    infile >> *this;
    infile.close();
    return(true);
}

istream& operator >>(istream &infile,input_map &obj) {
    string key,value;
    char firstchar;
    
    while ( (firstchar = infile.peek()) != EOF ) {
        if (firstchar == '#') {
            // SKIP COMMENT
            getline(infile,key,'\n');
            continue;
        }
        getline(infile,key,':');
        key = trim(key);
        getline(infile,value,'\n');
        value = trim(value);
        obj[key] = value;
    }
    
    return(infile);
}

ostream& operator <<(ostream &outfile,const input_map &mapout)  {
    map<string,string>::const_iterator mi;
    
    for (mi = mapout.begin(); mi != mapout.end(); ++mi) {
        outfile << mi->first << ":" << mi->second << endl;
    }
    return(outfile);
}

bool input_map::get(const std::string &keyword, double &vout) {
    
    std::istringstream data;
    std::map<std::string,std::string>::const_iterator mi;
    double value;
    
    mi = find(keyword);
    if (mi != end()) {
        data.str((*this)[keyword]);
        if (data >> vout && data.eof()) {
            if (echo) *log << echoprefix << keyword << ": " << vout << std::endl;
            return(true);
        }
        else {
            double (*my_erf)(double) = erf;
            double (*my_erfc)(double) = erfc;
            
            if (echo) *log << echoprefix << keyword << ": " << (*this)[keyword] << std::endl;
            /* TRY TO PARSE MATHEMATICAL EXPRESSION */
            mu::Parser P,P1;
            P.DefineFun("erf", my_erf, false);
            P.DefineFun("erfc", my_erfc, false);
            P1.DefineFun("erf", my_erf, false);
            P1.DefineFun("erfc", my_erfc, false);
            
            try {
                P.SetExpr((*this)[keyword]);
                mu::varmap_type variables = P.GetUsedVar();
                mu::varmap_type::const_iterator item = variables.begin();
                std::string echo_storage = echoprefix;
                echoprefix = echoprefix +'\t';
                for (; item!=variables.end(); ++item) {
                    if (!get(item->first,value)) {
                        // *log << echoprefix << "Trouble reading " << item->first << std::endl;
                        echoprefix = echo_storage;
                        return(false);
                    }
                    P1.DefineConst(item->first,value);
                }
                echoprefix = echo_storage;
                P1.SetExpr((*this)[keyword]);
                vout = P1.Eval();
                if (echo) *log << echoprefix << keyword << " = " << vout << std::endl;
                return(true);
            }
            catch (mu::Parser::exception_type &e) {
                std::cout << "Message:  " << e.GetMsg() << std::endl;
                std::cout << "Formula:  " << e.GetExpr() << std::endl;
                std::cout << "Token:     " << e.GetToken() << std::endl;
                std::cout << "Position: " << e.GetPos() << std::endl;
                std::cout << "Errc:      " << e.GetCode() << std::endl;
                return(false);
            }
            catch (...) { cout << "default exception"; }
        }
    }
    return(false);
}

bool input_map::get(const std::string &keyword, int &vout) {
    
    /* try to get a normal int */
    std::istringstream data;
    std::map<std::string,std::string>::const_iterator mi;
    
    mi = find(keyword);
    if (mi != end()) {
        data.str((*this)[keyword]);
        if (data >> vout && data.eof()) {
            if (echo) *log << echoprefix << keyword << ": " << vout << std::endl;
            return(true);
        }
        
        /* try to get a double and convert */
        double vdouble;
        if (get(keyword,vdouble)) {
            vout = vdouble;
            return(true);
        }
    }
    
    return(false);
}

bool input_map::get(const std::string &keyword, double *array, int nentry) {
    double value;
    std::string expression;
    std::istringstream data,data1;
    std::map<std::string,std::string>::const_iterator mi;
    
    mi = find(keyword);
    if (mi != end()) {
        data.str((*this)[keyword]);
        for (int n=0;n<nentry;++n) {
            data >> expression;
            data1.str(expression);
            if (!(data1 >> array[n])) {
                /* TRY TO PARSE MATHEMATICAL EXPRESSION */
                mu::Parser P,P1;
                try {
                    P.SetExpr(expression);
                }
                catch (mu::Parser::exception_type &e) {
                    std::cout << "Message:  " << e.GetMsg() << std::endl;
                    std::cout << "Formula:  " << e.GetExpr() << std::endl;
                    std::cout << "Token:     " << e.GetToken() << std::endl;
                    std::cout << "Position: " << e.GetPos() << std::endl;
                    std::cout << "Errc:      " << e.GetCode() << std::endl;
                    return(false);
                }
                mu::varmap_type variables = P.GetUsedVar();
                mu::varmap_type::const_iterator item = variables.begin();
                std::string echo_storage = echoprefix;
                echoprefix = echoprefix +'\t';
                for (; item!=variables.end(); ++item) {
                    if (!get(item->first,value)) {
                        echoprefix = echo_storage;
                        return(false);
                    }
                    P1.DefineConst(item->first,value);
                }
                echoprefix = echo_storage;
                try {
                    P1.SetExpr(expression);
                }
                catch (mu::Parser::exception_type &e) {
                    std::cout << "Message:  " << e.GetMsg() << std::endl;
                    std::cout << "Formula:  " << e.GetExpr() << std::endl;
                    std::cout << "Token:     " << e.GetToken() << std::endl;
                    std::cout << "Position: " << e.GetPos() << std::endl;
                    std::cout << "Errc:      " << e.GetCode() << std::endl;
                    return(false);
                }
                array[n] = P1.Eval();
            }
            data1.clear();
        }
        if (echo) {
            *log << echoprefix << keyword << ": ";
            for (int n=0;n<nentry;++n)
                *log << array[n] << ' ';
            *log << std::endl;
        }
        return(true);
    }
    return(false);
}

/* Function to get entire line as a string */
bool input_map::getline(const std::string &keyword, std::string& vout) {
    std::map<std::string,std::string>::const_iterator mi;
    
    mi = find(keyword);
    if (mi != end()) {
        vout = (*this)[keyword];
        if (echo) *log << echoprefix << keyword << ": " << vout << std::endl;
        return(true);
    }
    return(false);
}

/* Get w/default for getting an entire line */
void input_map::getlinewdefault(const std::string &keyword, std::string& vout, const std::string& dflt) {
    if (!getline(keyword,vout)) {
        vout = dflt;
        if (echo) *log << echoprefix << keyword << ": " << dflt << std::endl;
    }
    return;
}

int input_map::keys_with_ending(std::string const &ending,std::vector<std::string>& keywords) {
    int found = 0;
    for (auto mi = this->begin();mi != this->end();++mi) {
        if (mi->first.length() >= ending.length()) {
            if (!(mi->first.compare(mi->first.length() - ending.length(), ending.length(), ending))) {
                ++found;
                keywords.push_back(mi->first);
            }
        }
    }
    return(found);
}

int input_map::delete_entry(std::string exp) {
    auto mi = this->find(exp);
    if (mi != this->end()) {
        this->erase(mi);
        return(1);
    }
    return(0);
}

int input_map::delete_entries(std::string exp) {
    int found = 0;
    std::vector<string> to_delete;
    for (auto mi = this->begin();mi != this->end();++mi) {
        if (mi->first.find(exp) != std::string::npos) {
            to_delete.push_back(mi->first);
            ++found;
        }
    }
    
    for (auto key : to_delete) {
         this->erase(key);
     }
    return(found);
}


int input_map::rename_entries(std::string sfind, std::string sreplace) {
    int found = 0;
    std::vector<std::pair<std::string,std::string> > rename_list;
    for (auto mi = this->begin();mi != this->end();++mi) {
        auto pos = mi->first.find(sfind);
        if (pos != std::string::npos) {
            std::string newkey = mi->first;
            newkey.replace(pos, sfind.length(), sreplace);
            rename_list.push_back(std::pair<std::string,std::string>(mi->first,newkey));
            ++found;
        }
    }
    
    for (auto mypair : rename_list) {
        (*this)[mypair.second] = (*this)[mypair.first];
        this->erase(mypair.first);
     }
    
    return(found);
}





