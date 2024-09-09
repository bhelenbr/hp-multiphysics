/*
 *  input_map.h
 *  input map
 *
 *  Created by Brian Helenbrook on Fri Sep 06 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _input_map_h_
#define _input_map_h_
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <vector>

class input_map : public std::map<std::string,std::string> {
public:
    bool echo;
    std::string echoprefix;
    std::ostream *log;
    input_map() : echo(false), echoprefix(""), log(&std::cout) {
#ifdef BZ_DEBUG
    std::cerr << "#input_map: BZ_DEBUG is set\n";
#endif
#ifdef DEBUG
    std::cerr << "#input_map: Running in Xcode's DEBUG Mode\n";
#endif
    }
    bool input(const std::string &filename);
    friend std::istream& operator >>(std::istream &is,input_map &obj);
    friend std::ostream& operator <<(std::ostream &os,const input_map &obj);
    
    /* Generic get template routine */
    template<class T> bool get(const std::string &keyword, T &vout) {
        std::istringstream data;
        std::map<std::string,std::string>::const_iterator mi;
        mi = find(keyword);
        if (mi != end()) {
            data.str((*this)[keyword]);
            if (data >> vout) {
                if (echo) *log << echoprefix << keyword << ": " << vout << std::endl;
                return true;
            }
        }
        return(false);
    }
    
    /* Special case for doubles so can load from formula */
    bool get(const std::string &keyword, double &vout);
    
    /* Special case for ints so can load from formula for double and then convert */
    bool get(const std::string &keyword, int &vout);
    
    /* Special case for list of values on one line */
    template<class T> bool get(const std::string &keyword, T *array, int nentry) {
        std::istringstream data;
        std::map<std::string,std::string>::const_iterator mi;
        mi = find(keyword);
        if (mi != end()) {
            data.str((*this)[keyword]);
            for (int n=0;n<nentry;++n) {
                if (!(data >> array[n])) return(false);
            }
            if (echo) {
                *log << echoprefix << keyword << ": ";
                for (int n=0;n<nentry;++n)
                    *log << array[n] << ' ';
                *log << std::endl;
            }
            return true;
        }
        return(false);
    }
    
    /* Special case for doubles so can load from list of formulas */
    bool get(const std::string &keyword, double *array, int nentry);
    
    /* Get w/default function */
    template<class T> void getwdefault(const std::string &keyword, T &vout,const T dflt) {
        if (!get(keyword,vout)) {
            vout = dflt;
            if (echo) *log << echoprefix << keyword << ": " << dflt << std::endl;
        }
        return;
    }
    
    /* Get w/default function for a list on a single line */
    template<class T> void getwdefault(const std::string &keyword, T *array, int nentry, const T *dflt) {
        if (!get(keyword,array,nentry)) {
            for (int n=0;n<nentry;++n) array[n] = dflt[n];
            if (echo) {
                *log << echoprefix << keyword << ": ";
                for (int n=0;n<nentry;++n)
                    *log << array[n] << ' ';
                *log << std::endl;
            }
        }
        return;
    }
    
    /* Function to get entire line as a string */
    bool getline(const std::string &keyword, std::string& vout);
    
    /* Get w/default for getting an entire line */
    void getlinewdefault(const std::string &keyword, std::string& vout, const std::string& dflt);
    
    int keys_with_ending(std::string const &ending,std::vector<std::string>& keywords);
    int delete_entry(std::string key);
    int delete_entries(std::string exp);
    int rename_entries(std::string sfind, std::string sreplace);
};

#endif
