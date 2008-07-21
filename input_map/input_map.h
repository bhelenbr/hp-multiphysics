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

class input_map : public std::map<std::string,std::string> {
    public:
        bool echo; 
        std::string echoprefix;
        std::ostream *log;
        input_map() : echo(false), echoprefix(""), log(&std::cout) {}
        void input(const std::string &filename);
        friend std::istream& operator >>(std::istream &is,input_map &obj);
        friend std::ostream& operator <<(std::ostream &os,const input_map &obj);
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
        
        /* SPECIAL CASE FOR DOUBLES SO CAN LOAD FROM FORMULA */
        bool get(const std::string &keyword, double &vout);
         
        template<class T> void getwdefault(const std::string &keyword, T &vout,const T dflt) {
            if (!get(keyword,vout)) {
                vout = dflt;
                if (echo) *log << echoprefix << keyword << ": " << dflt << std::endl;
            }
            return;
        }
        
        /* LIST OF VALUES ON ONE LINE */
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
        
        /* SPECIAL CASE FOR DOUBLES SO CAN LOAD FROM LIST OF FORMULAS */
        bool get(const std::string &keyword, double *array, int nentry);
        
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
        
        bool getline(const std::string &keyword, std::string& vout) {
            std::map<std::string,std::string>::const_iterator mi;
            
            mi = find(keyword);
            if (mi != end()) {
                vout = (*this)[keyword];
                if (echo) *log << echoprefix << keyword << ": " << vout << std::endl;
                return(true);
            }
            return(false);
        }
        
        void getlinewdefault(const std::string &keyword, std::string& vout, const std::string& dflt) {            
            if (!getline(keyword,vout)) {
                vout = dflt;
                if (echo) *log << echoprefix << keyword << ": " << dflt << std::endl;
            }
            return;
        }
     
};

#endif
