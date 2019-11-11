/*
 *  main.cpp
 *  spline++
 *
 *  Created by Brian Helenbrook on 4/17/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include "spline.h"
#include <unistd.h>
#include <fstream>

// #define THIRD
const int ND = 2;
using namespace spline_functions2D;
using namespace blitz;


int main(int argc, char** argv) {
    TinyVector<double,ND> dx = 0.0, xfind, x, tan, curv, zero = 0.0;
    double angle = 0.0, size = 1.0, norm_dist = 0.0;
    int npts = 0;
    std::string coord_file;
    std::ifstream point_data;

    bool Find = false;
    bool UseFile = false;
    int opt;
    while ((opt = getopt (argc, argv, "m:s:r:o:i:hf")) != -1) {
        switch (opt) {
            case 'm': {
                istringstream input(optarg);
                for (int n=0;n<ND;++n) {
                    if (!(input >> dx(n))) {
                        std::cerr << "Offset should be a comma separated list of floats after -m argument (no spaces)\n";
                        return 1;
                    }
                    input.ignore(1); // skip comma
                }
                break;
            }
            case 's': {
                istringstream input(optarg);
                if (!(input >> size) || size <= 0) {
                    std::cerr << "Scale should be a positive float after -s argument\n";
                    return 1;
                }
                break;
            }
                
            case 'r': {
                istringstream input(optarg);
                if (!(input >> angle)) {
                    std::cerr << "Angle should be a float in degrees after -r argument\n";
                    return 1;
                }
                angle *= M_PI/180.0;
                break;
            }
            
            case 'o': {
                istringstream input(optarg);
                if (!(input >> norm_dist)) {
                    std::cerr << "Offset should be a float after -o argument\n";
                    return 1;
                }
                break;
            }
            case 'f': {
                Find = true;
                break;
            }
                
            case 'i': {
                UseFile = true;
                point_data.open(optarg);
                if (!point_data) {
                    std::cerr << "error: couldn't open file: " << optarg << std::endl;
                    return 1;
                }
                break;
            }

            case 'h': {
                std::cout << "Usage for finding points close to a spline: Spline -f -m <dx,dy> -s <scale> -r <angle> [-i <filename with list points to find>] <filename.spl> <x,y>" << std::endl;;
                std::cout << "Usage for interpolating on spline: Spline -m <dx,dy> -s <scale> -r <angle> -o <offset> [-i <filename with list of parametric coordinates] <filename.spl> <npoints> <curvature sensitivity>" << std::endl;
                std::cout << "No arguments are read beyond the spline filename if using -i option" << std::endl;
                std::cout << "curvature sensitivity defaults to 0 for uniform interpolation. Higher values causes insertion of points into highly curved segments" << std::endl;
                return 1;
            }
                
            case '?': {
                std::cerr << "Unknown option character " << optopt << std::endl;
                std::cerr << "Use spline -h for usage" << std::endl;
                return 1;
            }
            default: {
                std::cerr << "Use spline -h for usage" << std::endl;
                return 1;
            }
        }
    }
    
#ifdef THIRD
    spline3<ND> myspline;  // (3rd order spline)
#else
    spline<ND> myspline;  // (5th order spline)
#endif
    
    int index = optind;
    if (argc -index < 1) {
        std::cerr << "Missing spline filename?" << std::endl;
        return 1;
    }
    else {
        std::string file(argv[index++]);
        myspline.read(file);
    }
    
    if (Find && !UseFile) {
        if (argc-index == 1) {
            /* Read x,y coordinates of point to find */
            istringstream input(argv[index++]);
            for (int n=0;n<ND;++n) {
                if (!(input >> xfind(n))) {
                    std::cerr << "Point to find should be a comma separated list of floats after spline filename (no spaces)\n";
                    return 1;
                }
                input.ignore(1); // skip comma
            }
        }
        else {
            std::cerr << "Point to find should be a comma separated list of floats after spline filename (no spaces)\n";
            return 1;
        }
    }
    
    
    double alpha = 0.0; /* Curvature sensitivity */
    if (!Find && !UseFile) {
        if (argc-index >= 1) {
            /* Read number of points */
            istringstream input(argv[index++]);
            if (!(input >> npts)) {
                std::cerr << "Need to supply integer number of points to be interpolated\n";
                return 1;
            }

            if (argc-index == 1) {
                istringstream input(argv[index++]);
                if (!(input >> alpha)) {
                    std::cerr << "Error trying to read curvature sensitivity which should be a floating point number\n";
                    return 1;
                }
            }
        }
        else {
            std::cerr << "Need to supply integer number of points to be interpolated\n";
            return 1;
        }
    }
    
    double s;
    if (Find) {
        if (!UseFile) {
            transform2D(xfind,size,angle,dx);
            myspline.find(s,xfind);
            interpolate(myspline, s, size, angle, dx);
        }
        else {
            while (!point_data.eof()) {
                for (int n=0;n<ND;++n) {
                    point_data >> xfind(n);
                }
                myspline.find(s,xfind);
                interpolate(myspline, s, size, angle, dx);
            }
        }
    }
    else {
        if (!UseFile) {
            
            /* Uniform Spacing */
//            for (double s=myspline.start();s<=myspline.stop();s+=(myspline.stop()-myspline.start())/npts) {
//                interpolate(myspline, s, size, angle, dx, norm_dist);
//            }
            
            /* Spacing based on curvature and tangent */
            int maxsize = 2*npts;
            Array<double,1> svalues(maxsize);
            double ds = (myspline.stop()-myspline.start())/(npts-1);
            /* First guess is uniform spacing in parametric coordinate */
            svalues(0) = myspline.start();
            for (int i=0;i<npts;++i) {
                svalues(i) = myspline.start() +ds*i;
            }
            
            /* This gets uniform spacing in physical coordinate */
            double lambda = ds, lambdaold = 0.0;
            TinyVector<double,ND> x, tan, tan1, curv;
            int iter = 0;
            while(abs(lambda-lambdaold) > 0.01*ds && iter < 10) {
                for (int i=1;i<npts;++i) {
                    myspline.tangent(svalues(i), tan);
                    myspline.curvature(svalues(i), curv);
                    double t2 = dot(tan,tan);
                    curv = curv/t2 -dot(curv,tan)*tan/(t2*t2);
                    ds = lambda/sqrt(t2);
                    svalues(i) = svalues(i-1)+ds;
                }
                svalues(Range(2,npts-2)) = (svalues(Range(1,npts-3)) +svalues(Range(3,npts-1)))/2.0;
                double rescale = (myspline.stop()-myspline.start())/(svalues(npts-1)-svalues(0));
                lambdaold = lambda;
                lambda *= rescale;
                svalues = (svalues-myspline.start())*rescale +myspline.start();
                ++iter;
            }
            
            /* 2 level refinement of curved elements */
            for(iter = 0;iter<2;++iter) {
                myspline.tangent(svalues(0), tan1);
                for (int i=1;i<npts;++i) {
                    myspline.tangent(svalues(i), tan);
                    tan1 -= tan;
                    if (alpha*alpha*dot(tan1,tan1) > dot(tan,tan)) {
                        if (npts < maxsize) {
                            for(int j = npts-1;j>=i;--j) {
                                svalues(j+1) = svalues(j);
                            }
                            svalues(i) = 0.5*(svalues(i-1)+svalues(i));
                            ++npts;
                        }
                        else {
                            std::cerr << "allocated size wasn't large enough\n";
                            break;
                        }
                    }
                    tan1 = tan;
                }
            }
            
            for (int i=0;i<npts;++i) {
                interpolate(myspline, svalues(i), size, angle, dx, norm_dist);
            }
        }
        else {
            while (!point_data.eof()) {
                point_data >> s;
                interpolate(myspline, s, size, angle, dx, norm_dist);
            }
        }
        
        /* Normalized by curvature */
    }

	return(0);
}

