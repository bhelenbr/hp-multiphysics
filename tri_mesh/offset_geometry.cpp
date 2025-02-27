//
//  offset_geometry.cpp
//  tri_mesh
//
//  Created by Brian Helenbrook on 12/30/23.
//

#include "tri_mesh.h"
#include "tri_boundary.h"

//#define NOCOMM

void tri_mesh::offset_geometry(input_map& input) {
    FLT d;
    if (!input.get(gbl->idprefix +"_offset",d)) {
        *gbl->log << "no offset distance for " << gbl->idprefix << std::endl;
    }
    
    std::vector<bool> offset_flag(nebd,false);
    std::vector<bool> spline_flag(nebd,false);
    int noffsets = 0;
    int nonoffset_sides = 0;
    int maxeid = 0;
    for (int i=0;i<nebd;++i) {
        int inlist;
        *gbl->log << ebdry(i)->idprefix +"_offset" <<std::endl;
        if (input.get(ebdry(i)->idprefix +"_offset",inlist)) {
            input.delete_entry(ebdry(i)->idprefix +"_offset");
            ++noffsets;
            offset_flag[i] = true;
            if (ebdry(i)->mytype=="plain_spline")
                spline_flag[i] = true;
        }
        else {
            /* Count non-offset boundary sides so we can keep that the same */
            nonoffset_sides += ebdry(i)->nseg;
        }
        maxeid = max(ebdry(i)->idnum,maxeid);
        /* Fill in vertex information so we can move through the sides */
        tri_gbl->intwk(seg(ebdry(i)->seg(0)).pnt(0)) = i;
        pnt(seg(ebdry(i)->seg(ebdry(i)->nseg-1)).pnt(1)).info = i;
    }
    
    int maxvid = 0;
    for (int i = 0; i < nvbd; ++i) {
        maxvid = max(vbdry(i)->idnum,maxvid);
    }
    
    int newblock = 1;
//    for (int nprocessed = 0; nprocessed < noffsets;)
        /* Find a first boundary */
        
    for (int i=0;i<nebd;++i) {
        if (!offset_flag[i]) continue;
        int prev = pnt(seg(ebdry(i)->seg(0)).pnt(0)).info;
        if (prev < 0) {
            /* This is the beginning of an offset boundary */
            int vnum = ebdry(i)->vbdry(0);
            if (vnum > -1) {
                /* Check to see if there is an explicit designation for this vertex */
                std::string vtype;
                if (input.get(vbdry(vnum)->idprefix +"_offset_type",vtype)) {
                    /* Do something */
                }
            }
        }
        else if (prev == i) {
            /* Loop boundary */
            if (!spline_flag[i]) {
                *gbl->log << "Weird, Can't have a loop boundary that is not a spline in an offset" << std::endl;
                sim::abort(__LINE__,__FILE__,gbl->log);
            }
            else {
                spline_bdry<edge_bdry>* spbdry = dynamic_cast<spline_bdry<edge_bdry> *>(ebdry(i));
                TinyVector<FLT,ND> t1,t2,n1,n2;
                spbdry->my_spline.tangent(spbdry->smin,t1);
                spbdry->my_spline.tangent(spbdry->smax,t2);
                /* Spline is defined clockwise???? */
                t1 = t1/sqrt(t1[0]*t1[0] +t1[1]*t1[1]);
                n1(0) = -t1[1];
                n1(1) = t1[0];
                t2 = t2/sqrt(t2[0]*t2[0] +t2[1]*t2[1]);
                n2(0) = -t2[1];
                n2(1) = t2[0];
                int vnum = ebdry(i)->vbdry(0);
                if (vnum > -1) {
                    /* Check to see if there is an explicit designation for this vertex */
                    std::string vtype;
                    if (input.get(vbdry(vnum)->idprefix +"_offset_type",vtype)) {
                        input.delete_entry(vbdry(vnum)->idprefix +"_offset_type");
                        int vpnt = vbdry(vnum)->pnt;
                        FLT res = lngth(vpnt);
                        if (vtype == "polar") {
                            /* Calculate angles */
                            FLT thetaf = atan2(n2[1],n2[0]);
                            FLT dtheta = acos(n1[0]*n2[0]+n1[1]*n2[1]);
                            FLT theta0 = thetaf-dtheta;
                            
                            /* make new input file */
                            input["nblock"] = "3";
                            
                            /* Output offset domain stuff */
                            /* This is airfoil boundary layer domain */
                            ostringstream nstr;
                            nstr << "b" << newblock;
                            input[nstr.str()+"_type"] = "spline_mapped_mesh";
                            input[nstr.str()+"_mesh"] = nstr.str() +".d";
                            input[nstr.str()+"_spline"] = input[ebdry(i)->idprefix +"_filename"];
                            if (input.find(ebdry(i)->idprefix +"_theta") != input.end()) {
                                input[nstr.str()+"_theta"] = input[ebdry(i)->idprefix +"_theta"];
                            }
                            if (input.find(ebdry(i)->idprefix +"_center") != input.end()) {
                                input[nstr.str()+"_center"] = input[ebdry(i)->idprefix +"_center"];
                            }
                            if (input.find(ebdry(i)->idprefix +"_scale") != input.end()) {
                                input[nstr.str()+"_scale"] = input[ebdry(i)->idprefix +"_scale"];
                            }

                            /* Output polar domain stuff */
                            ostringstream nstr1, pntstring;
                            nstr1 << "b" << newblock+1; /* nstr1 is the polar domain */
                            input[nstr1.str()+"_type"] = "polar_mapped_mesh";
                            input[nstr1.str()+"_mesh"] = nstr1.str() +".d";
                            pntstring.setf(std::ios::scientific, std::ios::floatfield);
                            pntstring.precision(10);
                            pntstring << pnts(vbdry(vnum)->pnt)(0) << " " << pnts(vbdry(vnum)->pnt)(1);
                            input[nstr1.str()+"_pnt"] = pntstring.str();
                            input[nstr1.str()+"_theta_length"] = input[gbl->idprefix +"_offset"];
                            
                            /* Output outer domain stuff */
                            input[gbl->idprefix+"_mesh"] = gbl->idprefix +".d";
 
                            /* Output side definitions */
                            /* end boundary of b.l. domain */
                            /* start boundary of polar domain*/
                            ostringstream sstr;
                            sstr << maxeid+1;
                            input[nstr.str()+"_s" +sstr.str() +"_type"] = "prdc";
                            input[nstr1.str()+"_s" +sstr.str() +"_type"] = "prdc";
                            
                            /* top of offset domain*/
                            sstr.str("");
                            sstr.clear();
                            sstr << maxeid+5;
                            /* Rename boundary conditions for spline surface */
                            input.rename_entries(ebdry(i)->idprefix+"_",gbl->idprefix +"_s" +sstr.str() +"_");
#ifdef NOCOMM
                            input[nstr.str()+"_s" +sstr.str() +"_type"] = "plain";
                            input[gbl->idprefix +"_s" +sstr.str() +"_type"] = "spline";
#else
                            input[nstr.str()+"_s" +sstr.str() +"_type"] = "mapped_comm";
                            input[gbl->idprefix +"_s" +sstr.str() +"_type"] = "spline_comm";
#endif
                            input[gbl->idprefix +"_s" +sstr.str() +"_norm_dist"] = "-" +input[gbl->idprefix +"_offset"];

                            
                            /* start boundary of offset domain */
                            /* end boundary of polar domain */
                            sstr.str("");
                            sstr.clear();
                            sstr << maxeid+3;
                            input[nstr.str()+"_s" +sstr.str() +"_type"] = "prdc";
                            input[nstr1.str()+"_s" +sstr.str() +"_type"] = "prdc";
                            
                            /* airfoil surface (straight boundary) */
                            sstr.str("");
                            sstr.clear();
                            sstr << maxeid+6;
                            input[nstr.str()+"_s" +sstr.str() +"_type"] = "plain";
                            
                            /* Top of polar domain */
                            sstr.str("");
                            sstr.clear();
                            sstr << maxeid+2;
#ifdef NOCOMM
                            input[nstr1.str()+"_s" +sstr.str() +"_type"] = "plain";
                            input[gbl->idprefix +"_s" +sstr.str() +"_type"] = "circle";
#else
                            input[nstr1.str()+"_s" +sstr.str() +"_type"] = "mapped_comm";
                            input[gbl->idprefix +"_s" +sstr.str() +"_type"] = "circle_comm";
#endif
                            input[gbl->idprefix +"_s" +sstr.str() +"_radius"] = input[gbl->idprefix +"_offset"];
                            input[gbl->idprefix +"_s" +sstr.str() +"_center"] = pntstring.str();
                            
                            /* Output communication points */
                            sstr.str("");
                            sstr.clear();
                            sstr << maxvid+1;
#ifdef NOCOMM
                            input[nstr1.str()+"_v" +sstr.str() +"_type"] = "plain";
                            input[gbl->idprefix +"_v" +sstr.str() +"_type"] = "plain";
                            input[nstr.str()+"_v" +sstr.str() +"_type"] = "plain";
#else
                            input[nstr1.str()+"_v" +sstr.str() +"_type"] = "mapped_comm";
                            input[gbl->idprefix +"_v" +sstr.str() +"_type"] = "comm";
                            input[nstr.str()+"_v" +sstr.str() +"_type"] = "mapped_comm";
#endif
                            
                            /* Output communication points */
                            sstr.str("");
                            sstr.clear();
                            sstr << maxvid+2;
#ifdef NOCOMM
                            input[nstr1.str()+"_v" +sstr.str() +"_type"] = "plain";
                            input[gbl->idprefix +"_v" +sstr.str() +"_type"] = "plain";
                            input[nstr.str()+"_v" +sstr.str() +"_type"] = "plain";
#else
                            input[nstr1.str()+"_v" +sstr.str() +"_type"] = "mapped_comm";
                            input[gbl->idprefix +"_v" +sstr.str() +"_type"] = "comm";
                            input[nstr.str()+"_v" +sstr.str() +"_type"] = "mapped_comm";
#endif
                            
                            /* Output outer domain */
                            ofstream out;
                            out.open(gbl->idprefix + ".d");
                            out << nebd-1 +3 << std::endl;
                            TinyVector<FLT,ND> loc,tan,curv;
                            spline_functions2D::interpolate(loc, tan, curv, spbdry->my_spline, spbdry->my_spline.start(),spbdry->scale, spbdry->theta,spbdry->pos, -d);
                            out << "0: " << loc[0] << ' ' << loc[1] << ' ' << res << ' ' << maxvid+2 << std::endl;
                            
                            FLT smid = 0.5*(spbdry->my_spline.start()+spbdry->my_spline.stop());
                            spline_functions2D::interpolate(loc, tan, curv, spbdry->my_spline, smid,spbdry->scale, spbdry->theta,spbdry->pos, -d);
                            out << "1: " << loc[0] << ' ' << loc[1] << ' ' << res << ' ' << 0 << std::endl;

                            
                            spline_functions2D::interpolate(loc, tan, curv, spbdry->my_spline, spbdry->my_spline.stop(),spbdry->scale, spbdry->theta,spbdry->pos, -d);
                            out << "2: " << loc[0] << ' ' << loc[1] << ' ' << res << ' ' << maxvid+1 << std::endl;
                            
                            /* Assume outer boundary is also a loop */
                            int loop_start = (i+1) % nebd;
                            int loop_ind = loop_start;
                            std::vector<int> idnums;
                            int cnt = 3;
                            do {
                                int vnum = seg(ebdry(loop_ind)->seg(0)).pnt(0);
                                loc = pnts(vnum);
                                out << cnt << ": " << loc[0] << ' ' << loc[1] << ' ' << lngth(vnum) << ' ' << 0 << std::endl;
                                idnums.push_back(ebdry(loop_ind)->idnum);
                                vnum = seg(ebdry(loop_ind)->seg(ebdry(loop_ind)->nseg-1)).pnt(1);
                                loop_ind = tri_gbl->intwk(vnum);
                                ++cnt;
                                
                            } while (loop_ind != loop_start);
                            
                            out << nebd+2 << std::endl;
                            out << "0: 0 1 " << maxeid+5 << std::endl;
                            out << "1: 1 2 " << maxeid+5 << std::endl;
                            out << "2: 2 0 " << maxeid+2 << std::endl;
                            for (cnt = 3; cnt < nebd-2+3; ++cnt) {
                                out << cnt << ": " << cnt << ' ' << cnt +1 << ' ' << idnums[cnt-3] << std::endl;
                            }
                            out << cnt << ": " << cnt << ' ' << 3 << ' ' << idnums[cnt-3] << std::endl;
                            out.close();
                            
                            /* Output spline domain */
                            FLT dxds = sqrt(tan(0)*tan(0) +tan(1)*tan(1));
                            FLT res1 = res/dxds;
                            out.open(nstr.str() + ".d");
                            out << 5 << std::endl;
                            out << "0: " << spbdry->my_spline.start() << ' ' << 0.0  << ' ' << res1 << ' ' << 0 << std::endl;
                            out << "1: " << spbdry->my_spline.stop() << ' ' << 0.0 << ' ' << ' ' << res1 << ' ' << 0 << std::endl;
                            out << "2: " << spbdry->my_spline.stop() << ' ' << d << ' ' << res1 << ' ' << maxvid+1 << std::endl;
                            out << "3: " << 0.5*(spbdry->my_spline.start()+spbdry->my_spline.stop()) << ' ' << d << ' ' << res1 << ' ' << 0 << std::endl;
                            out << "4: " << spbdry->my_spline.start() << ' ' << d << ' ' << res1 << ' ' << maxvid+2 << std::endl;
                            out << 5 << std::endl;
                            out << "0: 0 1 " << maxeid+6 << std::endl;
                            out << "1: 1 2 " << maxeid+3 << std::endl;
                            out << "2: 2 3 " << maxeid+5 << std::endl;
                            out << "3: 3 4 " << maxeid+5 << std::endl;
                            out << "4: 4 0 " << maxeid+1 << std::endl;
                            out.close();
                            
                            /* Output polar domain */
                            out.open(nstr1.str() +".d");
                            out << 4 << std::endl;
                            out << "0: " << -thetaf*d << ' ' << 0.0 << ' ' << ' ' << res << ' ' << 0 << std::endl;
                            out << "1: " << -theta0*d << ' ' << 0.0  << ' ' << res << ' ' << 0 << std::endl;
                            out << "2: " << -theta0*d << ' ' << d << ' ' << res << ' ' << maxvid+2 << std::endl;
                            out << "3: " << -thetaf*d << ' ' << d << ' ' << res << ' ' << maxvid+1 << std::endl;
                            out << 4 << std::endl;
                            out << "0: 0 1 " << maxeid+4 << std::endl;
                            out << "1: 1 2 " << maxeid+1 << std::endl;
                            out << "2: 2 3 " << maxeid+2 << std::endl;
                            out << "3: 3 0 " << maxeid+3 << std::endl;
                            out.close();
                            

                            

                        }
                    }
                }
            }
        }
    }
    input.delete_entry("offset");
    input.delete_entry(gbl->idprefix +"_offset");
    input["ntstep"] = "2";
    
    ofstream out;
    out.open("offset.inpt");
    out << input;
    out.close();
 
}


//                    int tangent(const double spt, blitz::TinyVector<double,ND>& tan) const;
//                    int curvature(const double spt, blitz::TinyVector<double,ND>& curv) const;
//                    int find(double &spt, blitz::TinyVector<double,ND>& loc) const;
//                    double start() const {return(x(0));}
//                    double stop() const {return(x(npts-1));}
