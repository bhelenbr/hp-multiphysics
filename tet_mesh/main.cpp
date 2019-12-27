#include <blocks.h>
#include <string>
#include "tet_mesh.h"
#include <input_map.h>
#include <unistd.h>

#ifdef MPISRC
#include <mpi.h>
#endif

using namespace std;

int main(int argc, char *argv[]) {
	// tet_mesh zx;
	// zx.test();
	
#ifdef MPISRC
	int myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
#ifdef PTH
	// For debugging put interrupt here
	// On interrupt type this into gdb console: "handle SIGUSR1 nostop print pass"
	// for lldb: pro hand -p true -s false SIGUSR1
	// Then continue
	int rc = pth_init();
	if (!rc) {
		std::cerr << "couldn't start pth environment\n";
	}
#endif
	
  // parse args
    bool Append = false;
    bool Debugger = false;
    bool Coarsenby2 = false;
    bool Echo = false;
    bool Smooth = false;
    bool GMSHPartition = false;
    bool GMSHLabel = false;
    bool Coarsen_hp = false;
    bool Shift = false;
    bool Partition = false;
    bool Refineby2 = false;
    bool Scale = false;
    bool Vlngth = false;
    bool Format = false;
    bool Symmetrize = false;
    bool Cut = false;
    int informat = 4, outformat = 4;
    int p;
    int opt;
    while ((opt = getopt (argc, argv, "abcdefghi:lmo:p:rsvxyz")) != -1) {
        switch (opt) {
            case 'a': {
                Append = true;
                break;
            }
            case 'c': {
                Coarsenby2 = true;
                break;
            }
            case 'd': {
                Debugger = true;
                break;
            }
            case 'e': {
                Echo = true;
                break;
            }
            case 'f': {
                Smooth = true;
                break;
            }
            case 'g': {
                GMSHPartition = true;
                break;
            }
            case 'h': {
                std::cout << "tet_mesh utility" << std::endl;
                std::cout << "tet_mesh [-abcdefghlmpsvxyz] [-i input format] [-o output format] inputfile outputfile" << std::endl;
                std::cout << " -h prints usage information" << std::endl;
                std::cout << "-a append two meshes" << std::endl;
                std::cout << "-b create label of gmsh physical volumes" << std::endl;
                std::cout << "-c Coarsen mesh by factor of 2" << std::endl;
                std::cout << "-d Stop for Debugger" << std::endl;
                std::cout << "-e display basic mesh information" << std::endl;
                std::cout << "-f smooth mesh" << std::endl;
                std::cout << "-g partition gmsh physical volumes" << std::endl;
                std::cout << "-l coarsen substructured mesh" << std::endl;
                std::cout << "-m shift mesh position" << std::endl;
                std::cout << "-p # partition mesh into # parts" << std::endl;
                std::cout << "-r refine mesh by 2" << std::endl;
                std::cout << "-s scale mesh" << std::endl;
                std::cout << "-v Create a mesh resolution file" << std::endl;
                std::cout << "-x change format" << std::endl;
                std::cout << "-y Make mesh symmetric about y = 0" << std::endl;
                std::cout << "-z Cut mesh using indicator function" << std::endl;
                return 1;
            }
            case 'i': {
                std::istringstream input(optarg);
                if (!(input >> informat)) {
                    std::cerr << "Use tet_mesh -h for usage" << std::endl;
                    return 1;
                }
            }
            case 'l': {
                Coarsen_hp = true;
                break;
            }
            case 'm': {
                Shift = true;
                break;
            }
            case 'p': {
                Partition = true;
                std::istringstream input(optarg);
                if (!(input >> p)) {
                    std::cerr << "Unable to read number of partitions.  Use tet_mesh -h for usage" << std::endl;
                    return 1;
                }
                break;
            }
            case 'r': {
                Refineby2 = true;
                break;
            }
            case 'v': {
                Vlngth = true;
                break;
            }
            case 'x': {
                Format = true;
                break;
            }
            case 'y': {
                Symmetrize = true;
                break;
            }
            case 'z': {
                Cut = true;
                break;
            }
                
            case '?': {
                std::cerr << "Unknown option character " << optopt << std::endl;
                std::cerr << "Use tri_mesh -h for usage" << std::endl;
                return 1;
            }
            default: {
                std::cerr << "Use tri_mesh -h for usage" << std::endl;
                return 1;
            }
        }
    }
    int index = optind;
	
	tet_mesh::filetype in = static_cast<tet_mesh::filetype>(informat);
	tet_mesh::filetype out = static_cast<tet_mesh::filetype>(outformat);
	std::string bdry_nm(std::string(argv[1]) +"_bdry.inpt");
	ifstream intest;
	input_map inmap;
	intest.open(bdry_nm.c_str());
	if (intest) {
		intest.close();
		inmap.input(bdry_nm);
		inmap.echo = true;
		std::cout << "Using " << bdry_nm << std::endl;
	}
	//
	//    if (Cut) {
	//        class tet_mesh zx;
	//        zx.input(argv[1],in,1.0,inmap);
	//        for(int i=0;i<zx.npnt;++i)
	//            zx.gbl->fltwk(i) = zx.pnts(i)(0)*zx.pnts(i)(0) +zx.pnts(i)(1)*zx.pnts(i)(1) - 0.25;
	//
	//        zx.cut();
	//
	//        return(0);
	//    }
	//
	//
	//    /* TO SYMMETRIZE A MESH */
	//    if (Symmetrize) {
	//        class tet_mesh zx;
	//        zx.input(argv[1],in,8.0,inmap);
	//        zx.symmetrize();
	//        return 0;
	//    }
	//
	
	if (Vlngth) {
		class tet_mesh zx;
		
		zx.input(argv[1],in,8.0,inmap);
		std::string name;
		name = std::string(argv[1]) +".lngth";
		FILE *fp = fopen(name.c_str(),"w");
		for(int i=0;i<zx.npnt;++i) fprintf(fp,"%e\n",0.3); // 5.*zx.lngth(i));
		fclose(fp);
		return 0;
	}
	//
	//    if (Smooth) {
	//        class tet_mesh zx;
	//
	//        zx.input(argv[1],in,8.0,inmap);
	//        zx.smooth_cofa(2);
	//        zx.output(argv[2],out);
	//
	//        return 0;
	//    }
	//
	if (Refineby2) {
		class tet_mesh zx,zy;
		zx.input(argv[1],in,8.0,inmap);
		zy.refineby2(zx);
		zy.checkintegrity();
		zy.output(argv[2],out);
		return 0;
	}
	//
	//    if (Coarsen_hp) {
	//        class tet_mesh zx,zy;
	//
	//        int p;
	//        zx.input(argv[1],in,1.0,inmap);
	//        printf("input p\n");
	//        scanf("%d",&p);
	//        zy.coarsen_substructured(zx,p);
	//        zy.output(argv[2],out);
	//        return 0;
	//    }
	//
	if (Scale) {
		class tet_mesh zx;
		TinyVector<FLT,tet_mesh::ND> s;
		printf("Enter x y and z scaling\n");
		scanf("%le%le%le",&s(0),&s(1),&s(2));
		zx.input(argv[1],in,1.0,inmap);
		zx.scale(s);
		zx.output(argv[2],out);
		return 0;
	}
	
	if (Shift) {
		class tet_mesh zx;
		
		TinyVector<FLT,tet_mesh::ND> s;
		printf("Enter x y and z shift\n");
		scanf("%le%le%le",&s(0),&s(1),&s(2));
		zx.input(argv[1],in,1.0,inmap);
		zx.shift(s);
		zx.output(argv[2],out);
		return 0;
	}
	
	if (Format) {
		class tet_mesh zx;
		zx.input(argv[1],in,1.0,inmap);
		zx.output(argv[2],out);
		return(0);
	}

	
	if (GMSHPartition) {
		/* This separates different volumes in a GMSH mesh */
		/* The volumes must be numbered sequentially starting at 1 */
		class tet_mesh zx;
		std::string fname;
		ostringstream nstr;
		
		zx.input(argv[1],in,1.0,inmap);
		
		/* input calls setinfo but to make sure call it again because partition needs it to work */
		zx.tet_mesh::setinfo();
		
		for(int i=0;i<zx.ntet;++i)
			zx.tet(i).info = zx.gbl->fltwk(i)-1;
						
		Array<int,2> blist;
		Array<int,1> bnum;
		
		zx.setup_partition(p,blist,bnum);
		
		for(int i=0;i<p;++i) {
			nstr << "b" << i << std::flush;
			fname = "partition_" +nstr.str();
			std::cout << nstr.str() << "_mesh: " << fname << std::endl;
			nstr.str("");
			tet_mesh zpart;
			zpart.partition2(zx,i,p,blist,bnum);
			zpart.output(fname,tet_mesh::gmsh);
			zpart.output(fname);
		}
		return(0);
	}
	
	if (GMSHLabel) {
		/* This creates a file containing labels of the elements in different volumes in a GMSH mesh */
		class tet_mesh zx;
		zx.input(argv[1],in,1.0,inmap);

		string gridname(argv[1]), filename;
		size_t dotloc;
		dotloc = gridname.find_last_of('.');
		
		if (dotloc != string::npos) {
			/* Found and ending */
			filename = gridname.substr(0,dotloc) +".marks";
		}
		else {
			filename = gridname +".marks";
		}
		
		
		ofstream fout;
		fout.open(filename.c_str());
		
		for(int i=0;i<zx.ntet;++i)
			fout << i << ": " << static_cast<int>(zx.gbl->fltwk(i)-1) << '\n';
		
		fout.close();
		
		return(0);
	}
	
	if (Partition) {
#ifdef METIS
		class tet_mesh zx;
		std::string fname,mname;
		ostringstream nstr;
		zx.input(argv[1],in,1.0,inmap);
		
		/* input calls setinfo but to make sure call it again because partition needs it to work */
		zx.tet_mesh::setinfo();
		
		/* Uses metis to set tet(i).info to a partition number */
		zx.setpartition(p);
		
		//Array<int,2> blist;
		//Array<int,1> bnum;
		// zx.setup_partition(p,blist,bnum);
		zx.setup_partition2(p);
		
		/* The following is to load a marker file if it exists */
		ifstream fin;
		bool marks_flag = false;
		Array<int,1> marks;
		string gridname(argv[1]), filename;
		size_t dotloc;
		
		dotloc = gridname.find_last_of('.');
		if (dotloc != string::npos) {
			/* Found and ending */
			filename = gridname.substr(0,dotloc) +".marks";
		}
		else {
			filename = gridname +".marks";
		}
		fin.open(filename.c_str());
		if (fin) {
			marks.resize(zx.ntet);
			for (int i=0;i<zx.ntet;++i) {
				fin.ignore(80,':');
				fin >> marks(i);
				marks_flag = true;
			}
		}
		
			
		/* now partition mesh and marks if found */
		for(int i=0;i<p;++i) {
			tet_mesh zpart;
			nstr << "b" << i << std::flush;
			fname = "partition_" +nstr.str();
			mname = fname +".marks";
			// std::cout << nstr.str() << "_mesh: " << fname << std::endl;
			nstr.str("");
			// zpart(i).partition2(zx,i,p,blist,bnum);
			zpart.partition3(zx,i);
			if (marks_flag) {
				ofstream fout;
				fout.open(mname.c_str());
				
				/* tet.info gets wiped out so have to redo it */
				int ntet = 0;
				for(int tind=0;tind<zx.ntet;++tind) {
					if (zx.tet(tind).info == i) {
						fout << ntet << ": " << marks(tind) << '\n';
						++ntet;
					}
				}
				fout.close();
			}
			
			// zpart(i).output(fname,tet_mesh::gmsh);
			zpart.output(fname);
			zpart.checkintegrity();

		}
		
		//        for(int i=0;i<p;++i) {
		//			nstr << "b" << i << std::flush;
		//			fname = "partition_" +nstr.str();
		//			std::cout << nstr.str() << "_mesh: " << fname << std::endl;
		//			nstr.str("");
		//			zpart(i).partition(zx,i,p);
		//			zpart(i).checkintegrity();
		//			zpart(i).output(fname,out);
		//			zpart(i).output(fname,tet_mesh::gmsh);
		//
		//            //zpart(i).output(fname,tet_mesh::boundary);//temp fixme
		//        }
#else
		printf("Need metis package to partition\n");
#endif
		return(0);
	}
	
	//    if (Coarsen_Marks) {
	//        class tet_mesh zx;
	//
	//        zx.input(argv[1],in,1.0,inmap);
	//        FILE *fp = fopen(argv[3],"r");
	//
	//        for(int i=0;i<zx.npnt;++i) {
	//            fscanf(fp,"%d\n",&zx.pnt(i).info);
	//            zx.pnt(i).info = 1-zx.pnt(i).info;
	//        }
	//        zx.coarsen3();
	//        zx.output(argv[2],out);
	//        return(0);
	//    }
	
	
	if (argc == 2) {
		/* READ INPUT MAP FROM FILE */
		sim::blks.go(argv[1]);
	}
	else if (argc == 3) {
		/* READ INPUT MAP FROM FILE & OUTPUT TO FILE */
		sim::blks.go(argv[1],argv[2]);
	}
	
#ifdef PTH
	pth_kill();
#endif
#ifdef MPISRC
	MPI_Finalize();
#endif
	
	return(0);
}


multigrid_interface* block::getnewlevel(input_map& input) {
	int type;
	multigrid_interface *temp;
	
	if (!input.get(idprefix+"_type",type)) input.getwdefault("type",type,1);
	
	switch(type) {
		default: {
			temp = new tet_mesh();
			break;
		}
	} 
	
	return(temp);
}


