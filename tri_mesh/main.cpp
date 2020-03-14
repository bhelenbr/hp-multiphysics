#include "blocks.h"
#include <string>
#include "tri_mesh.h"
#include "math.h"
#include <input_map.h>
#include <symbolic_function.h>
#include <chrono>
#include <thread>
#include <unistd.h>

#ifdef MPISRC
#include <mpi.h>
#endif

using namespace std;

int main(int argc, char *argv[]) {
	clock_t cpu_time;

#ifdef MPISRC
	int myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
#ifdef PTH
	// For debugging put interrupt here
	// On interrupt type this into gdb console:
	// for gdb: handle SIGUSR1 nostop print pass
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
    bool Generate = false;
    bool Coarsen_hp = false;
    bool Shift = false;
    bool Partition = false;
    bool Refineby2 = false;
    bool Scale = false;
    bool Vlngth = false;
    bool Format = false;
    bool Symmetrize = false;
    bool Cut = false;
    int informat = 3, outformat = 3;
    int p;
    int opt;
    while ((opt = getopt (argc, argv, "acdefghi:lmo:p:rsvxyz")) != -1) {
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
                Generate = true;
                break;
            }
            case 'h': {
                std::cout << "tri_mesh utility" << std::endl;
                std::cout << "mod_map [-acdefghlmpsvxyz] [-i input format] [-o output format] inputfile outputfile" << std::endl;
                std::cout << " -h prints usage information" << std::endl;
                std::cout << "-a append two meshes" << std::endl;
                std::cout << "-c Coarsen mesh by factor of 2" << std::endl;
                std::cout << "-d Stop for Debugger" << std::endl;
                std::cout << "-e display basic mesh information" << std::endl;
                std::cout << "-f smooth mesh" << std::endl;
                std::cout << "-g generate mesh from .d file" << std::endl;
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
                    std::cerr << "Use tri_mesh -h for usage" << std::endl;
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
                    std::cerr << "Unable to read number of partitions.  Use tri_mesh -h for usage" << std::endl;
                    return 1;
                }
                break;
            }
            case 'r': {
                Refineby2 = true;
                break;
            }
            case 's': {
                Generate = true;
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

 
#if (defined(MPISRC) && !defined(petsc))
	if (Debugger) {
		int size;
		/*
		 we have to make sure that all processors have opened
		 connections to all other processors, otherwise once the
		 debugger has stated it is likely to receive a SIGUSR1
		 and kill the program.
		 */
		int ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);
		if (ierr != MPI_SUCCESS) {
			exit(1);
		}
		if (size > 2) {
			int dummy = 0;
			MPI_Status status;
			for (int i=0; i<size; i++) {
				if (myid != i) {
					ierr = MPI_Send(&dummy,1,MPI_INT,i,109,MPI_COMM_WORLD);
				}
			}
			for (int i=0; i<size; i++) {
				if (myid != i) {
					ierr = MPI_Recv(&dummy,1,MPI_INT,i,109,MPI_COMM_WORLD,&status);				}
			}
		}
		std::cout << "Waiting for debugger.  Process id is " << getpid() << std::endl;
#if __cplusplus >= 199711L
        std::this_thread::sleep_for(std::chrono::milliseconds(10000));
#endif
	}
#endif
	
	tri_mesh::filetype in = static_cast<tri_mesh::filetype>(informat);
	tri_mesh::filetype out = static_cast<tri_mesh::filetype>(outformat);
    
    /***********************************************************/
    /* Single input filename commands */
    /***********************************************************/
    std::string input_file;
    if (argc -index < 1) {
        std::cerr << "missing input filename" << std::endl;
        return(1);
    }
    else {
        input_file = argv[index++];
    }
	std::string bdry_nm(input_file +"_bdry.inpt");
	ifstream intest;
	input_map inmap;
	intest.open(bdry_nm.c_str());
	if (intest) {
		intest.close();
		inmap.input(bdry_nm);
		inmap.echo = true;
		std::cout << "Using " << bdry_nm << std::endl;
	}

	if (Echo) {
		class tri_mesh zx;
		zx.input(input_file,in,1.0,inmap);
		std::cout << "npnt: " << zx.npnt << " nseg: " << zx.nseg << " ntri: " << zx.ntri << std::endl;
		return(0);
	}
	

	if (Cut) {
		class tri_mesh zx;
		zx.input(input_file,in,1.0,inmap);
		for(int i=0;i<zx.npnt;++i)
			zx.gbl->fltwk(i) = zx.pnts(i)(0)*zx.pnts(i)(0) +zx.pnts(i)(1)*zx.pnts(i)(1) - 0.25;
		zx.cut();
		return(0);
	}
    
    /* TO SYMMETRIZE A MESH */
    if (Symmetrize) {
        class tri_mesh zx;
        zx.input(input_file,in,8.0,inmap);
        zx.symmetrize();
        return 0;
    }


    /***********************************************************/
    /* Two file input filename commands */
    /***********************************************************/
    std::string output_file;
    if (argc-index < 1) {
        output_file = "output";
    }
    else {
        output_file = argv[index++];
    }

	if (Vlngth) {
		class tri_mesh zx;
		zx.input(input_file,in,8.0,inmap);
		symbolic_function<2> length_modifier_function;
		input_map length_modifier_input;
		std::string length_modifier_string;
		std::cout << "Input function of (x0,x1) and t where t is the current mesh length" << std::endl;
		std::cin >> length_modifier_string;
		length_modifier_input["length_modifier"] = length_modifier_string;
		length_modifier_function.init(length_modifier_input,"length_modifier");
		for(int i=0;i<zx.npnt;++i) 
			zx.lngth(i) = length_modifier_function.Eval(zx.pnts(i),zx.lngth(i));
		zx.output(output_file,out);
		return 0; 
	}

	if (Smooth) {
		class tri_mesh zx;
		zx.input(input_file,in,8.0,inmap);
		zx.smooth_cofa(2);
		zx.output(output_file,out);

		return 0;
	}

	if (Refineby2) {
		class tri_mesh zx,zy;
		zx.input(input_file,in,8.0,inmap);
		zy.refineby2(zx);
		zy.checkintegrity();
		zy.output(output_file,out);
		return 0;
	}

	if (Coarsen_hp) {
		class tri_mesh zx,zy;

		int p;
		zx.input(input_file,in,1.0,inmap);
		printf("input p\n");
		scanf("%d",&p);
		zy.coarsen_substructured(zx,p);
		zy.output(output_file,out);
		return 0;
	}

	if (Scale) {
		class tri_mesh zx;
		TinyVector<FLT,2> s;
		printf("Enter x and y scaling\n");
		scanf("%le%le",&s(0),&s(1));
		zx.input(input_file,in,1.0,inmap);
		zx.scale(s);
		zx.output(output_file,out);
		return 0;
	}

	if (Shift) {
		class tri_mesh zx;

		TinyVector<FLT,2> s;
		printf("Enter x and y shift\n");
		scanf("%le %le",&s(0),&s(1));
		zx.input(input_file,in,1.0,inmap);
		zx.shift(s);
		zx.output(output_file,out);
		return 0;
	}

	if (Format) {
		class tri_mesh zx;
		zx.input(input_file,in,1.0,inmap);
		zx.output(output_file,out);
		return(0);
	}

	if (Partition) {
#ifdef METIS
		class tri_mesh zx;
		if (p > 1) {
			/* Load mesh */
			std::string fname;
			ostringstream nstr;
			zx.input(input_file,in,1.0,inmap);
			
			/* If there is a marks file then try to subdivide appended multiphysics mesh */
			ifstream marks_file;
			std::string marksfilename = input_file +".marks";
			marks_file.open(marksfilename.c_str());
			if (marks_file) {
				for(int i=0;i<zx.ntri;++i)
					marks_file >> zx.tri(i).info;
				zx.subpartition(p);
			}
			else {
				zx.setpartition(p);
			}

			/* Extract partitions */
			for(int i=0;i<p;++i) {
				tri_mesh zpart;
				nstr << "b" << i << std::flush;
				fname = "partition_" +nstr.str();
				std::cout << nstr.str() << "_mesh: " << fname << std::endl;
				zpart.partition(zx,i);
				zpart.checkintegrity();
				zpart.output("partition_"+nstr.str(),out);
				nstr.str("");
			}
		}
#else
		printf("Need metis package to partition\n");
#endif
		return(0);
	}
	
	if (Coarsenby2) {
		class tri_mesh zx,zy;

		zx.input(input_file,in,1.0,inmap);
		
		for(int i=0;i<zx.npnt;++i) {
			zx.lngth(i) *= 2.0;
		}
		clock();
		zy.coarsen(1.6,zx);
		cpu_time = clock();
		std::cout << "that took " << cpu_time << " cpu time" << std::endl;

		zy.output(output_file,out);
		return(0);
	}
	
//	if (Coarsen_Marks) {
//		class tri_mesh zx;
//
//		zx.input(input_file,in,1.0,inmap);
//		FILE *fp = fopen(argv[3],"r");
//
//		for(int i=0;i<zx.npnt;++i) {
//			fscanf(fp,"%d\n",&zx.pnt(i).info);
//			zx.pnt(i).info = 1-zx.pnt(i).info;
//		}
//		clock();
//		zx.coarsen3();
//		cpu_time = clock();
//		std::cout << "that took " << cpu_time << " cpu time" << std::endl;
//
//		zx.output(output_file,out);
//		return(0);
//	}
    
    if (Append) {
        class tri_mesh zx,zy;
        if (argc -index < 1) {
            std::cerr << "missing file to append or output filename" << std::endl;
            return(1);
        }
        zx.input(input_file,in,100.0,inmap);
        zy.input(output_file,in,1.0,inmap);
        zx.append(zy);
        zx.cleanup_after_adapt();
        zx.output(argv[index],out);
        /* Output a marks file so you know which element came from where */
        ofstream file;
        std::string filename;
        filename = std::string(argv[index]) +".marks";
        file.open(filename.c_str());
        for(int i=0;i<zx.ntri-zy.ntri;++i) {
            file << 0 << std::endl;
        }
        for (int i=0;i<zy.ntri;++i) {
            file << 1 << std::endl;
        }
        file.close();
        
        return(0);
    }

	if (argc == 2) {
		/* READ INPUT MAP FROM FILE */
		sim::blks.go(input_file);
	}
	else if (argc == 3) {
		/* READ INPUT MAP FROM FILE & OUTPUT TO FILE */
		sim::blks.go(input_file,output_file);
	}

#ifdef PTH
	pth_kill();
#endif
#ifdef MPISRC
	MPI_Finalize();
#endif

	return(0);
}
