#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"mesh.h"

int mesh::out_mesh(FLT (*vin)[ND], const char *filename, FILETYPE filetype = easymesh) const {
	char fnmapp[100];
	FILE *out;
	int i,j,tind,count;
   
   switch (filetype) {
   
      case (easymesh):
      /*	CREATE EASYMESH OUTPUT FILES */
         strcpy(fnmapp,filename);
         strcat(fnmapp,".n");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }
         fprintf(out,"%d\n",nvrtx);
         for(i=0;i<nvrtx;++i)
            fprintf(out,"\t%d: %17.10e %17.10e %d\n",i,vin[i][0],vin[i][1],vinfo[i]); 
            
         fclose(out);

      /*	SIDE FILE */		
         strcpy(fnmapp,filename);
         strcat(fnmapp,".s");
         out = fopen(fnmapp,"w");	
         fprintf(out,"%d\n",nside);
         for(i=0;i<nside;++i) {
            fprintf(out,"%d: %d %d %d %d %d\n",i,svrtx[i][0],svrtx[i][1],
            stri[i][0],stri[i][1],sinfo[i]);
         }	
         fclose(out);
      
         strcpy(fnmapp,filename);
         strcat(fnmapp,".e");
         out = fopen(fnmapp,"w");
         fprintf(out,"%d\n",ntri);
         for(i=0;i<ntri;++i) {
            fprintf(out,"%d: %d %d %d %d %d %d %d %d %d 0.0 0.0 %d\n",
            i,tvrtx[i][0],tvrtx[i][1],tvrtx[i][2],
            ttri[i][0],ttri[i][1],ttri[i][2],
            tside[i].side[0],tside[i].side[1],tside[i].side[2],tinfo[i]); 
         }	
      
         fclose(out);
         break;
      
      case (tecplot):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".dat");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }

         fprintf(out,"ZONE F=FEPOINT, ET=TRIANGLE, N=%d, E=%d\n",nvrtx,ntri);
         
         for(i=0;i<nvrtx;++i)
            fprintf(out,"%e  %e  \n",vin[i][0],vin[i][1]);

        	fprintf(out,"\n#CONNECTION DATA#\n");
         
         for(i=0;i<ntri;++i)
            fprintf(out,"%d %d %d\n",tvrtx[i][0]+1,tvrtx[i][1]+1,tvrtx[i][2]+1);
            
         fclose(out);
         break;

      case (gambit):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".FDNEUT");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }

         count = ntri;
         for(i=0;i<nsbd;++i)
            count += sbdry[i].num;
	
         fprintf(out,"** FIDAP NEUTRAL FILE\n");
         fprintf(out,"%s\n",filename);
         fprintf(out,"VERSION    8.01\n");
         fprintf(out,"29 Nov 1999    13:23:58\n");
         fprintf(out,"   NO. OF NODES   NO. ELEMENTS NO. ELT GROUPS          NDFCD          NDFVL\n");
         fprintf(out,"%15d%15d%15d%15d%15d\n",nvrtx,count,nsbd+1,2,2);
         fprintf(out,"   STEADY/TRANS     TURB. FLAG FREE SURF FLAG    COMPR. FLAG   RESULTS ONLY\n");
         fprintf(out,"%12d%12d%12d%12d%12d\n",0,0,0,0,0);
         fprintf(out,"TEMPERATURE/SPECIES FLAGS\n");
         fprintf(out," 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n");
         fprintf(out,"PRESSURE FLAGS - IDCTS, IPENY MPDF\n");
         fprintf(out,"         1         1         0\n");
         fprintf(out,"NODAL COORDINATES\n");
         
         for(i=0;i<nvrtx;++i) 
            fprintf(out,"%10d%20.10e%20.10e\n",i+1,vin[i][0],vin[i][1]);
            
         fprintf(out,"BOUNDARY CONDITIONS\n");
         fprintf(out,"         0         0         0     0.0\n");
         fprintf(out,"ELEMENT GROUPS\n");
         fprintf(out,"GROUP:%9d ELEMENTS:%10d NODES:%13d GEOMETRY:%5d TYPE:%4d\n",
         1,ntri,3,2,2);
         fprintf(out,"ENTITY NAME:   fluid\n");
         
         for(tind=0;tind<ntri;++tind) {
         	fprintf(out,"%8d%8d%8d%8d\n",tind+1,tvrtx[tind][0]+1,tvrtx[tind][1]+1,tvrtx[tind][2]+1);
         }
        
         count = ntri+1;
        
         for(i=0;i<nsbd;++i) {
            fprintf(out,"GROUP:%9d ELEMENTS:%10d NODES:%13d GEOMETRY:%5d TYPE:%4d\n",i+2,sbdry[i].num,2,0,sbdry[i].type);
            fprintf(out,"ENTITY NAME:   surface.1\n");
            for(j=0;j<sbdry[i].num;++j) 
                  fprintf(out,"%8d%8d%8d\n",count++,svrtx[sbdry[i].el[j]][0]+1,svrtx[sbdry[i].el[j]][1]+1);
         }
         fclose(out);
         break;

      case(text):
/*			JUST OUTPUT VERTEX POSITIONS FOR DEFORMING MESH */
         strcpy(fnmapp,filename);
         strcat(fnmapp,".txt");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }
         fprintf(out,"%d\n",nvrtx);
         for(i=0;i<nvrtx;++i)
            fprintf(out,"\t%d: %17.10e %17.10e\n",i,vin[i][0],vin[i][1]); 
            
         fclose(out);
         break;
         
      case(grid):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".grd");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }
         
/*			HEADER LINES */
         fprintf(out,"nvrtx: %d\t nside: %d\t ntri: %d\n",nvrtx,nside,ntri);

/*			VRTX INFO */                        
         for(i=0;i<nvrtx;++i)
            fprintf(out,"%17.10e %17.10e\n",vin[i][0],vin[i][1]);
                    
/*			SIDE INFO */
         for(i=0;i<nside;++i)
            fprintf(out,"%d %d\n",svrtx[i][0],svrtx[i][1]);

/*			THEN TRI INFO */
         for(i=0;i<ntri;++i)
            fprintf(out,"%d %d %d\n",tvrtx[i][0],tvrtx[i][1],tvrtx[i][2]);

/*			SIDE BOUNDARY INFO HEADER */
         fprintf(out,"nsbd: %d\n",nsbd);
         for(i=0;i<nsbd;++i)
            fprintf(out,"type: %d\t number %d\n",sbdry[i].type,sbdry[i].num);
         
         for(i=0;i<nsbd;++i)
            for(j=0;j<sbdry[i].num;++j)
               fprintf(out,"%d\n",sbdry[i].el[j]);
               
/*			VERTEX BOUNDARY INFO HEADER */
         fprintf(out,"nvbd: %d\n",nvbd);
         for(i=0;i<nvbd;++i)
            fprintf(out,"%d %d\n",vbdry[i].type,vbdry[i].num);
         
         for(i=0;i<nvbd;++i)
            for(j=0;j<vbdry[i].num;++j)
               fprintf(out,"%d\n",vbdry[i].el[j]);
         
         fclose(out);
         break;

      default:
         printf("That filetype is not supported yet\n");
         break;
   }

    return(1);
}

void mesh::setbcinfo() {
   int i,j;
   
/*	SET UP VRTX BC INFORMATION FOR OUTPUT */
   for(i=0;i<nvrtx;++i)
      vinfo[i] = 0;
   
   for(i=0;i<nvbd;++i)
      for(j=0;j<vbdry[i].num;++j)
         vinfo[vbdry[i].el[j]] = vbdry[i].type;

/*	SET UP SIDE BC INFORMATION FOR EASYMESH OUTPUT */
   for(i=0;i<nside;++i)
      sinfo[i] = 0;
   
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i].num;++j)
         sinfo[sbdry[i].el[j]] = sbdry[i].type;

/*	SET UP TRI INFO FOR EASYMESH OUTPUT */         
   for(i=0;i<ntri;++i)
      tinfo[i] = 0;
         
   return;
}
