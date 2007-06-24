/*
 *  box_impl.h
 *  quadtree
 *
 *  Created by Brian Helenbrook on Fri Sep 13 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include"quadtree.h"
#include<utilities.h>
#include<float.h>
#include<cstdio>
#include<assert.h>
#include<string.h>
#include<math.h>

/* STATIC VARIABLES FOR SEARCHING */ 
// class box<ND> **quadtree<ND>::srchlst;
// int quadtree<ND>::maxsrch = 0;
  
template<int ND> void quadtree<ND>::allocate(FLT (*v)[ND], int mxv) {
    int i;
    
    vrtx = v;
    maxvrtx = mxv;
    if (size != 0) {
        delete [] base;
        delete [] indx;
    }
    
/*	TAKES APPROXIMATELY 1:1 boxS TO VERTICES WITH 4 NODES / box */
/*	+10 IS FACTOR OF SAFETY FOR SMALL mxv */
    size = (int) 1.1*mxv +10; 
    base = new class box<ND>[size];
    indx = new class box<ND>*[maxvrtx];
    for(i=0;i<maxvrtx;++i)
        indx[i] = NULL;

    FLT x1[ND], x2[ND];
    for(i=0;i<ND;++i) {
        x1[i] = 0.0;
        x2[i] = 1.0;
    }
    base[0] = box<ND>(NULL,0,x1,x2);
    current = 1;
    
    if (maxsrch == 0) {
        srchlst = new class box<ND>*[size];
        maxsrch = size;
    }
    else {
        if (size > maxsrch) { 
            printf("#Warning: Better to allocate largest quadtree to smallest\n");
            delete []srchlst;
            srchlst = new class box<ND>*[size];
            maxsrch = size;
        }
    }
}

template<int ND> void quadtree<ND>::init(FLT x1[ND], FLT x2[ND]) {
    int i;
    
    for(i=0;i<ND;++i) {
        if(x1[i] >= x2[i]) {
         	printf("quadtree initialization error: zero domain area %d %f %f\n",i,x1[i],x2[i]);
            exit(1);
        }
    }
    
    

    base[0] = box<ND>(NULL,0,x1,x2);

    for(i = 0; i<maxvrtx; ++i)
        indx[i] = NULL;
        
    current = 1;
}

template<int ND> void quadtree<ND>::init() {

    base[0] = box<ND>(NULL,0,base[0].xmin,base[0].xmax);
    for(int i = 0; i<maxvrtx; ++i)
        indx[i] = NULL;
        
    current = 1;
}

template<int ND> void quadtree<ND>::copy(const class quadtree<ND>& tgt) {
    int i,j,n;
        
    vrtx = tgt.vrtx;
    if (size == 0) {
        maxvrtx = tgt.maxvrtx;
        size = tgt.size;
        base = new class box<ND>[size];
        indx = new class box<ND>*[maxvrtx];
    }
    else {
        assert(size >= tgt.current);
    }

    if (maxvrtx < tgt.maxvrtx) {
        for(i=0;i<maxvrtx;++i) {
            if (tgt.indx[i] == NULL) indx[i] = NULL;
            else indx[i] = base + (tgt.indx[i] -tgt.base);
        }
    
        for(i=maxvrtx;i<tgt.maxvrtx;++i)
            if (tgt.indx[i] != NULL) {
                printf("quadtree is too small for copy\n");
                exit(1);
        }
    }
    else {
        for(i=0;i<tgt.maxvrtx;++i) {
            if (tgt.indx[i] == NULL) indx[i] = NULL;
            else indx[i] = base + (tgt.indx[i] -tgt.base);
        }
        for(i=tgt.maxvrtx;i<maxvrtx;++i)
            indx[i] = NULL;
    }

    current = tgt.current;
    
    for(i=0;i<current;++i) {
        base[i].num = tgt.base[i].num;
        base[i].prnt = base +(tgt.base[i].prnt -tgt.base);
        base[i].pind = tgt.base[i].pind;
        for(n=0;n<ND;++n) {
            base[i].xmin[n] = tgt.base[i].xmin[n];
            base[i].xmax[n] = tgt.base[i].xmax[n];
        }
        if (base[i].num > 0) {
            for(j=0;j<(1<<ND);++j)
                base[i].node[j] = tgt.base[i].node[j];
        }
        else {
            for(j=0;j<(1<<ND);++j)
                base[i].dghtr[j] = base +(tgt.base[i].dghtr[j] -tgt.base);
        }
    }
    
    base[0].prnt = NULL;
    
    return;
}


template<int ND> void quadtree<ND>::reinit() {
    int i;
        
    assert(size != 0); 
    
    base[0].num = 0;
    current = 1;
    
    for(i=0;i<maxvrtx;++i)
        if (indx[i] != NULL) addpt(i);
    
    return;
}

template<int ND> void quadtree<ND>::addpt(int v0, class box<ND>* start) {
    int i,j,n;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
    class box<ND> *qpt;  // SAME
    int store[(1<<ND) +1];
    FLT avoid0,xshift,x1[ND],x2[ND],dx[ND];
    
    if (start == NULL) 
        qpt = base;
    else
        qpt = start;

/*	LOOP TO BOTTOM OF TREE */
    while (qpt->num < 0) {
        i = 0;
        for(n=0;n<ND;++n) {
            xshift = vrtx[v0][n] -0.5*(qpt->xmax[n] +qpt->xmin[n]);
            avoid0 = MAX(fabs(xshift), 10.*EPSILON);
            i = (i<<1) +((int) ((1+10.*EPSILON)*xshift/avoid0) +1)/2;
        }
        qpt = qpt->dghtr[i];
    }
    
    if (qpt->num < (1<<ND)) {
/*		box CAN ACCEPT NEW POINT */
        qpt->node[qpt->num++] = v0;
        indx[v0] = qpt;
        return;
    }

/*	box IS FULL SUBDIVIDE */
    if (current +(1<<ND) >= size) {
        printf("Need to allocate bigger quadtree %d\n",size);
        output("quad_error",quadtree::text);
        exit(1);
    }
    for (i=0;i<(1<<ND);++i)
        store[i] = qpt->node[i];
    store[(1<<ND)] = v0;
    
    qpt->num = -1;
    for(n=0;n<ND;++n)
        dx[n] = (qpt->xmax[n] -qpt->xmin[n])*0.5;

    for(i=0;i<(1<<ND);++i) {
        for(n=0;n<ND;++n) {
            j = (i>>(ND-n-1))&1;
            x1[n] = qpt->xmin[n] +dx[n]*j;
            x2[n] = x1[n] +dx[n];
        }
        base[current] = box<ND>(qpt,i,x1,x2);
        qpt->dghtr[i] = base +current;
        ++current;
    }
    
    for(i=0;i<(1<<ND)+1;++i)
        addpt(store[i],qpt);
        
    return;
}

/*	FIND CLOSEST POINT TO A POINT IN THE ARRAY */
template<int ND> FLT quadtree<ND>::nearpt(int const v0, int& pt) const {
    int i,n,nsrch,exclude;
    class box<ND> *topbox;
    FLT dist,dx[ND],xh[ND],xl[ND],mindist;
    class box<ND> *qpt;

    mindist = 0.0;
    for(n=0;n<ND;++n) {
        xl[n] = base[0].xmin[n];
        xh[n] = base[0].xmax[n];
        mindist += pow(xh[n]-xl[n],2);
    }
    nsrch = 0;
    
    if (indx[v0] == NULL) 
        srchlst[nsrch++] = base;
    else
        srchlst[nsrch++] = indx[v0];
        
    topbox = srchlst[0];	 

    for(;;) {
        while (nsrch > 0) {
            qpt = srchlst[--nsrch];
            for(n=0;n<ND;++n)
                if(qpt->xmin[n] > xh[n] || qpt->xmax[n] < xl[n]) goto NEXT;
            
            if (qpt->num < 0) {
/*				ADD DAUGHTERS TO LIST */
                for(i=0;i<(1<<ND);++i)
                    srchlst[nsrch++] = qpt->dghtr[i];
                continue;
            }
            else {
/*				CHECK POINTS */
                for(i=0;i<qpt->num;++i) {
                    if (qpt->node[i] == v0) continue;
                    dist = 0.0;
                    for(n=0;n<ND;++n) {
                        dx[n] = vrtx[v0][n] -vrtx[qpt->node[i]][n];
                        dist += pow(dx[n],2);
                    }
                    if (dist < mindist) {
                        pt = qpt->node[i];
                        mindist = dist;
                        dist = sqrt(dist);
                        for(n=0;n<ND;++n) {
                            xh[n] = vrtx[v0][n] +dist;
                            xl[n] = vrtx[v0][n] -dist;
                        }
                    }
                }
            }
            NEXT: continue;
        }
        for(n=0;n<ND;++n)
            if (topbox->xmin[n] > xl[n] || topbox->xmax[n] < xh[n]) goto MORECHECK;            
        return(mindist);

        MORECHECK:
        exclude = topbox->pind;
        topbox = topbox->prnt;

/*		CHECK IF WE ARE AT THE TOP LEVEL */
        if (topbox == NULL) return(mindist);
        
        for(i=0;i<(1<<ND);++i)
            if (i != exclude) srchlst[nsrch++] = topbox->dghtr[i];
    }

    return(mindist);
}

/*	FIND CLOSEST POINT TO A POINT IN THE ARRAY */
template<int ND> FLT quadtree<ND>::nearpt(FLT const x[ND], int& pt) const {
    int i,n,nsrch,exclude;
    class box<ND> *topbox;
    FLT dist,dx[ND],xh[ND],xl[ND],mindist;
    class box<ND> *qpt;

    mindist = 0.0;
    for(n=0;n<ND;++n) {
        xl[n] = base[0].xmin[n];
        xh[n] = base[0].xmax[n];
        mindist += pow(xh[n]-xl[n],2);
    }
    nsrch = 0;
    srchlst[nsrch++] = base;        
    topbox = srchlst[0];	 

    for(;;) {
        while (nsrch > 0) {
            qpt = srchlst[--nsrch];
            for(n=0;n<ND;++n)
                if(qpt->xmin[n] > xh[n] || qpt->xmax[n] < xl[n]) goto NEXT;
            
            if (qpt->num < 0) {
/*				ADD DAUGHTERS TO LIST */
                for(i=0;i<(1<<ND);++i)
                    srchlst[nsrch++] = qpt->dghtr[i];
                continue;
            }
            else {
/*				CHECK POINTS */
                for(i=0;i<qpt->num;++i) {
                    dist = 0.0;
                    for(n=0;n<ND;++n) {
                        dx[n] = x[n] -vrtx[qpt->node[i]][n];
                        dist += pow(dx[n],2);
                    }
                    if (dist < mindist) {
                        pt = qpt->node[i];
                        mindist = dist;
                        dist = sqrt(dist);
                        for(n=0;n<ND;++n) {
                            xh[n] = x[n] +dist;
                            xl[n] = x[n] -dist;
                        }
                    }
                }
            }
            NEXT: continue;
        }
        for(n=0;n<ND;++n)
            if (topbox->xmin[n] > xl[n] || topbox->xmax[n] < xh[n]) goto MORECHECK;            
        return(mindist);

        MORECHECK:
        exclude = topbox->pind;
        topbox = topbox->prnt;

/*		CHECK IF WE ARE AT THE TOP LEVEL */
        if (topbox == NULL) return(mindist);
        
        for(i=0;i<(1<<ND);++i)
            if (i != exclude) srchlst[nsrch++] = topbox->dghtr[i];
    }

    return(mindist);
}


template<int ND> void quadtree<ND>::dltpt(int v0) {
    int ind, i;
    class box<ND> *qpt;
    
    qpt = indx[v0];
    
    if (qpt == NULL) {
        printf("error: deleting point which isn't in tree %d\n",v0);
        exit(1);
    }
    
    indx[v0] = NULL;

    ind = -1;
    for(i=0;i<qpt->num;++i) {
        if (qpt->node[i] == v0) {
            ind = i;
            break;
        }
    }
    
    if (ind == -1) {
        printf("Couldn't find node to remove: %d\n",v0);
        exit(1);
    }
    
    for(i=ind+1;i<qpt->num;++i)
        qpt->node[i-1] = qpt->node[i];
    
    --qpt->num;
    
    return;
}

template<int ND> void quadtree<ND>::output(char *filename, FILETYPE type) {
    int i,j,n,nsrch;
    char fnmapp[100],etype[20],order[20];
    class box<ND> *qpt;
    FILE *out;
    
    strcpy(fnmapp,filename);
    
    switch (type) {
        case(text):
            strcat(fnmapp,".text");
            out = fopen(fnmapp,"w");
            if (out == NULL ) {
                printf("couldn't open output file %s for output\n",fnmapp);
                exit(1);
            }
            nsrch = 0;
            srchlst[nsrch++] = base;
            
            while (nsrch > 0) {
                qpt = srchlst[--nsrch];
                if (qpt->num < 0) {
        /*			ADD DAUGHTERS TO LIST */
                    for(i=0;i<(1<<ND);++i)
                        srchlst[nsrch++] = qpt->dghtr[i];
                    continue;
                }
                else {
        /*			PRINT POINTS */
                    fprintf(out,"box: \n");
                    for(n=0;n<ND;++n)
                        fprintf(out,"(%f %f) ",qpt->xmin[n],qpt->xmax[n]);
                    fprintf(out,"\n");
                    for(i=0;i<qpt->num;++i) {
                        fprintf(out,"%d ",i);
                        for(n=0;n<ND;++n)
                            fprintf(out,"%f ",vrtx[qpt->node[i]][n]);
                        fprintf(out,"\n");
                    }
                }
            }
            break;
        case(tecplot):
            if (ND == 2) {
                strcpy(etype,"QUADRILATERAL");
                strcpy(order,"1 3 4 2\n");
            }
            else if (ND == 3) {
                strcpy(etype,"BRICK");
                strcpy(order,"1 3 7 5 2 4 8 6");
            }
            strcat(fnmapp,".dat");
            out = fopen(fnmapp,"w");
            if (out == NULL ) {
                printf("couldn't open output file %s for output\n",fnmapp);
                exit(1);
            }
            nsrch = 0;
            srchlst[nsrch++] = base;
            
            while (nsrch > 0) {
                qpt = srchlst[--nsrch];
                if (qpt->num < 0) {
        /*			ADD DAUGHTERS TO LIST */
                    for(i=0;i<(1<<ND);++i)
                        srchlst[nsrch++] = qpt->dghtr[i];
                    continue;
                }
                else {
    /*				PRINT BOX & POINTS */
                    if (qpt->num > 0) {
                        fprintf(out,"ZONE N=%d, E=1, F=FEPOINT, ET=%s, C=RED\n",(1<<ND),etype);
                        for(i=0;i<(1<<ND);++i) {
                            for(n=0;n<ND;++n) {
                                j = (i>>(ND-n-1))&1;
                                if (j) fprintf(out,"%f  ",qpt->xmax[n]);
                                else fprintf(out,"%f  ",qpt->xmin[n]);
                            }
                            fprintf(out,"\n");
                        }
                        fprintf(out,"%s\n",order);
                        
                        fprintf(out,"ZONE I=%d, C=BLACK\n",qpt->num);
                        for(i=0;i<qpt->num;++i) {
                            for(n=0;n<ND;++n)
                                fprintf(out,"%f  ",vrtx[qpt->node[i]][n]);
                            fprintf(out,"\n");
                        }
                    }
                }
            }
            break;
    }

    return;
}

template<int ND> void quadtree<ND>::movept(int from, int to) {
    int i;
    class box<ND> *qpt;
    
    qpt = indx[from];
    
    if (qpt == NULL) {
        printf("error: moving point %d which isn't in tree to %d\n",from,to);
        exit(1);
    }
    
    if (indx[to] != NULL) {
        printf("error: moving %d to a point which exists %d\n",from,to);
        exit(1);
    }
    
    indx[to] = indx[from];
    indx[from] = NULL;

    for(i=0;i<qpt->num;++i) {
        if (qpt->node[i] == from) {
            qpt->node[i] = to;
            break;
        }
    }
    assert(i != qpt->num);

    return;
}

template<int ND> void quadtree<ND>::update(int bgn, int end) {
    int i,n;
    class box<ND> *qpt;
    
    for(i=bgn;i<end;++i) {
        qpt = indx[i];
        if (qpt == NULL) continue;
        for(n=0;n<ND;++n) {
            if (vrtx[i][n] < qpt->xmin[n] || vrtx[i][n] > qpt->xmax[n]) {
                update(i);
                break;
            }
        }
    }
    
    return;
}

template<int ND> void quadtree<ND>::update(int v0) {
    int i,n,nsrch,exclude;
    FLT x[ND];
    class box<ND> *topbox;
    class box<ND> *qpt;

    for(n=0;n<ND;++n)
        x[n] = vrtx[v0][n];
        
    nsrch = 0;
    srchlst[nsrch++] = indx[v0];
    topbox = srchlst[0];	

    for(;;) {
        while (nsrch > 0) {
            qpt = srchlst[--nsrch];
 
            if (qpt->num < 0) {
/*				ADD DAUGHTERS TO LIST */
                for(i=0;i<(1<<ND);++i)
                    srchlst[nsrch++] = qpt->dghtr[i];
                continue;
            }
            else {
                for(n=0;n<ND;++n)
                    if(qpt->xmin[n] > x[n] || qpt->xmax[n] < x[n]) goto UPDATE;
                continue;
                UPDATE: dltpt(v0);
                addpt(v0, qpt);
                return;
            }
        }
        
        exclude = topbox->pind;
        topbox = topbox->prnt;

/*		CHECK IF WE ARE AT THE TOP LEVEL */
        if (topbox == NULL) {
            printf("#Warning: point is outside of quadtree domain\n#");
            for(n=0;n<ND;++n)
                printf("(%f,%f): %f " ,base[0].xmin[n],base[0].xmax[n],x[n]);
            printf("\n");
            addpt(v0);
            return;
        }
        
        for(i=0;i<(1<<ND);++i)
            if (i != exclude) srchlst[nsrch++] = topbox->dghtr[i];
    }

    return;
}

