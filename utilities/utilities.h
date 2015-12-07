#include<stdlib.h>

#ifndef _utilities_h_
#define _utilities_h_

#ifndef SIGN
#define SIGN(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
#endif
#ifndef MIN
#define	MIN(a,b) (((a)<(b))?(a):(b))
#endif /* MIN */
#ifndef MAX
#define	MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifdef __cplusplus
extern "C" void *xmalloc(size_t a);	
extern "C" void number_str(char *out, const char *in, int n, int dgts);
extern "C" int RunningOnBigEndianMachine(void);
extern "C" void Swap2Bytes(unsigned char *data,int length);
extern "C" void Swap4Bytes(unsigned char *data,int length);
extern "C" void Swap8Bytes(unsigned char *data,int length);

#include<iostream>

class sharedmem {
   protected:
      bool corechkout;
      class core {
         public:
            int nrefs;
            bool locked;
            size_t size;
            char *p;
            
            core() : nrefs(1), locked(false), size(0)  {}
            core(size_t t) : nrefs(1), locked(false), size(0) {resize(t);}
            int resize(size_t t) {
               if (size > t) return(0);
               
               if (size) {
                  delete []p;
               }
               size = t;
               p = new char[t];
               return(0);
            }
            ~core() {
               if (locked) std::cerr << "locked memory is being freed" << std::endl;
               if (size) free(p);
            }
      }  *mycore;
      
   public:
      sharedmem() : corechkout(false) { mycore = new core; }
      sharedmem(size_t size) : corechkout(false) { mycore = new core(size);}
      void resize(size_t size) {
         if (corechkout || !mycore->locked) {
            mycore->resize(size);
            return;
         }
         std::cerr << "Can not change resize: core is locked" << std::endl;
      }
      int reference(const sharedmem &sin) {
         if (mycore == sin.mycore) return(0);
         --mycore->nrefs;
         if (!(mycore->nrefs)) delete mycore;
         mycore = sin.mycore;
         ++mycore->nrefs;
         return(0);
      }
      void checkout() {
         if (mycore->locked) {
            std::cerr << "memory already in use" << std::endl;
            return;
         }
         corechkout = true;
         mycore->locked=true;
         return;
      }
      void checkin() {corechkout = false; mycore->locked=false;}
      char * data() const {
         return(mycore->p);
      }
      size_t size() const {return(mycore->size);}
      bool is_mine() const {return(corechkout);}

      ~sharedmem() {
         --mycore->nrefs;
         if (!mycore->nrefs) delete mycore;
      }
};
#else
extern void *xmalloc(size_t a);	
extern void number_str(char *out, const char *in, int n, int dgts);
#endif

#define vect_alloc(a,b,type) \
{ /* printf("# " #a " %d bytes\n",(b)*sizeof(type)); */ \
a = (type *) malloc((b)*sizeof(type)); \
if (a == NULL) { printf(#a "couldn't allocate"); \
exit(1); } }

#define mat_alloc(a,i,j,type) {\
int k; \
/* printf("# " #a " %d bytes\n",(i)*sizeof(type *) + (i)*(j)*sizeof(type)); */ \
a=(type **) malloc((i)*sizeof(type *)); \
if (!a) { printf(#a "couldn't allocate"); \
exit(1); } \
a[0]=(type *) malloc((i)*(j)*sizeof(type)); \
if (!a[0]) { printf(#a "couldn't allocate"); \
exit(1); } \
for(k=1;k<i;++k) a[k]=a[k-1] +j; }

#define tens_alloc(a,i,j,k,type) {\
int l; \
/* printf("# " #a " %d bytes\n",(i)*sizeof(type **) + (i)*(j)*sizeof(type *) + \
(i)*(j)*(k)*sizeof(type)); */ \
a=(type ***) malloc((i)*sizeof(type **)); \
if (!a) { printf(#a "couldn't allocate"); \
exit(1); } \
a[0]=(type **) malloc((i)*(j)*sizeof(type *)); \
if (!a[0]) { printf(#a "couldn't allocate"); \
exit(1); } \
a[0][0] =(type *) malloc((i)*(j)*(k)*sizeof(type)); \
if (!a[0][0]) { printf(#a "couldn't allocate"); \
exit(1); } \
for(l=1;l<i;++l) {a[l] = a[l-1] + j;} \
for(l=1;l<(i)*(j);++l) {a[0][l] = a[0][l-1] + k; }\
}

#endif
