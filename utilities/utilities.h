#include<stdlib.h>

#define MAX(A,B) (A>B?A:B)
#define MIN(A,B) (A<B?A:B)
#define SIGN(A) (A > 0 ? 1 : -1)

#ifdef __cplusplus
extern "C" void *xmalloc(size_t a);	
extern "C" void number_str(char *out, const char *in, int n, int dgts);
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
