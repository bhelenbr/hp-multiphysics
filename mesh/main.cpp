#include<cstdio>
#include"blocks.h"
#include<iostream>

/* DEFINES/THINGS TO TIDDLE WITH */
/* r_mesh.h: FOURTH: SWITCH FROM LAPLACIAN TO BIHARMONIC */
/* blocks.cpp: GEOMETRIC: SWITCH FOR DETERMING K_IJ IN MULTIGRID (USE GEOMETRIC)*/
/* mesh.h: USEOLDBTYPE: remaps boundary identifies from my old types to new */
/* r_mesh.h: FIX?_MASK: DETERMINES DIRICHLET MESH MOVMENT BC'S */
/* r_mesh.h: FIX2?_MASK: DETERMINES DIRICHLET B.C.'S on Nabla^2 for Biharmonic */
/* mesh.h: ?DIR_MP: mask for communication boundaries */
/* r_mesh.cpp: different perturb functions for deformed mesh */
/* r_mesh.cpp: vnn/fadd - jacobi/multigrid parameters */

#include<cstring>
#include<math.h>
#include<utilities.h>
#include<time.h>

FLT center;
FLT amp;

int main(int argc, char *argv[]) {
   int i,step = 0;
   clock_t cpu_time;
   class blocks z;
   char outname[100];
   char *inname[2];
   char *out2[2];
   mesh x[10];
   

   /*
   x[0].in_mesh("/Volumes/work/helenbrk/Codes/spectral_hp/build/mesh134.0",grid);
   for (i=0;i<5;++i) {
      number_str(outname, "test", i, 1);
      x[i].setbcinfo();
      x[i].out_mesh(outname);
      x[i+1].coarsen(1.6,x[i]);
   }
   */ 
 
   /* START CLOCK TIMER */
   clock();
   
   /* THIS DEFORMS A MESH */
   /* CANONICAL TEST PROBLEM */
   inname[0] = "/Volumes/work/helenbrk/Codes/grids/BOAT/boat2";
   
   z.init(1,3,inname,easymesh,10.0);
   
   z.out_mesh("begin",tecplot);

   for(step = 1; step<=1;++step) {
      center += step/FLT(10);  // FOR THE MOVING CYLINDER PROBLEM
      amp = 0.25*step/10.;  // FOR THE DEFORMING SURFACE 
      z.ksrc();
      z.perturb();
      

      for(i=0;i<200;++i) {
         z.cycle(2);
         printf("%d ",i);
         z.maxres();
         printf("\n");
      }
      
      z.out_mesh("deformed",tecplot);
      
      z.restructure(0.66);
      number_str(outname, "test", step, 2);
      z.out_mesh(outname);
   }
   
   out2[0] = "bamp0.375";
   out2[1] = "tamp0.375";
   z.out_mesh(out2,easymesh);
   
   cpu_time = clock();
   printf("that took %ld cpu time\n",cpu_time);

   return(0);
   
   /* OUTPUT FROM TEST PROBLEM */
/*
#Boundaries
#0: type 8, number 50
#1: type 4, number 30
#2: type 1, number 42
#3: type 131080, number 50
#4: type 32, number 42
#5: type 65544, number 30
#
#
#COARSE MESH
#MAXVST 25416 VERTICES 764 SIDES 2173 ELEMENTS 1408
#MAX 1342 BDRY 0 TYPE 8 SIDES 25
#MAX 1342 BDRY 1 TYPE 4 SIDES 15
#MAX 1342 BDRY 2 TYPE 1 SIDES 21
#MAX 1342 BDRY 3 TYPE 131080 SIDES 25
#MAX 1342 BDRY 4 TYPE 32 SIDES 21
#MAX 1342 BDRY 5 TYPE 65544 SIDES 15
#
#
#COARSE MESH
#MAXVST 7261 VERTICES 221 SIDES 606 ELEMENTS 384
#MAX 671 BDRY 0 TYPE 8 SIDES 12
#MAX 671 BDRY 1 TYPE 4 SIDES 8
#MAX 671 BDRY 2 TYPE 1 SIDES 10
#MAX 671 BDRY 3 TYPE 131080 SIDES 12
#MAX 671 BDRY 4 TYPE 32 SIDES 10
#MAX 671 BDRY 5 TYPE 65544 SIDES 8
0 5.882e-01  2.566e-01  
1 1.564e-01  5.025e-02  
2 8.772e-02  2.808e-02  
3 6.307e-02  1.942e-02  
4 4.767e-02  1.285e-02  
5 3.718e-02  1.002e-02  
6 2.989e-02  8.135e-03  
7 2.424e-02  6.454e-03  
8 2.078e-02  5.061e-03  
9 1.824e-02  3.941e-03  
10 1.676e-02  3.372e-03  
11 1.558e-02  2.994e-03  
12 1.450e-02  2.686e-03  
13 1.352e-02  2.430e-03  
14 1.264e-02  2.214e-03  
15 1.184e-02  2.033e-03  
16 1.113e-02  1.928e-03  
17 1.055e-02  1.840e-03  
18 1.006e-02  1.759e-03  
19 9.618e-03  1.686e-03  
20 9.205e-03  1.618e-03  
21 8.820e-03  1.556e-03  
22 8.460e-03  1.499e-03  
23 8.123e-03  1.450e-03  
24 7.839e-03  1.407e-03  
25 7.583e-03  1.367e-03  
26 7.337e-03  1.330e-03  
27 7.101e-03  1.294e-03  
28 6.874e-03  1.261e-03  
29 6.658e-03  1.229e-03  
30 6.476e-03  1.199e-03  
31 6.299e-03  1.170e-03  
32 6.128e-03  1.146e-03  
33 5.961e-03  1.123e-03  
34 5.800e-03  1.101e-03  
35 5.643e-03  1.080e-03  
36 5.491e-03  1.059e-03  
37 5.343e-03  1.039e-03  
38 5.199e-03  1.020e-03  
39 5.059e-03  1.001e-03  
40 4.924e-03  9.848e-04  
41 4.793e-03  9.734e-04  
42 4.677e-03  9.623e-04  
43 4.580e-03  9.513e-04  
44 4.493e-03  9.405e-04  
45 4.407e-03  9.299e-04  
46 4.323e-03  9.194e-04  
47 4.240e-03  9.091e-04  
48 4.159e-03  8.989e-04  
49 4.080e-03  8.889e-04  
50 4.008e-03  8.790e-04  
51 3.936e-03  8.706e-04  
52 3.866e-03  8.624e-04  
53 3.797e-03  8.541e-04  
54 3.729e-03  8.460e-04  
55 3.663e-03  8.379e-04  
56 3.597e-03  8.299e-04  
57 3.533e-03  8.220e-04  
58 3.469e-03  8.141e-04  
59 3.407e-03  8.064e-04  
60 3.346e-03  7.986e-04  
61 3.286e-03  7.910e-04  
62 3.227e-03  7.835e-04  
63 3.170e-03  7.760e-04  
64 3.116e-03  7.686e-04  
65 3.066e-03  7.612e-04  
66 3.016e-03  7.540e-04  
67 2.966e-03  7.468e-04  
68 2.918e-03  7.397e-04  
69 2.871e-03  7.326e-04  
70 2.824e-03  7.257e-04  
71 2.779e-03  7.195e-04  
72 2.734e-03  7.136e-04  
73 2.690e-03  7.077e-04  
74 2.652e-03  7.018e-04  
75 2.616e-03  6.960e-04  
76 2.581e-03  6.903e-04  
77 2.547e-03  6.846e-04  
78 2.512e-03  6.789e-04  
79 2.479e-03  6.733e-04  
80 2.445e-03  6.677e-04  
81 2.412e-03  6.621e-04  
82 2.381e-03  6.566e-04  
83 2.351e-03  6.512e-04  
84 2.321e-03  6.458e-04  
85 2.292e-03  6.404e-04  
86 2.263e-03  6.351e-04  
87 2.234e-03  6.298e-04  
88 2.206e-03  6.246e-04  
89 2.177e-03  6.194e-04  
90 2.150e-03  6.142e-04  
91 2.122e-03  6.091e-04  
92 2.095e-03  6.041e-04  
93 2.068e-03  5.991e-04  
94 2.041e-03  5.941e-04  
95 2.015e-03  5.891e-04  
96 1.989e-03  5.843e-04  
97 1.963e-03  5.794e-04  
98 1.940e-03  5.746e-04  
99 1.917e-03  5.699e-04  
100 1.894e-03  5.651e-04  
101 1.874e-03  5.605e-04  
102 1.854e-03  5.558e-04  
103 1.835e-03  5.515e-04  
104 1.816e-03  5.474e-04  
105 1.797e-03  5.433e-04  
106 1.779e-03  5.392e-04  
107 1.760e-03  5.352e-04  
108 1.742e-03  5.312e-04  
109 1.725e-03  5.272e-04  
110 1.707e-03  5.232e-04  
111 1.689e-03  5.192e-04  
112 1.672e-03  5.153e-04  
113 1.656e-03  5.114e-04  
114 1.640e-03  5.076e-04  
115 1.624e-03  5.037e-04  
116 1.609e-03  4.999e-04  
117 1.593e-03  4.963e-04  
118 1.578e-03  4.928e-04  
119 1.564e-03  4.893e-04  
120 1.549e-03  4.859e-04  
121 1.534e-03  4.825e-04  
122 1.520e-03  4.791e-04  
123 1.507e-03  4.757e-04  
124 1.494e-03  4.724e-04  
125 1.481e-03  4.690e-04  
126 1.468e-03  4.657e-04  
127 1.456e-03  4.624e-04  
128 1.443e-03  4.591e-04  
129 1.431e-03  4.559e-04  
130 1.418e-03  4.527e-04  
131 1.406e-03  4.495e-04  
132 1.394e-03  4.465e-04  
133 1.382e-03  4.435e-04  
134 1.370e-03  4.406e-04  
135 1.358e-03  4.376e-04  
136 1.347e-03  4.347e-04  
137 1.335e-03  4.318e-04  
138 1.324e-03  4.289e-04  
139 1.312e-03  4.260e-04  
140 1.301e-03  4.231e-04  
141 1.290e-03  4.203e-04  
142 1.279e-03  4.174e-04  
143 1.268e-03  4.146e-04  
144 1.257e-03  4.118e-04  
145 1.246e-03  4.090e-04  
146 1.235e-03  4.062e-04  
147 1.225e-03  4.034e-04  
148 1.214e-03  4.007e-04  
149 1.204e-03  3.980e-04  
150 1.193e-03  3.952e-04  
151 1.183e-03  3.925e-04  
152 1.173e-03  3.898e-04  
153 1.163e-03  3.872e-04  
154 1.153e-03  3.845e-04  
155 1.143e-03  3.818e-04  
156 1.133e-03  3.792e-04  
157 1.123e-03  3.766e-04  
158 1.113e-03  3.740e-04  
159 1.104e-03  3.714e-04  
160 1.094e-03  3.688e-04  
161 1.085e-03  3.663e-04  
162 1.076e-03  3.637e-04  
163 1.066e-03  3.612e-04  
164 1.057e-03  3.587e-04  
165 1.048e-03  3.562e-04  
166 1.039e-03  3.537e-04  
167 1.030e-03  3.512e-04  
168 1.021e-03  3.488e-04  
169 1.012e-03  3.463e-04  
170 1.004e-03  3.439e-04  
171 9.951e-04  3.415e-04  
172 9.865e-04  3.391e-04  
173 9.781e-04  3.367e-04  
174 9.696e-04  3.344e-04  
175 9.613e-04  3.320e-04  
176 9.530e-04  3.297e-04  
177 9.449e-04  3.274e-04  
178 9.367e-04  3.251e-04  
179 9.287e-04  3.228e-04  
180 9.207e-04  3.205e-04  
181 9.128e-04  3.182e-04  
182 9.050e-04  3.160e-04  
183 8.972e-04  3.138e-04  
184 8.895e-04  3.116e-04  
185 8.819e-04  3.093e-04  
186 8.743e-04  3.072e-04  
187 8.668e-04  3.050e-04  
188 8.594e-04  3.028e-04  
189 8.520e-04  3.007e-04  
190 8.447e-04  2.985e-04  
191 8.374e-04  2.964e-04  
192 8.303e-04  2.943e-04  
193 8.232e-04  2.922e-04  
194 8.161e-04  2.902e-04  
195 8.091e-04  2.881e-04  
196 8.022e-04  2.860e-04  
197 7.953e-04  2.840e-04  
198 7.885e-04  2.820e-04  
199 7.818e-04  2.800e-04  
#Swap cycle finished: 938 sides swapped
#Swap cycle finished: 65 sides swapped
#Swap cycle finished: 0 sides swapped
#Yaber finished: 360 sides coarsened
#Rebay finished: new interior points 290, new boundary points 0
#
#
#COARSE MESH
#MAXVST 25416 VERTICES 782 SIDES 2227 ELEMENTS 1444
#MAX 1342 BDRY 0 TYPE 8 SIDES 25
#MAX 1342 BDRY 1 TYPE 4 SIDES 15
#MAX 1342 BDRY 2 TYPE 1 SIDES 21
#MAX 1342 BDRY 3 TYPE 131080 SIDES 25
#MAX 1342 BDRY 4 TYPE 32 SIDES 21
#MAX 1342 BDRY 5 TYPE 65544 SIDES 15
#
#
#COARSE MESH
#MAXVST 7261 VERTICES 251 SIDES 696 ELEMENTS 444
#MAX 671 BDRY 0 TYPE 8 SIDES 12
#MAX 671 BDRY 1 TYPE 4 SIDES 8
#MAX 671 BDRY 2 TYPE 1 SIDES 10
#MAX 671 BDRY 3 TYPE 131080 SIDES 12
#MAX 671 BDRY 4 TYPE 32 SIDES 10
#MAX 671 BDRY 5 TYPE 65544 SIDES 8
that took 692 cpu time

mblock has exited with status 0.
*/

}
