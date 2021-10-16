#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"

#define BILLION 1000000000L

double dgemm3(int n, MATRIX *x, MATRIX *y, MATRIX *z);

struct timespec start, end;


//WARNING: THIS HAS BUGS. I WAS NOT ABLE TO FULLY COMPLETE dgemm3. n=64 matrix is the only one that works correctly.
int main()
{    
    int n[6] = {64, 128, 256, 512, 1024, 2048};
    
    char * inputA[6] = {"./testFiles/inputA-64.dat", "./testFiles/inputA-128.dat", "./testFiles/inputA-256.dat", 
                            "./testFiles/inputA-512.dat", "./testFiles/inputA-1024.dat", "./testFiles/inputA-2048.dat"};
    
    char * inputB[6] = {"./testFiles/inputB-64.dat", "./testFiles/inputB-128.dat", "./testFiles/inputB-256.dat", 
                            "./testFiles/inputB-512.dat", "./testFiles/inputB-1024.dat", "./testFiles/inputB-2048.dat"};

    char * dgemm3OutputC[6] = {"./output/dgemm3/outputC-64.dat", "./output/dgemm3/outputC-128.dat", "./output/dgemm3/outputC-256.dat", 
                            "./output/dgemm3/outputC-512.dat", "./output/dgemm3/outputC-1024.dat", "./output/dgemm3/outputC-2048.dat"};
    int i, j;
    for (i=0; i<1; i++)
    {
        double dgemm3TimeAllRuns = 0.0;
        //double numOps = pow(n[i], 3);  
        for (j=0; j<1; j++)
        {
            printf("RUN %d: \n", j);
            MATRIX *x = newMATRIXFromFile(inputA[i], n[i]);
            MATRIX *y = newMATRIXFromFile(inputB[i], n[i]);
            MATRIX *z = newEmptyMATRIX(n[i]);   

            dgemm3TimeAllRuns += dgemm3(n[i], x, y, z);

            if (j==0) // we only need to write once. wastes time to overwrite the file many times
            {
                outputExistingMATRIXToFile(z, dgemm3OutputC[i]);
            }

            //freeMATRIX(x);
            //freeMATRIX(y);
            //freeMATRIX(z);
        }
        printf("\n");
        //double dgemm3GFLOPS = 36.0 * 1.0e-9 *numOps / (dgemm3TimeAllRuns / j);
        //printf("numOps: %f\n", numOps);
        //printf("TOTAL TIME FOR dgemm2: %f\n", dgemm3TimeAllRuns);
        //printf("AVERAGE GLOPS OF dgemm2 (RUNS: %d, n: %d) = %f\n", j, n[i], dgemm3GFLOPS);
        //printf("\n\n");
        
    }

    return 0;
}

/* Multiply n x n matrices x and y  */
double dgemm3(int n, MATRIX *x, MATRIX *y, MATRIX *z)
{
    double *a = x->data;
    double *b = y->data;
    double *c2 = z->data;
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k;
    for (i = 0; i < n; i+=3)
    {
        for (j = 0; j < n; j+=3)
        {    
            register int t = i*n+j; register int tt = t+n; register int ttt = t+n+n;
            register double c00 = c2[t]; register double c01 = c2[t+1];  register double c02 = c2[t+2];
            register double c10 = c2[tt]; register double c11 = c2[tt+1]; register double c12 = c2[tt+2];
            register double c20 = c2[ttt]; register double c21 = c2[ttt+1]; register double c22 = c2[ttt+2];

            for (k = 0; k < n; k+=3)
            {
                register int ta = i*n+k; register int tta = ta+n; register int ttta = ta+n+n;
                register int tb = k*n+j; register int ttb = tb+n; 
                register double a00 = a[ta]; register double a10 = a[tta]; register double a20 = a[ttta];
                register double b00 = b[tb]; register double b01 = b[tb+1]; register double b02 = b[tb+2];

                c00 += a00*b00 ; c01 += a00*b01 ; c02 += a00*b02; 
                c10 += a10*b00 ; c11 += a10*b01 ; c12 += a10*b02;
                c20 += a20*b00 ; c21 += a20*b01 ; c22 += a20*b02;

                a00 = a[ta+1]; a10 = a[tta+1]; a20 = a[ttta+2];
                b00 = b[ttb];  b01 = b[ttb+1]; b02 = b[ttb+2];

                c00 += a00*b00 ; c01 += a00*b01 ; c02 += a00*b02; 
                c10 += a10*b00 ; c11 += a10*b01 ; c12 += a10*b02;
                c20 += a20*b00 ; c21 += a20*b01 ; c22 += a20*b02;

            }
            
            c2[t] = c00;
            c2[t+1] = c01;
            c2[t+2] = c02;
            c2[tt] = c10;
            c2[tt+1] = c11;
            c2[tt+2] = c12;
            c2[ttt] = c20;
            c2[ttt+1] = c21;
            c2[ttt+2] = c22;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("dgemm3 (dim=%d) :: Total time in seconds is: %f\n", n, time_spent);
    return time_spent;
}

