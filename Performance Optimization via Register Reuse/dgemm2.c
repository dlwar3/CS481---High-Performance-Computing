#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"

#define BILLION 1000000000L

double dgemm2(int n, MATRIX *x, MATRIX *y, MATRIX *z);

struct timespec start, end;

int main()
{    
    int n[6] = {64, 128, 256, 512, 1024, 2048};
    
    char * inputA[6] = {"./testFiles/inputA-64.dat", "./testFiles/inputA-128.dat", "./testFiles/inputA-256.dat", 
                            "./testFiles/inputA-512.dat", "./testFiles/inputA-1024.dat", "./testFiles/inputA-2048.dat"};
    
    char * inputB[6] = {"./testFiles/inputB-64.dat", "./testFiles/inputB-128.dat", "./testFiles/inputB-256.dat", 
                            "./testFiles/inputB-512.dat", "./testFiles/inputB-1024.dat", "./testFiles/inputB-2048.dat"};

    char * dgemm2OutputC[6] = {"./output/dgemm2/outputC-64.dat", "./output/dgemm2/outputC-128.dat", "./output/dgemm2/outputC-256.dat", 
                            "./output/dgemm2/outputC-512.dat", "./output/dgemm2/outputC-1024.dat", "./output/dgemm2/outputC-2048.dat"};
    int i, j;
    for (i=0; i<6; i++)
    {
        double dgemm2TimeAllRuns = 0.0;
        double numOps = pow(n[i], 3);  
        for (j=0; j<10; j++)
        {
            printf("RUN %d: \n", j);
            MATRIX *x = newMATRIXFromFile(inputA[i], n[i]);
            MATRIX *y = newMATRIXFromFile(inputB[i], n[i]);
            MATRIX *z = newEmptyMATRIX(n[i]);   

            dgemm2TimeAllRuns += dgemm2(n[i], x, y, z);

            if (j==0) // we only need to write once. wastes time to overwrite the file many times
            {
                outputExistingMATRIXToFile(z, dgemm2OutputC[i]);
                printf("Max difference is %f\n", findMaxDifference(x, y));
            }

            freeMATRIX(x);
            freeMATRIX(y);
            freeMATRIX(z);
        }
        printf("\n");
        double dgemm2GFLOPS = 2 * 1.0e-9 *numOps / (dgemm2TimeAllRuns / j);
        printf("numOps: %f\n", numOps);
        printf("TOTAL TIME FOR dgemm2: %f\n", dgemm2TimeAllRuns);
        printf("AVERAGE GLOPS OF dgemm2 (RUNS: %d, n: %d) = %f\n", j, n[i], dgemm2GFLOPS);
        printf("\n\n");
    }

    return 0;
}

/* Multiply n x n matrices x and y  */
double dgemm2(int n, MATRIX *x, MATRIX *y, MATRIX *z)
{
    double *a = x->data;
    double *b = y->data;
    double *c2 = z->data;
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k;
    for (i = 0; i < n; i+=2)
    {
        for (j = 0; j < n; j+=2)
        {    
            register int t = i*n+j; register int tt = t+n; 
            register double c00 = c2[t]; register double c01 = c2[t+1];  register double c10 = c2[tt]; register double c11 = c2[tt+1];

            for (k = 0; k < n; k+=2)
            {
                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                //only need 4 registers if we switch out after multiplying quadrants
                register double a00 = a[ta]; register double a10 = a[tta]; register double b00 = b[tb]; register double b01 = b[tb+1]; 
                

                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
                a00 = a[ta+1]; a10 = a[tta+1]; b00 = b[ttb]; b01 = b[ttb+1]; // switch out a/b 00 and a10 and b01 for next product
                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
            }
            
            c2[t] = c00;
            c2[t+1] = c01;
            c2[tt] = c10;
            c2[tt+1] = c11;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("dgemm2 (dim=%d) :: Total time in seconds is: %f\n", n, time_spent);
    return time_spent;
}
