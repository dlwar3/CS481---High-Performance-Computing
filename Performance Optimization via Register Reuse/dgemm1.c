#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"

#define BILLION 1000000000L

double dgemm1(int n, MATRIX *x, MATRIX *y, MATRIX *z);

struct timespec start, end;

int main()
{    
    int n[6] = {64, 128, 256, 512, 1024, 2048};
    
    char * inputA[6] = {"./testFiles/inputA-64.dat", "./testFiles/inputA-128.dat", "./testFiles/inputA-256.dat", 
                            "./testFiles/inputA-512.dat", "./testFiles/inputA-1024.dat", "./testFiles/inputA-2048.dat"};
    
    char * inputB[6] = {"./testFiles/inputB-64.dat", "./testFiles/inputB-128.dat", "./testFiles/inputB-256.dat", 
                            "./testFiles/inputB-512.dat", "./testFiles/inputB-1024.dat", "./testFiles/inputB-2048.dat"};

    char * dgemm1OutputC[6] = {"./output/dgemm1/outputC-64.dat", "./output/dgemm1/outputC-128.dat", "./output/dgemm1/outputC-256.dat", 
                            "./output/dgemm1/outputC-512.dat", "./output/dgemm1/outputC-1024.dat", "./output/dgemm1/outputC-2048.dat"};

    int i, j;
    for (i=0; i<6; i++)
    {
        double dgemm1TimeAllRuns = 0.0;
        double numOps = pow(n[i], 3);  
        for (j=0; j<10; j++)
        {
            printf("RUN %d: \n", j);
            MATRIX *x = newMATRIXFromFile(inputA[i], n[i]);
            MATRIX *y = newMATRIXFromFile(inputB[i], n[i]);
            MATRIX *z = newEmptyMATRIX(n[i]);   

            dgemm1TimeAllRuns += dgemm1(n[i], x, y, z);

            if (j==0) // we only need to write once. wastes time to overwrite the file many times
            {
                outputExistingMATRIXToFile(z, dgemm1OutputC[i]);
                printf("Max difference is %f\n", findMaxDifference(x, y));
            }

            freeMATRIX(x);
            freeMATRIX(y);
            freeMATRIX(z);
        }

    
        printf("\n");
        double dgemm1GFLOPS = 2.0 * 1.0e-9 * numOps / (dgemm1TimeAllRuns / j);
        printf("numOps: %f\n", numOps);
        printf("TOTAL TIME FOR dgemm1: %f\n", dgemm1TimeAllRuns);
        printf("AVERAGE GLOPS OF dgemm1 (RUNS: %d, n: %d) = %f\n", j, n[i], dgemm1GFLOPS);
        printf("\n\n");
    }
    return 0;
}

/*dgemm1: simple ijk version triple loop algorithm with register reuse*/
double dgemm1(int n, MATRIX *x, MATRIX *y, MATRIX *z)
{
    double *a = x->data;
    double *b = y->data;
    double *c = z->data;
    int i,j,k;
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++) {
            register double r = c[i*n+j];
            for (k=0; k<n; k++)
            {
                r += a[i*n+k] * b[k*n+j];
            }
            c[i*n+j] = r;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("dgemm1 (dim=%d) :: Total time in seconds is: %f\n", n, time_spent);
    return time_spent;
}
