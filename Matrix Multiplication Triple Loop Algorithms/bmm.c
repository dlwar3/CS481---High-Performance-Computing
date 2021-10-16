/* single register reuse matrix multiplication*/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "matrix.h"


#define BILLION 1000000000L


struct timespec start, end;

// ijk blocking algo
double ijkBlock (double *A, double *B, double *C1, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += Bsize)
    {
        for (j = 0; j < n; j += Bsize)
        {
            for (k = 0; k < n; k += Bsize)
            {
                for (i1 = i; i1 < i+Bsize; i1++)
                {
                    for (j1 = j; j1 < j+Bsize; j1++)
                    {
                        register double r = C1[i1*n+j1];
                        for (k1 = k; k1 < k+Bsize; k1++)
                            r += A[i1*n+k1] * B[k1*n+j1];
                        C1[i1*n+j1] = r;
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

// jik blocking algo
double jikBlock (double *A, double *B, double *C2, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j += Bsize)
    {
        for (i = 0; i < n; i += Bsize)
        {
            for (k = 0; k < n; k += Bsize)
            {
                for (j1 = j; j1 < j+Bsize; j1++)
                {
                    for (i1 = i; i1 < i+Bsize; i1++)
                    {
                        register double r = C2[i1*n+j1];
                        for (k1 = k; k1 < k+Bsize; k1++)
                            r += A[i1*n+k1] * B[k1*n+j1];
                        C2[i1*n+j1] = r;
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}


// ikj blocking algo
double ikjBlock (double *A, double *B, double *C3, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += Bsize)
    {
        for (k = 0; k < n; k += Bsize)
        {
            for (j = 0; j < n; j += Bsize)
            {
                for (i1 = i; i1 < i+Bsize; i1++)
                {
                    for (k1 = k; k1 < k+Bsize; k1++)
                    {
                        register double r = A[i1*n+k1];
                        for (j1 = j; j1 < j+Bsize; j1++)
                            C3[i1*n+j1] += r * B[k1*n+j1];
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

// kij  blocking algo
double kijBlock (double *A, double *B, double *C4, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += Bsize)
    {
        for (i = 0; i < n; i += Bsize)
        {
            for (j = 0; j < n; j += Bsize)
            {
                for (k1 = k; k1 < k+Bsize; k1++)
                {
                    for (i1 = i; i1 < i+Bsize; i1++)
                    {
                        register double r = A[i1*n+k1];
                        for (j1 = j; j1 < j+Bsize; j1++)
                            C4[i1*n+j1] += r * B[k1*n+j1];
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}



// jki blocking algo
double jkiBlock (double *A, double *B, double *C5, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j += Bsize)
    {
        for (k = 0; k < n; k += Bsize)
        {
            for (i = 0; i < n; i += Bsize)
            {
                for (j1 = j; j1 < j+Bsize; j1++)
                {
                    for (k1 = k; k1 < k+Bsize; k1++)
                    {
                        register double r = B[k1*n+j1];
                        for (i1 = i; i1 < i+Bsize; i1++)
                            C5[i1*n+j1] += A[i1*n+k1] * r;
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

// kji  blocking algo
double kjiBlock (double *A, double *B, double *C6, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += Bsize)
    {
        for (j = 0; j < n; j += Bsize)
        {
            for (i = 0; i < n; i += Bsize)
            {
                for (k1 = k; k1 < k+Bsize; k1++)
                {
                    for (j1 = j; j1 < j+Bsize; j1++)
                    {
                        register double r = B[k1*n+j1];
                        for (i1 = i; i1 < i+Bsize; i1++)
                            C6[i1*n+j1] += A[i1*n+k1] * r;
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

void displayRunTimes(double time, int n, int B, char* algo)
{
    double numOps = 2 * pow(n,3);
    double GFLOPS = 2.0 * 1.0e-9 * numOps / time;
    printf("numOps: %f\n", numOps);
    printf("TOTAL TIME FOR bmm %s: %f\n", algo, time);
    printf("GFLOPS OF %s (n: %d, B:%d) = %f\n", algo, n, B, GFLOPS);
}

int main ()
{
    int n = 2048;
    int B[6] = {16, 32, 64, 128, 256, 512};

    char * bmmOutputC[6] = {"./output/bmm/outputC-ijk.dat", "./output/bmm/outputC-jik.dat", "./output/bmm/outputC-ikj.dat", 
                            "./output/bmm/outputC-kij.dat", "./output/bmm/outputC-jki.dat", "./output/bmm/outputC-kji.dat"};
  

    MATRIX *x = newMATRIXFromFile("./testFiles/inputA-2048.dat", n);
    MATRIX *y = newMATRIXFromFile("./testFiles/inputB-2048.dat", n);
    MATRIX *z1 = newEmptyMATRIX(n);
    MATRIX *z2 = newEmptyMATRIX(n);
    MATRIX *z3 = newEmptyMATRIX(n);
    MATRIX *z4 = newEmptyMATRIX(n);
    MATRIX *z5 = newEmptyMATRIX(n);
    MATRIX *z6 = newEmptyMATRIX(n);

    int i;
    for (i=0; i<6; i++)
    {
        double ijkTime = ijkBlock(x->data, y->data, z1->data, n, B[i]);
        outputExistingMATRIXToFile(z1, bmmOutputC[0]);
        displayRunTimes(ijkTime, n, B[i], "ijk");

        double jikTime =jikBlock(x->data, y->data, z2->data, n, B[i]);
        outputExistingMATRIXToFile(z2, bmmOutputC[1]);
        displayRunTimes(jikTime, n, B[i], "jik");

        double ikjTime =ikjBlock(x->data, y->data, z3->data, n, B[i]);
        outputExistingMATRIXToFile(z3, bmmOutputC[2]);
        displayRunTimes(ikjTime, n, B[i], "ikj");

        double kijTime =kijBlock(x->data, y->data, z4->data, n, B[i]);
        outputExistingMATRIXToFile(z4, bmmOutputC[3]);
        displayRunTimes(kijTime, n, B[i], "kij");

        double jkiTime =jkiBlock(x->data, y->data, z5->data, n, B[i]);
        outputExistingMATRIXToFile(z5, bmmOutputC[4]);
        displayRunTimes(jkiTime, n, B[i], "jki");

        double kjiTime =kjiBlock(x->data, y->data, z6->data, n, B[i]);
        outputExistingMATRIXToFile(z6, bmmOutputC[5]);
        displayRunTimes(kjiTime, n, B[i], "kji");
    
        printf("\nMax difference compare of ijk | jik: %f\n", findMaxDifference(z1, z2));
        printf("\nMax difference compare of ijk | ikj: %f\n", findMaxDifference(z1, z3));
        printf("\nMax difference compare of ijk | kij: %f\n", findMaxDifference(z1, z4));
        printf("\nMax difference compare of ijk | jki: %f\n", findMaxDifference(z1, z5));
        printf("\nMax difference compare of ijk | kji: %f\n", findMaxDifference(z1, z6));
    }

    

    freeMATRIX(x);
    freeMATRIX(y);
    freeMATRIX(z1);
    freeMATRIX(z2);
    freeMATRIX(z3);
    freeMATRIX(z4);
    freeMATRIX(z5);
    freeMATRIX(z6);

    return 0;
}

