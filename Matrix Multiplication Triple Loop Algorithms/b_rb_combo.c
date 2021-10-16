/* cache blocking AND register blocking combind*/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "matrix.h"


#define BILLION 1000000000L


struct timespec start, end;

// ijk combo blocking/block rr algo
double ijkCombo (double *A, double *B, double *C1, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += Bsize)
    {
        for (j = 0; j < n; j += Bsize)
        {
            for (k = 0; k < n; k += Bsize)
            {
                for (i1 = i; i1 < i+Bsize; i1+=2)
                {
                    for (j1 = j; j1 < j+Bsize; j1+=2)
                    {
                        register int tc = i1*n + j1;
                        register int ttc = tc + n;

                        register double C00 = C1[tc];
                        register double C01 = C1[tc+1];
                        register double C10 = C1[ttc];
                        register double C11 = C1[ttc+1];

                        for (k1 = k; k1 < k+Bsize; k1+=2){

                            register int ta = i1*n + k1;
                            register int tta = ta + n;

                            register int tb = k1*n + j1;
                            register int ttb = tb + n;

                            register double A00 = A[ta];
                            register double A01 = A[ta+1];
                            register double A10 = A[tta];
                            register double A11 = A[tta+1];

                            register double B00 = B[tb];
                            register double B01 = B[tb+1];
                            register double B10 = B[ttb];
                            register double B11 = B[ttb+1];

                            C00 += A00*B00 + A01 * B10;
                            C01 += A00*B01 + A01 * B11;
                            C10 += A10*B00 + A11 * B10;
                            C11 += A10*B01 + A11 * B11;
                        }

                        C1[tc] = C00;
                        C1[tc+1] = C01;
                        C1[ttc] = C10;
                        C1[ttc+1] = C11;
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

// jik combo blocking/block rr algo
double jikCombo (double *A, double *B, double *C2, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j += Bsize)
    {
        for (i = 0; i < n; i += Bsize)
        {
            for (k = 0; k < n; k += Bsize)
            {
                for (j1 = j; j1 < j+Bsize; j1+=2)
                {
                    for (i1 = i; i1 < i+Bsize; i1+=2)
                    {
                        register int tc = i1*n + j1;
                        register int ttc = tc + n;

                        register double C00 = C2[tc];
                        register double C01 = C2[tc+1];
                        register double C10 = C2[ttc];
                        register double C11 = C2[ttc+1];

                        for (k1 = k; k1 < k+Bsize; k1+=2){

                            register int ta = i1*n + k1;
                            register int tta = ta + n;

                            register int tb = k1*n + j1;
                            register int ttb = tb + n;

                            register double A00 = A[ta];
                            register double A01 = A[ta+1];
                            register double A10 = A[tta];
                            register double A11 = A[tta+1];

                            register double B00 = B[tb];
                            register double B01 = B[tb+1];
                            register double B10 = B[ttb];
                            register double B11 = B[ttb+1];

                            C00 += A00*B00 + A01 * B10;
                            C01 += A00*B01 + A01 * B11;
                            C10 += A10*B00 + A11 * B10;
                            C11 += A10*B01 + A11 * B11;
                        }

                        C2[tc] = C00;
                        C2[tc+1] = C01;
                        C2[ttc] = C10;
                        C2[ttc+1] = C11;
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}


// ikj combo blocking/block rr algo
double ikjCombo (double *A, double *B, double *C3, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i += Bsize)
    {
        for (k = 0; k < n; k += Bsize)
        {
            for (j = 0; j < n; j += Bsize)
            {
                for (i1 = i; i1 < i+Bsize; i1+=2)
                {
                    for (k1 = k; k1 < k+Bsize; k1+=2)
                    {
                        register int ta = i1*n + k1;
                        register int tta = ta + n;

                        register double A00 = A[ta];
                        register double A01 = A[ta+1];
                        register double A10 = A[tta];
                        register double A11 = A[tta+1];

                        for (j1 = j; j1 < j+Bsize; j1+=2){

                            register int tc = i1*n + j1;
                            register int ttc = tc + n;

                            register int tb = k1*n + j1;
                            register int ttb = tb + n;

                            register double B00 = B[tb];
                            register double B01 = B[tb+1];
                            register double B10 = B[ttb];
                            register double B11 = B[ttb+1];

                            C3[tc] += A00*B00 + A01 * B10;
                            C3[tc+1] += A00*B01 + A01 * B11;
                            C3[ttc] += A10*B00 + A11 * B10;
                            C3[ttc+1] += A10*B01 + A11 * B11;

                        }
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

// kij  combo blocking/block rr algo
double kijCombo (double *A, double *B, double *C4, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += Bsize)
    {
        for (i = 0; i < n; i += Bsize)
        {
            for (j = 0; j < n; j += Bsize)
            {
                for (k1 = k; k1 < k+Bsize; k1+=2)
                {
                    for (i1 = i; i1 < i+Bsize; i1+=2)
                    {
                        register int ta = i1*n + k1;
                        register int tta = ta + n;

                        register double A00 = A[ta];
                        register double A01 = A[ta+1];
                        register double A10 = A[tta];
                        register double A11 = A[tta+1];

                        for (j1 = j; j1 < j+Bsize; j1+=2) {

                            register int tc = i1*n + j1;
                            register int ttc = tc + n;

                            register int tb = k1*n + j1;
                            register int ttb = tb + n;

                            register double B00 = B[tb];
                            register double B01 = B[tb+1];
                            register double B10 = B[ttb];
                            register double B11 = B[ttb+1];

                            C4[tc] += A00*B00 + A01 * B10;
                            C4[tc+1] += A00*B01 + A01 * B11;
                            C4[ttc] += A10*B00 + A11 * B10;
                            C4[ttc+1] += A10*B01 + A11 * B11;
                        }

                    }
                }
            }
        }
    }    
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}



// jki combo blocking/block rr algo
double jkiCombo (double *A, double *B, double *C5, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j += Bsize)
    {
        for (k = 0; k < n; k += Bsize)
        {
            for (i = 0; i < n; i += Bsize)
            {
                for (j1 = j; j1 < j+Bsize; j1+=2)
                {
                    for (k1 = k; k1 < k+Bsize; k1+=2)
                    {

                        register int tb = k1*n + j1;
                        register int ttb = tb + n;

                        register double B00 = B[tb];
                        register double B01 = B[tb+1];
                        register double B10 = B[ttb];
                        register double B11 = B[ttb+1];

                        for (i1 = i; i1 < i+Bsize; i1+=2){

                            register int tc = i1*n + j1;
                            register int ttc = tc + n;

                            register int ta = i1*n + k1;
                            register int tta = ta + n;

                            register double A00 = A[ta];
                            register double A01 = A[ta+1];
                            register double A10 = A[tta];
                            register double A11 = A[tta+1];

                            C5[tc] += A00*B00 + A01 * B10;
                            C5[tc+1] += A00*B01 + A01 * B11;
                            C5[ttc] += A10*B00 + A11 * B10;
                            C5[ttc+1] += A10*B01 + A11 * B11;
                        }
                    }
                }
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_spent = ((end.tv_sec - start.tv_sec)*BILLION + end.tv_nsec - start.tv_nsec) * 1e-9;
    return time_spent;
}

// kji  combo blocking/block rr algo
double kjiCombo (double *A, double *B, double *C6, int n, int Bsize){
    double time_spent = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &start);
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k += Bsize)
    {
        for (j = 0; j < n; j += Bsize)
        {
            for (i = 0; i < n; i += Bsize)
            {
                for (k1 = k; k1 < k+Bsize; k1+=2)
                {
                    for (j1 = j; j1 < j+Bsize; j1+=2)
                    {
                        register int tb = k1*n+ j1;
                        register int ttb = tb + n;

                        register double B00 = B[tb];
                        register double B01 = B[tb+1];
                        register double B10 = B[ttb];
                        register double B11 = B[ttb+1];

                        for (i1 = i; i1 < i+Bsize; i1+=2){
                            register int tc = i1*n + j1;
                            register int ttc = tc + n;

                            register int ta = i1*n + k1;
                            register int tta = ta + n;

                            register double A00 = A[ta];
                            register double A01 = A[ta+1];
                            register double A10 = A[tta];
                            register double A11 = A[tta+1];

                            C6[tc] += A00*B00 + A01 * B10;
                            C6[tc+1] += A00*B01 + A01 * B11;
                            C6[ttc] += A10*B00 + A11 * B10;
                            C6[ttc+1] += A10*B01 + A11 * B11;
                        }
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
    printf("TOTAL TIME FOR combo %s: %f\n", algo, time);
    printf("GFLOPS OF %s (n: %d, B:%d) = %f\n", algo, n, B, GFLOPS);
}

int main ()
{
    int n = 2048;
    int B[6] = {16, 32, 64, 128, 256, 512};

    char * comboOutputC[6] = {"./output/combo/outputC-ijk.dat", "./output/combo/outputC-jik.dat", "./output/combo/outputC-ikj.dat", 
                            "./output/combo/outputC-kij.dat", "./output/combo/outputC-jki.dat", "./output/combo/outputC-kji.dat"};
  

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
        double ijkTime = ijkCombo(x->data, y->data, z1->data, n, B[i]);
        outputExistingMATRIXToFile(z1, comboOutputC[0]);
        displayRunTimes(ijkTime, n, B[i], "ijk");

        double jikTime =jikCombo(x->data, y->data, z2->data, n, B[i]);
        outputExistingMATRIXToFile(z2, comboOutputC[1]);
        displayRunTimes(jikTime, n, B[i], "jik");

        double ikjTime =ikjCombo(x->data, y->data, z3->data, n, B[i]);
        outputExistingMATRIXToFile(z3, comboOutputC[2]);
        displayRunTimes(ikjTime, n, B[i], "ikj");

        double kijTime =kijCombo(x->data, y->data, z4->data, n, B[i]);
        outputExistingMATRIXToFile(z4, comboOutputC[3]);
        displayRunTimes(kijTime, n, B[i], "kij");

        double jkiTime =jkiCombo(x->data, y->data, z5->data, n, B[i]);
        outputExistingMATRIXToFile(z5, comboOutputC[4]);
        displayRunTimes(jkiTime, n, B[i], "jki");

        double kjiTime =kjiCombo(x->data, y->data, z6->data, n, B[i]);
        outputExistingMATRIXToFile(z6, comboOutputC[5]);
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

