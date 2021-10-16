#ifndef matrix_h
#define matrix_h

typedef struct matrix
{
  int dim;
  float *data;
} MATRIX;

MATRIX* newEmptyMATRIX(int dim);
MATRIX* newMATRIXFromFile (char *file, int dim);
void generateNewMATRIXFile(char *file, int dim);
void populateMATRIX(MATRIX *m);
void outputExistingMATRIXToFile(MATRIX *m, char *file);
double findMaxDifference(MATRIX *a, MATRIX *b);
void freeMATRIX(MATRIX *m);

#endif