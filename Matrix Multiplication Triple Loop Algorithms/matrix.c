  #include <stdio.h>
  #include <stdlib.h>
  #include <assert.h>
  #include <math.h>
  #include <time.h>

  typedef struct matrix
  {
    int dim;
    double *data;
  } MATRIX;

  MATRIX* newEmptyMATRIX(int dim)
  {
      MATRIX *m = malloc(sizeof(MATRIX));
      assert(sizeof(m) > 0);
      m->dim = dim;
      m->data = (double*)calloc(sizeof(double), m->dim * m->dim);
      assert(sizeof(m->data) > 0);  //make sure we don't have empty data
      return m;
  }

  MATRIX* newMATRIXFromFile(char *file, int dim)
  {
    MATRIX *m = malloc(sizeof(MATRIX));
    FILE *fp = fopen(file, "r");
    assert(fp != NULL);
    m->dim = dim;
    //allocate matrix
    m->data = (double*)calloc(sizeof(double), m->dim * m->dim);
    assert(sizeof(m->data) > 0);  //make sure we don't have empty data

    //populate data array
    int i;
    for (i=0; i< dim * dim; i++)
    {
      fscanf(fp, "%lf", &m->data[i]);
    }
    fclose(fp);
    return m;
  }

  void generateNewMATRIXFile(char *file, int dim)
  {
      srand(time(NULL));
      FILE *fp = fopen(file, "w");
      assert(fp != NULL);
      int i;
      for (i = 0; i< dim * dim; i++)
      {
          fprintf(fp, "%f ", (double)rand());
      }
      fclose(fp);
  }

  void populateMATRIX(MATRIX *m)
  {
      srand(time(NULL));
      int i;
      for (i=0; i<m->dim * m->dim; i++)
      {
          m->data[i] = (double)rand();
      }

  }
  void outputExistingMATRIXToFile(MATRIX *m, char *file)
  {
      FILE *fp = fopen(file, "w");
      assert(fp != NULL);
      int i;
      for (i=0; i<m->dim * m->dim; i++)
      {
          fprintf(fp, "%f ", m->data[i]);
      }
      fclose(fp);
  }

  double findMaxDifference(MATRIX *a, MATRIX *b)
  {
    double maxDifference = 0.0;
    double temp = 0.0;
    int i;
    for (i=0; i<a->dim * a->dim; i++)
    {
      temp = fabs(a->data[i] - b->data[i]);
      if (temp > maxDifference)
      {
        maxDifference = temp;
      }
    }
    return maxDifference;
  }

  void freeMATRIX(MATRIX *m)
  {
    free(m->data);
    free(m);
  }

