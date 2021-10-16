# Makefile package and compute
# --- macros
OOPTS= -Wall -Wextra
CFLAGS= -lm -O0
OBJ= matrix.o
JOB0= 0
JOB1= 0
JOB2= 0

# --- targets
all: dgemm0 dgemm1 dgemm2 dgemm3

dgemm0: $(OBJ) dgemm0.o
	gcc $(OOPTS) $(OBJ) dgemm0.o -o dgemm0 $(CFLAGS)

dgemm1: $(OBJ) dgemm1.o
	gcc $(OOPTS) $(OBJ) dgemm1.o -o dgemm1 $(CFLAGS)

dgemm2: $(OBJ) dgemm2.o
	gcc $(OOPTS) $(OBJ) dgemm2.o -o dgemm2 $(CFLAGS)

dgemm3: $(OBJ) dgemm3.o
	gcc $(OOPTS) $(OBJ) dgemm3.o -o dgemm3 $(CFLAGS)

dgemm0.o: dgemm0.c matrix.h
	gcc $(OOPTS) -c dgemm0.c $(CFLAGS)

dgemm1.0: dgemm1.c matrix.h
	gcc $(OOPTS) -c dgemm1.c $(CFLAGS)

dgemm2.o: dgemm2.c matrix.h
	gcc $(OOPTS) -c dgemm2.c $(CFLAGS)

dgemm3.o: dgemm3.c matrix.h
	gcc $(OOPTS) -c dgemm3.c $(CFLAGS)

matrix.o: matrix.c matrix.h
	gcc $(OOPTS) -c matrix.c $(CFLAGS)

# --- remove binary and executable files
clean:
	rm -f $(OBJ) dgemm0.o dgemm1.o dgemm2.o dgemm3.o dgemm0 dgemm1 dgemm2 dgemm3

# --- submits dgemm-0.job dgemm-1.job dgemm-2.job. Each dimension has 10 iterations to average times.
submit:
	sbatch dgemm-0.job
	sbatch dgemm-1.job
	sbatch dgemm-2.job
	sbatch dgemm-3.job

# --- compare max differences of matrices. Must specify 
compareDifferences:
	grep "Max difference is" ./output/console/dgemm0/dgemm0.$(JOB0).out > temp1
	grep "Max difference is" ./output/console/dgemm1/dgemm1.$(JOB1).out > temp2
	grep "Max difference is" ./output/console/dgemm2/dgemm2.$(JOB2).out > temp3
	diff temp1 temp2
	diff temp2 temp3
	rm -f temp1 temp2 temp3
# --- compare all output files. 
compareFiles:
	diff ./output/dgemm0/outputC-64.dat ./output/dgemm1/outputC-64.dat
	diff ./output/dgemm0/outputC-128.dat ./output/dgemm1/outputC-128.dat
	diff ./output/dgemm0/outputC-256.dat ./output/dgemm1/outputC-256.dat
	diff ./output/dgemm0/outputC-512.dat ./output/dgemm1/outputC-512.dat
	diff ./output/dgemm0/outputC-1024.dat ./output/dgemm1/outputC-1024.dat
	diff ./output/dgemm0/outputC-2048.dat ./output/dgemm1/outputC-2048.dat
	diff ./output/dgemm2/outputC-64.dat ./output/dgemm1/outputC-64.dat
	diff ./output/dgemm2/outputC-128.dat ./output/dgemm1/outputC-128.dat
	diff ./output/dgemm2/outputC-256.dat ./output/dgemm1/outputC-256.dat
	diff ./output/dgemm2/outputC-512.dat ./output/dgemm1/outputC-512.dat
	diff ./output/dgemm2/outputC-1024.dat ./output/dgemm1/outputC-1024.dat
	diff ./output/dgemm2/outputC-2048.dat ./output/dgemm1/outputC-2048.dat