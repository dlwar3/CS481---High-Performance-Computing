# Makefile package and compute
# --- macros
CFLAGS= -lm -O0 -fopenmp

# --- targets
all: parallelMM_outer parallelMM_middle 3dheat 3dheatParallel

parallelMM_outer: $(OBJ) ./part1/parallelMM_outer.o
	gcc ./part1/parallelMM_outer.o -o ./part1/parallelMM_outer $(CFLAGS)

parallelMM_middle: $(OBJ) ./part1/parallelMM_middle.o
	gcc ./part1/parallelMM_middle.o -o ./part1/parallelMM_middle $(CFLAGS)

3dheat: $(OBJ) ./part2-3/3dheat.o
	gcc ./part2-3/3dheat.o -o ./part2-3/3dheat $(CFLAGS)

3dheatParallel: $(OBJ) ./part2-3/3dheatParallel.o
	gcc ./part2-3/3dheatParallel.o -o ./part2-3/3dheatParallel $(CFLAGS)

parallelMM_outer.o: ./part1/parallelMM_outer.c 
	gcc -c ./part1/parallelMM_outer.c $(CFLAGS)

parallelMM_middle.o: ./part1/parallelMM_middle.c 
	gcc -c ./part1/parallelMM_middle.c $(CFLAGS)

3dheat.o: 3dheat.c
	gcc -c ./part2-3/3dheat.c $(CFLAGS)

3dheatParallel.o: 3dheatParallel.c
	gcc -c ./part2-3/3dheatParallel.c $(CFLAGS)

# --- remove binary and executable files
clean:
	rm -f $(OBJ) ./part1/parallelMM_outer ./part1/parallelMM_middle ./part2-3/3dheat ./part2-3/3dheatParallel ./part1/parallelMM_outer.o ./part1/parallelMM_middle.o ./part2-3/3dheat.o ./part2-3/3dheatParallel.o

# --- submits outer, middle, sequential 3dheat, and parallel version of 3dheat 1,4,8,16,32 threads
submit:
	sbatch mmOuter.job
	sbatch mmMiddle.job
	sbatch 3dHeat.job
	sbatch 3dheatParallel.job