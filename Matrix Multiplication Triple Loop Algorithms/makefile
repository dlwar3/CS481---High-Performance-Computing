# Makefile package and compute
# --- macros
OOPTS= -Wall -Wextra
CFLAGS= -lm -O0
OBJ= matrix.o

# --- targets
all: srrmm bmm combo

srrmm: $(OBJ) srrmm.o
	gcc $(OOPTS) $(OBJ) srrmm.o -o srrmm $(CFLAGS)

bmm: $(OBJ) bmm.o
	gcc $(OOPTS) $(OBJ) bmm.o -o bmm $(CFLAGS)

combo: $(OBJ) b_rb_combo.o
	gcc $(OOPTS) $(OBJ) b_rb_combo.o -o combo $(CFLAGS)

srrmm.o: srrmm.c matrix.h
	gcc $(OOPTS) -c srrmm.c $(CFLAGS)

bmm.o: bmm.c matrix.h
	gcc $(OOPTS) -c bmm.c $(CFLAGS)

b_rb_combo.o: b_rb_combo.c matrix.h
	gcc $(OOPTS) -c b_rb_combo.c $(CFLAGS)

matrix.o: matrix.c matrix.h
	gcc $(OOPTS) -c matrix.c $(CFLAGS)

# --- remove binary and executable files
clean:
	rm -f $(OBJ) srrmm.o srrmm bmm.o bmm b_rb_combo.o combo

# --- submits dgemm-0.job dgemm-1.job dgemm-2.job. Each dimension has 10 iterations to average times.
submit:
	sbatch srrmm.job
	sbatch bmm.job
	sbatch combo.job