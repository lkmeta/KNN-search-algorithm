SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc -O3
MPICC = mpicc
CFLAGS=-O3
RM = rm -f

# ==========================================
# TARGETS


EXECUTABLES = v1 arxikh v2

all: $(EXECUTABLES)

v0: v0.c
	$(CC) $< -o $@ -lopenblas -lpthread -lm

v0alt: v0alt.c
	$(CC) $< -o $@ -lopenblas -lpthread -lm -fopenmp

v1: v0alt.c v1.c
	$(MPICC) $(CFLAGS) -o v1 v1.c -lopenblas -lpthread -lm -fopenmp

arxikh: v0alt.c arxikh.c
	$(MPICC) $(CFLAGS) -o arxikh arxikh.c -lopenblas -lpthread -lm -fopenmp

v2: v2.c
	$(CC) $< -o $@ -lm 

test:
	#@printf "\n** Testing v0\n"
	#./v0
	#@printf "\n** Testing v0alt with 4 threads by default\n"
	#./v0alt 4
	@printf "\n** Testing MPI v1 with 3 mpi processes \n"
	mpirun -np 3 v1

test1:
	@printf "\n** Testing MPI arxikh with 3 mpi processes \n"
	mpirun -np 3 arxikh

test2:
	@printf "\n** Testing v2 \n"
	./v2

clean:
	$(RM) *.o *~ $(EXECUTABLES)

