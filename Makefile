SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc -O3
MPICC = mpicc
CFLAGS=-O3
RM = rm -f

# ==========================================
# TARGETS


EXECUTABLES = v1 v2_c

all: $(EXECUTABLES)

v0: v0.c
	$(CC) $< -o $@ -lopenblas -lpthread -lm

v0alt: v0alt.c
	$(CC) $< -o $@ -lopenblas -lpthread -lm -fopenmp

v1: v0alt.c v1.c
	$(MPICC) $(CFLAGS) -o v1 v1.c -lopenblas -lpthread -lm -fopenmp

#v2: v2.c
#	$(CC) $< -o $@ -lm 
	
v2_:v2.c v2_.c
	$(MPICC) $(CFLAGS) -o v2_ v2_.c -lm 

test1:
	#@printf "\n** Testing v0\n"
	#./v0
	#@printf "\n** Testing v0alt with 4 threads by default\n"
	#./v0alt 4
	@printf "\n** Testing MPI v1 with 3 mpi processes \n"
	mpirun -np 3 v1

	
test2:
	@printf "\n** Testing MPI v2 with 3 mpi processes \n"
	mpirun -np 3 v2_

clean:
	$(RM) *.o *~ $(EXECUTABLES)

