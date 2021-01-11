SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc -O3
MPICC = mpicc
CFLAGS=-O3
RM = rm -f

# ==========================================
# TARGETS

EXECUTABLES = v0 v1 v2

default: all

all: $(EXECUTABLES)

v0: v0.c
	$(CC) $< -o $@ -lopenblas -lpthread -lm -fopenmp

v1: v1.c
	$(MPICC) $(CFLAGS) -o v1 v1.c -lopenblas -lpthread -lm -fopenmp
	
v2: v2.c v2.c
	$(MPICC) $(CFLAGS) -o v2 v2.c -lm 

.PHONY: clean

# ==========================================
# TESTS

testv0: 
	@printf "\n** Testing v0 with random array\n"
	./v0 10000 5 20

testv1: 
	@printf "\n** Testing v1 using MPI with random array\n"
	mpirun -np 4 ./v1 10000 5 20

testv2: 
	@printf "\n** Testing v2 using MPI with random array\n"
	mpirun -np 4 ./v2 10000 5 20

# ==========================================
# CLEAN 

clean:
	$(RM) *.o *~ $(EXECUTABLES)
