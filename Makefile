SHELL := /bin/bash


# ============================================
# COMMANDS

CC = gcc -O3
RM = rm -f

# ==========================================
# TARGETS


EXECUTABLES = v0

all: $(EXECUTABLES)

v0: v0.c
	$(CC) $< -o $@ -lopenblas -lpthread -lm

# hello_openblas: hello_openblas.c
# 	$(CC) $< -o $@ -lopenblas -lpthread

test:
	@printf "\n** Testing v0\n"
	./v0

clean:
	$(RM) *.o *~ $(EXECUTABLES)
