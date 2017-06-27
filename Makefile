# Choose compiler
CC = gcc
# setenv GSL_RNG_TYPE r250 -- to set environment variable
# Optimization flags: cc is generic, icc is Intel's compiler, gcc is the standard gnu compiler
CFLAGS_gcc = -g -Wall -DHAVE_INLINE #Usually gcc is needed for debugging
CFLAGS_icc = -O3 -xP  -DHAVE_INLINE

# Set linker flags: i.e. mpe libraries or options like -static

#LFLAGS_gcc = -static
#LFLAGS_icc = -static

LFLAGS_gcc = 
LFLAGS_icc = 

# Compilation and linking flags 
CFLAGS = $(CFLAGS_$(CC)) -c
LFLAGS = $(LFLAGS_$(CC)) -lm -lgsl -lgslcblas

C_SRC     = main.c initmt.c polymerize.c misc.c interface.c rod.c forces.c rod_init.c propagator.c rod_misc.c eulerstep.c
OBJ_flex  = main.o initmt.o polymerize.o misc.o interface.o rod.o forces.o rod_init.o propagator.o rod_misc.o eulerstep.o

.c.o:
	$(CC)  $(CFLAGS) $*.c

flex: $(OBJ_flex)
	$(CC)  $(LFLAGS) -o mt $(OBJ_flex)

clean:
	rm -f *.o mt
