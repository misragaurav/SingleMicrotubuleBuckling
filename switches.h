// The top most header - contains switches and constants
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define Pi   acos(-1.0)

#define real double

// Cell shape and size
//#define circle
#define square
#define Dia 40  // Cell radius or square side
#define thickness 40

// clamped or hinged MTs

int iprint, iter, t;
char datadir[64] ;
