#include "switches.h"
#include "rod_struct.h"

// Functions needed across files
void init(struct rod *rod, struct rodparams *rparams);                     // function in initrod.c file
void euler(struct rod *rod, struct rodparams *rparams, real dt);
void force(struct rod *rod, struct rodparams *rparams);           // function in forces.c file
void propagator(struct rod *rod, struct rodparams *rparams, real dt);

void rotmatrix_node(struct rod *rod);
void rotmatrix_elem(struct rod *rod);

void file_name (char *name, char *work_dir, int task_number);               // function in misc.c file

void com(struct rod *rod, struct rodparams *rparams);
void energy(struct rod *rod, struct rodparams *rparams, int t);

void friction(struct rod *rod, real dt, real factor);
