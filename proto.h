// Function prototypes
#include "struct.h"
#include "rod_struct.h"

void initmt(struct params *params, struct centro *centro, struct mt *mt);
void initrod(struct params *params, struct rodparams *rparams, struct rod *rod, real dt);
void polymerize(struct mt *mt, struct params *params, real dt);
void randomorient(struct mt *mt);
void addseg(struct mt *mt);
void delseg(struct mt *mt);
void rotmatrix(struct seg *ps);
int checkboundary(real rx, real ry, real rz);

void interfacein (struct mt *mt, struct rod *rod, struct centro *centro);
void rodflex(struct rodparams *rparams, struct rod *rod, int t);
void interfaceout(struct mt *mt, struct rod *rod);

void tempinterfacein (struct mt *mt, struct rod *rod, struct centro *centro, struct rodparams *rparams);
void tempinterfaceout(struct mt *mt, struct rod *rod);

void eulerupdate(struct mt *mt, real dt);
void buildmt(struct mt *mt);
void movecentro(struct params *params, struct centro *centro, struct mt *mt, real dt);

void Qinv(real Q[][8]);
