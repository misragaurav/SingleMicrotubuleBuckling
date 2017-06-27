// MT data structures
#include "params.h"

struct seg {				     
  real rx, ry, rz;             // Coordinates
  real vx, vy, vz;             // Velocity (for node 0)
  real th, ph, si;             // orientation
  real w1, w2, w3;
  real R[3][3];                // rotation matrix - space to body
  int index;
};

struct mt {
  struct seg *seg;           // Nodal data structure
  int N, status;             // Number of segments
  real dL;
  real fcx, fcy, fcz;        // force between mt and centro in space frame
};

struct centro {
  real rx, ry, rz;
  real vx, vy, vz;
  real fx, fy, fz;
};

// Global variables
gsl_rng * rng;

int Nmt;   // number of microtubules
int Nmax;
real ds;   // segment size

