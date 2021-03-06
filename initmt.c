#include "proto.h"

void initmt(struct params *params, struct centro *centro, struct mt *mt)
{
  struct seg *ps;
  int n = 0;

  // Init params
  params->YI          = 25; // in pN micron^2
  params->Vpoly       = 1e-4*100; //micron per ms
  params->Vdepoly     = 3e-4; //micron per ms

  params->k_catast    = 5e-5; // per ms
  params->k_recov     = 2e-4; // per ms
  params->Vmax        = 1e-3; // micron per ms
  params->Fmax 	      = 10;   // pN
  params->koff 	      = 1e-3; // per ms
  params->kappa       = 1000; // pN/micron

  params->rho_motor   = 10;    // per micron

  params->F0          = params->rho_motor*params->Fmax/( params->Fmax*params->koff/(params->kappa*params->Vmax) + 1 );
  params->fric_para   = params->F0/params->Vmax;
  params->fric_perp   = params->rho_motor*params->kappa/params->koff;
  //printf("F0: %f \t fric_para: %f \t fric_perp: %f \n", params->F0, params->fric_para, params->fric_perp );

  params->fric_centro = params->fric_para*100;

  params->k           = 100; // pN/micron - centrosome MT spring

  // Init Centro
  centro->rx = 0;
  centro->ry = 0;
  centro->rz = 0;

  centro->vx = 0;
  centro->vy = 0;
  centro->vz = 0;

  centro->fx = 0;
  centro->fy = 0;
  centro->fz = 0;

  // Init MT

  for(n=0; n<Nmt; n++)
    {
      mt[n].status = 0 ; // Initially growing status for all MTs
      mt[n].dL     = 0 ;
      mt[n].N      = 1 ; // Must start with 1 for the current code

      mt[n].seg = (struct seg *) calloc(Nmax, sizeof(struct seg) ); //Nmax used

      ps = &mt[n].seg[mt[n].N-1]; // placing all MTs at centrosome
      ps->rx = centro->rx;
      ps->ry = centro->ry;
      ps->rz = centro->rz;      

      randomorient(&mt[n]);
    }
}
