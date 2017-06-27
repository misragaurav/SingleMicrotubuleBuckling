#include "proto.h"
#define boundary_depoly

void polymerize(struct mt *mt, struct params *params, real dt)
{
  int i =0;

  // switching between poly and depoly
  /*
  if (mt->status == 0 && gsl_rng_uniform(rng) < params->k_catast*dt)
    mt->status = 1;
  if (mt->status == 1 && gsl_rng_uniform(rng) < params->k_recov*dt)
    mt->status = 0;
  */

    mt->status = 0;

  // growing and shrinking
  if (mt->status == 0) // grow
    {
      mt->dL = mt->dL + params->Vpoly*dt;
      if (mt->dL > ds)
	{
	  for(i=0; i<(int)(mt->dL/ds); i++)
	    {
	      if (mt->N >= Nmax)
		{
		  mt->dL = 0;  // reset dL if N=Nmax and add no more
		}
	      else
		{
		  addseg(mt);
		}
	    }
	}
    }
  
  else  //shrink
    {
      mt->dL = mt->dL - params->Vdepoly*dt;
      if (mt->dL <0)
	{
	  for(i=0; i>(int)(mt->dL/ds); i--)
	    {	  
	      if (mt->N <= 1)
		{
		  randomorient(mt);
		  mt->dL = 0;   // reset dL if N=1
		}
	      else
		{
		  delseg(mt);
		}
	    }
	}
    }
}

void addseg(struct mt *mt)
{
  struct seg *ps1, *ps2;
  real rx, ry, rz, drx, dry, drz, ds1, ds2, ds3;

  ds1=0; 
  ds2=0; 
  ds3=ds;

  ps1   = &mt->seg[mt->N-1];

  drx = ps1->R[0][0]*ds1 + ps1->R[1][0]*ds2 + ps1->R[2][0]*ds3;
  dry = ps1->R[0][1]*ds1 + ps1->R[1][1]*ds2 + ps1->R[2][1]*ds3;
  drz = ps1->R[0][2]*ds1 + ps1->R[1][2]*ds2 + ps1->R[2][2]*ds3;

  rx = ps1->rx + drx;
  ry = ps1->ry + dry;
  rz = ps1->rz + drz;
 
  // check boundary
  if(checkboundary(rx, ry, rz) == 1)
    {
      mt->dL = mt->dL - ds;
      mt->N  = mt->N + 1;
      //      printf("Nmax: %i \t N: %i \n",Nmax, mt->N);
      ps2    = &mt->seg[mt->N-1];
      ps2->rx = rx;
      ps2->ry = ry;
      ps2->rz = rz;
      ps2->th = ps1->th;
      ps2->ph = ps1->ph;
      ps2->si = ps1->si;
      ps2->index = ps1->index;
      rotmatrix(ps2);
    }
  else
    {
      //      printf("growing out of bounds \n");
      mt->dL = 0; // reset dL if it hits the boundary
#ifdef boundary_depoly
      mt->status = 1;
#endif
    }
}

void delseg(struct mt *mt)
{
  struct seg *ps;

  ps   = &mt->seg[mt->N-1];

  ps->rx = 0;
  ps->ry = 0;
  ps->rz = 0;
  ps->th = 0;
  ps->ph = 0;
  ps->si = 0;
  ps->index = 0; // Left out the rotation matrix

  mt->N =  mt->N - 1;
  mt->dL = mt->dL + ds;
}
