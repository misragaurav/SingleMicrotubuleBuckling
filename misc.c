#include "proto.h"

void randomorient(struct mt *mt)
{
  struct seg *ps;
  real x=1, y=1, z=1;

  while(x*x+y*y+z*z>=1)
    {
      // generate points in a cube
      x = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
      y = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
      z = 2*gsl_rng_uniform(rng)-1; // [-1 +1)
      // accept those which lie in the enclosed sphere
    }

  ps =&mt->seg[0];  //// Always 0th node because randomorient() is called only when there is just 1 node
 
  /*
  ps->th = acos(z);
  
  if(ps->th < Pi*0.2 || ps->th > Pi*0.8)
     ps->th = Pi/2.0;
  */
  ps->th = Pi/2.0;
  if(atan2(y, x) < 0)
    ps->ph = atan2(y, x) + 2*Pi; // to convert from (-pi pi) to (0 2pi)
  else 
    ps->ph = atan2(y, x);
  
  ps->si = 0;

  ps->index = 0;

  //////////////////////////
  /*
  if(2*gsl_rng_uniform(rng)-1 > 0)
    ps->ph = 3*Pi/2.0;
  else
    ps->ph = Pi/2.0;

  printf("%f \n", ps->ph);

  ps->th = Pi/2.0;  
  ps->si = 0;
  */
  //////////

  rotmatrix(ps);
}

void rotmatrix(struct seg *ps)
{
  real th, ph, si;
  th = ps->th;
  ph = ps->ph;
  si = ps->si;
  // following the x convention as used in Goldstein
  ps->R[0][0] = cos(si)*cos(ph)-cos(th)*sin(ph)*sin(si); ps->R[0][1] = cos(si)*sin(ph)+cos(th)*cos(ph)*sin(si); ps->R[0][2] = sin(si)*sin(th);
  ps->R[1][0] =-sin(si)*cos(ph)-cos(th)*sin(ph)*cos(si); ps->R[1][1] =-sin(si)*sin(ph)+cos(th)*cos(ph)*cos(si); ps->R[1][2] = cos(si)*sin(th);
  ps->R[2][0] = sin(th)*sin(ph);                         ps->R[2][1] =-sin(th)*cos(ph);                         ps->R[2][2] = cos(th);
}
/*
void eulerupdateangles(struct mt *mt, real dt)
{
  struct seg *ps;
  real thdot, phdot, sidot;
  real thold, phold, siold;
  real thini, phini, siini;
  real rx, ry, rz;
  real tol= 10e-12;
  real err = 1;
  int i;
  
  for(i=0;i<mt->N;i++)
    {
      ps = &mt->seg[i];

      thini = ps->th;
      phini = ps->ph;
      siini = ps->si;

      err = 1;
      
      while(err>tol)
	{
	  thold = ps->th;
	  phold = ps->ph;
	  siold = ps->si;

	  // printf("%i th: \t %.16e \t %.16e \t %.16e \n", i, ps->th, ps->ph, ps->si);

	  sidot = ps->w3/(1+cos(ps->th));
	  phdot = ( ps->w1*sin(ps->si) + ps->w2*cos(ps->si) )/sin(ps->th);
	  thdot =   ps->w1*cos(ps->si) - ps->w2*sin(ps->si);

	  //	  	  printf("%i \t %.16e \t %.16e \t %.16e \n", i, ps->w1 ,ps->w2 ,ps->w3 );

	  ps->th = thini + sidot*dt; //Euler update
	  ps->ph = phini + phdot*dt;
	  ps->si = siini + thdot*dt;

	  //	   printf("%i th: \t %.16e \t %.16e \t %.16e \n", i, ps->th, ps->ph, ps->si);

	  err = fabs( fabs(ps->th)+fabs(ps->ph)+fabs(ps->si)-fabs(thold)-fabs(phold)-fabs(siold) );
	  //	  printf("euler update. err = %.16e \t tol= %.16e \n", err, tol);

	}
      
      rotmatrix(ps);
    }

  ps = &mt->seg[0];
  rx = ps->rx + ps->vx*dt; //Euler update
  ry = ps->ry + ps->vy*dt;
  rz = ps->rz + ps->vz*dt;
  
  if(checkboundary(rx, ry, rz) == 1)
    {
      ps->rx = rx;
      ps->ry = ry;
      ps->rz = rz;
    }
}
*/
void positionupdate(struct rod *rod, struct rodparams *rparams, real dt)
{
  struct node *pn;
  int  m;

  for (m = 0; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];
      
      pn->rx = pn->rx + pn->px ;
      pn->ry = pn->ry + pn->py ;
      pn->rz = pn->rz + pn->pz ;
    }
}

void buildmt(struct mt *mt)
{
  struct seg *ps, *ps2;
  real rx, ry, rz, drx, dry, drz, ds1, ds2, ds3;
  int i;

  ds1=0; ds2=0; ds3=ds;

  for(i=1;i<mt->N;i++)
    {
      ps  = &mt->seg[i-1];
      ps2 = &mt->seg[i];

      drx = ps->R[0][0]*ds1 + ps->R[1][0]*ds2 + ps->R[2][0]*ds3;
      dry = ps->R[0][1]*ds1 + ps->R[1][1]*ds2 + ps->R[2][1]*ds3;
      drz = ps->R[0][2]*ds1 + ps->R[1][2]*ds2 + ps->R[2][2]*ds3;

      rx = ps->rx + drx;
      ry = ps->ry + dry;
      rz = ps->rz + drz;
 
      if(checkboundary(rx, ry, rz) == 1)
	{
	  ps2->rx = rx;
	  ps2->ry = ry;
	  ps2->rz = rz;
	}
    }
}

void movecentro(struct params *params, struct centro *centro, struct mt *mt, real dt)
{
  struct seg *ps;
  real rx, ry, rz; 
  int n;

  centro->fx = 0;
  centro->fy = 0;
  centro->fz = 0;

  for(n=0; n<Nmt; n++)
    { 
      ps = &mt[n].seg[0];
      mt[n].fcx = params->k*(ps->rx - centro->rx); // individual MT force
      mt[n].fcy = params->k*(ps->ry - centro->ry);
      mt[n].fcz = params->k*(ps->rz - centro->rz);
    }

  for(n=0; n<Nmt; n++)
    {
      centro->fx = centro->fx + mt[n].fcx;  //sum of all MT forces
      centro->fy = centro->fy + mt[n].fcy;
      centro->fz = centro->fz + mt[n].fcz;
    }

  rx = centro->rx + (centro->fx/params->fric_centro)*dt;
  ry = centro->ry + (centro->fy/params->fric_centro)*dt;
  rz = centro->rz + (centro->fz/params->fric_centro)*dt;

  if(checkboundary(rx, ry, rz) == 1)
    {
      centro->rx = rx;
      centro->ry = ry;
      centro->rz = rz;
    }
}
