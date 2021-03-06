#include "struct.h"
#include "rod_proto.h"

#define convertdist 1e3   //micro to nm
#define convertt    1e3   //ms to micro s
#define convertYI   1e6   //pN micron^2 to fN nm^2
#define convertk    1e-3  //spring constant - pN/micron to pN/nm
#define convertfric 1e-3  //pN.ms/micron^2 to N.s/m^2

void Qinv (real Q[4][8])  /* Invert Q matrix */
{
  real ratio;
  int    i, j, k;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      Q[i][j+4] = 0.0;
  for (i = 0; i < 4; i++)
    Q[i][i+4] = 1.0;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      {
        if (i == j)  continue;
        ratio = Q[j][i]/Q[i][i];
        for (k = 0; k < 8; k++)
          Q[j][k] -= ratio*Q[i][k];
      }
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      Q[i][j+4] /= Q[i][i];
}


int checkboundary(real rx, real ry, real rz)
{
#ifdef circle
  {
    if(rx*rx + ry*ry > Dia*Dia/4.0 || rz > thickness/2.0 || rz < -thickness/2.0)
      { return 0; }
    else
      { return 1; }
  }
#else

#ifdef square
  {
    if(rx > Dia/2.0 || rx < -Dia/2.0 || ry > Dia/2.0 || ry < -Dia/2.0 || rz > thickness/2.0 || rz < -thickness/2.0)
      { return 0; }
    else
      { return 1; }
  }
#else
  { printf("\n shape misspelt! bailing out \n"); exit(1); }
#endif
#endif
}

void initrod(struct params *params, struct rodparams *rparams, struct rod *rod, real dt)
{
  int n;

  // time step for position update
  rparams->dt  = convertt*dt;

  //parameters from the main code
  rparams->le  = convertdist*ds;                // microns to nm
  rparams->YI  = convertYI*params->YI;
  rparams->k   = convertk*params->k;           // force per unit length - pN/micron to fN/nm 
  rparams->F0  = convertk*params->F0;          // force per unit length - pN/micron to fN/nm
  
  rparams->zt1 = convertfric*params->fric_perp; // 
  rparams->zt2 = convertfric*params->fric_perp;	// fric/length converted to fric/mass in momenta.c
  rparams->zt3 = convertfric*params->fric_para;					       

  // Independent parameters - in units of rod code
  rparams->dia   = 25;    //nm
  rparams->rho   = 1e-9; //fg/nm^3
  rparams->sigma = 0.33;  //dimensionless

  // Derived parameters
  rparams->A = Pi*rparams->dia*rparams->dia/4.0; //Area
  rparams->I = rparams->A*rparams->A/(4.0*Pi);   //Area moment
  rparams->Y = rparams->YI/rparams->I;
  rparams->G = rparams->Y/(2+2*rparams->sigma);

  rparams->xm =   rparams->A*rparams->rho; // mass per unit length
            
  rparams->i1 =   rparams->I*rparams->rho; // Moment of inertia per unit length
  rparams->i2 =   rparams->I*rparams->rho;
  rparams->i3 = 2*rparams->I*rparams->rho;

  rparams->cs1 =   rparams->A*rparams->G;  // Moments of inertia per unit length
  rparams->cs2 =   rparams->A*rparams->G;       
  rparams->cs3 =   rparams->A*rparams->Y;
  rparams->cb1 =   rparams->I*rparams->Y;
  rparams->cb2 =   rparams->I*rparams->Y;
  rparams->cb3 = 2*rparams->I*rparams->G;

  // Reference strains
  rparams->s01 = 0;
  rparams->s02 = 0;					       
  rparams->s03 = 1;					       
  rparams->b01 = 0; 				       
  rparams->b02 = 0;                                        
  rparams->b03 = 0;

  // Rolling fric
  rparams->zr1 = 1e-4*rparams->zt1; //need critical damping for the higher rotational modes 
  rparams->zr2 = 1e-4*rparams->zt1;                                        
  rparams->zr3 = 1e-4*rparams->zt1;

  //Allocate memory to nodes and segments
  for(n=0; n<Nmt; n++)
    {
      rod[n].nodes = (struct node *) calloc(Nmax, sizeof(struct node) ); //Nmax used
      rod[n].elems = (struct elem *) calloc(Nmax, sizeof(struct elem) ); //Nmax used
    }
}

void tempinterfacein(struct mt *mt, struct rod *rod, struct centro *centro, struct rodparams *rparams)
{
  struct seg *ps;
  struct node *pn;
  struct elem *pe;

  real Q[4][8], Qplus[4][4];
  real q0, qx, qy, qz, q;
  real dq0, dqx, dqy, dqz;
  real TH =0, PH=0, PS=0;
  real b1 = 0.00005, b2 = 0, b3 =0;
  real le = rparams->le;
  int k, l, m;
  int N = 20;

  rod->n_nodes = N;
  rod->n_elems = N+1;

  q0 = cos(TH/2.0)*cos((PS+PH)/2.0); qx = sin(TH/2.0)*cos((PH-PS)/2.0);
  qy = sin(TH/2.0)*sin((PH-PS)/2.0); qz = cos(TH/2.0)*sin((PS+PH)/2.0);

  pn = &rod->nodes[0];
  pn->q0 = q0;  pn->qx = qx;  pn->qy = qy; pn->qz = qz;

  for (m = 1; m < rod->n_nodes; m++)
    {
      pn = &rod->nodes[m];

      if(m<=rod->n_nodes/4)
	b1=0;
      if(m>rod->n_nodes/4 && m<=3*rod->n_nodes/4)
	b1= 0.0005*cos(2*Pi*(m-rod->n_nodes/4)/(2*rod->n_nodes/4));
      if(m>3*rod->n_nodes/4)
	b1=0;

      Q[0][0] =  0.0; Q[0][1] = -b1;  Q[0][2] = -b2;  Q[0][3] = -b3;
      Q[1][0] =  b1;  Q[1][1] = 0.0;  Q[1][2] =  b3;  Q[1][3] = -b2;
      Q[2][0] =  b2;  Q[2][1] = -b3;  Q[2][2] = 0.0;  Q[2][3] = -b1;
      Q[3][0] =  b3;  Q[3][1] =  b2;  Q[3][2] = -b1;  Q[3][3] = 0.0;

      for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	  {
	    Q[k][l] *= 0.25*le;
	    Qplus[k][l] = Q[k][l];
	  }

      for (k = 0; k < 4; k++)
	for (l = 0; l < 4; l++)
	  Q[k][l] *= -1.0;
      for (k = 0; k < 4; k++)
	{
	  Q[k][k] += 1.0;
	  Qplus[k][k] += 1.0;
	}
      Qinv (Q);

      dq0 = Qplus[0][0]*q0 + Qplus[0][1]*qx + Qplus[0][2]*qy + Qplus[0][3]*qz;
      dqx = Qplus[1][0]*q0 + Qplus[1][1]*qx + Qplus[1][2]*qy + Qplus[1][3]*qz;
      dqy = Qplus[2][0]*q0 + Qplus[2][1]*qx + Qplus[2][2]*qy + Qplus[2][3]*qz;
      dqz = Qplus[3][0]*q0 + Qplus[3][1]*qx + Qplus[3][2]*qy + Qplus[3][3]*qz;	  

      q0 =  Q[0][4]*dq0 + Q[0][5]*dqx + Q[0][6]*dqy + Q[0][7]*dqz;
      qx =  Q[1][4]*dq0 + Q[1][5]*dqx + Q[1][6]*dqy + Q[1][7]*dqz;
      qy =  Q[2][4]*dq0 + Q[2][5]*dqx + Q[2][6]*dqy + Q[2][7]*dqz;
      qz =  Q[3][4]*dq0 + Q[3][5]*dqx + Q[3][6]*dqy + Q[3][7]*dqz;

      q  = sqrt( q0*q0 + qx*qx + qy*qy + qz*qz);
      q0 /= q; qx /= q; qy /= q;  qz /= q; 
      pn->q0 = q0; pn->qx = qx; pn->qy = qy; pn->qz = qz; 
    }

  rotmatrix_node(rod);
  rotmatrix_elem(rod);

  ps = &mt->seg[0];
  pn = &rod->nodes[0]; // Get the first node

  pn->rx = convertdist*ps->rx;
  pn->ry = convertdist*ps->ry;
  pn->rz = convertdist*ps->rz;

  rod->centrorx = convertdist*centro->rx;
  rod->centrory = convertdist*centro->ry;
  rod->centrorz = convertdist*centro->rz;
 
    // Init strain
  for (m = 0; m < rod->n_elems; m++)
    {
      pe = &rod->elems[m];
      pe->s1  = 0;   
      pe->s2  = 0;   
      pe->s3  = 1;
      pe->b1  = 0;   
      pe->b2  = 0;   
      pe->b3  = 0;
    }
  
  // Initialize positions of nodes other than the first one
  real rx, ry, rz;
  pn = &rod->nodes[0];
  rx = pn->rx;  
  ry = pn->ry;  
  rz = pn->rz;
  
  pn->px = 0.0; pn->py = 0.0; pn->pz = 0.0;  //init p and l
  pn->l1 = 0.0; pn->l2 = 0.0; pn->l3 = 0.0;

  for (m = 1; m < rod->n_nodes; m++)
    {
      pe = &rod->elems[m];
      pn = &rod->nodes[m];

      pn->px = 0.0; pn->py = 0.0; pn->pz = 0.0;  //init p and l
      pn->l1 = 0.0; pn->l2 = 0.0; pn->l3 = 0.0;

      rx += rparams->le*(pe->d[0][0]*pe->s1 + pe->d[1][0]*pe->s2 + pe->d[2][0]*pe->s3);
      ry += rparams->le*(pe->d[0][1]*pe->s1 + pe->d[1][1]*pe->s2 + pe->d[2][1]*pe->s3);
      rz += rparams->le*(pe->d[0][2]*pe->s1 + pe->d[1][2]*pe->s2 + pe->d[2][2]*pe->s3);

      pn->rx = rx;  pn->ry = ry;  pn->rz = rz;
    }
}

void tempinterfaceout(struct mt *mt, struct rod *rod)
{
  struct seg  *ps;
  struct node *pn;
  int i;

  //will need the angles back from the rod code later on when I have to polymerize.

  mt->N = rod->n_nodes;

  for(i=0;i<rod->n_nodes;i++)
    {  
      pn = &rod->nodes[i];
      ps = &mt->seg[i];

      ps->rx = pn->rx/convertdist;
      ps->ry = pn->ry/convertdist;
      ps->rz = pn->rz/convertdist;
    }
}

/*
  void interfacein(struct mt *mt, struct rod *rod, struct centro *centro)
  {
  struct seg *ps;
  struct node *pn;
  real q2, q;
  int i;

  rod->n_nodes = mt->N;
  rod->n_elems = mt->N+1;

  // Populate the quaternions in rod nodes
  for(i=0;i<rod->n_nodes;i++)
  {  
  ps = &mt->seg[i];
  pn = &rod->nodes[i];

  pn->q0 = cos(ps->th/2.0)*cos((ps->ph+ps->si)/2.0);
  pn->qx = sin(ps->th/2.0)*cos((ps->ph-ps->si)/2.0);
  pn->qy = sin(ps->th/2.0)*sin((ps->ph-ps->si)/2.0);
  pn->qz = cos(ps->th/2.0)*sin((ps->ph+ps->si)/2.0);
      
  q2 = pn->q0*pn->q0 + pn->qx*pn->qx + pn->qy*pn->qy + pn->qz*pn->qz;
  q  = sqrt(q2);

  pn->q0 = pn->q0/q;
  pn->qx = pn->qx/q;
  pn->qy = pn->qy/q;
  pn->qz = pn->qz/q;
  }
      
  ps = &mt->seg[0];
  pn = &rod->nodes[0]; // Get the first node position each time

  pn->rx = convertdist*ps->rx;
  pn->ry = convertdist*ps->ry;
  pn->rz = convertdist*ps->rz;

  rod->centrorx = convertdist*centro->rx;
  rod->centrory = convertdist*centro->ry;
  rod->centrorz = convertdist*centro->rz;
  }

  void interfaceout(struct mt *mt, struct rod *rod)
  {
  struct seg  *ps;
  struct node *pn;
  int i;

  for(i=0;i<rod->n_nodes;i++)
  {  
  ps = &mt->seg[i];
  pn = &rod->nodes[i];

  ps->w1 = pn->w1*convertt;
  ps->w2 = pn->w2*convertt;
  ps->w3 = pn->w3*convertt;
  }

  pn = &rod->nodes[0];
  ps = &mt->seg[0];

  ps->vx = (convertt/convertdist)*pn->vx;
  ps->vy = (convertt/convertdist)*pn->vy;
  ps->vz = (convertt/convertdist)*pn->vz;

  }
*/
