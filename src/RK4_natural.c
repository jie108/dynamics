#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/////////   c code for implenmenting  fourth order Runge-Kutta method for solving an ODE
////////    using cubic B-spline as basis functions
////////    date: 5-18-2008
////////    date: 5-21-08: add hessian

/////////  accessary functions

void NCBasis(double *t, double *knot, int Lt, int Mk, double *result);
void NPCBasis(double *t, double *knot, int Lt, int Mk, double *result);
void NHCBasis(double *t, double *knot, int Lt, int Mk, double *result);
void gOneN(double *t, double *beta, double *knot, int Lt, int Mk, double *result);
void gTwoN(double *t, double *beta, double *knot, int Lt, int Mk, double *result);
void gThreeN(double *t, double *beta, double *knot, int Lt, int Mk, double *result);
double Pos(double x);


//void xPathN(int *NN, double *hh, double *x0, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);
//void xDerivAN(int *NN, double *hh, double *y, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);
//void xDerivThetaN(int *NN, double *hh, double *y, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);
//void xDerivBetaN(int *NN, double *hh, double *y, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);
//void xHessThetaN(int *NN, double *hh, double *y, double *v, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);
//void xHessAN(int *NN, double *hh, double *y, double *z, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);
//void xHessBetaN(int *NN, double *hh, double *y, double *w, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result);

///// RK solution of the  sample path for the non-linear dynamical model
void xPathN(int *NN, double *hh, double *x0, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int N, Lx, Mk;
double h;
int i,j;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;

ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));


for(i=0;i<Lx;i++)
 {
	result[i*(N+1)]=x0[i];
	ycur[i]=x0[i];
 }

for(j=0;j<N;j++)
 {
    gOneN(ycur, beta, knot, Lx, Mk, kcur1);

    for (i=0;i<Lx;i++)
    {
		kcur1[i]=kcur1[i]*exp(theta[i]);
		ytemp[i]=ycur[i]+h*kcur1[i]/2;

	}
    gOneN(ytemp, beta, knot, Lx, Mk, kcur2);


    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=kcur2[i]*exp(theta[i]);
			ytemp[i]=ycur[i]+h*kcur2[i]/2;

		}
     gOneN(ytemp, beta, knot, Lx, Mk, kcur3);


    for (i=0;i<Lx;i++)
	    {
			kcur3[i]=kcur3[i]*exp(theta[i]);
			ytemp[i]=ycur[i]+h*kcur3[i];

		}

    gOneN(ytemp, beta, knot, Lx, Mk, kcur4);


    for (i=0;i<Lx;i++)
	    {
			kcur4[i]=kcur4[i]*exp(theta[i]);
			ycur[i]=ycur[i]+h*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
            result[i*(N+1)+(j+1)]=ycur[i];
		}


 }

free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);

}

/// RK solution for derivatives given the sample path
void xDerivAN(int *NN, double *hh, double *y, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int Nh;
int i,j;
double h2;
double  *g2,*gtemp;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;
int N, Lx, Mk;
double h;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;

h2=2*h;
Nh=N/2;

g2=(double *) malloc (Lx*(N+1)*sizeof(double));
gtemp=(double *) malloc (Lx*sizeof(double));
ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));

 for (j=0;j<(N+1);j++)
 {
	 for(i=0;i<Lx;i++)
	 {
		 ycur[i]=y[i*(N+1)+j];

	 }

	 gTwoN(ycur, beta, knot, Lx, Mk, gtemp);


	for(i=0;i<Lx;i++)
	 {
	 	g2[i*(N+1)+j]=gtemp[i];

	 }

 }


for(i=0;i<Lx;i++)
 {
	result[i*(Nh+1)]=1;
	ycur[i]=1;
 }


for(j=0;j<Nh;j++)
 {

    for (i=0;i<Lx;i++)
    {
		kcur1[i]=exp(theta[i])*ycur[i]*g2[i*(N+1)+(2*j)];
		ytemp[i]=ycur[i]+h2*kcur1[i]/2;

	}

    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=exp(theta[i])*ytemp[i]*g2[i*(N+1)+(2*j+1)];
			ytemp[i]=ycur[i]+h2*kcur2[i]/2;

		}

    for (i=0;i<Lx;i++)
	    {
			kcur3[i]=exp(theta[i])*ytemp[i]*g2[i*(N+1)+(2*j+1)];
			ytemp[i]=ycur[i]+h2*kcur3[i];

		}


    for (i=0;i<Lx;i++)
	    {
			kcur4[i]=exp(theta[i])*ytemp[i]*g2[i*(N+1)+(2*j+2)];
			ycur[i]=ycur[i]+h2*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
            result[i*(Nh+1)+(j+1)]=ycur[i];
		}


 }

free(g2);
free(gtemp);
free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);

}



void xDerivThetaN(int *NN, double *hh, double *y, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int Nh;
int i,j;
double h2;
double  *g1, *g2,*gtemp;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;
int N, Lx, Mk;
double h;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;

h2=2*h;
Nh=N/2;

g1=(double *) malloc (Lx*(N+1)*sizeof(double));
g2=(double *) malloc (Lx*(N+1)*sizeof(double));
gtemp=(double *) malloc (Lx*sizeof(double));
ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));

 for (j=0;j<(N+1);j++)
 {
	 for(i=0;i<Lx;i++)
	 {
		 ycur[i]=y[i*(N+1)+j];

	 }

	 gTwoN(ycur, beta, knot, Lx, Mk, gtemp);


	for(i=0;i<Lx;i++)
	 {
	 	g2[i*(N+1)+j]=gtemp[i];

	 }

     gOneN(ycur, beta, knot, Lx, Mk, gtemp);


		for(i=0;i<Lx;i++)
		 {
		 	g1[i*(N+1)+j]=gtemp[i];


	 }
 }


for(i=0;i<Lx;i++)
 {
	result[i*(Nh+1)]=0;
	ycur[i]=0;
 }


for(j=0;j<Nh;j++)
 {
  for (i=0;i<Lx;i++)
    {
		kcur1[i]=exp(theta[i])*(ycur[i]*g2[i*(N+1)+(2*j)]+g1[i*(N+1)+(2*j)]);
		ytemp[i]=ycur[i]+h2*kcur1[i]/2;

	}

    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(2*j+1)]+g1[i*(N+1)+(2*j+1)]);
			ytemp[i]=ycur[i]+h2*kcur2[i]/2;

		}

    for (i=0;i<Lx;i++)
	    {
			kcur3[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(2*j+1)]+g1[i*(N+1)+(2*j+1)]);
			ytemp[i]=ycur[i]+h2*kcur3[i];

		}


    for (i=0;i<Lx;i++)
	    {
			kcur4[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(2*j+2)]+g1[i*(N+1)+(2*j+2)]);
			ycur[i]=ycur[i]+h2*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
            result[i*(Nh+1)+(j+1)]=ycur[i];
		}


 }

free(g1);
free(g2);
free(gtemp);
free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);

}




void xDerivBetaN(int *NN, double *hh, double *y, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int Nh;
int i,j,k;
double h2;
double  *g1, *g2,*gtemp;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;
double  *betacur;
int N, Lx, Mk;
double h;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;


h2=2*h;
Nh=N/2;

g1=(double *) malloc (Lx*(N+1)*sizeof(double));
g2=(double *) malloc (Lx*(N+1)*sizeof(double));
gtemp=(double *) malloc (Lx*sizeof(double));
ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));
betacur=(double *) malloc ((Mk+3)*sizeof(double));

 for (j=0;j<(N+1);j++)
 {
	 for(i=0;i<Lx;i++)
	 {
		 ycur[i]=y[i*(N+1)+j];

	 }

	 gTwoN(ycur, beta, knot, Lx, Mk, gtemp);


	for(i=0;i<Lx;i++)
	 {
	 	g2[i*(N+1)+j]=gtemp[i];

	 }
 }

 for(k=0;k<(Mk+3);k++)
 {
	 betacur[k]=0;

 }


for (k=0;k<(Mk+3);k++)
{

 betacur[k]=1;


  for (j=0;j<(N+1);j++)
  {
 	 for(i=0;i<Lx;i++)
 	 {
 		 ycur[i]=y[i*(N+1)+j];

 	 }

 	 gOneN(ycur, betacur, knot, Lx, Mk, gtemp);


 	for(i=0;i<Lx;i++)
 	 {
 	 	g1[i*(N+1)+j]=gtemp[i];

 	 }
 }

  betacur[k]=0;

 for(i=0;i<Lx;i++)
 {
	result[k*(Nh+1)*Lx+i]=0;
	ycur[i]=0;
 }


for(j=0;j<Nh;j++)
 {

  for (i=0;i<Lx;i++)
    {
		kcur1[i]=exp(theta[i])*(ycur[i]*g2[i*(N+1)+(2*j)]+g1[i*(N+1)+(2*j)]);
		ytemp[i]=ycur[i]+h2*kcur1[i]/2;

	}

    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(2*j+1)]+g1[i*(N+1)+(2*j+1)]);
			ytemp[i]=ycur[i]+h2*kcur2[i]/2;

		}

    for (i=0;i<Lx;i++)
	    {
			kcur3[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(2*j+1)]+g1[i*(N+1)+(2*j+1)]);
			ytemp[i]=ycur[i]+h2*kcur3[i];

		}



    for (i=0;i<Lx;i++)
	    {
			kcur4[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(2*j+2)]+g1[i*(N+1)+(2*j+2)]);
			ycur[i]=ycur[i]+h2*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
			//result[i*(Nh+1)*Mk+(j+1)*Mk+k]=ycur[i];
            result[k*(Nh+1)*Lx+(j+1)*Lx+i]=ycur[i];
		}


 }


}

free(g1);
free(g2);
free(gtemp);
free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);
free(betacur);

}


//// RK solution for Hessians given the derivatives and sample path
void xHessThetaN(int *NN, double *hh, double *y, double *v, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int Nh, Nh2;
int i,j;
double h4;
double  *g1, *g2,*g3, *gtemp;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;
int N, Lx, Mk;
double h;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;

h4=4*h;
Nh=N/4;
Nh2=N/2;


g1=(double *) malloc (Lx*(N+1)*sizeof(double));
g2=(double *) malloc (Lx*(N+1)*sizeof(double));
g3=(double *) malloc (Lx*(N+1)*sizeof(double));
gtemp=(double *) malloc (Lx*sizeof(double));
ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));

 for (j=0;j<(N+1);j++)
 {
	 for(i=0;i<Lx;i++)
	 {
		 ycur[i]=y[i*(N+1)+j];

	 }


	 gTwoN(ycur, beta, knot, Lx, Mk, gtemp);
	for(i=0;i<Lx;i++)
	 {
	 	g2[i*(N+1)+j]=gtemp[i];

	 }


     gOneN(ycur, beta, knot, Lx, Mk, gtemp);
		for(i=0;i<Lx;i++)
		 {
		 	g1[i*(N+1)+j]=gtemp[i];
	 }



	  gThreeN(ycur, beta, knot, Lx, Mk, gtemp);
	 		for(i=0;i<Lx;i++)
	 		 {
	 		 	g3[i*(N+1)+j]=gtemp[i];
	 }

}


for(i=0;i<Lx;i++)
 {
	result[i*(Nh+1)]=0;
	ycur[i]=0;
 }


for(j=0;j<Nh;j++)
 {


  for (i=0;i<Lx;i++)
    {
		kcur1[i]=exp(theta[i])*(g1[i*(N+1)+(4*j)]+(ycur[i]+2*v[i*(Nh2+1)+(2*j)])*g2[i*(N+1)+(4*j)]+v[i*(Nh2+1)+(2*j)]*v[i*(Nh2+1)+(2*j)]*g3[i*(N+1)+(4*j)]);
		ytemp[i]=ycur[i]+h4*kcur1[i]/2;

	}


    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=exp(theta[i])*(g1[i*(N+1)+(4*j+2)]+(ytemp[i]+2*v[i*(Nh2+1)+(2*j+1)])*g2[i*(N+1)+(4*j+2)]+v[i*(Nh2+1)+(2*j+1)]*v[i*(Nh2+1)+(2*j+1)]*g3[i*(N+1)+(4*j+2)]);
			ytemp[i]=ycur[i]+h4*kcur2[i]/2;

		}

   for (i=0;i<Lx;i++)
	    {
			kcur3[i]=exp(theta[i])*(g1[i*(N+1)+(4*j+2)]+(ytemp[i]+2*v[i*(Nh2+1)+(2*j+1)])*g2[i*(N+1)+(4*j+2)]+v[i*(Nh2+1)+(2*j+1)]*v[i*(Nh2+1)+(2*j+1)]*g3[i*(N+1)+(4*j+2)]);
			ytemp[i]=ycur[i]+h4*kcur3[i];
		}
   for (i=0;i<Lx;i++)
	    {
			kcur4[i]=exp(theta[i])*(g1[i*(N+1)+(4*j+4)]+(ytemp[i]+2*v[i*(Nh2+1)+(2*j+2)])*g2[i*(N+1)+(4*j+4)]+v[i*(Nh2+1)+(2*j+2)]*v[i*(Nh2+1)+(2*j+2)]*g3[i*(N+1)+(4*j+4)]);
			ycur[i]=ycur[i]+h4*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
            result[i*(Nh+1)+(j+1)]=ycur[i];
		}


 }

free(g1);
free(g2);
free(g3);
free(gtemp);
free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);

}

void xHessAN(int *NN, double *hh, double *y, double *z, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int Nh, Nh2;
int i,j;
double h4;
double  *g2,*g3, *gtemp;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;
int N, Lx, Mk;
double h;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;

h4=4*h;
Nh=N/4;
Nh2=N/2;



g2=(double *) malloc (Lx*(N+1)*sizeof(double));
g3=(double *) malloc (Lx*(N+1)*sizeof(double));
gtemp=(double *) malloc (Lx*sizeof(double));
ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));

 for (j=0;j<(N+1);j++)
 {
	 for(i=0;i<Lx;i++)
	 {
		 ycur[i]=y[i*(N+1)+j];

	 }


	 gTwoN(ycur, beta, knot, Lx, Mk, gtemp);
	for(i=0;i<Lx;i++)
	 {
	 	g2[i*(N+1)+j]=gtemp[i];

	 }




	  gThreeN(ycur, beta, knot, Lx, Mk, gtemp);
	 		for(i=0;i<Lx;i++)
	 		 {
	 		 	g3[i*(N+1)+j]=gtemp[i];
	 }

}


for(i=0;i<Lx;i++)
 {
	result[i*(Nh+1)]=0;
	ycur[i]=0;
 }


for(j=0;j<Nh;j++)
 {


  for (i=0;i<Lx;i++)
    {
		kcur1[i]=exp(theta[i])*(ycur[i]*g2[i*(N+1)+(4*j)]+z[i*(Nh2+1)+(2*j)]*z[i*(Nh2+1)+(2*j)]*g3[i*(N+1)+(4*j)]);
		ytemp[i]=ycur[i]+h4*kcur1[i]/2;

	}


    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(4*j+2)]+z[i*(Nh2+1)+(2*j+1)]*z[i*(Nh2+1)+(2*j+1)]*g3[i*(N+1)+(4*j+2)]);
			ytemp[i]=ycur[i]+h4*kcur2[i]/2;

		}

   for (i=0;i<Lx;i++)
	    {
			kcur3[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(4*j+2)]+z[i*(Nh2+1)+(2*j+1)]*z[i*(Nh2+1)+(2*j+1)]*g3[i*(N+1)+(4*j+2)]);
			ytemp[i]=ycur[i]+h4*kcur3[i];
		}
   for (i=0;i<Lx;i++)
	    {
			kcur4[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(4*j+4)]+z[i*(Nh2+1)+(2*j+2)]*z[i*(Nh2+1)+(2*j+2)]*g3[i*(N+1)+(4*j+4)]);
			ycur[i]=ycur[i]+h4*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
            result[i*(Nh+1)+(j+1)]=ycur[i];
		}


 }


free(g2);
free(g3);
free(gtemp);
free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);

}




void xHessBetaN(int *NN, double *hh, double *y, double *w, double *theta, double *beta, double *knot, int *Lxx, int *Mkk, double *result)
{
int Nh, Nh2;
int i,j,k, k1,k2;
double h4;
double  *g2,*g3, *g2l, *gtemp;
double  *ycur, *ytemp;
double  *kcur1, *kcur2, *kcur3, *kcur4;
double *betal;
int N, Lx, Mk;
double h;

 N=*NN;
 Lx=*Lxx;
 Mk=*Mkk;
 h=*hh;

h4=4*h;
Nh=N/4;
Nh2=N/2;



g2=(double *) malloc (Lx*(N+1)*sizeof(double));
g3=(double *) malloc (Lx*(N+1)*sizeof(double));
g2l=(double *) malloc ((Mk+3)*Lx*(N+1)*sizeof(double));
gtemp=(double *) malloc (Lx*sizeof(double));
ycur=(double *) malloc (Lx*sizeof(double));
ytemp=(double *) malloc (Lx*sizeof(double));
kcur1=(double *) malloc (Lx*sizeof(double));
kcur2=(double *) malloc (Lx*sizeof(double));
kcur3=(double *) malloc (Lx*sizeof(double));
kcur4=(double *) malloc (Lx*sizeof(double));
betal=(double *) malloc ((Mk+3)*sizeof(double));


 for (j=0;j<(N+1);j++)
 {
	 for(i=0;i<Lx;i++)
	 {
		 ycur[i]=y[i*(N+1)+j];

	 }



	 gTwoN(ycur, beta, knot, Lx, Mk, gtemp);
	for(i=0;i<Lx;i++)
	 {
	 	g2[i*(N+1)+j]=gtemp[i];

	 }




	  gThreeN(ycur, beta, knot, Lx, Mk, gtemp);
	 		for(i=0;i<Lx;i++)
	 		 {
	 		 	g3[i*(N+1)+j]=gtemp[i];
	 }

}

 for(k=0;k<Mk+3;k++)
 {

	 betal[k]=0;

 }

 for(k=0;k<Mk+3;k++)
 {
	betal[k]=1;

	 for (j=0;j<(N+1);j++)
	 {
		 for(i=0;i<Lx;i++)
		 {
			 ycur[i]=y[i*(N+1)+j];

		 }



		 gTwoN(ycur, betal, knot, Lx, Mk, gtemp);
		for(i=0;i<Lx;i++)
		 {
		 	g2l[k*Lx*(N+1)+i*(N+1)+j]=gtemp[i];

		 }

    }

	betal[k]=0;
 }




for(k1=0;k1<Mk+3;k1++)
{

  for(k2=0;k2<Mk+3;k2++)
  {
      for(i=0;i<Lx;i++)
      {
     	result[k1*(Mk+3)*Lx*(Nh+1)+k2*Lx*(Nh+1)+i]=0;
    	ycur[i]=0;
      }



for(j=0;j<Nh;j++)
 {

 for (i=0;i<Lx;i++)
    {

		kcur1[i]=exp(theta[i])*(ycur[i]*g2[i*(N+1)+(4*j)]+w[k1*(Nh2+1)*Lx+(2*j)*Lx+i]*g2l[k2*Lx*(N+1)+i*(N+1)+(4*j)]+w[k2*(Nh2+1)*Lx+(2*j)*Lx+i]*g2l[k1*Lx*(N+1)+i*(N+1)+(4*j)]+w[k1*(Nh2+1)*Lx+(2*j)*Lx+i]*w[k2*(Nh2+1)*Lx+(2*j)*Lx+i]*g3[i*(N+1)+(4*j)]);
		ytemp[i]=ycur[i]+h4*kcur1[i]/2;

	}


    for (i=0;i<Lx;i++)
	    {
			kcur2[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(4*j+2)]+w[k1*(Nh2+1)*Lx+(2*j+1)*Lx+i]*g2l[k2*Lx*(N+1)+i*(N+1)+(4*j+2)]+w[k2*(Nh2+1)*Lx+(2*j+1)*Lx+i]*g2l[k1*Lx*(N+1)+i*(N+1)+(4*j+2)]+w[k1*(Nh2+1)*Lx+(2*j+1)*Lx+i]*w[k2*(Nh2+1)*Lx+(2*j+1)*Lx+i]*g3[i*(N+1)+(4*j+2)]);
		    ytemp[i]=ycur[i]+h4*kcur2[i]/2;

		}

   for (i=0;i<Lx;i++)
	    {
			kcur3[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(4*j+2)]+w[k1*(Nh2+1)*Lx+(2*j+1)*Lx+i]*g2l[k2*Lx*(N+1)+i*(N+1)+(4*j+2)]+w[k2*(Nh2+1)*Lx+(2*j+1)*Lx+i]*g2l[k1*Lx*(N+1)+i*(N+1)+(4*j+2)]+w[k1*(Nh2+1)*Lx+(2*j+1)*Lx+i]*w[k2*(Nh2+1)*Lx+(2*j+1)*Lx+i]*g3[i*(N+1)+(4*j+2)]);
		    ytemp[i]=ycur[i]+h4*kcur3[i];

	    }


   for (i=0;i<Lx;i++)
	    {
			kcur4[i]=exp(theta[i])*(ytemp[i]*g2[i*(N+1)+(4*j+4)]+w[k1*(Nh2+1)*Lx+(2*j+2)*Lx+i]*g2l[k2*Lx*(N+1)+i*(N+1)+(4*j+4)]+w[k2*(Nh2+1)*Lx+(2*j+2)*Lx+i]*g2l[k1*Lx*(N+1)+i*(N+1)+(4*j+4)]+w[k1*(Nh2+1)*Lx+(2*j+2)*Lx+i]*w[k2*(Nh2+1)*Lx+(2*j+2)*Lx+i]*g3[i*(N+1)+(4*j+4)]);
			ycur[i]=ycur[i]+h4*(kcur1[i]+2*kcur2[i]+2*kcur3[i]+kcur4[i])/6;
            result[k1*(Mk+3)*Lx*(Nh+1)+k2*Lx*(Nh+1)+(j+1)*Lx+i]=ycur[i];
		}
  }

 }
}

free(g2);
free(g3);
free(g2l);
free(gtemp);
free(ycur);
free(ytemp);
free(kcur1);
free(kcur2);
free(kcur3);
free(kcur4);
free(betal);
}




/////////////////////////////////
/// natural-spline

void NCBasis(double *t, double *knot, int Lt, int Mk, double *result)
{

  int i,j;
  int Mkk;
  double temp;

  Mkk=Mk+3;

  for(i=0;i<Lt;i++)
  {
	  temp=t[i];
	  result[i*Mkk]=temp;
	  result[i*Mkk+1]=temp*temp;
	  result[i*Mkk+2]=temp*temp*temp;

	 for(j=0;j<Mk;j++)
	 {
	   temp=t[i]-knot[j];
	   temp=Pos(temp);
	   result[i*Mkk+(j+3)]=temp*temp*temp;

	 }
  }


}


///natural-spline derivative



void NPCBasis(double *t, double *knot, int Lt, int Mk, double *result)
{

  int i,j;
  int Mkk;
  double temp;

  Mkk=Mk+3;

  for(i=0;i<Lt;i++)
  {
	  temp=t[i];
	  result[i*Mkk]=1;
	  result[i*Mkk+1]=2*temp;
	  result[i*Mkk+2]=3*temp*temp;

	 for(j=0;j<Mk;j++)
	 {
	   temp=t[i]-knot[j];
	   temp=Pos(temp);
	   result[i*Mkk+(j+3)]=3*temp*temp;

	 }
  }

}



///natural-spline second derivative

void NHCBasis(double *t, double *knot, int Lt, int Mk, double *result)
{

 int i,j;
  int Mkk;
  double temp;

  Mkk=Mk+3;

  for(i=0;i<Lt;i++)
  {
	  temp=t[i];
	  result[i*Mkk]=0;
	  result[i*Mkk+1]=2;
	  result[i*Mkk+2]=6*temp;

	 for(j=0;j<Mk;j++)
	 {
	   temp=t[i]-knot[j];
	   temp=Pos(temp);
	   result[i*Mkk+(j+3)]=6*temp;

	 }
  }

}

//// link function g and its derivative
void gOneN(double *t, double *beta, double *knot, int Lt, int Mk, double *result)
{

double *eval,temp;
int i,j;
int Mkk;

Mkk=Mk+3;
eval=(double *) malloc (Lt*Mkk*sizeof(double));

NCBasis(t, knot, Lt, Mk, eval);

 for (i=0;i<Lt;i++)
 {
	 temp=0;
   for (j=0;j<Mkk;j++)
   {
	  temp=temp+eval[i*Mkk+j]*beta[j];

   }
  result[i]=temp;
  }

 free(eval);
}


void gTwoN(double *t, double *beta, double *knot, int Lt, int Mk, double *result)
{

double *eval,temp;
int i,j;
int Mkk;

Mkk=Mk+3;
eval=(double *) malloc (Lt*Mkk*sizeof(double));

NPCBasis(t, knot, Lt, Mk, eval);

 for (i=0;i<Lt;i++)
 {
	 temp=0;
   for (j=0;j<Mkk;j++)
   {
	  temp=temp+eval[i*Mkk+j]*beta[j];

   }
  result[i]=temp;
  }

  free(eval);
}



void gThreeN(double *t, double *beta, double *knot, int Lt, int Mk, double *result)
{

double *eval,temp;
int i,j;
int Mkk;

Mkk=Mk+3;
eval=(double *) malloc (Lt*Mkk*sizeof(double));

NHCBasis(t, knot, Lt, Mk, eval);

 for (i=0;i<Lt;i++)
 {
	 temp=0;
   for (j=0;j<Mkk;j++)
   {
	  temp=temp+eval[i*Mkk+j]*beta[j];

   }
  result[i]=temp;
  }

  free(eval);
}



////auxiliary functions
double Pos(double x)
{
	double result;

	result=x;
	if(x<0)
	{
		result=0;
	}

	return(result);

}



