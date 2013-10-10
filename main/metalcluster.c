#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double h,Umax,R0,pi,rho0,mu,uu[2000],u[100][2000],Vcoul[2000];
int n=1500;
double rmax=18.;

double num1(double E, double L)
     /* risolve
        du/dr=z
	dz/dr=(veff(r)-E)*u con condizioni al contorno:
        l' autovalore viene fissato risolvendo l' equazione
        F(E)=0 dove F=U(x_n)-U_0(x_n) 
     */
{
  int i,loop,iwrite=1;
  int good_enough;
  double u2,u0,u1,k0,k1,k2,r,veff,v1,aux;
  Umax=0;
  u0=0.;
  u1=pow(h,L+1);
  uu[1]=0;
  uu[2]=pow(h,L+1);
  k0=0.;
  r=h;
  if(r<R0) veff=2.*pi*rho0*(r*r/3.-R0*R0)+L*(L+1)/(2.*r*r)+Vcoul[1];
  if(r>=R0) veff=2.*pi*rho0*(-2.*R0*R0*R0/(r*3.))+L*(L+1)/(2.*r*r)+Vcoul[1];
  k1=2.*(E-veff);
  for (i=2;i <= n;i++)
    {
      r=h*i;
  	if(r<R0) veff=2.*pi*rho0*(r*r/3.-R0*R0)+L*(L+1)/(2.*r*r)+Vcoul[i];
 	 if(r>=R0) veff=2.*pi*rho0*(-2.*R0*R0*R0/(r*3.))+L*(L+1)/(2.*r*r)+Vcoul[i];
      k2=2.*(E-veff);
      /* Numerov algorythm */
      aux=1./(1.+h*h/12.*k2);
      u2=2.*(1-5.*h*h/12.*k1)*u1-(1.+h*h/12.*k0)*u0;
      u2=u2*aux;
      Umax=(Umax>u2?Umax:u2);
      uu[i]=u2;
      u0=u1;
      u1=u2;
      k0=k1;
      k1=k2;
    };
 return u2;
}


main(void)
{
 double temp,temp1,Eig,Eigd,Eigu,Ud,Uu,U,U0,LL,unorm,r;
 double  rho[2000];
 double deg,Eig0,rhoint,pnum; 
 double rs,nions,v,rr,vsum;
 double diff,de;
 double E[100];
 int i,j,jj,n,nmax,L,Lmax,nn,npr,nfill;
 int LLL[100],N[100],idx[100];
 FILE *out,*out1;

 printf(" COMPUTES THE SINGLE PARTICLE SPECTRUM FOR A JELLIUM SPHERE  \n\n");

 printf ("Enter maximum angular momentum: ");
 scanf("%i",&Lmax);
 printf ("Enter maximum principal number: ");
 scanf("%i",&nmax);
 printf ("Enter the Wigner-Seiz parameter: ");
 scanf("%lf",&rs);
 printf ("Enter the charge: ");
 scanf("%lf",&nions);
 printf ("Enter the number of filled levels: ");
 scanf("%i",&nfill);
 out=fopen("u.out","w+");
 U0=0.;
 n=1500;
 /* setup constants */
 pi=4.*atan(1.);
 rho0=3./(4.*pi*rs*rs*rs);
 R0=pow(3.*nions/(4.*pi*rho0),1./3.);
 Eig0=-2.*pi*rho0*R0*R0*0.99;
 h=rmax/n;
 printf("Sharp radius: %10.5f \n",R0);
 printf("Jellium density: %10.5f \n",rho0);

 de=0.001;
 for(i=0;i<n;i++)
 {
	 r=h*i;
	 rho[i]=rho0/(1+exp((r-R0)/0.3)); 
 };
/* Use the density for computing Coulomb potential */
   for(i=1;i<n;i++)
   {
	   r=h*i;
	   vsum =0.;
	   for(j=0;j<i;j++)
	   {
		   rr=h*j;
		   vsum = vsum + rr*rr*rho[j];
	   };
	   vsum = -vsum/r; 
	   for(j=i;j<n;j++)
	   {
		   rr=h*j;
		   vsum=vsum+rr*rho[i];
	   };
	   vsum = vsum*h;
	   Vcoul[i] = 4.*pi*vsum;
   };
		   
 j=0;
   for(L=0;L<Lmax+1;L++)
   {
     Eig=Eig0;
     npr=0;
 	for (nn=0;nn<nmax;nn++)
       	{
	 /* try to bracket the solution first */
	 /* get a first solution and check the sign of the difference */
	 LL=(double)L;
	 U=num1(Eig,LL);
	 temp1=copysign(1.,U-U0);
	 temp=temp1;
	 Eigd=Eig;
	 Ud=U-U0;

	 while(temp*temp1>0.)
	   {
	     Eig=Eig+de;
	     U=num1(Eig,LL);
	     temp=copysign(1.,U-U0);
	   };
	 Uu=U-U0;
	 Eigu=Eig;

	 /* Now use Newton-Raphson to find the solution */
	 i=0;
	 while(fabs(U-U0)/Umax>1.e-5 && i<10000)
	   {
	     diff = (Uu-Ud)/(Eigu-Eigd);
	     Eig = Eigu - Uu/diff;
	     U = num1(Eig,LL);
	     if(U-U0>0) 
	       {
		 Uu=U;
		 Eigu=Eig;
	       };
	     if(U-U0<0)
	       {
		 Ud=U;
		 Eigd=Eig;
	       };
	     i++;
	   };
	 if(i==10000) printf("WARNING: level not converged  %10i  \n",j);
	 E[j]=Eig;
	 N[j]=nn+1;
	 LLL[j]=L;
/* Normalizzazione */ 	 
	 unorm=0.;
     for(i=0;i<n;i++)
      	{
	 unorm=unorm+uu[i]*uu[i];
   	};
         unorm=12.56637061*unorm*h;
	 unorm=sqrt(unorm);
     for(i=1;i<n;i++)
      	{
	 u[j][i]=uu[i]/unorm;
	};
        j++;
        Eig=Eig*0.99;
       };
   };

 
 for(i=0;i<100;i++) idx[i]=i;
 for(i=nmax*(Lmax+1);i<100;i++) E[i]=10.;
   sort(E,idx); 
     rhoint=0.;
     pnum=0.;
     for(i=1;i<n;i++) rho[i]=0.; 
     printf("%5i \n",nfill);
     for(i=0;i<n;i++) rho[i]=0.;
     for(j=0;j<nfill;j++)
     {
     jj=idx[j];
     for(i=1;i<n;i++)
      	{
	 r=h*i;
	 deg=4.*LLL[jj]+2.;
	 rho[i]=rho[i]+deg*u[jj][i]*u[jj][i]/r/r;
	 fprintf(out,"%10.5f    %10.5f   \n",r,u[jj][i]/r);
     	};
	 pnum=pnum+deg;
	 fprintf(out,"\n");
      };
     for(i=1;i<n;i++)
     {
	     r=h*i;
	     rhoint=rhoint+rho[i]*r*r;
     };
/* Use the density for computing Coulomb potential */
   for(i=1;i<n;i++)
   {
	   r=h*i;
	   vsum =0.;
	   for(j=0;j<i;j++)
	   {
		   rr=h*j;
		   vsum = vsum + rr*rr*rho[j];
	   };
	   vsum = -vsum/r; 
	   for(j=i;j<n;j++)
	   {
		   rr=h*j;
		   vsum=vsum+rr*rho[i];
	   };
	   vsum = vsum*h;
	   Vcoul[i] = 4.*pi*vsum;
   }
 fclose(out);
 for(j=0;j<nmax*(Lmax+1);j++) 
   {
     printf("N: %3i   L:  %3i    E: %10.5f \n",N[idx[j]],LLL[idx[j]],E[j]);
   };
 printf("Number of particles from integral: %10.5f   Exact: %10.5f \n",12.56637061*rhoint*h,pnum);
 out1=fopen("d.out","w+");
     for(i=1;i<n;i++)
      	{
	 r=h*i;
  	if(r<R0) v=2.*pi*rho0*(r*r/3.-R0*R0);
 	 if(r>=R0) v=2.*pi*rho0*(-2.*R0*R0*R0/(r*3.));
	 fprintf(out,"%10.5f    %10.5f     %10.5f     %10.5f \n",r,rho[i],v,Vcoul[i]);
   	};
}

sort(double tab[100],int idx[100])
{
  int i,j,iaux,iswap;
  double aux;
   
 

  iswap=1;
  j=0;
  while(iswap!=0)
    {
      iswap=0;
      for(i=0;i<99;i++)
	{
	  if(tab[i]>tab[i+1])
	    {
	      aux=tab[i];
	      tab[i]=tab[i+1];
	      tab[i+1]=aux;
	      iaux=idx[i];
	      idx[i]=idx[i+1];
	      idx[i+1]=iaux;
	      iswap=1;
	    };
	};
    };
}
