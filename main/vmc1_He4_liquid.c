/* Esempio 5:                         
                               
 VMC per 4^He liquido:
                       - potenziale di LJ
                       - funzione d' onda Jastrow-McMillan                                
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


const int npmax=500;
const int nest=5;
const double KB=1.;  /* Energie in Kelvin, lunghezze in Angstrom */
const double hbo2m=6.0599278; /* hbar^2/2m (dipende dalla massa!) */

double delta,rho,sigma_lj,eps_lj,l,li,l2;
double Pi,acc,vbox,b,tpb,tjf;
int npart,nstack,ncx,ncy,ncz;


struct particle
{ 
  double x,y,z;
};

struct table
{
  double point[1000];
};



double total_potential(),pair_wf(),pair_pot(),vlj(),vmod();
double kinetic_energy(),vtail(),rgauss(),u();
struct table grnow(),addtable();  

struct table addtable(struct table table1, struct table table2)
{
  int i;
  struct table temp;
  for(i=0;i<nstack;i++)
    {
      temp.point[i]=table1.point[i]+table2.point[i];
    };
  return temp;
}

copy_points(struct particle a[npmax],struct particle b[npmax])
{
  int i;
  for (i=0;i<npmax;i++)
    {
      b[i].x=a[i].x;
      b[i].y=a[i].y;
      b[i].z=a[i].z;

    };
}

    
main() 
{
  struct particle walker[npmax];    
  struct table gr,vcorr;
  double etot,vout;
  double dr,l22;
  double dest[nest],estave[nest],estave2[nest];
  double estsum[nest],estsum2[nest],estnow[nest],estvar[nest];  /* stimatori */
  int idest[nest];
  double vbulk,accum;
  double gr_norm,r,norm,rout,rin,tau,cvv;
  float delta_in,rho_in,b_in;
  int nstep1,nstep2;
  int i,j,k;
  FILE *out;

  out=fopen("gr.out","w+");
  Pi=4.*atan(1.);

  printf("Fluido di Lennard-Jones in 3D: \n");


  sigma_lj=2.556;
  eps_lj=10.4;
  printf("Numero di particelle: ");
  scanf("%i",&npart);
   printf("Numero di celle fcc nelle tre direzioni: ");
  scanf("%i %i %i",&ncx,&ncy,&ncz);
  printf("Densita:");
  scanf("%f",&rho_in);
  rho=rho_in;
  printf("Parametro variazionale b: ");
  scanf("%f",&b_in);
  b=b_in;
  printf("Time-step: ");
  scanf("%f",&delta_in);
  delta=delta_in;
  printf("Numero di passi di equilibratura: ");
  scanf("%i",&nstep1);  
  printf("Numero di passi: ");
  scanf("%i",&nstep2);  
  printf("Numero di punti per g(r): ");
  scanf("%i",&nstack);  



  l=pow((double)npart/rho,1./3.);
  li=1./l;
  l2=0.5*l;


  printf ("Lato della scatola di simulazione %10.8f \n",l);
  vbox=0.;
  l22=l2*l2-1.e-10;
  vbox=vlj(l22);
  vout=vtail(l);
  vbulk=Pi*vbox*(npart-1)/12.; 
  printf ("Correzione di coda: %10.5e\n",vout);
  printf ("Potenziale a L/2: %10.5e %10.5e\n",vbox,vlj(l2*l2));
  printf ("Correzione interna: %10.5e\n",vbulk);

  srand(1);
  get_initial_positions(walker);
  accum=0.;
  for(i=0;i<nstep1;i++)    /* loop principale */
    {
      advance(walker);
      accum=accum+acc;
      if((i+1)%100==0) printf("Step: %6i  Acc: %10.5f \n",i+1,accum/(npart*(i+1)));
    };
  printf("\n Fine della fase di equilibratura \n");
  for(j=0;j<nest;j++)
    {
      estsum[j]=0.;
      estsum2[j]=0.;
    };
  printf("  Step      Vnow      Vave       Verr       Tave        Terr       Eave       Evar\n");
  for(i=0;i<201;i++) 
  {
	  gr.point[i]=0;
	  vcorr.point[i]=0.;
  };
  for(i=0;i<nstep2;i++)   
    {
      advance(walker);
      tpsi(walker);
      estnow[0]=total_potential(walker)/(double)npart+vout+vbulk;
      estnow[1]=estnow[0]+tpb/(double)npart;
      estnow[3]=tpb/(double)npart;
      estnow[4]=tjf/(double)npart;
      estnow[2]=acc/(double)npart;
      gr=addtable(gr,grnow(walker));
      for(j=0;j<nest;j++)
	{
	  estsum[j]=estsum[j]+estnow[j];
	  estsum2[j]=estsum2[j]+estnow[j]*estnow[j];
	  estave[j]=estsum[j]/(i+1);
	  estave2[j]=estsum2[j]/(i+1);
	  estvar[j]=fabs(estave[j]*estave[j]-estave2[j]);
	  dest[j]=sqrt(estvar[j]/(double)(i+1));
	  idest[j]=dest[j]*1.e5;
	};
      if((i+1)%100==0)\
       printf("%6i %10.5f %10.5f (%5i) %10.5f %10.5f %10.5f %10.5f\n",\
               (i+1),estnow[1],estave[1],idest[1],estave[3],estave[4],estave[0],estave[2]);
    };
  for(j=0;j<nest;j++)
	{
	  estave[j]=estsum[j]/nstep2;
	  estave2[j]=estsum2[j]/nstep2;
	  dest[j]=sqrt(fabs(estave[j]*estave[j]-estave2[j])/(double)nstep2);
	};
  printf("Energia potenziale media: %10.8f +-  %10.8f\n",\
                 estave[0],dest[0]);
  printf("Energia cinetica PB media: %10.8f +-  %10.8f\n",estave[3],dest[3]);
  printf("Energia cinetica JF media: %10.8f +-  %10.8f\n",estave[4],dest[4]);
  printf("Energia totale media: %10.8f +-  %10.8f\n",estave[1],dest[1]);
  printf("Varianza energia totale: %10.8f\n",estvar[4]);

  for(i=0;i<nstack;i++) 
   {
     dr=l*0.5/(double)(nstack+1); /* g(r) */
     r=(i+0.5)*dr;
     rin=i*dr;
     rout=(i+1)*dr;
     norm=4./3.*Pi*(rout*rout*rout-rin*rin*rin)*rho;
     gr_norm=(double)gr.point[i]/((double)nstep2*norm*npart);
     fprintf(out,"%10.5f %10.5f\n",r,gr_norm);
   };
}

get_initial_positions(struct particle walker[npmax])
{
  struct particle trial_position;

/*  int ncx=4,ncy=4,ncz=4; */

  
  double basisx[4]={0.,0.5,0.0,0.5};
  double basisy[4]={0.,0.5,0.5,0.0};
  double basisz[4]={0.,0.0,0.5,0.5};


  int i,j,k,ll;
  double side,xd,yd,zd;
  int index; 
  

  side=pow(4./rho,1./3.);

  printf("elementary cell side: %10.5f\n",side);
  
  index=0;
  for(i=0;i<ncx;i++)
    {
      for(j=0;j<ncy;j++)
	{
	  for(k=0;k<ncz;k++)
	    {
	      for(ll=0;ll<4;ll++)
		{
		  walker[index].x=((double)i+basisx[ll])*side;
		  walker[index].y=((double)j+basisy[ll])*side;
	          walker[index].z=((double)k+basisz[ll])*side;
                  index++;
		};
	    };
	};
    };
  if(index!=npart)
    {
      printf("ERROR: INCONSISTENT NUMBER OF PARTICLES  %5i  %5i",index,npart);
      abort();
    };
  for(j=0;j<npart;j++)
    {      
      xd=walker[j].x;               
      walker[j].x=xd-l*rint(xd*li);
      yd=walker[j].y;
      walker[j].y=yd-l*rint(yd*li);
      zd=walker[j].z;
      walker[j].z=zd-l*rint(zd*li);
      printf("%10.5f  %10.5f   %10.5f\n",walker[j].x,walker[j].y,walker[j].z);
    };      
}

advance(struct particle walker[npmax])

{
  struct particle new;
  int i,j;
  double dx,dy,dz,arg;
  double uo,un,p,csi;
  acc=0.;
  for(j=0;j<npart;j++)
    {
      uo=pair_wf(j,walker[j].x,walker[j].y,walker[j].z,walker);
      dx=delta*(0.5-(double)rand()/(double)RAND_MAX);
      new.x=walker[j].x+dx;
      new.x=new.x-l*rint(new.x*li);
      dy=delta*(0.5-(double)rand()/(double)RAND_MAX);
      new.y=walker[j].y+dy;
      new.y=new.y-l*rint(new.y*li);
      dz=delta*(0.5-(double)rand()/(double)RAND_MAX);
      new.z=walker[j].z+dz;
      new.z=new.z-l*rint(new.z*li);
      un=pair_wf(j,new.x,new.y,new.z,walker);
      arg=un-uo;
      p=exp(-arg);
      csi=(double)rand()/(double)RAND_MAX;

      if(p>csi)
	{
	  acc=acc+1.;
          walker[j].x=new.x;
          walker[j].y=new.y;
          walker[j].z=new.z;
	};    
      /*      printf("arg: %10.5f p: %10.5f acc: %10.5f \n",arg,p,acc); */ 
    };
}



double total_potential(struct particle walker[npmax])
{
  int i,j;
  double dx,dy,dz,rr,v;
  v=0.;
  for(i=0;i<npart;i++)
    {
      for(j=0;j<i;j++)
	{
	  dx=walker[i].x-walker[j].x;
	  dx=dx-l*rint(dx*li);
      	  dy=walker[i].y-walker[j].y;
	  dy=dy-l*rint(dy*li);
	  dz=walker[i].z-walker[j].z;
	  dz=dz-l*rint(dz*li);
	  rr=dx*dx+dy*dy+dz*dz;
          v=v+vlj(rr);
	};
    };
  return v;
}

double pair_wf(int i,double x, double y, double z, struct particle walker[npmax])
{
  int j;
  double dx,dy,dz,rr,v;
  
  v=0.;
      for(j=0;j<i;j++)
	{
	  dx=x-walker[j].x;
	  dx=dx-l*rint(dx*li);
      	  dy=y-walker[j].y;
	  dy=dy-l*rint(dy*li);
	  dz=z-walker[j].z;
          dz=dz-l*rint(dz*li);
	  rr=dx*dx+dy*dy+dz*dz;
	  v=v+u(rr);
	};
      for(j=i+1;j<npart;j++)
	{
	  dx=x-walker[j].x;
          dx=dx-l*rint(dx*li);
      	  dy=y-walker[j].y;
	  dy=dy-l*rint(dy*li);
	  dz=z-walker[j].z;
	  dz=dz-l*rint(dz*li);
	  rr=dx*dx+dy*dy+dz*dz;
	  v=v+u(rr);
	};
  return v;
}
	      
double vlj(double r2)
{
  double r2i,r6i,r12i,v,l22;
  double rcore;

  l22=l2*l2;
  if(r2>l22)
    {
      v=0;
      return v;
    };
  
  rcore=0.3*sigma_lj;    

  if(r2<rcore*rcore) r2=rcore*rcore;
    r2i=sigma_lj*sigma_lj/r2;
    r6i=r2i*r2i*r2i;
    r12i=r6i*r6i;
    v=4.*eps_lj*(r12i-r6i)-vbox;
    return v;
}
  
double u(double r2)
{
  double ri,r2i,r5i,r,v,l22;
  double rcore;

  l22=l2*l2;
  if(r2>l22)
    {
      v=0;
      return v;
    };
  
  rcore=0.3*sigma_lj;    

  if(r2<rcore*rcore) r2=rcore*rcore;
    ri=b/sqrt(r2);
    r2i=b*b/r2;
    r5i=r2i*r2i*ri;
    v=r5i;
    return v;
}  

tpsi(struct particle walker[npmax])
{
  struct particle dpsi[npmax];
  double d2psi;
  int i,j;
  double dx,dy,dz,r2,v,rcore;
  double l22,dux,duy,duz,sr2i,r2i,r6i,ri;
  v=0.;
  for(i=0;i<npart;i++)
    {
      dpsi[i].x=0.;
      dpsi[i].y=0.;
      dpsi[i].z=0.;
    };
  d2psi=0.;
  for(i=0;i<npart;i++)  /* calcola grad(Psi)/Psi */
    {
      for(j=0;j<i;j++)
	{
	  dx=walker[i].x-walker[j].x;
	  dx=dx-l*rint(dx*li);
	  dy=walker[i].y-walker[j].y;
	  dy=dy-l*rint(dy*li);
	  dz=walker[i].z-walker[j].z;
	  dz=dz-l*rint(dz*li);
	  r2=dx*dx+dy*dy+dz*dz;
	  l22=l2*l2;
	  if(r2<l22)
	    {
	      rcore=0.3*sigma_lj;    
	      if(r2<rcore*rcore) r2=rcore*rcore;
	      r2i=b*b/r2;
              ri=1./(b*sqrt(r2));
	      r6i=r2i*r2i*r2i;
	      v=-5.*r6i*ri;  /* calcola u'(r)/r */
	      dux=v*dx;          /* calcola \vec(f_ij)=-u'(r_ij)/r_ij*\vec(r_ij) */
	      duy=v*dy;
	      duz=v*dz;
	      dpsi[i].x=dpsi[i].x+dux; 
	      dpsi[i].y=dpsi[i].y+duy;
	      dpsi[i].z=dpsi[i].z+duz;
	      dpsi[j].x=dpsi[j].x-dux; /* usiamo grad(u)_ij = -gradu_ji */
	      dpsi[j].y=dpsi[j].y-duy;
	      dpsi[j].z=dpsi[j].z-duz;
	      d2psi=d2psi+20.*r6i*ri;
	    }
	   
	};
    };
  tpb=0.;
  tjf=0.;
  for(i=0;i<npart;i++)
    {
      tpb=tpb+0.25*(dpsi[i].x*dpsi[i].x);
      tpb=tpb+0.25*(dpsi[i].y*dpsi[i].y);
      tpb=tpb+0.25*(dpsi[i].z*dpsi[i].z);
    };
  tpb=hbo2m*(d2psi-tpb);
  tjf=0.5*hbo2m*d2psi;
}  

struct table grnow(struct particle walker[npmax])
{
  int i,j;
  double dx,dy,dz,rr;
  double r,dh;
  int idx;
  struct table gg;
  dh=l2/(double)(nstack+1);

  for(i=0;i<201;i++) gg.point[i]=0;
  for(i=0;i<npart;i++)
    {
      for(j=0;j<i;j++)
	{
	  dx=walker[i].x-walker[j].x;
	  dx=dx-l*rint(dx*li);
      	  dy=walker[i].y-walker[j].y;
	  dy=dy-l*rint(dy*li);
	  dz=walker[i].z-walker[j].z;
	  dz=dz-l*rint(dz*li);
	  rr=dx*dx+dy*dy+dz*dz;
          r=sqrt(rr);
          if(r<l2)
	    {
	      idx=(int)(r/dh);
	      gg.point[idx]=gg.point[idx]+2;
	    };
	};
    };
  return gg;
}

double vtail(double l)
{
  double aux,aux3,aux9,v;
  
  aux=2.*sigma_lj*li;
  aux3=aux*aux*aux;
  aux9=aux3*aux3*aux3;
  v=8./9.*Pi*eps_lj*rho*(aux9-3.*aux3);
  return v;
}


double rgauss()
{
  double eta1,eta2,x1;   /* inizializzazione */
    
  eta1=(double)rand()/(double)RAND_MAX;
  eta2=(double)rand()/(double)RAND_MAX;
  x1=sqrt(-2.*log(eta1))*cos(2.*Pi*eta2);
  return x1;
};


