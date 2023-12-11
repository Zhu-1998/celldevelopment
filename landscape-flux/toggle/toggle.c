
#include <math.h>
#include <stdlib.h>
#include "stdio.h"

#define TINY 1.0e-20;
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define Dim 2
#define La 200
#define Lb 200
#define XMAX 6.0
#define DC 0.1

FILE *fp1,*fp2;
char st[20],st1[20],st2[20];
main()
{
	int i,j,m,k,index1,index2,A,B,filenum;
	int grid[Dim],p_num[La*Lb]={0},p[La][Lb]={0};
	double s=0,eps=1.0e-4,f_cyc[La][Lb]={0},f_cdk1[La][Lb]={0};
	double fvec[Dim]={0},x[Dim]={0.0,0.0},xc[2],xi[2],xf[2];
	void rk2(double h,double x[],double fvec[],double xmin[],double xmax[],long *point);
	double x_min[Dim],x_max[Dim],tau,tau0;
	double h=0.001,x_old[Dim];
	double iter=1.0,Tn=2.0e7,Tni=1.0e5,Tnf,n_write=1.0e4;
	long   r=1,*point;
	Tnf = Tn+Tni;
	point=&r;
	for (i=0;i<Dim;i++) {
		x_max[i]=XMAX;
		x_min[i]=0.0;
		x[i]=1.0;
		}
for(filenum=1;filenum<2;filenum++)
{
	for(A=0;A<La;A++)
	{
		for(B=0;B<Lb;B++)
		{
			p[A][B]= 0.0;
			f_cyc[A][B] = 0.0;
			f_cdk1[A][B] = 0.0;
		}
	}
	//a1 = 0.5 + (filenum-1.0)*0.5;
	index1 = 1;
	index2 = 2;
	printf("filenum:%d\n",filenum);
	//sprintf(st1,"bistable_a%0.1f.txt",a);
	sprintf(st1,"bistable_Xm%0.1f.txt",XMAX);
	if((fp1=fopen(st1,"w+"))==NULL){
		printf("Cannot open file. \n");exit(0);}
	tau = 0;
	tau0 = 0;
	j = 1;
	for(iter=1.0;iter<=Tnf;iter=iter+1)
	{
		//for (i=0;i<Dim;i++) x_old[i]=x[i]; 
		if(iter==Tni+1)
		{
			fprintf(fp1,"%0.0f	%f",iter-Tni-1,tau0);
			//printf("%e %f",iter-Tni-1,tau0);
			for (i=0;i<Dim;i++){
				fprintf(fp1,"	%f",x[i]);
				//printf("	%f",x[i]);
			}
			fprintf(fp1,"\n");
			//printf("\n");	
		}
		rk2(h,x,fvec,x_min,x_max,point);	
		if(iter>Tni&iter<=Tni+n_write)
		{
			tau0=tau0+h;
			fprintf(fp1,"%0.0f	%f",iter-Tni,tau0);
			//printf("%e %f",iter-Tni,tau0);
			for (i=0;i<Dim;i++){
				fprintf(fp1,"	%f",x[i]);
				//printf("	%f",x[i]);
			}
			fprintf(fp1,"\n");
			//printf("\n");
			//getchar();
		}
		if(j%1000000==0){
			printf("iter=%e\n",iter);
			j=0;
		}
		j++;
		xc[0]=x[index1-1]; xi[0]=x_min[index1-1]; xf[0]=x_max[index1-1];
		xc[1]=x[index2-1]; xi[1]=x_min[index2-1]; xf[1]=x_max[index2-1];
		if(iter>=Tni) 
		{
			A = (int)((xc[0]-xi[0])*La/(xf[0]-xi[0]));
			B = (int)((xc[1]-xi[1])*Lb/(xf[1]-xi[1]));
			p[A][B] = p[A][B] + 1;
			f_cyc[A][B] = f_cyc[A][B] + fvec[index1-1];
			f_cdk1[A][B] = f_cdk1[A][B] + fvec[index2-1];
		}
		tau=tau+h;
	}
		//sprintf(st2,"p1_a%0.1f.txt",a);
		sprintf(st2,"landscape.txt",XMAX);
		if ((fp2=fopen(st2,"w+"))==NULL){
			printf("Cannot open file. \n");exit(0);}
		for(A=0;A<La;A++){
			for(B=0;B<La;B++){
					if(p[A][B]==0) fprintf(fp2,"%d	%d	0	0.0	0.0\n",A,B);
					else fprintf(fp2,"%d	%d	%d	%e	%e\n",A,B,p[A][B],f_cyc[A][B]/p[A][B],f_cdk1[A][B]/p[A][B]);
			}
		}
}
	fclose(fp1);
	fclose(fp2);

}

void force(double xp[],double fvec[])
{

	int i;
	double x[Dim];

	for (i=0;i<Dim;i++){
		x[i]=xp[i];		}
	fvec[0]=5/(1+pow(x[1],2))-x[0];
	fvec[1]=5/(1+pow(x[0],2))-x[1];

}

void rk2(double h,double x[],double fvec[],double xmin[],double xmax[],long *point)
{
	int i,j,k,indx[Dim];
	//long   r=1,*point;
	double xh[Dim],fxh[Dim],sqrt_h;
	double diff[Dim],gx[Dim],k1[Dim],xn[Dim],noise[Dim];
	void force(double x[],double fvec[]);
	double gasdev(long *idum);
	//point=&r;
	for (i=0;i<Dim;i++) {
		diff[i]=DC;
		gx[i]=sqrt(2*diff[i]);
		}	
	//while(1)
	//{
		sqrt_h = sqrt(h);
		force(x,fvec);
		for(i=0;i<Dim;i++)
		{
			noise[i]=sqrt_h*gx[i]*gasdev(point);
			k1[i]=h*fvec[i]+noise[i];
			xh[i]=x[i]+k1[i];
		}
		force(xh,fxh);
		for(i=0;i<Dim;i++){
			xn[i]=x[i]+0.5*(k1[i]+h*fxh[i]+noise[i]);// Heun method
			if(xn[i]<xmin[i]) x[i]=2*xmin[i]-xn[i];//reflecting boundary condition
			else if(xn[i]>xmax[i]) x[i]=2*xmax[i]-xn[i]; 
			else x[i]=xn[i];
		}	
	//}
}

double gasdev(long *idum)
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
