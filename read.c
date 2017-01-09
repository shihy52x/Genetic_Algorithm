#include "stdio.h"
#include "stdlib.h"
#include "string.h"

void read_coord(FILE *fp1,int cnatm,double *x) 
{

 int i,j,k;
j=0;
 int skip,skip2;
 int error;
 char s[100];
       

        	fscanf(fp1,"%d\n",&cnatm);
        	fscanf(fp1,"%s\n",s);

	for (j=0;j<cnatm;j++)
		{
   		error=fscanf(fp1," %s %lf %lf  %lf\n",       s,    x+j*3+0,   x+j*3+1,    x+j*3+2)  ;
      //             printf("j=%d, error=%d  x=%lf \n",j,error,*(x+j*3));
                   if (error==-1)
                   {printf("\n\nNot enough data in pt.coord!!!\n\n");exit (0);}
		}
            // printf("%lf %lf %lf\n",*x,*(x+1),*(x+2));


}



void write_coord(FILE *fp2,int natm, int cnatm,double *x, double *y, double energy) 
{

 int i,j,k;
j=0;
 int skip,skip2;
 int nskip=12;
 char s[100];

               fprintf(fp2,"%d\n",natm+cnatm);
               fprintf(fp2,"Atoms.Energy:%lf\n",energy);
	for (j=0;j<natm;j++)
		{
                fprintf(fp2,"      %6s %12.4lf %12.4lf %12.4lf \n", "Pt",*( x+j*3+0),*(x+j*3+1), *(x+j*3+2) );
		}
	for (j=0;j<cnatm;j++)
		{
                fprintf(fp2,"      %6s %12.4lf %12.4lf %12.4lf \n", "C",*( y+j*3+0),*(y+j*3+1), *(y+j*3+2) );
		}


}





void write_his(FILE *fp2,int natm,double *x) 
{
 int i,j,k;
j=0;
 int skip,skip2;
 int nskip=12;
 char s[100];

               fprintf(fp2,"%d\n",natm);
               fprintf(fp2,"Atoms. Timestep: 0\n");
	for (j=0;j<natm;j++)
		{
                fprintf(fp2,"      %d %lf %lf %lf \n", 1,*( x+j*3+0),*(x+j*3+1), *(x+j*3+2) );
		}


}
