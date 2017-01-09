#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define TARGET          27
#define T               300
#define pi              3.1415926
#define Max_Atom        200
#define CMax_Atom        500

typedef struct
{
 double fitness;
 double  gen[Max_Atom][3];
 double  cgen[CMax_Atom][3];
 double ep;
}ga_struct;


void rand_init(int natm,double *x, double boxl)

{
double RR=0.0;
int i=0;
int ii=0;
int j=0;
  for ( i=0;i<natm;i++)
    {
      for ( j=0;j<3;j++)
        {
          *(x+i*3+j)=(rand()/(RAND_MAX*1.0)*2-1)*boxl/2.0;
        }

      for( ii=0;ii<i;ii++)
        {
          RR=0.0;
          RR=   pow(*(x+i*3+0)-*(x+(ii)*3+0),2);
          RR=RR+pow(*(x+i*3+1)-*(x+(ii)*3+1),2);
          RR=RR+pow(*(x+i*3+2)-*(x+(ii)*3+2),2);
//          printf ("ii=%d,i=%d,RR=%f\n",ii,i,RR);
          if (RR<4)
            {i=i-1;
//          printf ("k=%d,i=%d,RR=%f\n",k,i,RR);
            }
         }
   }
}




int cmp(const void *a,const void *b)
{
      return ((double *)b)[0] >((double*)a)[0] ? 1:-1;
}


//Function: center: centerize the geometryand rotate around the center

void
center(int natm ,double a[Max_Atom][3])
{
 double x0=0,y0=0,z0=0;
 double xx=0,yy=0,zz=0,xxold=0;
// random rotation angle generator, rotate about y axis and z axis
//
 double  dtheta1=rand()%360;
 double  dtheta2=rand()%360;
         dtheta1=dtheta1/360.0*2*pi;
         dtheta2=dtheta2/360.0*2*pi;




 int i=0;
// calculation of center coordination x0,y0,z0
  for(i=0;i<natm;i++)
    {
      x0=a[i][0]+x0;
      y0=a[i][1]+y0;
      z0=a[i][2]+z0;
    }
  x0=x0/natm;
  y0=y0/natm;
  z0=z0/natm;

// translation of center to 0,0,0 and rotation
//
  for(i=0;i<natm;i++)
    {   
      xx=a[i][0]-x0;
      yy=a[i][1]-y0;
      zz=a[i][2]-z0;
      xxold=xx;
      xx=xx*cos(dtheta1)-yy*sin(dtheta1);
      yy=xxold*sin(dtheta1)+yy*cos(dtheta1);
      xxold=xx;
      xx=xx*cos(dtheta2)-zz*sin(dtheta2);
      zz=xxold*sin(dtheta2)+zz*cos(dtheta2);
      a[i][0]=xx;
      a[i][1]=yy;
      a[i][2]=zz;
//    printf("after  rotation 1  %12.4f %12.4f %12.4f\n",a[i][0],a[i][1],a[i][2]);
    } 

//Function: mate:

}


void
shift(int natm ,double a[Max_Atom][3],double ddptc)
{
 double zmin=0;
 int i=0;
//initialize geometry center 

  for(i=0;i<natm;i++)
   {
   
   if(zmin>a[i][2]) {zmin=a[i][2];} 
//   printf("before center 1  %12.4f %12.4f %12.4f\n",a[i][0],a[i][1],a[i][2]);
   }
 
double  dgap=ddptc-(zmin-0); //0 is the z coord of graphne sheet 
   
// translation of cluster 
//
  for(i=0;i<natm;i++)
   {
    a[i][2]=a[i][2]+dgap;//dgap is the atuall distance. ddptc is required distance
   }

}


void
init_carbon_i(int cnatm,int POPSIZE, ga_struct *beta_population,int i)
{

double *cx = (double *) malloc(3*cnatm*sizeof(double));
int j=0;
FILE  *fp_input=fopen("c.coord","r");
read_coord(fp_input,cnatm,cx);
  for (j=0;j<cnatm;j++)
    {
      beta_population[i].cgen[j][0]=(double)(cx[j*3+0]);
      beta_population[i].cgen[j][1]=(double)(cx[j*3+1]);
      beta_population[i].cgen[j][2]=(double)(cx[j*3+2]);
    }
free (cx);
fclose(fp_input);
  
}


void
init_carbon(int cnatm,int POPSIZE, ga_struct *beta_population)
{
double *cx = (double *) malloc(3*cnatm*sizeof(double));
int i=0,j=0;
FILE  *fp_input=fopen("c.coord","r");
read_coord(fp_input,cnatm,cx);
for(i=0;i<POPSIZE;i++)
  {
    for (j=0;j<cnatm;j++)
      {
        beta_population[i].cgen[j][0]=(double)(cx[j*3+0]);
        beta_population[i].cgen[j][1]=(double)(cx[j*3+1]);
        beta_population[i].cgen[j][2]=(double)(cx[j*3+2]);
  
      }
  }
free(cx);
fclose(fp_input);
  
}



void
init_population(int ptnatm, int cnatm, int POPSIZE,ga_struct *population, ga_struct *beta_population,double boxl,int flag_res)
{
 double *x = (double *) malloc(3*ptnatm*sizeof(double));
 double *cx = (double *) malloc(3*cnatm*sizeof(double));
 double efinal=0;
 FILE  *fp_coordc=fopen("c.coord","r");
 FILE  *fp_coordpt=fopen("pt.coord","r");

 if (flag_res==1&&fp_coordpt==NULL)
    {printf("ERROR, can not open pt.coord!!\n");exit(0);}
//read coordinates of graphene from c.coord

 read_coord(fp_coordc,cnatm,cx);
 int i=0;
 int j=0;


 for(i=0;i<POPSIZE;i++)
   {
     population[i].fitness=0;
     beta_population[i].fitness=0;

//reading from first generation
//
//read_coord(fp_input2,natm,x);
     if(flag_res==0)
        {rand_init(ptnatm,x,boxl);}
     else if(flag_res==1)
        {
          read_coord(fp_coordpt,ptnatm,x);
        }


     for (j=0;j<ptnatm;j++)
       {
//                printf("%d %lf %lf %lf \n", skip,x,y,z);
		
         population[i].gen[j][0]=(double)(x[j*3+0]);
         population[i].gen[j][1]=(double)(x[j*3+1]);
         population[i].gen[j][2]=(double)(x[j*3+2]);
         beta_population[i].gen[j][0]=0;
         beta_population[i].gen[j][1]=0;
         beta_population[i].gen[j][2]=0;

        }

      for (j=0;j<cnatm;j++)
        {
          population[i].cgen[j][0]=(double)(cx[j*3+0]);
          population[i].cgen[j][1]=(double)(cx[j*3+1]);
          population[i].cgen[j][2]=(double)(cx[j*3+2]);
          beta_population[i].cgen[j][0]=(double)(cx[j*3+0]);
          beta_population[i].cgen[j][1]=(double)(cx[j*3+1]);
          beta_population[i].cgen[j][2]=(double)(cx[j*3+2]);
         }

//       center(natm,population[i].gen);  
          printf("%d %lf %lf %lf \n",i,population[i].cgen[cnatm-1][0],population[i].cgen[cnatm-1][1],population[i].cgen[cnatm-1][2]);
          printf("%d %lf %lf %lf \n",i,population[i].gen[0][0],population[i].gen[ptnatm-1][1],population[i].gen[ptnatm-1][2]);
     }
    free (x);
    free (cx);
    fclose(fp_coordc);
    if(fp_coordpt!=NULL){ fclose(fp_coordpt);}

}



void
cal_pop_energy(int POPSIZE,ga_struct *population,int ptnatom, int cnatom)
{
  double efinal=0;
  int i=0;

  for(i=0;i<POPSIZE;i++)
    { 
       cal_energy(&efinal,population[i].gen,population[i].cgen,0,0,ptnatom,  cnatom, 100);
       population[i].ep=efinal;
    }
 
}





void
normal_fitness(int POPSIZE,ga_struct *population) 
{      

double  mine=population[0].ep;
double  maxe=population[POPSIZE-1].ep;
int     i=0;

    for (i=0;i<POPSIZE;i++)
      {        
        if(mine==maxe)
          {
         population[i].fitness=1.0;
          } 
        else
          {
          population[i].fitness=1-0.7*(population[i].ep-mine)/(maxe-mine);
        
          }
      }

}
 

int
sort_func (const void *e1,const void *e2)

{
	return (((ga_struct *) e1)->ep > ((ga_struct *)e2)->ep ? 1:-1);
}



// function count: count the number of atoms below or above xx plane
int
count_x(int natm, double a[Max_Atom][3],int xx)
{

int  i=0;

  for(i=0;i<natm;i++)
   {
	    if( a[i][0] >0)
		 {xx=xx+1;   
 		 }
   }
 return xx;
}

//Function: mate:


double min(double a,double b)
{
if(a>b){return b;}else{return a;}
}

void 
elitism(int esize,int ptnatm, int cnatm ,ga_struct *population,ga_struct *beta_population)
{
 int i=0,j=0;
 for (i=0;i<esize;i++)
   {
//      center(natm,population[i].gen);  
      beta_population[i].ep=population[i].ep;
      for(j=0;j<ptnatm;j++)
        {
          beta_population[i].gen[j][0]=population[i].gen[j][0];
          beta_population[i].gen[j][1]=population[i].gen[j][1];
          beta_population[i].gen[j][2]=population[i].gen[j][2];
        }
      for(j=0;j<cnatm;j++)
        {
          beta_population[i].cgen[j][0]=population[i].cgen[j][0];
          beta_population[i].cgen[j][1]=population[i].cgen[j][1];
          beta_population[i].cgen[j][2]=population[i].cgen[j][2];
         }
   }
}

void 

mutate ( int natm, double Mu,int POPSIZE,ga_struct *population)
{

 int rand_pop=0;
 int rand_atom=0;
 int rand_direct=0;
 int rand_N=0;
 double rand_l=0.0;
 int i=0,j=0;

 for(i=0;i<POPSIZE*Mu;i++)
   {

     rand_pop=rand()%POPSIZE;
     rand_atom=rand()%natm;
     rand_N=rand()%50;
     for(j=0;j<rand_N;j++)
       {
         rand_l=(rand()/(RAND_MAX*1.0)-0.5)*4;
         rand_direct=rand()%3;
  //    printf("randpop= %d, rand_atom= %d rand_N= %d, rand_direct= %d rand_l= %lf\n", rand_pop,rand_atom,rand_N,rand_direct,rand_l);
          population[rand_pop].gen[rand_atom][rand_direct]=population[rand_pop].gen[rand_atom][rand_direct]+rand_l;

        }  
    }

}

void 
mate( int ptnatm, int cnatm, int esize,int POPSIZE,ga_struct *population,ga_struct *beta_population, int Temp, double e_mate, double ddptc,int min_step)
{
 int    rand_index=0;
 int    randa=0;
 int    randb=0;
 int    min_index=0;
 double rand_p=0.0;
 double rand_ddptc=0.0;
 
 int    Na=0,Nb=0;
 double plandelta=0;
 double bulkd=1.0; //the distance between two bulk after join;
 double parentrate=0.2;
 double Ne_mate=e_mate*ptnatm;

 FILE  *fbad=fopen("bad.xyz","a+");
// all graphene in this generation will be reset to c.coord.


 int    i=0,j=0;
 for(i=esize;i<POPSIZE;i++)
   {
   
// select two parents 
     for (j=0;j<2;j++)
       {	
         do
           {
             rand_index=rand()%POPSIZE;
             rand_p=rand()/(RAND_MAX*1.0);
           } while(population[rand_index].fitness<rand_p);
         if(j==0){randa=rand_index;}
         if(j==1){randb=rand_index;}
              //   if(randa==randb){j=j-1;}
        }  
     center(ptnatm,population[randa].gen);  
     center(ptnatm,population[randb].gen);

     qsort(population[randb].gen,ptnatm,sizeof(population[randb].gen[0]),cmp);   
     qsort(population[randa].gen,ptnatm,sizeof(population[randa].gen[0]),cmp);     

     Na=count_x(ptnatm,population[randa].gen,0);
     Nb=ptnatm-Na;

     plandelta=population[randa].gen[Na-1][0]-population[randb].gen[Na][0];

     for (j=0;j<ptnatm;j++)
       {
         if(j<Na)
           {
             beta_population[i].gen[j][0]=population[randa].gen[j][0];
             beta_population[i].gen[j][1]=population[randa].gen[j][1];
             beta_population[i].gen[j][2]=population[randa].gen[j][2];
           } 
         else
           {
             beta_population[i].gen[j][0]=population[randb].gen[j][0]+plandelta-bulkd*plandelta/fabs(plandelta);
             beta_population[i].gen[j][1]=population[randb].gen[j][1];
             beta_population[i].gen[j][2]=population[randb].gen[j][2];
            }
        } 
     rand_ddptc=2-rand()/(RAND_MAX*1.0)*4;
     init_carbon_i(cnatm,POPSIZE, beta_population,i);
     shift(ptnatm,beta_population[i].gen,rand_ddptc);
//   shift(cnatm,beta_population[i].cgen,0);
//   printf("before calcualtion %12.4f %12.4f %12.4f\n",beta_population[i].cgen[0][0],beta_population[i].cgen[0][1],beta_population[i].cgen[0][2]); 
//   printf("before calcualtion %12.4f %12.4f %12.4f\n",beta_population[i].gen[0][0],beta_population[i].gen[0][1],beta_population[i].gen[0][2]); 
     cal_energy(&beta_population[i].ep,beta_population[i].gen,beta_population[i].cgen,1,Temp,ptnatm,cnatm,min_step);
//   write_coord(fbad,ptnatm,cnatm,beta_population[i].gen,beta_population[i].cgen);
//   printf("after calculation %12.4f %12.4f %12.4f\n",beta_population[i].cgen[0][0],beta_population[i].cgen[0][1],beta_population[i].cgen[0][2]); 
//   printf("after calcualtion %12.4f %12.4f %12.4f\n",beta_population[i].gen[0][0],beta_population[i].gen[0][1],beta_population[i].gen[0][2]); 
 
    
     if(population[randa].ep<population[randb].ep)
       {min_index=randa;}
     else
       {min_index=randb;}
          	    
  
     if(beta_population[i].ep-population[0].ep<ptnatm*0.1)
       {
         write_coord(fbad,ptnatm,cnatm,beta_population[i].gen,beta_population[i].cgen,beta_population[i].ep);
       }
//                if(beta_population[i].ep-0.0000001<population[POPSIZE-1].ep)
//                beta is the newpopulation  
     if(beta_population[i].ep-Ne_mate<population[min_index].ep)
       {
         printf(" Eold= %12.4lf %d %d SUCCESSFULLY MATED !!! ", population[min_index].ep,randa,randb);
//                   printf("old parant randa= %d randb= %d eranda= %lf, erandb= %lf newe=%lf\n",randa,randb,population[randa].ep,population[randb].ep,beta_population[i].ep );
       }
      else  
        {
          printf(" Eold= %12.4lf %d %d NOT          MATED !!!" ,population[min_index].ep,randa,randb);
//                   printf("old parant randa= %d randb= %d eranda= %lf, erandb= %lf newe=%lf\n",randa,randb,population[randa].ep,population[randb].ep,beta_population[i].ep );
          i=i-1;
   
        }
 }
 fclose(fbad);
}


 






void
swap ( ga_struct ** p1, ga_struct* * p2)
{
  ga_struct *tmp = *p1;
  *p1 = *p2;
  *p2 = tmp;
}


 void filter(int natm,int POPSIZE, int *  NEWPOPSIZE, double delte,ga_struct *population)
{
  int i,j,k;
 *NEWPOPSIZE=POPSIZE;
  for (i=1;i<POPSIZE;i++)
      {
       for (j=i+1;j<POPSIZE;j++)
        {
         if(fabs(population[j].ep-population[i].ep)<delte)
           
           {
             population[j].ep=population[0].ep;
               for (k=0;k<natm;k++)
                {
                population[j].gen[k][0]=population[0].gen[k][0];
                population[j].gen[k][1]=population[0].gen[k][1];
                population[j].gen[k][2]=population[0].gen[k][2];
                }
                  
//              population[j].ep=1000;
//           *NEWPOPSIZE=*NEWPOPSIZE-1;
           
           }
       } 
      }

   qsort (population, POPSIZE, sizeof(ga_struct),sort_func);

}


int  main()
{
 FILE     *fpenergy=fopen("./energy","w"); // Print the energy evolution
 FILE     *fprestart=fopen("./restart","w"); //print the generation and popsize coordinates
 FILE     *fpoptim=fopen("./optim.xyz","wb");
 FILE *fp=fopen("data.txt","r");
 char skip[10];
// FILE *fp_input = fopen("first_population","rb");
 FILE     *fp_input2 = fopen("ga_input2.5","r");
 int       NSTEP=0;
 int       step=0;
 double    ee_mate=0 ;       // the energy cretiria for accepting children cluster. the larger, the less restrict. can be 0.1 0 or -0.1 
 double    ELITRATE=0.2;
 int       POPSIZE=0;
 int       min_step=0;
 int       NEWPOPSIZE=POPSIZE;
 double    delte=0.00001 ; //0.00001;
 int       ptnatm,cnatm;
 double    boxl=0;
 double    Mu=0.2;
 double    dptc;
 int       glob=0,globconvg=20;
 double    globe[30];
 int       Temp; 
 srand(time(NULL));
 time_t    current_time;
 char*     c_time_string,c_time_final;
 double    seconds,start,end,seconds_new,seconds_old,seconds_total;
 struct    timespec now,tmstart;
 int       flag_res=0;
  //print time with noraml format
 clock_gettime(CLOCK_REALTIME, &tmstart);
 start = (double)clock() /(double) CLOCKS_PER_SEC;
 seconds_old=  (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9);
//initilization
 
 current_time = time(NULL);
 c_time_string = ctime(&current_time);
 printf("Current time is %s\n", c_time_string);





   
 fscanf(fp_input2,"%d %s\n", &ptnatm,skip);
 fscanf(fp_input2,"%d %s\n", &cnatm,skip);
 fscanf(fp_input2,"%d %s\n ", &POPSIZE,skip);
 fscanf(fp_input2,"%d %s\n", &NSTEP,skip);
 fscanf(fp_input2,"%d %s\n", &globconvg,skip);
 fscanf(fp_input2,"%lf %s\n",&ELITRATE,skip);
 fscanf(fp_input2,"%lf %s\n",&delte,skip);
 fscanf(fp_input2,"%d %s\n",&Temp,skip);
 fscanf(fp_input2,"%d %s\n",&min_step,skip);
 fscanf(fp_input2,"%lf %s\n",&ee_mate,skip);
 fscanf(fp_input2,"%lf %s\n",&dptc,skip);
 fscanf(fp_input2,"%lf %s\n",&boxl,skip);
 fscanf(fp_input2,"%d %s\n",&flag_res,skip);
 printf("********************JOB started*****************************\n");
 printf("********************JOB started*****************************\n");
 printf("********************JOB started*****************************\n\n\n");
 printf("Number of atoms %d\n", ptnatm);
 printf("Number of catoms %d\n", cnatm);
 printf("POPSIZE         %d\n",POPSIZE);
 printf("NSTEP           %d\n",NSTEP);
 printf("globconvg       %d\n",globconvg);
 printf("ELITRATE        %lf\n",ELITRATE);
 printf("delte           %lf\n",delte);
 printf("Temprature      %d\n",Temp);
 printf("minimiz   step  %d\n",min_step);
 printf("ee_mate         %lf\n",ee_mate);
 printf("dptc            %lf\n",dptc);
 printf("intial boxl     %lf\n",boxl);
 printf("reading from pt_coord?  %d\n",flag_res);
 printf("\n\n\n******** end reading input information ***********************\n\n\n");

 ga_struct *population = malloc(sizeof(ga_struct)*POPSIZE);
 ga_struct *beta_population = malloc(sizeof(ga_struct)*POPSIZE);

 init_population(ptnatm,cnatm,POPSIZE,population,beta_population,boxl,flag_res);
 cal_pop_energy(POPSIZE,population,ptnatm,cnatm);

int  i=0,j=0;
int  esize=POPSIZE*ELITRATE;
   for (i=0;i<POPSIZE;i++)
   {
   center(ptnatm,population[i].gen);  
   shift(ptnatm,population[i].gen,dptc);  
//   shift(cnatm,population[i].cgen,0);  
//    write_coord(fpoptim,natm,cnatm,population[i].gen,population[i].cgen);
//   write_coord(fpoptim,ptnatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
   } 

 for (step=0;step<NSTEP;step++)
 { 

      
   cal_pop_energy(POPSIZE,population,ptnatm,cnatm);
   
 printf("\n\n\n\n***********************************************\n");
 printf(  "***********************************************\n");
 printf(  "***********************************************\n");
   printf("Gen= %d starting optimization..................\n\n",step); 

   qsort (population, POPSIZE, sizeof(ga_struct),sort_func);

   normal_fitness(POPSIZE,population);

//             if(step>3)              
//             { 
//             filter(natm,POPSIZE,&NEWPOPSIZE,delte,population);
//             }

// print the coordinates and energy information of current generation into files 
          //"energy file" records the energy and fitness information

//          fprintf(fphistory,"GENERATION=  %d\n",step);
///               write_coord(fpoptim,natm,cnatm,population[0].gen,population[0].cgen);
 	  
   for(i=0;i<POPSIZE;i++)
   {
     fprintf(fpenergy,"%d num %d   %lf\n",step, i,population[i].ep);
     fflush(fpenergy);
//	       fprintf(fphistory,"POPULATION  %d energy %f\n",i,population[i].ep);
//	             for (j=0;j<natm;j++)
//			  { 
//			  fprintf(fphistory,"%d %12.4f %12.4f %12.4f\n",j,population[i].gen[j][0],population[i].gen[j][1],population[i].gen[j][2]); 
//			  }

//	       printf(" after  sort %lf\n",population[i].ep/10000);
   } 


// Preserve the first esize parentes from the previous generation. copy from population to beta_population
   printf("fabs %lf\n",fabs(population[0].ep-population[POPSIZE-1].ep));
   if(fabs(population[0].ep-population[POPSIZE-1].ep)<delte)
   {
    
     fprintf(fpenergy,"%d %lf  %lf global minimum \n",step,population[0].ep,population[0].fitness);
     globe[glob]=population[POPSIZE-1].ep;
     glob=glob+1;           
//       break;
   }

   if(glob>20)
   {
     if(fabs(globe[glob]-globe[glob-20])<delte)
     {
       fprintf(fpenergy,"%d %lf  %lf final global minimum \n",step,population[0].ep,population[0].fitness);
       break;
     }
   }


   esize=POPSIZE*ELITRATE;
   if (esize<1){esize=1;}
   elitism(esize,ptnatm, cnatm,population,beta_population);


// Generate the rest part of beta_generation by mating process

   mate(ptnatm,cnatm,esize,POPSIZE,population,beta_population,Temp,ee_mate,dptc,min_step);
              
//               mutate ( natm, Mu, POPSIZE, beta_population);

   swap(&population,&beta_population);
//   printf("%12.4f %12.4f %12.4f\n",population[0].cgen[0][0],population[0].cgen[0][1],population[0].cgen[0][2]); 
//   printf("%12.4f %12.4f %12.4f\n",population[1].cgen[0][0],population[1].cgen[0][1],population[1].cgen[0][2]); 
//   printf("%12.4f %12.4f %12.4f\n",beta_population[1].gen[0][0],beta_population[1].gen[0][1],beta_population[1].gen[0][2]); 
   fprestart=fopen("./restart","w");
   for (i=0;i<POPSIZE;i++)
   {
//    write_coord(fpoptim,natm,cnatm,population[i].gen,population[i].cgen);
   
      write_coord(fpoptim,ptnatm,cnatm,population[i].gen,population[i].cgen,population[i].ep);
      write_coord(fprestart,ptnatm,0,population[i].gen,population[i].cgen,population[i].ep);
     fflush(fpoptim);
     fflush(fprestart);
   } 
 fclose(fprestart);
    
   
   clock_gettime(CLOCK_REALTIME, &now);
   seconds_new = (double)((now.tv_sec+now.tv_nsec*1e-9));
   seconds=seconds_new-seconds_old;
   seconds_total = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
   printf("\nWall time for this generation is %lf s\n",seconds);
   printf("\nWall time totally  %lf s\n",seconds_total);
   seconds_old=seconds_new;
 }
 fclose(fpenergy);
 fclose(fpoptim);
//   fclose(fp_input );





 printf("\n********************JOB FINISHED*****************************\n");

// time information
 current_time = time(NULL);
 c_time_string = ctime(&current_time);
 printf("Current time is %s\n", c_time_string);
 
  // measure elapsed wall time
 clock_gettime(CLOCK_REALTIME, &now);
 seconds = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));
 printf("wall time %fs\n", seconds);

  // measure CPU  time
 end = (double)clock() / (double) CLOCKS_PER_SEC;
 printf("cpu time %fs\n", end - start);
 printf("\n********************JOB FINISHED*****************************\n");
 free(population);
 free(beta_population);

}
