/* cppsayhello.cpp */
#include <iostream>
#include<iostream>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include <sstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string> // this should be already included in <sstream>


#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
#include "min.h"
#include "minimize.h"
#include "output.h"
#include "universe.h"
#include "update.h"
#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "library.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "input.h"
#include "lammps.h"
#include "mpi.h"
#include "min.h"
#include "minimize.h"
#include "output.h"
#include "universe.h"
#include "update.h"


using namespace LAMMPS_NS;

extern "C" void get_coord(void *ptr,double *coords);
extern "C" double get_energy(void *ptr,double *efinal);
extern "C" void scatter_natm(void *ptr, int natm);
extern "C" double  cal_energy( double *efinal, double *xx, double *yy ,int flagdym ,int Temp,int ptnatom, int cnatom);
//extern void lammps_get_unwrapped (void*, double*);
//extern void lammps_put_wrapped   (void*, double*);
//extern void lammps_get_force     (void*, double*);
//extern void lammps_get_type      (void*, int*);



double  cal_energy( double *efinal, double *xx, double *yy, int flagdym, int Temp,int ptnatom,int cnatom )
//flagdym==1 then do dynamics after and before optimization
{
 
  /* setup MPI and various communicators
     driver runs on all procs in MPI_COMM_WORLD
     comm_lammps only has 1st P procs (could be all or any subset) */

  char **arg;
  arg =  new char*[5];
  arg[0]=new char [20];
  arg[1]=new char [20];
  arg[2]=new char [20];
  int narg=3;
  strcpy(arg[1],"1");
  strcpy(arg[2],"in_Pt13.lmp");


  MPI_Init(&narg,&arg);


  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  int nprocs_lammps = atoi(arg[1]);
  if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int lammps;
  int flagpause;
  if (me < nprocs_lammps) lammps = 1;
  else lammps = MPI_UNDEFINED;
  MPI_Comm comm_lammps;
  MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
  LAMMPS *lmp;
  lmp = new LAMMPS (0, NULL, MPI_COMM_WORLD);

///**************************** LAMMPS, be quiet!*******************
    lmp->screen=0;
//  lmp->logfile=0;
//   lmp->universe->uscreen=0;
//   lmp->universe->ulogfile=0;
//*****************************************************************




/**********************************************************************
     run the input script thru LAMMPS one line at a time until end-of-file
     driver proc 0 reads a line, Bcasts it to all procs
     (could just send it to proc 0 of comm_lammps and let it Bcast)
     all LAMMPS procs call lammps_command() on the line */

   FILE *fp;
   if (me == 0) 
   {
     fp = fopen("in_Pt13.lmp","r");
    if (fp == NULL) 
    {
      printf("ERROR1: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD,1);
     }
   }

  int n;
  char line[1024];

    int natom;
  
         flagpause=0;
  double *xxx = (double *) malloc(3*(ptnatom+cnatom)*sizeof(double));

  while (1) 
  {

//   printf("me=%d\n",me);
    if (me == 0) 
    {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
//      if (n == 0) fclose(fp);
    }
//   printf("line=%s\n",line);
//     fp = fopen("thanks","wb");
//     fclose(fp);

    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);

    if(line[0]=='#'&&line[1]=='P'&&line[2]=='A'&&line[3]=='U'&&line[4]=='S') 
      {
          flagpause=1;//this command will not be read into lammps
          natom = lammps_get_natoms(lmp);
          if((ptnatom+cnatom)!=natom)
             {
               printf("waning atom number error!!\n");
               printf("input N=%d, Now N=%d\n",ptnatom+cnatom,natom);
               exit(1);
             }
           
          for(int  i=0;i<ptnatom;i++)
          {
            *( xxx+(i)*3+0)=*( xx+i*3+0) ; 
            *( xxx+(i)*3+1)=*( xx+i*3+1)  ;
            *( xxx+(i)*3+2)=*( xx+i*3+2)  ;
          }
          for(int  i=0;i<cnatom;i++)
          {
            *( xxx+(i+ptnatom)*3+0)=*( yy+i*3+0) ; 
            *( xxx+(i+ptnatom)*3+1)=*( yy+i*3+1)  ;
            *( xxx+(i+ptnatom)*3+2)=*( yy+i*3+2)  ;
          }
          lammps_scatter_atoms(lmp,"x",1,3,xxx);
    
          char nvt_command[100];
      
          sprintf(nvt_command,"velocity all create %d 4928459 rot yes dist gaussian",Temp);
          lmp->input->one(nvt_command);
    
          lmp->input->one("fix nve1 all nve ");
          for (int ii=Temp;ii>0;ii--)
           {
             sprintf(nvt_command,"fix scale1 all temp/rescale 1 %d %d 1 1",ii+1,ii);
             lmp->input->one(nvt_command);
             lmp->input->one("run 1 ");
             lmp->input->one("unfix scale1 ");
           }

          lmp->input->one("unfix nve1 ");
          
            
//         lmp->input->one("minimize        0.0 1.e-5 10 10");
     }
     lmp->input->one(line);
 
 
   } 
   fclose(fp);


   if(flagpause==0){printf("ERROR, 'PAUSE HERE information missing'\n "); MPI_Abort(MPI_COMM_WORLD,1);}

//*********************end of reading lammps input file************************


//*****************************************************************************
//******************* add extra command into lammps input file*****************

  
    *efinal=lmp->update->minimize->efinal;
    printf(" E=  %8.4f  ",*efinal);
    lammps_gather_atoms(lmp,"x",1,3,xxx);
    
        for(int  i=0;i<ptnatom;i++)
    {
      *( xx+i*3+0) =*( xxx+(i)*3+0);
      *( xx+i*3+1) =*( xxx+(i)*3+1) ;
      *( xx+i*3+2) =*( xxx+(i)*3+2) ;
    }
//    yy=xxx+3*ptnatom;
        for(int  i=0;i<cnatom;i++)
    {
      *( yy+i*3+0) =*( xxx+(i+ptnatom)*3+0);
      *( yy+i*3+1) =*( xxx+(i+ptnatom)*3+1) ;
      *( yy+i*3+2) =*( xxx+(i+ptnatom)*3+2) ;
    }
//    int ii=0;
//       printf("in optimization %5d,%lf,%lf,%lf\n",ii,*( xx+ii*3+0),*( xx+ii*3+1),*( xx+ii*3+2));
//       printf("in optimization %5d,%lf,%lf,%lf\n",ii,*( yy+ii*3+0),*( yy+ii*3+1),*( yy+ii*3+2));

// it is very wired that if using yy=xxx+ptatom*3 inside this function, yy value is changed, but yy can not be transfrred back to outsie.    
//summary. if define int **yy=lmp->atom->x. then yy[0][0]=3 will change internal value. if define int *x, then must use scatter function to change the value of internal.

    

     

//***************************************************************************
//******************** add extra command into lammps input file***************

    
    delete lmp;
//    free(arg);
    free(xxx);
//
    MPI_Finalize();
}
