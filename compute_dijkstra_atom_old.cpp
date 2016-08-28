/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_dijkstra_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "group.h"

#define INFTY  9999
#define DOOR_X 20.0 
#define DOOR_Y 10.0
#define SKIN   0.001
 
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDijkstraAtom::ComputeDijkstraAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal compute dijkstra/atom command");

   lower = force->numeric(FLERR,arg[3]);
   upper = force->numeric(FLERR,arg[4]);
   rad   = force->numeric(FLERR,arg[5]);
   label = force->numeric(FLERR,arg[6]);

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;

  nmax = 0;
  dijkstraID = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDijkstraAtom::~ComputeDijkstraAtom()
{
  memory->destroy(dijkstraID);
}

/* ---------------------------------------------------------------------- */

void ComputeDijkstraAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute dijkstra/atom unless atoms have IDs");
  if (force->pair == NULL)
    error->all(FLERR,"Compute dijkstra/atom requires a pair style be defined");

  // need an occasional full neighbor list
  // full required so that pair of atoms on 2 procs both set their dijkstraID

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"dijkstra/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute dijkstra/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeDijkstraAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeDijkstraAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow dijkstraID array if necessary

  if (atom->nlocal+atom->nghost > nmax) {
    memory->destroy(dijkstraID);
    nmax = atom->nmax;
    memory->create(dijkstraID,nmax,"dijkstra/atom:dijkstraID");
    vector_atom = dijkstraID;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // if group is dynamic, insure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = 0;
    comm->forward_comm_compute(this);
  }

  // every atom starts in its own cluster, with dijkstraID = atomID

  tagint *tag = atom->tag;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
   dijkstraID[i] = 0;
    }
  // loop until no more changes on any proc:
  // acquire dijkstraIDs of ghost atoms
  // loop over my atoms, checking distance to neighbors
  // if both atoms are in cluster, assign lowest dijkstraID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = 1;
  double **x = atom->x;

  int n = atom->nlocal;

  int    d,m,min,s,v,vv,change,anychange,h,dist[n];
  int    cost[n][n],prev[n],selected[n],path[n];
  double room_1,room_2,d1,d2,d1min,d2min,diam2;

  v=0;
  vv=0;

  room_1=lower;
  room_2=upper;
  diam2=4.0*rad*rad;

  d1min=DOOR_X*DOOR_X;
  d2min=DOOR_X*DOOR_X;

  for (ii=0;ii<n;ii++){
    for (jj=0;jj<n;jj++){
          cost[ii][jj]=INFTY;
    }
  }
  
  if (1) {
    comm->forward_comm_compute(this);

    change = 0;
 

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
     
    

      if (xtmp<DOOR_X && xtmp+rad>DOOR_X-SKIN ) {

          

        if (ytmp-rad<lower) d1=(xtmp-DOOR_X)*(xtmp-DOOR_X)+(ytmp-lower)*(ytmp-lower);
        if (ytmp+rad>upper) d2=(xtmp-DOOR_X)*(xtmp-DOOR_X)+(ytmp-upper)*(ytmp-upper);

        if (d1<d1min) {v=tag[i] ; d1min=d1; }
        if (d2<d2min) {vv=tag[i]; d2min=d2;  }
      }

   
  
     jlist = firstneigh[i];
     jnum = numneigh[i];

     for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (!(mask[j] & groupbit)) continue;       

        

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < diam2 && xtmp<DOOR_X && (xtmp-delx)<DOOR_X) {
          cost[tag[i]][tag[j]]=1; 
          cost[tag[j]][tag[i]]=1;
        }
      }
    }
  
    j=0;

    if (v*vv)
     {
       // inicializo las variables
    
       for(i=0;i<n;i++)
         {
           
           *(dist+i)=INFTY;
           *(prev+i)=-1;
           *(selected+i)=0;
           *(path+i)=-1;
        
         }

       s=v-1;
       *(selected+s)=1;
       *(dist+s)=0;

       // recorro la red hasta que alcance a vv (target)
       // calculo la distancia tentativa d (de los no visitados).
       // si es minima la acepto.

       while(*(selected+vv-1)==0)
         {
           m=n+1;
           min=INFTY;
        
          for(i=0;i<n;i++)
            {
              d=(*(dist+s))+cost[s][i]; 
              if(d<(*(dist+i)) && *(selected+i)==0)
              {
                  *(dist+i)=d;
                  *(prev+i)=s;
              }
              if(*(dist+i)<min && *(selected+i)==0)
                {
                  min=*(dist+i);
                  m=i;
                }
            }

          if(m<n+1) {s=m; *(selected+s)=1;} 
          else {vv=0; break;}
         }

       // almaceno el camino mas corto en path

     
       s=vv-1;
       while(s!=-1)
         {
           *(path+j)=s;
           s=*(prev+s);
           j++;
         }     
        
     }


    // ordeno los nodos de menor a mayor

    for(h=0;h<j;h++)
      for(i=0;i<j-1;i++)
        if(*(path+i)>(*(path+i+1))) 
          {
            m=*(path+i);
            *(path+i)=(*(path+i+1));
            *(path+i+1)=m;
          }
  
        for(i=0;i<n;i++){
          dijkstraID[i]=0;
        }
    

      for(i=0;i<j;i++)
        {
           h=*(path+i);
         printf("%d\t",h );
           dijkstraID[h]=(int)label;
        }   

        printf("\n\n\n");
    
    // stop if all procs are done

    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    //if (!anychange) break;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeDijkstraAtom::pack_forward_comm(int n, int *list, double *buf,
                                          int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  if (commflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = dijkstraID[j];
    }
  } else {
    int *mask = atom->mask;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(mask[j]).d;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeDijkstraAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag)
    for (i = first; i < last; i++) dijkstraID[i] = buf[m++];
  else {
    int *mask = atom->mask;
    for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDijkstraAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
