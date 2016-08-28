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

int    dijkstra(double *people,int *path,int n);

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDijkstraAtom::ComputeDijkstraAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute Dijkstra/atom command");

  double cutoff = force->numeric(FLERR,arg[3]);
  cutsq = cutoff*cutoff;

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;

  nmax = 0;
  DijkstraID = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDijkstraAtom::~ComputeDijkstraAtom()
{
  memory->destroy(DijkstraID);
}

/* ---------------------------------------------------------------------- */

void ComputeDijkstraAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute Dijkstra/atom unless atoms have IDs");
  if (force->pair == NULL)
    error->all(FLERR,"Compute Dijkstra/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,
               "Compute Dijkstra/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list
  // full required so that pair of atoms on 2 procs both set their DijkstraID

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"Dijkstra/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute Dijkstra/atom");
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

  // grow DijkstraID array if necessary

  if (atom->nlocal+atom->nghost > nmax) {
    memory->destroy(DijkstraID);
    nmax = atom->nmax;
    memory->create(DijkstraID,nmax,"Dijkstra/atom:DijkstraID");
    vector_atom = DijkstraID;
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

  // every atom starts in its own Dijkstra, with DijkstraID = atomID

  tagint *tag = atom->tag;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    DijkstraID[i] = 0;
  }

  // loop until no more changes on any proc:
  // acquire DijkstraIDs of ghost atoms
  // loop over my atoms, checking distance to neighbors
  // if both atoms are in Dijkstra, assign lowest DijkstraID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = 1;
  double **x = atom->x;

  int change,done,anychange;

  while (1) {
    comm->forward_comm_compute(this);

    change = 0;
    while (1) {
      done = 1;
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (!(mask[j] & groupbit)) continue;
          if (DijkstraID[i] == DijkstraID[j]) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            DijkstraID[i] = DijkstraID[j] = MIN(DijkstraID[i],DijkstraID[j]);
            done = 0;
          }
        }
      }
      if (!done) change = 1;
      if (done) break;
    }

    // stop if all procs are done

    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;
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
      buf[m++] = DijkstraID[j];
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
    for (i = first; i < last; i++) DijkstraID[i] = buf[m++];
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

/* ----------------------------------------------------------------------
   dijkstra from panico_basico.c
------------------------------------------------------------------------- */

int dijkstra(double *people,int *path,int n)


{
  /*
 int    d,h,i,j,m,min,s,v,vv;
 int    cost[n*n],dist[n],prev[n],selected[n];        
 double xi,yi,ri,xj,yj,rj,dd,rr,d1,d2,d1min,d2min;

 v=0;
 vv=0;

 d1min=DOOR_X*DOOR_X;
 d2min=DOOR_X*DOOR_X;

 for(i=0;i<n-1;i++)
   {
     *(cost+n*i+i)=INFTY;
     xi=*(people+3*i+0);
     yi=*(people+3*i+1);
     ri=*(people+3*i+2)-RCR;

     if((xi<DOOR_X) && (xi+ri>DOOR_X))
       {
         d1=DOOR_X*DOOR_X;
         d2=DOOR_X*DOOR_X;

         if(yi-ri<ROOM_1) d1=(xi-DOOR_X)*(xi-DOOR_X)+(yi-DOOR_Y)*(yi-DOOR_Y);
         if(yi+ri>ROOM_2) d2=(xi-DOOR_X)*(xi-DOOR_X)+(yi-DOOR_Y)*(yi-DOOR_Y);

         if(d1<d1min) {v=i+1;  d1min=d1;}
         if(d2<d2min) {vv=i+1; d2min=d2;}
       }

     for(j=i+1;j<n;j++)
       {
         xj=*(people+3*j+0);
         yj=*(people+3*j+1);
         rj=*(people+3*j+2)-RCR;

         dd=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj);
         rr=(ri+rj)*(ri+rj);

         if((xi<DOOR_X) && (xj<DOOR_X) && (dd<rr)) {*(cost+n*i+j)=1; *(cost+n*j+i)=1;}
         else                                      {*(cost+n*i+j)=INFTY; *(cost+n*j+i)=INFTY;}
       }
  }

  *(cost+n*i+i)=INFTY;

 j=0;

 if(v*vv)
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
  
     // recorro la red hasta que alcance a vv (traget)
     // calculo la distancia tentativa d (de los no visitados).
     // si es minima la acepto.

     while(*(selected+vv-1)==0)
       {
         m=n+1;;
         min=INFTY;
        
        for(i=0;i<n;i++)
          {
            d=(*(dist+s))+(*(cost+n*s+i)); 
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
  
   //for(i=0;i<j;i++) printf("%d\t",*(path+i)); printf("\n");

  return j; 
  */
}



