// B. Joachimi
// 14.05.2015 
// calculates tomographic shear/clustering power spectrum covariance (Gaussian)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "bjutils.h"

#define NL 12
#define NZ 3

int pscov();

FILE *bj_fileopen(char *rw,char *name);
double *bj_alloc(int dim);
double ***bj_alloc_3D(int dim1,int dim2,int dim3);
void bj_free_3D(double ***arr,int dim1,int dim2,int dim3);

int main(int argc, char *argv[])
{
  int i,j,k,status;
  FILE *dat;
  char name[200];

  // set covariance parameters
  double ngal[NZ]={10./3.,10./3.,10./3.};
  double sig[NZ]={0.35,0.35,0.35};
  double area = 5000.0; //Survey area in deg^2

  // allocate arrays
  double *ell=bj_alloc(NL);
  double ***ps=bj_alloc_3D(4,NZ*NZ,NL);
  const int ND=NZ*(2*NZ+1);
  double ***cov=bj_alloc_3D(ND,ND,NL);  // [NZ*NZ][NZ*NZ][NL] suffices for individual WL or clustering covariances

  // read power spectra from file
  sprintf(name,"ps_test.dat");
  dat=bj_fileopen("r",name);

  for (i=0;i<NL;i++) { 
   fscanf(dat,"%lf",&ell[i]);
    for (j=0;j<NZ;j++) { 
      for (k=j;k<NZ;k++) { 
	fscanf(dat,"%lf",&ps[0][k+j*NZ][i]);
	ps[0][j+k*NZ][i]=ps[0][k+j*NZ][i];  // symmetrise
      }
    }
  }
  fclose(dat);

  // calculate covariance
  status=pscov(ell,ps,cov,area,ngal,sig,NZ,NL,1);
  printf("%i\n",status);

  // clean up
  bj_free_3D(ps,4,NZ*NZ,NL);
  bj_free_3D(cov,ND,ND,NL);
  free(ell);
  return 0;
}




/*** Notes for pscov routine ***/
/*
returns status flag [0: all fine; !=0: something went wrong]

 ARGUMENTS:
 corr_type: flag that selects type of power spectra used [0: combined WL and clustering; 1: WL only; 2: clustering only]
 nz: no. of tomographic redshift bins
 nl: no. of angular frequency bins
 ell[nl]: angular frequencies
 ps[4][nz*nz][nl]: power spectra (modified in routine!)
   1st: 0: ee(WL; for corr_type=1 only this used); 1: en; 2: ne; 3: nn(clustering; for corr_type=2 only this used)   
   2nd: tomographic bins in order 11,..,1nz,21,..,2nz,nz1,..,nznz
   3rd: angular frequencies
 cov[nd][nd][nl]: contains covariance on return
   1st: nd=nz^2 for corr_type=1,2; nd=nz*(2*nz+1) for corr_type=0 (here double-counting of data is avoided; structure as in Eq. (38) of Joachimi & Bridle (2010))
   2nd: see 1st
   3rd: angular frequencies
 survey_area: survey area in [deg^2]
 ndens[nz]: array of (effective) galaxy number density for each tomographic bin in units of arcmin^{-2} [Note: assumed to be the same for clustering and weak lensing]
 sigma_eps[nz]: array of total (as opposed to per component) ellipticity dispersion values for each tomographic bin [ignored if corr_type=1]
*/

#define PI 3.14159
#define deg2_to_sterad 3.046174198e-4   // converts deg^2 to sterad
int pscov(double *ell,double ***ps,double ***cov,double survey_area,double *ndens,double *sigma_eps,int nz,int nl,int corr_type)
{
  int i,j,k,i2,j2,m,n,p,q,p2,q2,index,index2;
  double noise,dlnl;
  const double prefacconst=2.*PI/(deg2_to_sterad*survey_area);
  double *delta=calloc(nl,sizeof(double));


  // determine bin width
  if (nl==1) {
    printf("Error: cannot determine bin width if only one angular frequency bin given.\n");
    return(-1);
  }
  else {
    dlnl=(log(ell[nl-1])-log(ell[0]))/(1.*nl-1.);
  }
  
  for(k=0;k<nl;k++) {
    delta[k]=ell[k]*(exp(dlnl/2.)-exp(-dlnl/2.));
  }


  // add shot noise
  for(i=0;i<nz;i++) {   
    if((corr_type==0)||(corr_type==1)) {
      noise=sigma_eps[i]*sigma_eps[i]/(2.*ndens[i]/(deg2_to_sterad/3600.));   
      for(k=0;k<nl;k++) {
	ps[0][i+nz*i][k]+=noise;
      }
    }

    if((corr_type==0)||(corr_type==2)) {
      noise=1./(ndens[i]/(deg2_to_sterad/3600.));
      for(k=0;k<nl;k++) {
	ps[3][i+nz*i][k]+=noise;
      }
    }
  }
  

  // calculate covariance
  if(corr_type>=1) {
    q=(corr_type-1)*3;  //q=0: WL; q=3: clustering

    for(k=0;k<nl;k++) {
     for(i=0;i<nz;i++) {
      for(j=0;j<nz;j++) {
       for(i2=0;i2<nz;i2++) {
	for(j2=0;j2<nz;j2++) {
	 cov[j+nz*i][j2+nz*i2][k]=prefacconst/(ell[k]*delta[k])*(ps[q][i2+nz*i][k]*ps[q][j2+nz*j][k]+ps[q][j2+nz*i][k]*ps[q][i2+nz*j][k]);
	}
       }
      }
     }
    }
  }
  else {  // joint WL and clustering
    for(k=0;k<nl;k++) {
     index=0;
     for(m=0;m<4;m++) {
      p=m/2;
      q=m%2;
      for(i=0;i<nz; i++) {
       for(j=i;j<nz; j++) {

	index2=0;
	for(n=0;n<4;n++) {
	 p2=n/2;
	 q2=n%2;
	 for(i2=0;i2<nz; i2++) {
	  for(j2=i2;j2<nz; j2++) {

	   if (index2>=index) {
	    if (!(((q+2*p==2)&&(i==j))||((q2+2*p2==2)&&(i2==j2)))) {

	     cov[index][index2][k]=prefacconst/(ell[k]*delta[k])*(ps[p2+2*p][i2+nz*i][k]*ps[q2+2*q][j2+nz*j][k]+ps[q2+2*p][j2+nz*i][k]*ps[p2+2*q][i2+nz*j][k]);

	    } // excludes auto-corr for ne -> identical to en
	   }  // ensures that no values below diagonal are written

	   index2++;
	  }
	 }
	}

	index++;
       }
      }
     }

    }
  }


  // clean up
  free(delta);
  return 0;
}
#undef deg2_to_sterad
#undef PI





FILE *bj_fileopen(char *rw,char *name)
{
  FILE *dat;
  if ((dat=fopen(name,rw))==NULL) {
    printf("Could not open file %s\n",name);
    exit(-1);
  }
  return dat;
}


double *bj_alloc(int dim)
{
  double *arr=calloc(dim,sizeof(double));
  return arr;
}


double ***bj_alloc_3D(int dim1,int dim2,int dim3)
{
  int i,j;
  double ***arr=calloc(dim1,sizeof(double *));
  for (i=0;i<dim1;i++) {
    arr[i]=calloc(dim2,sizeof(double *));
    for (j=0;j<dim2;j++) {
      arr[i][j]=calloc(dim3,sizeof(double));
    }
  }
  return arr;
}

void bj_free_3D(double ***arr,int dim1,int dim2,int dim3)
{
  int i,j;
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(arr[i][j]);
    }
    free(arr[i]);
  }
  free(arr);
  return;
}
