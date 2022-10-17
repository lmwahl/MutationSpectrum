/* simulated adaptive walks with variable K in NK landscape */
/* see tstv.c for a carefully commented version of the code  */

// currently set up to write out fben for Ts and Tv separately to f01 data files //

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define MAXN 251     // max possible genome length
#define MAXSTEPS 500 // max number of accepted steps in walk
#define MAXMUTES 5001 // max number of mutants in DFEs
#define MAXWALKS 501  // max number of walks (replicates)
#define MAXWALKLENGTHS 20  // max number of different walk lengths to try
#define MAXSPECS 20  // max number of different bias values to test
#define K 1  // K in the NK fitness landscape, number of neighbours

/* code can be edited to switch from ts:tv to at:gc mutations */
void transition(int seq[], int n, int newseq[]);
void transversion(int seq[], int n, int newseq[]);
void mutation(int seq[],int newseq[],double ffixed, double fts); // at:gc

int N=200;  // actual genome length, change this to change N


void main()  //(int argc, char **argv)
{
int walklengths[18] = {5e4}, nlengths = 1;
int i,j,k,n, ilength, currentwalklength, nlookups, lookind;
int nsteps, maxwalklength = 6e4;
double ftstep=1;
double ffixed = .1;
int nwalks = 500;
double fts, ftsmu, fwt[MAXWALKS][MAXWALKLENGTHS], dfe;
int nstepsall[MAXSTEPS][MAXWALKLENGTHS];
int ift, iwalk, ispec, itry, pos;
double oneoverNe = .00001;
 int nmutes = 5000;  // should not be greater than 3*N 
int neighs[K][MAXN], anc[MAXN], tmpseq[MAXN];
double w0, wseq[MAXSTEPS], wtmp, s;
char filename[100];
FILE *fpout;
double sumpos, sumneg;
//double fben[MAXWALKS][MAXWALKLENGTHS], fdel[MAXWALKS][MAXWALKLENGTHS];
//double sbarben[MAXWALKS][MAXWALKLENGTHS],sbardel[MAXWALKS][MAXWALKLENGTHS];
//double  fbenmu[MAXWALKS][MAXSPECS][MAXWALKLENGTHS], fdelmu[MAXWALKS][MAXSPECS][MAXWALKLENGTHS];
//double sbarbenmu[MAXWALKS][MAXSPECS][MAXWALKLENGTHS],sbardelmu[MAXWALKS][MAXSPECS][MAXWALKLENGTHS];
int npos, nneg, imute;
double sums[20], dfemu;

 nlookups = 4;  for (i=1;i<=K;i++) nlookups = nlookups*4;
 double **wis = (double **)malloc(nlookups * sizeof(double *)); 
 for (i=0; i<nlookups; i++) 
         wis[i] = (double *)malloc(MAXN * sizeof(double));
 double **fben = (double **)malloc(nwalks * sizeof(double *)); 
 for (i=0; i<nwalks; i++) 
         fben[i] = (double *)malloc(nlengths * sizeof(double));
 double **fdel = (double **)malloc(nwalks * sizeof(double *)); 
 for (i=0; i<nwalks; i++) 
         fdel[i] = (double *)malloc(nlengths * sizeof(double));
 double **sbarben = (double **)malloc(nwalks * sizeof(double *)); 
 for (i=0; i<nwalks; i++) 
         sbarben[i] = (double *)malloc(nlengths * sizeof(double));
 double **sbardel = (double **)malloc(nwalks * sizeof(double *)); 
 for (i=0; i<nwalks; i++) 
         sbardel[i] = (double *)malloc(nlengths * sizeof(double));
 double ***fbenmu = (double ***)malloc(nwalks * sizeof(double **)); 
 for (i=0; i<nwalks; i++) {
         fbenmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            fbenmu[i][j] = (double *)malloc(nlengths * sizeof(double));}
 double ***fdelmu = (double ***)malloc(nwalks * sizeof(double **)); 
 for (i=0; i<nwalks; i++) {
         fdelmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            fdelmu[i][j] = (double *)malloc(nlengths * sizeof(double));}
 double ***sbarbenmu = (double ***)malloc(nwalks * sizeof(double **)); 
 for (i=0; i<nwalks; i++) {
         sbarbenmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            sbarbenmu[i][j] = (double *)malloc(nlengths * sizeof(double));}
 double ***sbardelmu = (double ***)malloc(nwalks * sizeof(double **)); 
 for (i=0; i<nwalks; i++) {
         sbardelmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            sbardelmu[i][j] = (double *)malloc(nlengths * sizeof(double)); }


 srand(1);
 fprintf(stdout,"Running Ts:Tv, K =%d\n",K);

 // wt bias is set here; mutant steps through ftstep in a later loop (ispec) //
 //     for (fts=ftstep;fts<=1-ffixed;fts+=ftstep) {
 for (fts=1.0/3.0; fts<=1.0/3.0; fts+=ftstep) {

 for (iwalk=0;iwalk<nwalks;iwalk++) {
   for (i=0;i<nlookups;i++)
      for (n=0;n<N;n++) wis[i][n] = (rand()/(double)RAND_MAX);
   for (k=0;k<K;k++) for (n=0;n<N;n++) neighs[k][n] = (int)(N*(rand()/(double)RAND_MAX));
   for (n=0;n<N;n++) anc[n] = (int)(4*(rand()/(double)RAND_MAX));
   w0 = 0;
   for (n=0;n<N;n++) {
     lookind = anc[n];
     for (k=0;k<K;k++) lookind += anc[neighs[k][n]]*pow(4.0,(double)(k+1));
     w0 += wis[lookind][n];
   }
   wseq[0] = w0;
   nsteps = 0;
   currentwalklength = 0;
   for (ilength=0;ilength<nlengths;ilength++) {
     maxwalklength = walklengths[ilength];
     // test opening the file before doing any more computation
     sprintf(filename, "kany_data_f01_%d/N%d_Sw_walk%d_fts%1.2f.out",K,N,maxwalklength,fts);
     fpout = fopen(filename,"w");
     if (fpout==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}  
     fclose(fpout);
     
   for (itry = 1; itry<=maxwalklength-currentwalklength; itry++) {
     /* mutation(anc,tmpseq,ffixed,fts); */
     pos = (int)(N*(rand()/(double)RAND_MAX));
     if ((rand()/(double)RAND_MAX)<fts)
       transition(anc,pos,tmpseq);
     else
     transversion(anc,pos,tmpseq);
     wtmp = 0;
     for (n=0;n<N;n++) {
       lookind = tmpseq[n];
       for (k=0;k<K;k++) lookind += tmpseq[neighs[k][n]]*pow(4.0,(double)(k+1));
       wtmp += wis[lookind][n];
     }
    s = wtmp/wseq[nsteps]-1;
    if (s>0) {
      if ((rand()/(double)RAND_MAX)<2*s) {
        nsteps = nsteps+1;
        for (k=0;k<N;k++) anc[k] = tmpseq[k];
        wseq[nsteps] = wtmp;
      }
    }  else { 
      if ((rand()/(double)RAND_MAX)<oneoverNe*exp(2*s)) {
       nsteps = nsteps+1;
       for (k=0;k<N;k++) anc[k] = tmpseq[k];
       wseq[nsteps] = wtmp;
      }
   }
   }  // loop on itry
   
   fwt[iwalk][ilength] = wseq[nsteps];
   npos=0; nneg = 0; sumpos=0; sumneg=0;
   for (imute=0;imute<nmutes;imute++) {
     /*     mutation(anc,tmpseq,ffixed,fts); */
     pos = (int)(N*(rand()/(double)RAND_MAX));
     if ((rand()/(double)RAND_MAX)<fts)
       transition(anc,pos,tmpseq);
     else
     transversion(anc,pos,tmpseq); 
     wtmp = 0;
     for (n=0;n<N;n++) {
       lookind = tmpseq[n];
       for (k=0;k<K;k++) lookind += tmpseq[neighs[k][n]]*pow(4.0,(double)(k+1));
       wtmp += wis[lookind][n];
     }
     dfe = wtmp/fwt[iwalk][ilength] - 1;
     if (dfe>0) {npos++; sumpos+= dfe;}
     if (dfe<=0) {nneg++; sumneg+= dfe;}
   }

   nstepsall[iwalk][ilength] = nsteps;
   fben[iwalk][ilength] = (double)npos/nmutes;
   fdel[iwalk][ilength] = (double)nneg/nmutes;
   if (npos>0) sbarben[iwalk][ilength] = sumpos/npos;
   else sbarben[iwalk][ilength]=0;
   if (nneg>0) sbardel[iwalk][ilength] = sumneg/nneg;
   else sbardel[iwalk][ilength]=0;
   
  ispec = 0;
   for (ftsmu=0;ftsmu<=1;ftsmu+=1.0) {
     // for (ftsmu=ftstep;ftsmu<=ftstep;ftsmu+=ftstep) {
     //   for (ftsmu=ftstep;ftsmu<=1-ftstep;ftsmu+=ftstep) {
    npos=0; nneg = 0; sumpos=0; sumneg=0;
    for (imute=0;imute<nmutes;imute++) {
      /* mutation(anc,tmpseq,ffixed,ftsmu); */
     pos = (int)(N*(rand()/(double)RAND_MAX));
     if ((rand()/(double)RAND_MAX)<ftsmu)
       transition(anc,pos,tmpseq);
     else
     transversion(anc,pos,tmpseq);
     wtmp = 0;
     for (n=0;n<N;n++) {
       lookind = tmpseq[n];
       for (k=0;k<K;k++) lookind += tmpseq[neighs[k][n]]*pow(4.0,(double)(k+1));
       wtmp += wis[lookind][n];
     }
     dfemu = wtmp/fwt[iwalk][ilength] - 1;
     if (dfemu>0) {npos++; sumpos+= dfemu;}
     if (dfemu<0) {nneg++; sumneg+= dfemu;}
    }
   fbenmu[iwalk][ispec][ilength] = (double)npos/nmutes;
   fdelmu[iwalk][ispec][ilength] = (double)nneg/nmutes;
   if (npos>0) sbarbenmu[iwalk][ispec][ilength] = sumpos/npos;
   else sbarbenmu[iwalk][ispec][ilength]=0;
   if (nneg>0) sbardelmu[iwalk][ispec][ilength] = sumneg/nneg;
   else sbardelmu[iwalk][ispec][ilength] = 0;
   ispec++;
   }  // loop on ispec (ftsmu)
   }  // loop on ilength
 }  // loop on iwalk

 // all the walks are now finished, save the results to files
   for (ilength=0;ilength<nlengths;ilength++) {
     //     sprintf(filename, "kany_data_f01_%d/N%d_Sw_walk%d_fts%1.2f_TsTv.out",K,N,walklengths[ilength],fts);
     //     fpout = fopen(filename,"w");
     //     for (iwalk=0;iwalk<nwalks;iwalk++)
     //       fprintf(fpout,"%f %f\n",fbenmu[iwalk][0][ilength],fbenmu[iwalk][1][ilength]);
     //     fclose(fpout);
     sprintf(filename, "kany_data_f01_%d/N%d_Sw_walk%d_fts%1.2f.out",K,N,walklengths[ilength],fts);
     fpout = fopen(filename,"w");
     if (fpout==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}  
     
    sums[1] =0; sums[2]=0; sums[3]= 0; sums[4] = 0; sums[5]=0; sums[6]=0;
    int nzb=0; int nzd=0;
  for (iwalk=0;iwalk<nwalks;iwalk++) {
   sums[1]+= fben[iwalk][ilength]; sums[2]+= sbarben[iwalk][ilength];
   sums[3]+= fdel[iwalk][ilength]; sums[4]+= sbardel[iwalk][ilength];
   sums[5]+= fwt[iwalk][ilength];  sums[6]+= nstepsall[iwalk][ilength];
   if (sbarben[iwalk][ilength]>0) nzb++;
   if (sbardel[iwalk][ilength]<0) nzd++;
  }
   fprintf(fpout,"%f %f %f %f\n",fts,ftstep,sums[5]/nwalks,sums[6]/nwalks);
   fprintf(fpout,"%f %f %f %f\n",sums[1]/nwalks,sums[2]/nzb,sums[3]/nwalks,sums[4]/nzd);

  for (i=0;i<ispec;i++) {
   sums[1] =0; sums[2]=0; sums[3]= 0; sums[4] = 0; nzb=0; nzd=0; 
   for (iwalk=0;iwalk<nwalks;iwalk++) {
    sums[1]+= fbenmu[iwalk][i][ilength]; sums[2]+= sbarbenmu[iwalk][i][ilength];
    sums[3]+= fdelmu[iwalk][i][ilength]; sums[4]+= sbardelmu[iwalk][i][ilength];
    if (sbarbenmu[iwalk][i][ilength]>0) nzb++;
    if (sbardelmu[iwalk][i][ilength]<0) nzd++;
   }
   fprintf(fpout,"%f %f %f %f\n",sums[1]/nwalks,sums[2]/nzb,sums[3]/nwalks,sums[4]/nzd);
   }

  fclose(fpout);

  } //loop on ilength for saving to file
/*
 fprintf(stdout,"Note: printing all data to output files\n");
 fpout = fopen("wt_data.txt","w");
 for (iwalk=0;iwalk<nwalks;iwalk++)
   fprintf(fpout,"%f %f %f %f %f %d\n",fben[iwalk],sbarben[iwalk],fdel[iwalk],sbardel[iwalk],fwt[iwalk],nstepsall[iwalk]);
 fclose(fpout);
 fpout = fopen("mt_data.txt","w");
 for  (iwalk=0;iwalk<nwalks;iwalk++)
   for (i=0;i<ispec;i++)
     fprintf(fpout," %f %f %f %f\n",fbenmu[iwalk][i],sbarbenmu[iwalk][i],fdelmu[iwalk][i],sbardelmu[iwalk][i]);
 fclose(fpout);
*/
}  // loop on fts of wildtype

}  // end of main

void transition(int seq[MAXN],int pos,int newseq[MAXN])
// using standard codon model order
// 0 = T, 1 = C, 2 = A, 3 = G
{
for (int i=0;i<N;i++) newseq[i]=seq[i];
switch(seq[pos]) {
  case 0  : newseq[pos] = 1; break;
  case 1  : newseq[pos] = 0; break;
  case 2  : newseq[pos] = 3; break;
  case 3  : newseq[pos] = 2; break;
}
}

void transversion(int seq[MAXN],int pos,int newseq[MAXN])
{
for (int i=0;i<N;i++) newseq[i]=seq[i];
switch(seq[pos]) {
  case 0  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 2; 
  else newseq[pos] = 3;
   break;
  case 1  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 2;
  else newseq[pos] = 3;
  break;
  case 2  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
  else newseq[pos] = 1;
  break;
  case 3  :
  if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
  else newseq[pos] = 1;
  break;
}
}

void mutation(int seq[MAXN],int newseq[MAXN],double ffixed, double fts)
{
  int i, pos;
  double muprob;

/* This subroutine mutates the sequence. The mutations will be AT to AT
or GC to GC with probability ffixed.  With probability 1-ffixed, the mutn
will be GC to AT or vice versa.  Among the latter mutations,
proportion fts are GC to AT, and proportion 1-fts are AT to GC.
Note that this proportion (fts) will be maintained independent of the
GC content of the genome.  So fts corresponds to the actually realized
proportion of GC to AT mutations.  This should be kept in mind when
interpreting results if the genome GC content changes.  */
  
for (i=0;i<N;i++) newseq[i]=seq[i];
 muprob = rand()/(double)RAND_MAX;
 if (muprob<ffixed) {  // we want AT to AT or GC to GC
  pos = (int)(N*(rand()/(double)RAND_MAX));
  switch(seq[pos]) {
    case 0: newseq[pos] = 2; break;
    case 2: newseq[pos] = 0; break;
    case 1: newseq[pos] = 3; break;
    case 3: newseq[pos] = 1; break;
   }
  }
 else {
  if (muprob<(ffixed+fts)) {   // we want GC to AT
   pos = (int)(N*(rand()/(double)RAND_MAX));
   while ((seq[pos]==2||seq[pos]==0)  pos = (int)(N*(rand()/(double)RAND_MAX));
   if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
   else newseq[pos] = 2;
  }
  else {// we want AT to GC
   pos = (int)(N*(rand()/(double)RAND_MAX));
   while ((seq[pos]==1)||seq[pos]==3)  pos = (int)(N*(rand()/(double)RAND_MAX));
   if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 1;
   else newseq[pos] = 3;
  }
 }
}
