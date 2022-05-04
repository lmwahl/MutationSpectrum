/* revised version of code with variable K in NK landscape 
 using a population size NN and ~Wright-Fisher sampling such that
 multiple mutations may segregate.
 Instead of walks of different lengths, we will evolve for a 
 specified number of generations
 In this 'aa.c' version, we have NC codons and N = NC/3 amino acids.
 The fitness landscape depends only on the amino acids.  Therefore,
 synonymous changes have no effect on fitness. */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define MAXN 251  /* max length of genome */  
#define MAXTYPES 2000    /* number of distinct genotypes in the popn */
#define MAXSPECS 20     /* max number of ts:tv biases to shift to when making dfes */
#define K 3

void transition(int seq[], int n, int newseq[]);
void transversion(int seq[], int n, int newseq[]);
void mutation(int seq[],int newseq[],double ffixed, double fts);

int NB=99;        /* sequence length, must be a multiple of three */
int N;           /* number of amino acids, length of sequence in fitness landscape model */
int ntypesave = 50;    /* number of genotypes to start */
int Noftypesave = 20000;
int carryingcap = 100000;

void main()  //(int argc, char **argv)
{
int nreps = 10000;
double mu = 1e-4;  // mutation rate per GENOME per generation, not per base, 1e-4, -5 and -6 are good
int ngens = 5001;
int dfetimes[9] = {0, 100, 200, 500, 1000, 1500, 2e3, 4e3, 5e3}; //, 10e3, 20e3, 25e3};	    
int ntimes = 9;
double ftsave = 0.55;  // initial fts, freq of ts
double biases[10] = { .1, .2, .3, .4, .5, .55, .6, .7, .8, .9};  // ts fractions when making dfes
int nspecs = 10, ispec;
int  ndfe = 500;  // number of fitness measurements in each dfe
int ctoa[4][4][4]; // lookup table that translates codons to amino acids
 double pftmu = 0;  // probability that the mutation affects freq of ts
double lethalfrac = 0; // fraction of fitness contributions that are lethal
double lethalwi = -1e6; // fitness contribution if lethal
int itime = 0;  // index of the next dfetime we will do 
int i,j,k,n, igen, itype, jtype, irep, nlookups, lookind;
long seed=-1;
double fts, w0, sumpos, sumneg, bias;
 int npos, nneg, pos, ntotal, ntypes, aaneigh;
int neighs[K][MAXN]; // epistatic neighbours in NK model
int seqs[MAXTYPES][MAXN]; //sequences of genotypes that make up the popn
int dfeseq[MAXN];  // temporary sequence of potential indl contributing to dfe
int Noftype[MAXTYPES]; // number of individuals of each distinct genotype
double ws[MAXTYPES], ftseq[MAXTYPES]; // fitness and freq of tv for each genotype
double wbar, ftbar, wtmp, nweights[MAXTYPES];
char filename[100], filenamet[100], dfefilename[100];
FILE *fpout, *fpdfe, *fpoutt;
float gasdev(long *), poidev(float, long *), ran1(long *);

 N = NB/3;
 nlookups = 21;  for (i=1;i<=K;i++) nlookups = nlookups*21;
 double **Nnwbarsum = (double **)malloc(ngens * sizeof(double *));
 for (i=0; i<ngens; i++)
   Nnwbarsum[i] = (double *)malloc(3 * sizeof(double *));
 double **wis = (double **)malloc(nlookups * sizeof(double *)); 
 for (i=0; i<nlookups; i++) 
         wis[i] = (double *)malloc(MAXN * sizeof(double));
 double ***fbenmu = (double ***)malloc(nreps * sizeof(double **)); 
 for (i=0; i<nreps; i++) {
         fbenmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            fbenmu[i][j] = (double *)malloc(ntimes * sizeof(double));}
 double ***fdelmu = (double ***)malloc(nreps * sizeof(double **)); 
 for (i=0; i<nreps; i++) {
         fdelmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            fdelmu[i][j] = (double *)malloc(ntimes * sizeof(double));}
 double ***sbarbenmu = (double ***)malloc(nreps * sizeof(double **)); 
 for (i=0; i<nreps; i++) {
         sbarbenmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            sbarbenmu[i][j] = (double *)malloc(ntimes * sizeof(double));}
 double ***sbardelmu = (double ***)malloc(nreps * sizeof(double **)); 
 for (i=0; i<nreps; i++) {
         sbardelmu[i] = (double **)malloc(MAXSPECS * sizeof(double *));
         for (j=0; j<MAXSPECS; j++)
            sbardelmu[i][j] = (double *)malloc(ntimes * sizeof(double)); }
 for (i=0;i<ngens;i++) { Nnwbarsum[i][0] = 0; Nnwbarsum[i][1] =0; Nnwbarsum[i][2] = 0; }
 initialize_ctoa(ctoa);
 srand(1);
 
 fprintf(stdout,"Running Ts:Tv, N = %d, K =%d, ntypes = %d, ngens = %d, ts = %f, mu=%e\n",N,K,ntypesave,ngens,ftsave,mu);
     // test opening the files before doing any more computation
 sprintf(filename, "data/N%d_K%d_mu%4.0e_ts%d_N_Ntypes_wbar.out",N,K,mu,(int)(100.0*ftsave));
 fpout = fopen(filename,"w");
 if (fpout==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}
 sprintf(filenamet, "data/N%d_K%d_ts%d_times.out",N,K,(int)(100.0*ftsave));
 fpoutt = fopen(filenamet,"w");
 if (fpoutt==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}
 for (i=0;i<ntimes;i++) fprintf(fpoutt,"%d\n",dfetimes[i]);
 fclose(fpoutt);

 for (irep=0;irep<nreps;irep++) {  // replicate the whole process nreps times
   if ((nreps<=20) || (irep/10.0 == floor(irep/10.0))) { fprintf(stdout,"%d ",irep); fflush(stdout); }
   ntypes = ntypesave; itime = 0;
   // make the fitness landscape first
   for (i=0;i<nlookups;i++)
      for (n=0;n<N;n++) {
        if (ran1(&seed)>lethalfrac)
	  wis[i][n] = (rand()/(double)RAND_MAX);
	else wis[i][n] = lethalwi;
      }
   for (k=0;k<K;k++) for (n=0;n<N;n++) neighs[k][n] = (int)(N*(rand()/(double)RAND_MAX));
   //  then make the initial population of random sequences
   for (itype=0;itype<ntypes;itype++) {
   for (n=0;n<NB;n++) seqs[itype][n] = (int)(4*(rand()/(double)RAND_MAX));
   w0 = 0;
   for (n=0;n<N;n++) {
     lookind = ctoa[seqs[itype][3*n]][seqs[itype][3*n+1]][seqs[itype][3*n+2]];
     for (k=0;k<K;k++) {
       aaneigh = ctoa[seqs[itype][3*neighs[k][n]]][seqs[itype][3*neighs[k][n]+1]][seqs[itype][3*neighs[k][n]+2]];
       lookind += aaneigh*pow(21.0,(double)(k+1));
     }  
     w0 += wis[lookind][n];
   }
   ws[itype] = fmax(w0,0.0); Noftype[itype] = Noftypesave;
   ftseq[itype] = ftsave;
   } // end loop, ntypes genotypes have been created, with 1 individual each
     // here we give them all the same starting freq of tv
     
   for (igen = 0; igen<ngens; igen++) {   // loop over ngens generations
   // first compute total number of individuals and mean fitness
   double wwsum = 0; ntotal = 0;
   for (itype=0;itype<ntypes;itype++) {
       wwsum += ws[itype]*Noftype[itype];
       ntotal += Noftype[itype];
   }
   wbar = wwsum/ntotal;
   if (igen == dfetimes[itime]) {
     // make a dfe and print it out for this generation
     // also make dfes that we would get if the ts:tv bias changed
     // currently just making a dfe of the most common genotype in the popn
     // could also pick individuals at random to construct a popn-wide dfe
      // find the most common genotype
     int tmpN = Noftype[0]; int bigtype = 0;
     for (itype=1; itype<ntypes;itype++)
       if (Noftype[itype] > tmpN) { tmpN = Noftype[itype]; bigtype = itype; }
     // find fitness of most common genotype
     double wbig = 0;
     for (n=0;n<N;n++) {
        lookind = ctoa[seqs[bigtype][3*n]][seqs[bigtype][3*n+1]][seqs[bigtype][3*n+2]];
        for (k=0;k<K;k++) {
	  aaneigh = ctoa[seqs[bigtype][3*neighs[k][n]]][seqs[bigtype][3*neighs[k][n]+1]][seqs[bigtype][3*neighs[k][n]+2]];
          lookind += aaneigh*pow(21.0,(double)(k+1));
        }  
        wbig += wis[lookind][n];
     }
     // loop through possible fs ratios (bias) and compute dfe entries
     for (ispec=0; ispec<nspecs; ispec++) {
       bias = biases[ispec];
       npos=0; nneg = 0; sumpos=0; sumneg=0;
       //sprintf(dfefilename, "data/dfe_N%d_K%d_ts%d_%d_%d_%d.out",N,K,(int)(100.0*ftsave),irep,igen,(int)round(100.0*bias));
       //fpdfe = fopen(dfefilename,"w");
       //if (fpdfe==NULL) { fprintf(stdout,"Unable to open dfe file\n"); exit(1);}
       for (j = 0; j<ndfe; j++) {
         pos = (int)(NB*(rand()/(double)RAND_MAX));
         if ((rand()/(double)RAND_MAX)<bias)
	   transition(seqs[bigtype],pos,dfeseq);
         else
           transversion(seqs[bigtype],pos,dfeseq);
         wtmp = 0;
         for (n=0;n<N;n++) {
           lookind = ctoa[dfeseq[3*n]][dfeseq[3*n+1]][dfeseq[3*n+2]];
           for (k=0;k<K;k++) {
             aaneigh = ctoa[dfeseq[3*neighs[k][n]]][dfeseq[3*neighs[k][n]+1]][dfeseq[3*neighs[k][n]+2]];
             lookind += aaneigh*pow(21.0,(double)(k+1));
           }  
           wtmp += wis[lookind][n];
         }
         wtmp = fmax(wtmp,0)/wbig-1.0;
	 if (wtmp>0) {npos++; sumpos+= wtmp;}
         if (wtmp<0) {nneg++; sumneg+= wtmp;}
         //fprintf(fpdfe,"%f\n",(float)fmax(wtmp,0));
       }  // end loop on j entries of the dfe
       //fclose(fpdfe);
       fbenmu[irep][ispec][itime] = (double)npos/ndfe;
       fdelmu[irep][ispec][itime] = (double)nneg/ndfe;
       if (npos>0) sbarbenmu[irep][ispec][itime] = sumpos/npos;
       else sbarbenmu[irep][ispec][itime]=0;
       if (nneg>0) sbardelmu[irep][ispec][itime] = sumneg/nneg;
       else sbardelmu[irep][ispec][itime] = 0;
     } // end loop across all ts:tv biases for the dfes
     itime++;
   }  // end of if statement for igen being a time to write out a dfe
   // each individual in the next generation is chosen according to
   // the freq of each type in the current generation, weighted by fitness and Ricker factor
   double ricker = exp(1.0 - (double)ntotal/(double)carryingcap);
   //fprintf(stdout, "ricker is %f\n",(float)ricker);
   ntotal = 0;
   for (itype=0;itype<ntypes;itype++) {
     Noftype[itype] = poidev((float)(Noftype[itype]*ricker*ws[itype]/wbar),&seed);
     ntotal += Noftype[itype];
   }
   //fprintf(stdout,"ntotal is %d\n", ntotal);
   // After making the new generation, we allow mutation and recompute fitness
   nweights[0] = (float)Noftype[0]/ntotal;
   for (itype=1;itype<ntypes;itype++)
     nweights[itype] = nweights[itype-1] + (float)Noftype[itype]/ntotal;
   nweights[ntypes-1]=1;  // fix any rounding errors
   // figure out how many new mutations
   double nmutes;
   nmutes = mu*ntotal;
   if (ran1(&seed)<(nmutes -floor(nmutes))) nmutes = ceil(nmutes);
   else nmutes = floor(nmutes);
   for (int imutes = 0;imutes<nmutes;imutes++) {
     float r = ran1(&seed);
     int mutype = 0;
     while (r > nweights[mutype]) mutype++;
     Noftype[mutype]--;  //subtract mutant from its type
     ntypes++;  // add a new type
     if (ntypes>MAXTYPES) { fprintf(stdout,"Too many types!\n");  exit(2);};
     Noftype[ntypes-1]=1;  // add mutant to the new type
     //  now we still need to fill in ws, ftseq and seqs of ntypes-1 
     if (ran1(&seed)<pftmu) {  // mutation affects freq of ts
       //fts = ftseq[mutype] + sdfts*gasdev(&seed);  
         fts = ran1(&seed);  
         if (fts>1.0) fts=1;
         if (fts<0.0) fts=0;
	 ws[ntypes-1] = ws[mutype];
	 ftseq[ntypes-1] = fts;
         for (i=0;i<NB;i++) seqs[ntypes-1][i] = seqs[mutype][i];	 
     }
     else { // mutation occurs in the sequence
      /* mutation(anc,tmpseq,ffixed,fts); */
     pos = (int)(NB*(rand()/(double)RAND_MAX));
     if ((rand()/(double)RAND_MAX)<ftseq[mutype])
       transition(seqs[mutype],pos,seqs[ntypes-1]);
     else
       transversion(seqs[mutype],pos,seqs[ntypes-1]);
     wtmp = 0;
     for (n=0;n<N;n++) {
       lookind = ctoa[seqs[ntypes-1][3*n]][seqs[ntypes-1][3*n+1]][seqs[ntypes-1][3*n+2]];
       for (k=0;k<K;k++) {
         aaneigh = ctoa[seqs[ntypes-1][3*neighs[k][n]]][seqs[ntypes-1][3*neighs[k][n]+1]][seqs[ntypes-1][3*neighs[k][n]+2]];
         lookind += aaneigh*pow(21.0,(double)(k+1));
       } 
       wtmp += wis[lookind][n];
     }
     ws[ntypes-1] = fmax(wtmp,0.0);
     ftseq[ntypes-1] = ftseq[mutype];
     }  // mutation occurred in the sequence
   }  // finish loop on imutes, each mutation
   // remove any types that have zero individuals
   for (itype=0;itype<ntypes;itype++) {
    if (Noftype[itype]==0) {
       jtype = ntypes-1;
       while ( (Noftype[jtype]==0) && jtype>itype) {
         jtype--;
	 ntypes--;
       }
       if (jtype == itype) ntypes--;
       else {
         Noftype[itype] = Noftype[jtype];
         ntypes--;
         ws[itype] = ws[jtype];
         ftseq[itype] = ftseq[jtype];
         for (i=0;i<NB;i++) seqs[itype][i] = seqs[jtype][i];
       }
     }
   }  // loop on itype, removing zero entries
   // compute some summary statistics for the output file, written each generation
   ntotal = 0;  ftbar = 0; wbar = 0;
   for (itype=0;itype<ntypes;itype++) {
     ntotal += Noftype[itype];
     ftbar += ftseq[itype]*Noftype[itype];
     wbar += ws[itype]*Noftype[itype];
   }
   Nnwbarsum[igen][0] += ntotal/carryingcap;
   Nnwbarsum[igen][1] += ntypes;
   Nnwbarsum[igen][2] += wbar/ntotal;
   //fprintf(fpout,"%d %d %f %f\n",ntotal, ntypes, (float)ftbar/ntotal, (float)wbar/ntotal);
   if (ntotal==0) { fprintf(stdout,"Population went extinct at generation %d\n",igen);
                    igen = ngens +1;
		  }
   }  // loop on igen, a generation has been completed
   
   // fclose(fpout);
  }  // loop on irep

 // all the reps are now finished, save the results to files
   for (igen=0;igen<ngens;igen++) fprintf(fpout,"%f %f %f\n",(float)carryingcap*Nnwbarsum[igen][0]/nreps,(float)Nnwbarsum[igen][1]/nreps,(float)Nnwbarsum[igen][2]/nreps);
   fclose(fpout);
   for (itime=0;itime<ntimes;itime++) {
     sprintf(dfefilename, "data/dfe_N%d_K%d_mu%4.0e_ts%d_%d.out",N,K,mu,(int)round(100*ftsave),dfetimes[itime]);
     fpdfe = fopen(dfefilename,"w");
     if (fpdfe==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}  
     int nzb=0; int nzd=0; double sums[5];
     for (ispec=0;ispec<nspecs;ispec++) {
      bias = biases[ispec];
      sums[1] =0; sums[2]=0; sums[3]= 0; sums[4] = 0; nzb=0; nzd=0; 
      for (irep=0;irep<nreps;irep++) {
       sums[1]+= fbenmu[irep][ispec][itime]; sums[2]+= sbarbenmu[irep][ispec][itime];
       sums[3]+= fdelmu[irep][ispec][itime]; sums[4]+= sbardelmu[irep][ispec][itime];
       if (sbarbenmu[irep][ispec][itime]>0) nzb++;
       if (sbardelmu[irep][ispec][itime]<0) nzd++;
      }
      fprintf(fpdfe,"%f %f %f %f %f\n",(float)bias,(float)sums[1]/nreps,(float)sums[2]/nzb,(float)sums[3]/nreps,(float)sums[4]/nzd);
     }  // loop on bias for writing files
     fclose(fpdfe);
   } // loop on itime for writing files
}  // end of main

void transition(int seq[MAXN],int pos,int newseq[MAXN])
{
for (int i=0;i<NB;i++) newseq[i]=seq[i];
switch(seq[pos]) {
  case 0  : newseq[pos] = 1; break;
  case 1  : newseq[pos] = 0; break;
  case 2  : newseq[pos] = 3; break;
  case 3  : newseq[pos] = 2; break;
}
}

void transversion(int seq[MAXN],int pos,int newseq[MAXN])
{
for (int i=0;i<NB;i++) newseq[i]=seq[i];
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
  
for (i=0;i<NB;i++) newseq[i]=seq[i];
 muprob = rand()/(double)RAND_MAX;
 if (muprob<ffixed) {  // we want AT to AT or GC to GC
  pos = (int)(NB*(rand()/(double)RAND_MAX));
  switch(seq[pos]) {
    case 0: newseq[pos] = 1; break;
    case 1: newseq[pos] = 0; break;
    case 2: newseq[pos] = 3; break;
    case 3: newseq[pos] = 2; break;
   }
  }
 else {
  if (muprob<(ffixed+fts)) {   // we want GC to AT
   pos = (int)(NB*(rand()/(double)RAND_MAX));
   while (seq[pos]<2)  pos = (int)(NB*(rand()/(double)RAND_MAX));
   if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 0;
   else newseq[pos] = 1;
  }
  else {// we want AT to GC
   pos = (int)(NB*(rand()/(double)RAND_MAX));
   while (seq[pos]>1)  pos = (int)(NB*(rand()/(double)RAND_MAX));
   if ((rand()/(double)RAND_MAX)<0.5) newseq[pos] = 2;
   else newseq[pos] = 3;
  }
 }
}

void initialize_ctoa(int ctoa[4][4][4])
{
  ctoa[0][0][0] = 0;
  ctoa[0][0][1] = 0;
  ctoa[0][0][2] = 1;
  ctoa[0][0][3] = 1;
  ctoa[0][1][0] = 2;
  ctoa[0][1][1] = 2;
  ctoa[0][1][2] = 2;
  ctoa[0][1][3] = 2;
  ctoa[0][2][0] = 3;
  ctoa[0][2][1] = 3;
  ctoa[0][2][2] = 4;
  ctoa[0][2][3] = 4;
  ctoa[0][3][0] = 5;
  ctoa[0][3][1] = 5;
  ctoa[0][3][2] = 4;
  ctoa[0][3][3] = 6;

  ctoa[1][0][0] = 1;
  ctoa[1][0][1] = 1;
  ctoa[1][0][2] = 1;
  ctoa[1][0][3] = 1;
  ctoa[1][1][0] = 7;
  ctoa[1][1][1] = 7;
  ctoa[1][1][2] = 7;
  ctoa[1][1][3] = 7;
  ctoa[1][2][0] = 8;
  ctoa[1][2][1] = 8;
  ctoa[1][2][2] = 9;
  ctoa[1][2][3] = 9;
  ctoa[1][3][0] = 10;
  ctoa[1][3][1] = 10;
  ctoa[1][3][2] = 10;
  ctoa[1][3][3] = 10;

  ctoa[2][0][0] = 11;
  ctoa[2][0][1] = 11;
  ctoa[2][0][2] = 11;
  ctoa[2][0][3] = 12;
  ctoa[2][1][0] = 13;
  ctoa[2][1][1] = 13;
  ctoa[2][1][2] = 13;
  ctoa[2][1][3] = 13;
  ctoa[2][2][0] = 14;
  ctoa[2][2][1] = 14;
  ctoa[2][2][2] = 15;
  ctoa[2][2][3] = 15;
  ctoa[2][3][0] = 2;
  ctoa[2][3][1] = 2;
  ctoa[2][3][2] = 10;
  ctoa[2][3][3] = 10;

  ctoa[3][0][0] = 16;
  ctoa[3][0][1] = 16;
  ctoa[3][0][2] = 16;
  ctoa[3][0][3] = 16;
  ctoa[3][1][0] = 17;
  ctoa[3][1][1] = 17;
  ctoa[3][1][2] = 17;
  ctoa[3][1][3] = 17;
  ctoa[3][2][0] = 18;
  ctoa[3][2][1] = 18;
  ctoa[3][2][2] = 19;
  ctoa[3][2][3] = 19;
  ctoa[3][3][0] = 20;
  ctoa[3][3][1] = 20;
  ctoa[3][3][2] = 20;
  ctoa[3][3][3] = 20;
}
