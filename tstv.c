#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define MAXN 251        //max genome length
#define MAXSTEPS 500    //max number of steps in adaptive walk
#define MAXMUTES 500    //max number of mutations in DFE
#define MAXWALKS 600    //max number of adaptive walks 
#define MAXSPECS 100    //max number of different ts:tv values for output files

void transition(int seq[], int n, int newseq[]);
void transversion(int seq[], int n, int newseq[]);
void mutation(int seq[],int newseq[],double ffixed, double fts);

int N=200;  // genome length, global variable

void main() 
{

int gc_flag = 1;  // if this flag = 0, run ts:tv case; if 1, run at:gc case
int nsteps;  // nsteps is number of accepted steps in adaptive walk
int maxwalklength = 1e4;  // maxwalklength is how many possible new mutations
                           // are generated and tested during the adaptive walk.
                           // These are only accepted if s>0 and if a random
                           // number is less than 2s (mimicking fixation).
			   // Thus "walklengths" are much longer than "steps".
double ftstep=.1;  // increment for fraction of transitions
double ffixed = .1; // fixed fraction of GC-GC or AT-AT mutations in AT:GC case
int nwalks = 500; // number of adaptive walks
double fts, ftsmu;  // Ts fraction in wt and mutant
double fwt[MAXWALKS];  // final wt fitness at the end of the walk, for each walk
double dfe[MAXWALKS][MAXMUTES];  // s values (DFE) for all mutants for each walk 
int nstepsall[MAXWALKS];  // number of steps taken in each walk
int i,j,k, ift, iwalk, ispec, itry, pos;  // temporary counters 
double oneoverNe = .00001;  // deleterious mutations can also fix with
  // probability exp(2s)*oneoverNe i.e. exp(2s)/Ne for Ne = 1e5 here
  // This happened so rarely results were indistinguishable whether this was
  // allowed or commented out.
double wis[4][4][MAXN];  // w_i values for each locus, given base at that locus
                        //  and base at the epistatically coupled locus
int nmutes = 300;   // number of fitnesses to use to calculate DFE
int neighs[MAXN];   // locus of epistatic neighbour for each locus in sequence
int anc[MAXN], tmpseq[MAXN];  // ancestor sequence, temporary sequence
double w0, wseq[MAXSTEPS];  // initial fitness, fitness at each step;
char filename[100];  // output file name
FILE *fpout;   // output file pointer
double fben[MAXWALKS], fdel[MAXWALKS];  // fractions of ben and del mutns in wt DFE
double sbarben[MAXWALKS],sbardel[MAXWALKS]; // mean s_ben and s_del in wt DFE
double fbenmu[MAXWALKS][MAXSPECS], fdelmu[MAXWALKS][MAXSPECS];  // as above for mutants
double sbarbenmu[MAXWALKS][MAXSPECS],sbardelmu[MAXWALKS][MAXSPECS]; // ditto
// note that for mutant variables we keep track for mutants of "MAXSPECS" possible
// ts fractions
double sumpos, sumneg;  // housekeeping, counting, arithmetic
int npos, nneg, imute, ilength;
double sums[20], dfemu;
double wtmp, s;

// Edit the line below (several examples shown) to set the walklengths required.
// Use only one walklengths and nlengths=1 for waterfall plots, for example.
// Note that this code is extremely inefficient but runs a completely new independent
// set of walks for each length.  So if you have walklengths = {1e4, 2e4} for example,
// the code will generate the required number of walks to length 1e4, then start over
// and generate NEW walks from time 0 to time 2e4, rather than continue the walks used
// for 1e4.  This could be easily changed.

// int walklengths[1] = {2.5e4}, nlengths = 1;
 int walklengths[6] = {1e4, 2e4, 5e4, 1e5, 2e5, 5e5}, nlengths = 6;
// int walklengths[12] = {3e4, 4e4, 6e4, 7e4, 8e4, 9e4, 12e4, 14e4, 16e4, 18e4, 3e5, 4e5}, nlengths = 12;

srand(1);  // initialize random number generator

if (gc_flag)  fprintf(stdout,"NOTE: ts = bias GC -> AT (0,1), tv = bias AC->GC (2,3), ffixed = bias AT->AT or GC->GC fixed to %f\n",ffixed); 

for (ilength=0;ilength<nlengths;ilength++) {  // loop over walk lengths
   maxwalklength = walklengths[ilength];  // how many adaptive steps will be accepted

    // uncomment only one of the two lines below, depending on whether we want
    // to investigate a range of wt fts (waterfall plot), or a single fixed wt fts
// for (fts=ftstep;fts<=1-ffixed;fts+=ftstep) {
for (fts=0.35; fts<=0.35; fts+=ftstep) {

// put your output filename here, in this example "Sw" means s-weighted walk
 sprintf(filename, "gc_data/N%d_Sw_walk%d_fts%1.2f.out",N,maxwalklength,fts);
 fpout = fopen(filename,"w");
  
 for (iwalk=0;iwalk<nwalks;iwalk++) {  // loop over walks
   // first, create a new landscape for this walk, by setting the w_i and neighbours
   for (i=0;i<=3;i++) for (j=0;j<=3;j++)
      for (k=0;k<N;k++) wis[i][j][k] = (rand()/(double)RAND_MAX);
   for (k=0;k<N;k++) neighs[k] = (int)(N*(rand()/(double)RAND_MAX));
   // pick a new random ancestor sequence
   for (k=0;k<N;k++) anc[k] = (int)(4*(rand()/(double)RAND_MAX));
   // compute the fitness of the ancestor 
   w0 = 0;
   for (k=0;k<N;k++) 
     w0 += wis[anc[k]][anc[neighs[k]]][k];
   // step 0.  fitness at step 0 is w0, and nsteps = 0.
   wseq[0] = w0;
   nsteps = 0;
   // now try to make a new step in the walk, up to maxwalklength number of attempts
   for (itry = 1; itry<=maxwalklength; itry++) {
     // use either "mutation" or "transition/transversion" depending on which case.
     // Since the gc case is more complicated, we move it all to the subroutine.
     if (gc_flag) mutation(anc,tmpseq,ffixed,fts);
     else {
       pos = (int)(N*(rand()/(double)RAND_MAX));  // find position to mutate
       if ((rand()/(double)RAND_MAX)<fts)  // transition with probability fts
         transition(anc,pos,tmpseq);
       else
       transversion(anc,pos,tmpseq);
     }
     // compute fitness of the mutated sequence,  tmpseq
     wtmp = 0;
     for (k=0;k<N;k++) 
       wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
     // and then compute s, relative to the current wt fitness
     s = wtmp/wseq[nsteps]-1;
    if (s>0) {  // if beneficial mutation
      if ((rand()/(double)RAND_MAX)<2*s) {  // if it survives drift (prob 2s)
        nsteps = nsteps+1;  // we accept this new sequence as a step in the walk
        for (k=0;k<N;k++) anc[k] = tmpseq[k]; // new sequence replaces the ancestor
        wseq[nsteps] = wtmp;  // new fitness gets added to the list of fitnesses on this walk
      }
    }  else {  // if the mutation is deleterious or neutral, we might still accept it
    // with probability exp(2s)/Ne
      if ((rand()/(double)RAND_MAX)<oneoverNe*exp(2*s)) {
       nsteps = nsteps+1;
       for (k=0;k<N;k++) anc[k] = tmpseq[k];
       wseq[nsteps] = wtmp;
      }
   }
   }  // finish loop on itry, we are now finished this adaptive walk
   
   fwt[iwalk] = wseq[nsteps];  // save the final fitness of the wt at end of walk
   // now make the DFE for the wt
   npos=0; nneg = 0; sumpos=0; sumneg=0;
   for (imute=0;imute<nmutes;imute++) {  // loop over new mutations to add to wt DFE
   // first, make a single substitution in the wt sequence, as described above
     if (gc_flag) mutation(anc,tmpseq,ffixed,fts);
     else {
       pos = (int)(N*(rand()/(double)RAND_MAX));
       if ((rand()/(double)RAND_MAX)<fts)
        transition(anc,pos,tmpseq);
       else
         transversion(anc,pos,tmpseq);
     }
     // then compute the fitness and then s for the mutated sequence, tmpseq
     wtmp = 0;
     for (k=0;k<N;k++) 
       wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
       // add this s value to the dfe for the wt for this walk
     dfe[iwalk][imute] = wtmp/fwt[iwalk] - 1;
     // keep track of how many beneficial and deleterious mutations in this dfe
     if (dfe[iwalk][imute]>0) {npos++; sumpos+= dfe[iwalk][imute];}
     if (dfe[iwalk][imute]<0) {nneg++; sumneg+= dfe[iwalk][imute];}
   }

   nstepsall[iwalk] = nsteps;  // keep track of how many accepted steps in this walk
   fben[iwalk] = (double)npos/nmutes;  // fraction beneficial in wt DFE
   fdel[iwalk] = (double)nneg/nmutes;  // fraction deleterious in wt DFE
   if (npos>0) sbarben[iwalk] = sumpos/npos; // mean s of beneficial
   else sbarben[iwalk]=0;
   if (nneg>0) sbardel[iwalk] = sumneg/nneg;  // mean s of deleterious
   else sbardel[iwalk]=0;

  // now compute dfe for a mutator with a different ts fraction from the wt
  // The mutator has the same sequence as the wt, and differs only in ts:tv bias
  // (or gc:at bias).
  ispec = 0;
  // loop over ts fraction of mutator, edit the line below for gc:at case as needed
  //  for (ftsmu=ftstep;ftsmu<=1-ftstep;ftsmu+=ftstep) { 
  for (ftsmu=0.05; ftsmu<=1-ffixed; ftsmu+=ftstep) {  // example for gc:at case
    // make a dfe for this mutator, exactly as described above for the wt
    npos=0; nneg = 0; sumpos=0; sumneg=0;
    for (imute=0;imute<nmutes;imute++) { 
      if (gc_flag) mutation(anc,tmpseq,ffixed,ftsmu);
      else {
        pos = (int)(N*(rand()/(double)RAND_MAX));
        if ((rand()/(double)RAND_MAX)<ftsmu)
          transition(anc,pos,tmpseq);
        else
          transversion(anc,pos,tmpseq);
      }
     wtmp = 0;
     for (k=0;k<N;k++) 
       wtmp += wis[tmpseq[k]][tmpseq[neighs[k]]][k];
     dfemu = wtmp/fwt[iwalk] - 1;
     if (dfemu>0) {npos++; sumpos+= dfemu;}
     if (dfemu<0) {nneg++; sumneg+= dfemu;}
    }
   fbenmu[iwalk][ispec] = (double)npos/nmutes;
   fdelmu[iwalk][ispec] = (double)nneg/nmutes;
   if (npos>0) sbarbenmu[iwalk][ispec] = sumpos/npos;
   else sbarbenmu[iwalk][ispec]=0;
   if (nneg>0) sbardelmu[iwalk][ispec] = sumneg/nneg;
   else sbardelmu[iwalk][ispec] = 0;
   ispec++;
   }  // end of loop on ispec (stepping through mutator ts values)
   // at this point, we have completed an adaptive walk across this landscape
   // We have created a dfe for the wt
   // We have created a dfe for a set of mutators that have the same sequence
   // as the wt, but differ in their ts fraction.
   // Now we repeat this be creating a new landscape, new ancestor, and walking again.
}  // end of loop on iwalks

// keep track of the average fraction beneficial, etc, across all walks
sums[1] =0; sums[2]=0; sums[3]= 0; sums[4] = 0; sums[5]=0; sums[6]=0;
int nzb=0; int nzd=0;
for (iwalk=0;iwalk<nwalks;iwalk++) {
  sums[1]+= fben[iwalk]; sums[2]+= sbarben[iwalk];
  sums[3]+= fdel[iwalk]; sums[4]+= sbardel[iwalk];
  sums[5]+= fwt[iwalk];  sums[6]+= nstepsall[iwalk];
  if (sbarben[iwalk]>0) nzb++;  // number of beneficial mean s values that are non-zero
  if (sbardel[iwalk]<0) nzd++; // ditto, deleterious
}
// print results to output file, first for wt and then for mutants
// note first row of output file contains: fts, ftsstep, fwt, nsteps
// following rows contain fben, sbarben, fdel, sbardel for wt and then all mutators
fprintf(fpout,"%f %f %f %f\n",fts,ftstep,sums[5]/nwalks,sums[6]/nwalks);
fprintf(fpout,"%f %f %f %f\n",sums[1]/nwalks,sums[2]/nzb,sums[3]/nwalks,sums[4]/nzd);

for (i=0;i<ispec;i++) {
sums[1] =0; sums[2]=0; sums[3]= 0; sums[4] = 0; nzb=0; nzd=0; 
for (iwalk=0;iwalk<nwalks;iwalk++) {
  sums[1]+= fbenmu[iwalk][i]; sums[2]+= sbarbenmu[iwalk][i];
  sums[3]+= fdelmu[iwalk][i]; sums[4]+= sbardelmu[iwalk][i];
  if (sbarbenmu[iwalk][i]>0) nzb++;
  if (sbardelmu[iwalk][i]<0) nzd++;
}
fprintf(fpout,"%f %f %f %f\n",sums[1]/nwalks,sums[2]/nzb,sums[3]/nwalks,sums[4]/nzd);
}

fclose(fpout);
// uncomment below for more detailed output as required
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
}  // loop on fts
 }  // loop on walklength

}

void transition(int seq[MAXN],int pos,int newseq[MAXN])
// using standard codon model ordering:
// T = 0, C = 1
// A = 2, G = 3
{
for (int i=0;i<N;i++) newseq[i]=seq[i]; // first make a copy of seq into newseq
// then make a transition at position pos 
switch(seq[pos]) {
  case 0  : newseq[pos] = 1; break;
  case 1  : newseq[pos] = 0; break;
  case 2  : newseq[pos] = 3; break;
  case 3  : newseq[pos] = 2; break;
}
}

void transversion(int seq[MAXN],int pos,int newseq[MAXN])
// give sequence seq, make a random transversion at position pos to form newseq
{
for (int i=0;i<N;i++) newseq[i]=seq[i];  // first copy the seq to newseq
// then make a transversion at position pos (depends on base at position pos)
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
// given some fixed AT->AT and GC->GC fraction "ffixed",
// and GC->AT fraction "fts"
// make a mutation to seq to form newseq
{
  int i, pos;
  double muprob;
  
for (i=0;i<N;i++) newseq[i]=seq[i];  // first copy seq to newseq
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
   while ((seq[pos]==0)||seq[pos]==2)  pos = (int)(N*(rand()/(double)RAND_MAX));
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
