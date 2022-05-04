# MutationSpectrum

Contains the C code to run the adaptation simulations in Sane et al. 2022

The files are:

tstv.c   -- This is the basic adaptive walk simulation where a wildtype individual has a fixed mutation bias, 
and the DFE is created at different points along the walk, for the WT and also for possible mutants that might
have a different bias from the WT. The bias can be either transition:transversion bias or GC->AT bias.
The fitness landscape is an NK model with variable N (sequence length) and K=1.  This version of the code has 
extensive comments.  Start here.

tstv_kany.c -- An extension of tstv.c that allows K in the NK model to be specified.

mm.c -- An extension of tstv_kany.c.  In an adaptive walk, the population is assumed to be genetically 
homogeneous; one beneficial mutation occurs and fixes (or is lost) before the next occurs.  
This assumption that is relaxed here.  Now a *population* of individuals adapt, such
that multiple mutations are simultaneously segregating.

aa.c -- An extention of mm.c in which mutations occur at the base pair level (as they did in the previous
code) but in aa.c the bases are translated to amino acids before computing fitness.  Fitness depends only 
on the amino acid sequence.  Thus synonymous mutations have no effect on fitness.  The fitness landscape is
still an NK fitness model but in aa.c, N refers to the length of the amino acid sequence.  The length of the base pair
sequence is 3*N.

program_name.makefile -- each of the C programs also has a makefile which can be used to compile it.

Note that some of the programs above (as noted in their makefiles) depend on subroutines from Numerical Recipes in C (ANSI-C edition).
The code from Numerical Recipes is proprietory and is not reproduced here.
