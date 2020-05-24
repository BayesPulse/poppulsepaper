/*******************************************************************/
/*************************pop_mcmc.c *****************************/
/*******************************************************************/

#include "pop_deconvolution_main.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 mmm: Order statistic used for distribution of pulse locations.
      This is inputted by the user and is typically 3
 adelta: counter that increments when a halflife estimate is accetped in the MH draw
 ndelta: counter that increments for each loop of the halflife MH draw
 atime: counter that increments when a pulse location estimate is accetped in the MH draw
 ntime: counter that increments for each loop of the pulse location MH draw
 arem: counter that increments when an individual pulse mass estimate is accetped in the MH draw
 nrem: counter that increments for each loop of the individual pulse mass MH draw
 arew: counter that increments when an individual pulse width estimate is accepted in the MH draw
 nrew: counter that increments for each loop of the individual pulse width MH draw
 afepmv: counter that increments when an overall mean pulse mass standard deviation (sma) is accetped in the MH draw
 nfepmv: counter that increments for each loop of the overall mean pulse mass standard deviation (sma) MH draw
 afepwv: counter that increments when an overall mean pulse width standard deviation (smw) is accetped in the MH draw
 nfepwv: counter that increments for each loop of the overall mean pulse width standard deviation (smw) MH draw
 afemv: counter that increments when an overall pulse mass standard deviation (sa) is accetped in the MH draw
 nfemv: counter that increments for each loop of the overall pulse mass standard deviation (sa) MH draw
 afewv: counter that increments when an overall pulse width standard deviation (sw) is accetped in the MH draw
 nfewv: counter that increments for each loop of the overall pulse width standard deviation (sw) MH draw
 afebv: counter that increments when a baseline standard deviation (sb) is accetped in the MH draw
 nfebv: counter that increments for each loop of the baseline standard deviation (sb) MH draw
 afehv: counter that increments when a half-life standard deviation (sh) is accetped in the MH draw
 nfehv: counter that increments for each loop of the half-life standard deviation (sh) MH draw
 fitstart: The first time in hours that a pulse may occur
 fitend: The last time in hours that a pulse may occur

*********************************************************************/

extern int mmm;
long adelta = 0;
long ndelta = 0;
long atime = 0;
long ntime = 0;
long arem = 0;
long nrem = 0;
long arew = 0;
long nrew = 0;
long afepmv = 0;
long nfepmv = 0;
long afepwv = 0;
long nfepwv = 0;
long afemv = 0;
long nfemv = 0;
long afewv = 0;
long nfewv = 0;
long afebv = 0;
long nfebv = 0;
long afehv = 0;
long nfehv = 0;
extern double fitstart;
extern double fitend;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 mcmc 
 draw_fe_priors
 draw_fe_priors2
 draw_bh_mean
 draw_bh_var
 draw_fixed_effects
 draw_fe_precision
 draw_bh
 draw_times
 draw_random_effects
 error_squared
 adjust_acceptance
 adjust2_acceptance
 
**********************************************************************/


/*********************************************************************/
                    /*START OF MCMC SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mcmc: this runs the BDMCMC process
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects that exist;
               Common_parms *parms; the current values of the common parameters;
               double **ts; this is the matrix of observed data (a column of 
                  times and S columns of log(concentration);
               long iter; the number of iterations to run;
               int N; the number of observations in **ts;
               Priors *priors; the parameters of the prior distributions;
               unsigned long *seed; seed values needed for the randon number generator;
               char *file1; the output file name for common parameters
               double propvar[]; vector of proposal variances input by user
               Hyper_priors *hyper; hyperpriors input into algorithm by user
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k,l: generic counters
 num_node: counter of num_node for loops
 num_node2: number of pulses
 NN=50: Output every NNth iteration to files
 NNN=5000: Output every NNNth interation to screen or output file
 vtime: proposal variance for individual pulse locations
 vrem: proposal variance for individual masses
 vrew: proposal variance for individual widths
 vfepmv: proposal variance for standard deviation of overall mean mass, sma
 vfepwv: proposal variance for standard deviation of overall mean width, smw
 vfemv: proposal variance for standard deviation of pulse mass, sa
 vfewv: proposal variance for standard deviation of pulse width, sw
 vfebv: proposal variance for standard deviation of baseline, sb
 vfehv: proposal variance for standard deviation of half-life, sh
 ssq: sum of squared differences between log(concentration) and expected value
 **pmd_var: proposal variance-covariance matrix for baseline and halflife
 **pmd_vch: cholesky decomposed proposal variance-covariance matrix for
    baseline and halflife
 *new_node: index used to go through nodes in file output
 *subject: index used to go through subjects in file output
 *common: name given to inputted file name for common parameters
 *parm: name given to inputted file name for pulse parameters

 SUBROUTINES USED
 inverse_gamma: found in randgen.c; draws from the inverse gamma distribution
 error_squared: found in this file; calculates sum of squared error
 draw_bh: found in this file; draws baseline and half-life
 draw_times: found in this file; draws individual pulse locations
 birth_death: found in birthdeath.c; runs birthdeath process
 draw_fixed_effects: found in this file; draws subject specific mean mass and width
 draw_fe_precision: found in this file; draws st dev of pulse mass and width
 draw_random_effects: found in this file; draws individual masses and widths
 draw_bh_mean: found in this file; draws mean of baseline and half-life
 draw_bh_var: found in this file; draws st dev of baseline and half-life
 draw_fe_priors: found in this file; draws mean of mean pulse mass and width
 draw_fe_priors2: found in this file; draws st dev of mean pulse mass and width
 adjust_acceptance: found in this file; adjusts proposal variance
 adjust_acceptance2: found in this file; adjusts 2-d proposal variance*/
 /***********************************************************************/

void mcmc(Subject_type *sublist,Common_parms *parms,double **ts,long iter,int N,
          Priors *priors,unsigned long *seed,char *file1,double propvar[],Hyper_priors *hyper)
{
 int i,j,k,l,num_node,num_node2,NN=50,NNN=5000;
 double vfepmv,vfepwv,vrem,vrew,vfemv,vfewv,vfebv,vfehv,vtime,ssq,**pmd_var,**pmd_vch;
 Node_type *new_node;
 Subject_type *subject;
 FILE *common,*parm;

 double inverse_gamma(double alpha,double beta,unsigned long *seed);
 double error_squared(double **ts,Subject_type *list,Common_parms *parms,int N);
 void draw_bh(Subject_type *,Common_parms *,Priors *,double **,
           int,unsigned long *,double **);
 void draw_times(Subject_type *,Common_parms *,double **,
                  int,unsigned long *,double);
 void birth_death(Subject_type *list,double **ts,Common_parms *parms,int N,
                  unsigned long *seed,int);
 void draw_fixed_effects(Subject_type *,Priors *,Common_parms *,unsigned long *);
 void draw_fe_precision(Subject_type *,Priors *,Common_parms *,double,double,unsigned long *);
 void draw_random_effects(double **ts,Subject_type *,Common_parms *,int,double,double,unsigned long *);
 void draw_tscale(Node_type *, Common_parms *, double, double,
              long *, long * , long *, long *);
 void draw_bh_mean(Subject_type *,Priors *,Common_parms *,unsigned long *,Hyper_priors *);
 void draw_bh_var(Subject_type *,Priors *,double,double,unsigned long *,Hyper_priors *);
 void draw_fe_priors(Subject_type *,Priors *,Common_parms *,unsigned long *,Hyper_priors *);
 void draw_fe_priors2(Subject_type *,Priors *,double,double,unsigned long *,Hyper_priors *);
 void adjust_acceptance(double,double *);
 void adjust2_acceptance(double,double **,double);
    int cholesky_decomp(double **,int);

 strcat(file1,".out");
 common = fopen(file1,"w");

 subject=sublist->succ;
 while (subject != NULL){
     subject->csub = fopen(subject->common,"w");
     subject->psub = fopen(subject->pulse,"w");
     subject=subject->succ;
    }

 /*proposal variances for variances of subject mean mass and width (sigma-ma and sigma-mw)*/
 vfepmv = propvar[0];
 vfepwv = propvar[3];

 /*proposal variances for variances of overall mean mass and width (sigma-a and sigma-w)*/
 vfemv = propvar[1];
 vfewv = propvar[4];

 /*proposal variances for individual masses and widths*/
 vrem = propvar[2];
 vrew = propvar[5];

 /*proposal variances for variances of baseline and halflife*/
 vfebv = propvar[6];
 vfehv = propvar[8];

 /*proposal variances for subject baselines and halflifes*/
 pmd_var = (double **)calloc(2,sizeof(double *));
 for (i=0;i<2;i++)
    pmd_var[i] = (double *)calloc(2,sizeof(double));

 pmd_var[0][0] = propvar[7];
 pmd_var[1][1] = propvar[9];
 pmd_var[0][1] = pmd_var[1][0] = -0.90*sqrt(pmd_var[0][0])*sqrt(pmd_var[1][1]);

 /*proposal variance for pulse locations*/
 vtime = propvar[10];

/***cholesky decompose the proposal var-covar matrix for b and hl***/
    pmd_vch = (double **)calloc(2,sizeof(double *));
    for (i=0;i<2;i++)
      pmd_vch[i] = (double *)calloc(2,sizeof(double));
      
for (i=0;i<2;i++)
    for (j=0;j<2;j++)
        pmd_vch[i][j] = pmd_var[i][j];

    if (!cholesky_decomp(pmd_vch,2)){
      printf("not PSD matrix A\n");
      exit(0);
    }
       
 for (i=0;i<iter;i++) {

/* run the birth-death algorithm to select pulses */
   
   parms->iter=i;

   birth_death(sublist,ts,parms,N,seed,i);

   /* draw fixed effects priors; mua and muw*/

   draw_fe_priors(sublist,priors,parms,seed,hyper);

   /* draw fixed effects prior variances; sma and smw*/

   draw_fe_priors2(sublist,priors,vfepmv,vfepwv,seed,hyper);

   /* draw mean baseline and halflife; mub and muh*/

   draw_bh_mean(sublist,priors,parms,seed,hyper);

   /* draw variance for baseline and halflife; sb and sh*/

   draw_bh_var(sublist,priors,vfebv,vfehv,seed,hyper);

   /* draw subject specific means; muak and muwk*/

   draw_fixed_effects(sublist,priors,parms,seed);

   /* draw mass and width variances; sa and sw*/

   draw_fe_precision(sublist,priors,parms,vfemv,vfewv,seed);

   /* draw subject specific baselines and halflifes; Bk and Hk*/

   draw_bh(sublist,parms,priors,ts,N,seed,pmd_vch);

   /* draw pulse locations; tauki*/

   draw_times(sublist,parms,ts,N,seed,vtime);

   /* draw pulse mass and width; Aki and s2pki*/

   draw_random_effects(ts,sublist,parms,N,vrem,vrew,seed);

   /* draw model error; s2e*/
   ssq = error_squared(ts,sublist,parms,N);
   parms->sigma = inverse_gamma(priors->alpha + (parms->numsub*N)/2,priors->beta + 0.5*ssq,seed);
   parms->lsigma = log(parms->sigma);

   fflush(stdout);

/***************************************************************/

   /* Every 50th iteration, print estimates to files*/
   if (!(i%NN)) {
     subject=sublist->succ;
     while (subject != NULL){
         num_node = 0;
         new_node = subject->list->succ;
         while (new_node !=NULL){
             num_node++;
             new_node=new_node->succ;
         }
         /*Print subject specific parameters to c1sk.out files*/
         fprintf(subject->csub,"%d %d %lf %lf %lf %lf\n",i,num_node,subject->theta[0],subject->theta[1],subject->basehalf[0],subject->basehalf[1]);
         fflush(subject->csub);
         new_node = subject->list->succ;
         num_node2=0;
         /*Print pulse specific parameters to p1sk.out files*/
         while (new_node != NULL) {
            fprintf(subject->psub,"%d %d %d %lf %lf %lf\n",i,num_node,num_node2,new_node->theta[0],new_node->theta[1],new_node->time);
            num_node2++;
            new_node = new_node->succ;
            }
         
         fflush(subject->psub);
         subject=subject->succ;

     }
     /*Print common parameters to c1.out file*/
     fprintf(common,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,priors->fe_mean[0],priors->fe_mean[1],priors->fe_precision[0],priors->fe_precision[1],parms->re_precision[0],parms->re_precision[1],priors->meanbh[0],priors->meanbh[1],priors->varbh[0],priors->varbh[1],parms->sigma);
     fflush(common);
     fflush(stdout);
   }

   /*Print acceptance rates*/
   if (!(i%NNN)) {
     printf("\n\n");
     printf("iter = %d \n",i);
     printf("mu_b %.2lf mu_h %.4lf mu_a %.2lf mu_w %.2lf  v %.4le\n",
            priors->meanbh[0],priors->meanbh[1],priors->fe_mean[0],priors->fe_mean[1],parms->sigma);
     subject=sublist->succ;
     printf("pct s_ma = %.2lf pct s_a = %.2lf pct A_ki = %.2lf \n",
            (double)afepmv/(double)nfepmv,(double)afemv/(double)nfemv,(double)arem/(double)nrem);
     printf("pct s_mw = %.2lf pct s_w = %.2lf pct s2p_ki = %.2lf \n",
            (double)afepwv/(double)nfepwv,(double)afewv/(double)nfewv,(double)arew/(double)nrew);
     printf("pct s_b = %.2lf pct s_h = %.2lf pct B-HL = %.2lf pct time = %.2lf \n",
            (double)afebv/(double)nfebv,(double)afehv/(double)nfehv,(double)adelta/(double)ndelta,(double)atime/(double)ntime);

     fflush(stdout);
   }

   /*Every 500th iteration until 25000, update proposal variances*/

   if (!(i%500) && i<25000 && i >0) {

      adjust2_acceptance((double)adelta/(double)ndelta,pmd_var,-0.90);
      adjust_acceptance((double)atime/(double)ntime,&vtime);
      adjust_acceptance((double)arem/(double)nrem,&vrem);
      adjust_acceptance((double)arew/(double)nrew,&vrew);
      adjust_acceptance((double)afepwv/(double)nfepwv,&vfepwv);
      adjust_acceptance((double)afepmv/(double)nfepmv,&vfepmv);
      adjust_acceptance((double)afemv/(double)nfemv,&vfemv);
      adjust_acceptance((double)afewv/(double)nfewv,&vfewv);
      adjust_acceptance((double)afebv/(double)nfebv,&vfebv);
      adjust_acceptance((double)afehv/(double)nfehv,&vfehv);

/*After adjusting proposal variances, reset the counters for acceptance rates*/
      adelta = ndelta = 0;
      atime = ntime = 0;
      arem = nrem = 0;
      arew = nrew = 0;
      afepmv = nfepmv = 0;
      afepwv = nfepwv = 0;
      afemv = nfemv = 0;
      afewv = nfewv = 0;
      afebv = nfebv = 0;
      afehv = nfehv = 0;

      /***cholesky decompose the proposal var-covar matrix***/
 for (k=0;k<2;k++)
    for (l=0;l<2;l++)
        pmd_vch[k][l] = pmd_var[k][l];

    if (!cholesky_decomp(pmd_vch,2)){
      printf("pmd not PSD matrix\n");
      exit(0);
    }
           
    } /*End of adjusting acceptance if-then statement*/

 } /*End of loop through iterations*/

 fclose(common);
 fclose(parm);
} /*End of MCMC */

/*********************************************************************/
               /*START OF draw_fe_priors SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_priors: this runs the M-H sampler draw for overall mean mass and width;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Priors *priors; the current values of the prior parameters;
               Common_parms *parms; the current values of the common parameters;
               unsigned long *seed; seed values needed for the randon number generator;
               Hyper_priors *hyper; parameters of the prior distributions;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 fesum: sum of subject specific mean amplitudes or widths;
        part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics

SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution  */
/**************************************************************************/

void draw_fe_priors(Subject_type *sublist,Priors *priors,Common_parms *parms,unsigned long *seed,Hyper_priors *hyper)
{
/* declare variables */
    int j,t;
    double *tmp1,*new_mean,prop_new,prop_old,psum,pcomp,prior_old,prior_new;
    double prior_ratio,prop_ratio,temp,alpha;
    Subject_type *subject;
    double kiss(unsigned long *);
    double rnorm(double,double,unsigned long *);

    /*Allocate Memory*/
    new_mean = (double *)calloc(2,sizeof(double));
    tmp1 = (double *)calloc(2,sizeof(double));

    /*Add 1 to the counters for acceptance rates of pop mean pulse mass and width*/
    nfepmv++;
    nfepwv++;

    /*Assign current acceptance counts to temporary vector*/
    tmp1[0] = afepmv;
    tmp1[1] = afepwv;
    
/*OLD MAY NOT NEED*/
   int j;
   double fesum,gmean,gvar;
   Subject_type *subject;
/***********************************/
    
/*declare functions */
   double rnorm(double,double,unsigned long *);
   double kiss(unsigned long *);

    for (j=0;j<2;j++){

        fesum=0;
        subject=sublist->succ;

        /*This while loop goes through the pulses to get the sum involved*/
        while(subject != NULL){
            fesum+=subject->theta[j];
            subject=subject->succ;
            } /*end of loop through pulses*/

        /*Gibbs sampler mean*/
        gmean = (hyper->hvar[j]*fesum) + (priors->fe_precision[j]*priors->fe_precision[j]*hyper->hmean[j]);
        gmean /= (parms->numsub*hyper->hvar[j])+(priors->fe_precision[j]*priors->fe_precision[j]);

        /*Gibbs sampler variance*/
        gvar = hyper->hvar[j]*priors->fe_precision[j]*priors->fe_precision[j];
        gvar /= (parms->numsub*hyper->hvar[j])+(priors->fe_precision[j]*priors->fe_precision[j]);

        /*This draws the new value for mua/muw using a Gibbs sampler*/
        priors->fe_mean[j] = rnorm(gmean,sqrt(gvar),seed);
    }

}

/*********************************************************************/
               /*START OF draw_fe_priors2 SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_priors2: this runs the M-H draw for standard deviations of subject
                   level means of mass and width;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Priors *priors; the current values of the prior parameters;
               double v1; the proposal variance for overall st-dev of mean mass;
               double v2; the proposal variance for overall st-dev of mean width;
               unsigned long *seed; seed values needed for the randon number generator;
               Hyper_priors *hyper; parameters of the prior distributions;
    RETURNS: None; all updates are made internally; this does update the global
             variables regarding acceptance rate: afepmv, afepwv, nfepmv, and nfepwv  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j,t: generic counters
 *tmp1: current acceptance counters are saved to this vector so we can increment
    inside a loop
 *new_var: proposed values of st. dev. of mass and width
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
    proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subject: current list of subjects and their qualities

SUBROUTINES USED
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution  */
/**************************************************************************/

void draw_fe_priors2(Subject_type *sublist,Priors *priors,double v1,double v2,unsigned long *seed,Hyper_priors *hyper)
{
  int j,t;
  double *tmp1,*new_var,prop_new,prop_old,psum,pcomp,prior_old,prior_new;
  double prior_ratio,prop_ratio,temp,alpha;
  Subject_type *subject;
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);

  /*Allocate Memory*/
  new_var = (double *)calloc(2,sizeof(double));
  tmp1 = (double *)calloc(2,sizeof(double));

  /*Add 1 to the counters for acceptance rates of sigma_ma and sigma_mw*/
  nfepmv++;
  nfepwv++;

  /*Assign current acceptance counts to temporary vector*/
  tmp1[0] = afepmv;
  tmp1[1] = afepwv;

  /*Draw proposed values for sigma_ma and sigma_mw*/
  new_var[0]=rnorm(priors->fe_precision[0],v1,seed);
  new_var[1]=rnorm(priors->fe_precision[1],v2,seed);

  /*Accept or Reject sigma_ma, then accept or reject for sigma_mw*/
  for (t=0; t<2; t++){
    /*t and j are both used to test what happens if we reverse the
    order of sma and smw*/
      j=1-t;
  /*Calculate ratio of Cauchy priors*/
/*      prior_new = new_var[j]/hyper->prec[j];
      prior_old = priors->fe_precision[j]/hyper->prec[j];
      prior_new *= prior_new;
      prior_old *= prior_old;
      prior_new += 1.0;
      prior_old += 1.0;
      prior_ratio = prior_old/prior_new;
*/
      /*Compute the sum included in the "likelihood"*/
      /*Also compute the old value divided by the new value raised to num_node*/
      /*This sum is the same assuming the current value of sigma_j in both*/
      /*the numerator and denominator of rho*/
      psum = 0;
      subject = sublist->succ;
      prior_ratio = 1;
      while (subject != NULL){
        prior_ratio *= priors->fe_precision[j]/new_var[j];
        pcomp = subject->theta[j] - priors->fe_mean[j];
        pcomp *= pcomp;
        psum += pcomp;
        subject = subject->succ;
        }

      /*Log prior ratio value*/
      prior_ratio = log(prior_ratio);

      /*Complete the "likelihood" portion of rho*/
      prop_old = 1.0/priors->fe_precision[j];
      prop_old /= priors->fe_precision[j];
      prop_new = 1.0/new_var[j];
      prop_new /= new_var[j];
      prop_ratio = psum*0.5*(prop_old - prop_new);

      /*We only can accept the proposed value if it is positive*/
      if(new_var[j]>0 && new_var[j]<hyper->prec[j]){

          /*Compute log rho, and set alpha equal to min(log rho,0)*/
          alpha = (0<(temp=(prior_ratio+prop_ratio)))?0:temp;

          /*If log(U) < log rho, accept the proposed value*/
          /*Increase acceptance count by 1*/
          if(log(kiss(seed))<alpha){
              tmp1[j]++;
              priors->fe_precision[j]=new_var[j];
          }
      }
   }

  /*Set acceptance count equal to temp vector components*/
  afepmv = tmp1[0];
  afepwv = tmp1[1];

  free(new_var);
  free(tmp1);

}

/*********************************************************************/
               /*START OF draw_bh_mean SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh_mean: this runs the Gibbs sampler draw for mean baseline and half life;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Priors *priors; the current values of the prior parameters;
               Common_parms *parms; the current values of the common parameters;
               unsigned long *seed; seed values needed for the randon number generator;
               Hyper_priors *hyper; parameters of the prior distributions;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 fesum: sum of subject specific mean baselines and half-lifes;
        part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics

SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution  */
/**************************************************************************/

void draw_bh_mean(Subject_type *sublist,Priors *priors,Common_parms *parms,unsigned long *seed,Hyper_priors *hyper)
{
/* declare variables */
   int j;
   double fesum,gmean,gvar;
   Subject_type *subject;

/*declare functions */
   double rnorm(double,double,unsigned long *);
   double kiss(unsigned long *);

    for (j=0;j<2;j++){

        fesum=0;
        subject=sublist->succ;

        /*This while loop goes through the pulses to get the sum involved*/
        while(subject != NULL){
            fesum+=subject->basehalf[j];
            subject=subject->succ;
            } /*end of loop through pulses*/

        gmean = (hyper->meanvarbh[j]*fesum) + (priors->varbh[j]*priors->varbh[j]*hyper->meanmeanbh[j]);
        gmean /= (parms->numsub*hyper->meanvarbh[j])+(priors->varbh[j]*priors->varbh[j]);

        gvar = hyper->meanvarbh[j]*priors->varbh[j]*priors->varbh[j];
        gvar /= (parms->numsub*hyper->meanvarbh[j])+(priors->varbh[j]*priors->varbh[j]);

        priors->meanbh[j] = rnorm(gmean,sqrt(gvar),seed);
    }

}

/*********************************************************************/
               /*START OF draw_bh_var SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh_var: this runs the M-H draw for standard deviations of subject
                   level baselines and half-lifes;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Priors *priors; the current values of the prior parameters;
               double v1; the proposal variance for overall st-dev of baseline;
               double v2; the proposal variance for overall st-dev of half-life;
               unsigned long *seed; seed values needed for the randon number generator;
               Hyper_priors *hyper; parameters of the prior distributions;
    RETURNS: None; all updates are made internally; this does update the global
             variables regarding acceptance rate: afebv, afehv, nfebv, and nfehv  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 *tmp1: current acceptance counters are saved to this vector so we can increment
    inside a loop
 *new_var: proposed values of st. dev. of baseline and half-life
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
    proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subject: current list of subjects and their qualities

SUBROUTINES USED
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution  */
/**************************************************************************/

void draw_bh_var(Subject_type *sublist,Priors *priors,double v1,double v2,unsigned long *seed,Hyper_priors *hyper)
{
  int j;
  double *tmp1,*new_var,prop_new,prop_old,psum,pcomp,prior_old,prior_new;
  double prior_ratio,prop_ratio,temp,alpha;
  Subject_type *subject;
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);

  /*Allocate Memory*/
  new_var = (double *)calloc(2,sizeof(double));
  tmp1 = (double *)calloc(2,sizeof(double));

  /*Add 1 to the counters for acceptance rates of sigma_b and sigma_h*/
  nfebv++;
  nfehv++;

  /*Assign current acceptance counts to temporary vector*/
  tmp1[0] = afebv;
  tmp1[1] = afehv;

  /*Draw proposed values for sigma_b and sigma_h*/
  new_var[0]=rnorm(priors->varbh[0],v1,seed);
  new_var[1]=rnorm(priors->varbh[1],v2,seed);

  /*Accept or Reject sigma_b, then accept or reject for sigma_h*/
  for (j=0; j<2; j++){
        /*Calculate ratio of Cauchy priors*/
      /*prior_new = new_var[j]/hyper->varbh[j];
      prior_old = priors->varbh[j]/hyper->varbh[j];
      prior_new *= prior_new;
      prior_old *= prior_old;
      prior_new += 1.0;
      prior_old += 1.0;
      prior_ratio = prior_old/prior_new;
*/
      /*Compute the sum included in the "likelihood"*/
      /*Also compute the old value divided by the new value raised to num_node*/
      /*This sum is the same assuming the current value of sigma_j in both*/
      /*the numerator and denominator of rho*/
      psum = 0;
      subject = sublist->succ;
      prior_ratio = 1;
      while (subject != NULL){
        prior_ratio *= priors->varbh[j]/new_var[j];
        pcomp = subject->basehalf[j] - priors->meanbh[j];
        pcomp *= pcomp;
        psum += pcomp;
        subject = subject->succ;
        }

      /*Log prior ratio value*/
      prior_ratio = log(prior_ratio);

      /*Complete the "likelihood" portion of rho*/
      prop_old = 1.0/priors->varbh[j];
      prop_old /= priors->varbh[j];
      prop_new = 1.0/new_var[j];
      prop_new /= new_var[j];
      prop_ratio = psum*0.5*(prop_old - prop_new);

      /*We only can accept the proposed value if it is positive*/
      /*With a uniform(0,2) prior, skip this step if we have a sd less than 2*/
      if(new_var[j]>0 && new_var[j]<hyper->varbh[j]){

          /*Compute log rho, and set alpha equal to min(log rho,0)*/
          /*With a uniform prior, get rid of prior_ratio*/
          alpha = (0<(temp=(prior_ratio+prop_ratio)))?0:temp;

          /*If log(U) < log rho, accept the proposed value*/
          /*Increase acceptance count by 1*/
          if(log(kiss(seed))<alpha){
              tmp1[j]++;
              priors->varbh[j]=new_var[j];
          }
      }
  }

  /*Set acceptance count equal to temp vector components*/
  afebv = tmp1[0];
  afehv = tmp1[1];

  free(new_var);
  free(tmp1);

}

/*********************************************************************/
               /*START OF draw_fixed_effects SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fixed_effects: this runs the Gibbs sampler draw for subject-specific
                      mean masses and widths;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Priors *priors; the current values of the prior parameters;
               Common_parms *parms; the current values of the common parameters;
               unsigned long *seed; seed values needed for the randon number generator;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 numnode: counter of pulses
 resum: sum of pulse amplitudes or widths for a given subject;
        part of evaluation of gmean and gvar
 gmean: mean of the distribution used in Gibbs sampler
 gvar: variance of distribution used in Gibbs sampler
 *subject: current list of subjects and their characteristics
 *subnode: for a given subject, list of pulses and their characteristics

SUBROUTINES USED
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution  */
/**************************************************************************/

void draw_fixed_effects(Subject_type *sublist,Priors *priors,Common_parms *parms,unsigned long *seed)
{
  int j,numnode;
  double resum, gmean, gvar;
  Subject_type *subject;
  Node_type *subnode;
  void mean_contribution(Node_type *,double **,Common_parms *,int,double);
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);

  /*Go to start of subject list*/
  subject = sublist->succ;

  /*Go through each existing subject*/
  while (subject != NULL) {

      for(j=0;j<2;j++){

      /*Go to start of pulse list for this subject*/
      subnode = subject->list->succ;
      resum = 0.0;
      numnode=0;

      /*Go through each existing pulse and add up log amplitude/width
       Also count pulses*/
      while(subnode != NULL){
          resum+=log(subnode->theta[j]);
          numnode++;
          subnode = subnode->succ;
        }

      /*Calculate mean of distribution used for Gibbs sampler*/
      gmean = priors->fe_precision[j]*priors->fe_precision[j];
      gmean *= resum;
      gmean += parms->re_precision[j]*parms->re_precision[j]*priors->fe_mean[j];
      gmean /= (numnode*priors->fe_precision[j]*priors->fe_precision[j])+(parms->re_precision[j]*parms->re_precision[j]);

      /*Calculate variance of distribution used for Gibbs sampler*/
      gvar = (priors->fe_precision[j]*priors->fe_precision[j])*(parms->re_precision[j]*parms->re_precision[j]);
      gvar /= (numnode*priors->fe_precision[j]*priors->fe_precision[j])+(parms->re_precision[j]*parms->re_precision[j]);

      /*Use Gibbs sampler to draw new value*/
      subject->theta[j]=rnorm(gmean,sqrt(gvar),seed);

        } /*end of loop through mass & width */

    /*Advance to next subject*/
    subject = subject->succ;

  } /*end of loop through subjects*/

}

/*********************************************************************/
               /*START OF draw_fe_precision SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fe_precision: this runs the M-H draw for standard deviations of pulse
                     masses and widths;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Priors *priors; the current values of the prior parameters;
               Common_parms *parms; the current values of the common parameters;
               double v1; the proposal variance for overall st-dev of pulse mass;
               double v2; the proposal variance for overall st-dev of pulse width;
               unsigned long *seed; seed values needed for the randon number generator;
    RETURNS: None; all updates are made internally; this does update the global
             variables regarding acceptance rate: afemv, afewv, nfemv, and nfewv  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 *tmp1: current acceptance counters are saved to this vector so we can increment
    inside a loop
 *new_var: proposed values of st. dev. of pulse mass and width
 prop_new: part of evaluation of prop_ratio portion of log(rho)
 prop_old: same as prop_new but with current (not proposed) value
 psum: part of evaluation of prop_ratio portion of log(rho)
 pcomp: part of evaluation of prop_ratio portion of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 prop_ratio: comparison between random effects distributions under current and
    proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *subnode: for a given subject, list of pulses and their characteristics
 *subject: current list of subjects and their qualities

SUBROUTINES USED
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution  */
/**************************************************************************/

void draw_fe_precision(Subject_type *sublist,Priors *priors,Common_parms *parms,double v1,double v2,unsigned long *seed)
{
  int j;
  double *tmp1,*new_var,prop_new,prop_old,psum,pcomp,prior_old,prior_new;
  double prior_ratio,prop_ratio,temp,alpha;
  Node_type *subnode;
  Subject_type *subject;
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);

  /*Allocate Memory*/
  new_var = (double *)calloc(2,sizeof(double));
  tmp1 = (double *)calloc(2,sizeof(double));

  /*Add 1 to the counters for acceptance rates of sigma_a and sigma_w*/
  nfemv++;
  nfewv++;

  /*Assign current acceptance counts to temporary vector*/
  tmp1[0] = afemv;
  tmp1[1] = afewv;

  /*Draw proposed values for sigma_a and sigma_w*/
  new_var[0]=rnorm(parms->re_precision[0],v1,seed);
  new_var[1]=rnorm(parms->re_precision[1],v2,seed);

  /*Accept or Reject sigma_a, then accept or reject for sigma_w*/
  for (j=0; j<2; j++){

      /*Calculate ratio of Cauchy priors*/
/*      prior_new = new_var[j]/priors->re_var[j];
      prior_old = parms->re_precision[j]/priors->re_var[j];
      prior_new *= prior_new;
      prior_old *= prior_old;
      prior_new += 1.0;
      prior_old += 1.0;
      prior_ratio = prior_old/prior_new;
*/
      /*Compute the sum included in the "likelihood"*/
      /*Also compute the old value divided by the new value raised to num_node*/
      /*This sum is the same assuming the current value of sigma_j in both*/
      /*the numerator and denominator of rho*/
      psum = 0;
      subject = sublist->succ;
      prior_ratio = 1;
      while (subject != NULL){
          subnode = subject->list->succ;
          while (subnode !=NULL){
              prior_ratio *= parms->re_precision[j]/new_var[j];
              pcomp = log(subnode->theta[j]) - subject->theta[j];
              pcomp *= pcomp;
              psum += pcomp;
              subnode = subnode->succ;
              }
          subject=subject->succ;
          }

      /*Log prior ratio value*/
      prior_ratio = log(prior_ratio);

      /*Complete the "likelihood" portion of rho*/
      prop_old = 1.0/parms->re_precision[j];
      prop_old /= parms->re_precision[j];
      prop_new = 1.0/new_var[j];
      prop_new /= new_var[j];
      prop_ratio = psum*0.5*(prop_old - prop_new);

      /*We only can accept the proposed value if it is positive*/
      if(new_var[j]>0.1 && new_var[j]<priors->re_var[j]){

          /*Compute log rho, and set alpha equal to min(log rho,0)*/
          alpha = (0<(temp=(prior_ratio+prop_ratio)))?0:temp;
          /*If log(U) < log rho, accept the proposed value*/
          /*Increase acceptance count by 1*/
          if(log(kiss(seed))<alpha){
              tmp1[j]++;
              parms->re_precision[j]=new_var[j];
          }
      }
  } /*end of loop through s_a and s_w */

  /*Set acceptance count equal to temp vector components*/
  afemv = tmp1[0];
  afewv = tmp1[1];

  free(new_var);
  free(tmp1);

}

/*********************************************************************/
                    /*START OF draw_bh SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_bh: this runs the M-H draw for baseline and halflife (drawn together);
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Common_parms *parms; the current values of the common parameters;
               Priors *priors; the current values of the prior parameters;
               double **ts; this is the matrix of observed data (a column of
                  times and S columns of log(concentration);
               int N; the number of observations in each column of **ts;
               unsigned long *seed; seed values needed for the randon number generator;
               double **var; the proposal variance-covariance matrix for
                  baseline and halflife
    RETURNS: None; all updates are made internally; this does update the global
             variables regarding acceptance rate: adelta and ndelta  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 num_node: number of pulses
 *subnode: for a given subject, current list of pulses and their qualities
 *new_node: counter for number of pulses
 *subject: current list of subjects and their characteristics
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 *pmd: vector of proposed baseline and halflife
 like_ratio: difference between current and proposed log-likelihoods
 priorb_old: part of evaluation of prior_ratio portion of log(rho)
 priorb_new: same as priorb_old, but with proposed baseline
 priorh_old: part of evaluation of prior_ratio portion of log(rho)
 priorh_new: same as priorb_old, but with proposed halflife
 prior_ratio: comparison between prior distributions under current and proposed values
 currentmd[2]: current baseline and halflife saved to this variable name so we can
    evaluate likelihood under proposed baseline and halflife
 **currentmc: current mean_contrib for each pulse is saved to this matrix
    so we can evaluate likelihood under proposed pulses' mean_contrib
 temp: value used in evaluation of log(rho)
 current_decay: current decay rate is saved to this variable name so we can
    evaluate likelihood under proposed decay rate
 current_like: current likelihood for a given subject; saved so we can compare

SUBROUTINES USED
 rmvnorm: found in randgen.c; draws from the multivariate normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
    pulse with specific parameters
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
    given subject under current parameters*/
/**************************************************************************/

void draw_bh(Subject_type *sublist,Common_parms *parms,Priors *priors,double **ts,
           int N,unsigned long *seed,double **var)
{
  int i,j,k,num_node;
  Node_type *subnode, *new_node;
  Subject_type *subject;
  double alpha,plikelihood,*pmd,like_ratio,priorb_old,priorb_new,priorh_old;
  double priorh_new,prior_ratio,currentmd[2],**currentmc,temp,current_decay;
  double current_like;
  int rmvnorm(double *,double **,int,double *,unsigned long *,int);
  double kiss(unsigned long *);
  void mean_contribution(Node_type *,double **,Common_parms *,int,double);
  double likelihood2(Node_type *,double **,Common_parms *,int,Node_type *,double);

  /*Go to start of list of subjects*/
  subject = sublist->succ;

  /*this tells likelihood2 which subject's likelihood to calculate*/
  parms->subindex = 0;

  /*Loop thru all subjects*/
  while (subject != NULL){

   /*Calculate the current likelihood for this subject*/
   current_like = likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);

   /*Count number of pulses for this subject*/
   num_node= 0;
    new_node = subject->list->succ;
    while (new_node != NULL) {
        new_node = new_node->succ;
        num_node++;
    }

    /*Allocate memory*/
    pmd = (double *)calloc(2,sizeof(double));
    currentmc = (double **)calloc(num_node,sizeof(double *));
    for (i=0;i<num_node;i++)
        currentmc[i] = (double *)calloc(N,sizeof(double));

   /*Increase denominator of acceptance rate for b and hl*/
   ndelta++;

  /*Draw proposal values for b and hl*/
  rmvnorm(pmd,var,2,subject->basehalf,seed,1);

  /*Only proceed if we draw reasonable values*/
  if (pmd[0] > 0 && pmd[1] > 3) {

    /*Compute ratio of prior densities*/
    priorb_old = subject->basehalf[0] - priors->meanbh[0];
    priorb_old *= 0.5*priorb_old/(priors->varbh[0]*priors->varbh[0]);
    priorb_new = pmd[0]-priors->meanbh[0];
    priorb_new *= 0.5*priorb_new/(priors->varbh[0]*priors->varbh[0]);
    priorh_old = subject->basehalf[1] - priors->meanbh[1];
    priorh_old *= 0.5*priorh_old/(priors->varbh[1]*priors->varbh[1]);
    priorh_new = pmd[1]-priors->meanbh[1];
    priorh_new *= 0.5*priorh_new/(priors->varbh[1]*priors->varbh[1]);

    prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

      /*Save current values of b and hl and set new current values equal to
      proposal values */
    for (k=0;k<2;k++) {
      currentmd[k] = subject->basehalf[k];
      subject->basehalf[k] = pmd[k];
    }

    /*Save current decay rate; calculate new current decay rate based on
      proposal value of hl */
    current_decay = subject->decay;
    subject->decay = log(2)/subject->basehalf[1];

    /*Save current mean_contrib for each pulse; calculate new current
      mean_contrib for each pulse based on proposed values of b and hl*/
    i = 0;
    subnode = subject->list->succ;
    while (subnode != NULL) {
      for (j=0;j<N;j++)
	currentmc[i][j] = subnode->mean_contrib[j];
        mean_contribution(subnode,ts,parms,N,subject->basehalf[1]);
      i++;
      subnode = subnode->succ;
    }

    /*Calculate proposed likelihood and then calculate likelihood ratio*/
    plikelihood = likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);
    like_ratio = plikelihood - current_like;

    /*Calculate log rho; set alpha equal to min(0,log rho) */
    alpha = (0 < (temp = (prior_ratio+like_ratio))) ? 0:temp;

    /*If log U < log rho, increase acceptance rate by 1  */
    if (log(kiss(seed)) < alpha) {
      adelta++;
    }

    /*Otherwise, we need to revert values back to their current state */
    else {

      /*Set b and hl back equal to current values*/
      subject->basehalf[0] = currentmd[0];
      subject->basehalf[1] = currentmd[1];

      /*Set mean_contrib back equal to current values for each pulse*/
      i = 0;
      subnode = subject->list->succ;
      while (subnode != NULL) {
	for (j=0;j<N;j++) {
	  subnode->mean_contrib[j] = currentmc[i][j];
        }
	i++;
	subnode = subnode->succ;
      }

      /*Set decay rate back equal to current value*/
      subject->decay = current_decay;

    } /*End of if else statement*/

  }

  /*Go to next subject*/
  subject = subject->succ;
  parms->subindex++;

  /*Free memory */
  for (i=0;i<num_node;i++)
    free(currentmc[i]);
  free(currentmc);
  free(pmd);

  } /*End of loop through subjects*/

}

/*********************************************************************/
                    /*START OF draw_times SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_times: this runs the M-H draws for each individual pulse location;
    ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
               Common_parms *parms; the current values of the common parameters;
               double **ts; this is the matrix of observed data (a column of
                  times and S columns of log(concentration);
               int N; the number of observations in **ts;
               unsigned long *seed; seed values needed for the randon number generator;
               double v; the proposal variance for pulse location
    RETURNS: None; all updates are made internally; this does update the global
             variables regarding acceptance rate: atime and ntime  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,k: generic counters
 *subject: current list of subjects and their qualities
 *subnode: for a given subject, current list of pulses and thier qualities
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 ptime: proposed pulse location
 like_ratio: difference between current and proposed log-likelihoods
 prior_ratio: comparison between prior distributions under current and proposed values
 tmp: value used in evaluation of log(rho)
 current_time: current pulse location saved to this variable name so we can
    evaluate likelihood under proposed pulse location
 *tmp_mean_contrib: current mean_contrib for a pulse is saved to this vector
    so we can evaluate likelihood under proposed pulse's mean_contrib
 current_like: current likelihood for a given subject; used for comparison
 time_diff1: part of evaluation of prior_ratio portion of log(rho)
 time_diff1_new: same as time_diff1, but with proposed pulse location
 time_diff2: part of evaluation of prior_ratio portion of log(rho)
 time_diff2_new: same as time_diff1, but with proposed pulse location

 SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
    pulse with specific parameters
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
    given subject under current parameters*/
/**************************************************************************/

void draw_times(Subject_type *sublist,Common_parms *parms,double **ts,
             int N,unsigned long *seed,double v) {
 int i,k;
 Subject_type *subject;
 Node_type *subnode;
 double alpha,plikelihood,ptime,like_ratio,prior_ratio,tmp,current_time,*tmp_mean_contrib;
 double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new;
 double rnorm(double,double,unsigned long *);
 double kiss(unsigned long *);
 void mean_contribution(Node_type *,double **,Common_parms *,int,double);
 double likelihood2(Node_type *,double **,Common_parms *,int,Node_type *,double);

 /*Allocate Memory*/
 tmp_mean_contrib = (double *)calloc(N,sizeof(double));

 /*Go to start of subject list*/
 subject = sublist->succ;

 /*Index used for likelihood2*/
 parms->subindex=0;

 /*Loop thru subjects*/
 while (subject != NULL){
    k=0;
    /*Go to first pulse for this subject*/
    subnode = subject->list->succ;

    /*Loop thru pulses for this subject*/
    while (subnode!=NULL) {

        /*Increase denominator of acceptance rate for time*/
        ntime++;

        /*Find current likelihood for this subject*/
        current_like = likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);

        /*Compute proposal time*/
        ptime = rnorm(subnode->time,v,seed);

        /*Only proceed if our proposed time is reasonable */
        if (ptime <=fitend && ptime > fitstart) {

            /*If we're at the first pulse, use these as the first numerator/
            denominator of the prior ratio*/
            if (subnode->pred != NULL) {
                time_diff1_new = ptime - subnode->pred->time;
                time_diff1 = subnode->time - subnode->pred->time;
            }

            /*Otherwise, use these as the first numerator/denominator of
            the prior ratio*/
            else {
                time_diff1_new = ptime - fitstart;
                time_diff1 = subnode->time - fitstart;
            }

            /*If we're at the last pulse, use these as the second numerator/
            denominator of the prior ratio*/
            if (subnode->succ != NULL) {
                time_diff2_new = subnode->succ->time - ptime;
                time_diff2 = subnode->succ->time - subnode->time;
            }

            /*Otherwise, use these as the second numerator/denominator*/
            else {
                time_diff2_new = fitend - ptime;
                time_diff2_new = fitend - subnode->time;
            }

            /*Combine it all for the prior ratio*/
            prior_ratio = (mmm-1)*log((time_diff1_new*time_diff2_new)/(time_diff1 * time_diff2));

            /*Save current time*/
            current_time = subnode->time;

            /*Set pulse's time to the proposal value*/
            subnode->time = ptime;

            /*Save mean_contrib of this pulse*/
            for (i=0;i<N;i++)
                tmp_mean_contrib[i] = subnode->mean_contrib[i];

            /*Recalculate mean_contrib of this pulse assuming proposed value*/
            mean_contribution(subnode,ts,parms,N,subject->basehalf[1]);

            /*Calculate the likelihood for this subject under the proposed value*/
            plikelihood = likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);

            /*Calculate the likelihood ratio*/
            like_ratio = plikelihood - current_like;

            /*Calculate log rho; set alpha equal to min(0,log rho)*/
            alpha = (0 < (tmp = (prior_ratio + like_ratio))) ? 0:tmp;

            /*If log U < log rho, we accept proposed value. Increase acceptance
            count by one*/
            if (log(kiss(seed)) < alpha) {
                atime++;
            }

            /*Otherwise, we reject, and we have to revert back to current values*/
            else {

                /*Set pulse's time back to current time*/
                subnode->time = current_time;

                /*Set mean_contrib of this pulse back to current value*/
                for (i=0;i<N;i++)
                    subnode->mean_contrib[i] = tmp_mean_contrib[i];
            }
        } /*End of if time is feasible statement*/

    /*Advance one pulse*/
    subnode = subnode->succ;
    k++;
    }/*End of loop through pulses*/

    /*Go to next subject*/
   subject=subject->succ;
   parms->subindex++;

    }/*End of loop through subjects*/

 /*Free Memory*/
 free(tmp_mean_contrib);

}

/*********************************************************************/
              /*START OF draw_random_effects SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_random_effects: this runs the M-H draw for individual pulse masses and
                      widths;
    ARGUMENTS: double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               Subject_type *list; this is the current list of subjects that exist;
               Common_parms *parms; the current values of the common parameters;
               int N; the number of observations in **ts;
               double v1; the proposal variance for individual pulse masses;
               double v2; the proposal variance for individual pulse widths;
               unsigned long *seed; seed values needed for the randon number generator;
    RETURNS: None; all updates are made internally; this does update the global
             variables regarding acceptance rate: arem, arew, nrem, nrew  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 temp: value used in evaluation of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 old_val: current value of the random effect is saved to this variable so we
    can evaluate likelihood under the proposed value
 *tmp1: current acceptance counters are saved to this vector so we can increment
    inside a loop
 *pRE: proposed values of individual mass and width
 prior_ratio: comparison between prior distributions under current and proposed values
 like_ratio: difference between current and proposed log-likelihoods
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 plikelihood: likelihood under proposed value
 *old_contrib: current mean_contrib for a pulse is saved to this vector
    so we can evaluate likelihood under proposed pulse's mean_contrib
 current_like: current likelihood for a given subject; used for comparison
 *subnode: for a given subject, current list of pulses and their qualities
 *subject: current list of subjects and their qualities

SUBROUTINES USED
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
    pulse with specific parameters
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution
 likelihood2: found in birthdeath.c; calculates and returns likelihood for a
    given subject under current parameters  */
/**************************************************************************/

void draw_random_effects(double **ts,Subject_type *sublist,Common_parms *parms,int N,double v1, double v2,unsigned long *seed)
{
  int i,j,k;
  double temp,prior_old,prior_new,old_val,*tmp1,*pRE,prior_ratio,like_ratio;
  double alpha,plikelihood,*old_contrib,current_like;
  Node_type *subnode;
  Subject_type *subject;
  void mean_contribution(Node_type *,double **,Common_parms *,int,double);
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);
  double likelihood2(Node_type *,double **,Common_parms *,int,Node_type *,double);

  /*Set acceptance counts equal to temporary vector*/
  tmp1 = (double *)calloc(2,sizeof(double));
  tmp1[0] = arem;
  tmp1[1] = arew;

  /*Go to first subject*/
  subject = sublist->succ;

  /*Index of subjects used for likelihood2*/
  parms->subindex=0;

  /*Loop thru subjects*/
  while (subject !=NULL){

      /*Go to start of pulses for this subject*/
      subnode = subject->list->succ;
      k=0;

      /*Go through each existing pulse for this subject*/
      while (subnode != NULL) {

        /*Allocate Memory*/
        old_contrib = (double *)calloc(N,sizeof(double));
        pRE = (double *)calloc(2,sizeof(double));

        /*Increase the denominators of the acceptance rates*/
          nrem++;
          nrew++;

          /*Draw proposed values of current pulse's mass and width*/
          pRE[0]=rnorm(log(subnode->theta[0]),v1,seed);
          pRE[1]=rnorm(log(subnode->theta[1]),v2,seed);

          /*Determine if we accept or reject proposed pulse mass then determine*/
          /*if we accept or reject proposed pulse width*/
          for(j=0;j<2;j++){
            current_like = likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);

            /*Compute the log of the ratio of the priors*/
            prior_old = log(subnode->theta[j]) - subject->theta[j];
            prior_old *= 0.5*prior_old;
            prior_new = pRE[j] - subject->theta[j];
            prior_new *= 0.5*prior_new;
            prior_ratio = prior_old - prior_new;
            prior_ratio /= parms -> re_precision[j];
            prior_ratio /= parms -> re_precision[j];

            /*Save the current value of mass/width*/
            old_val = subnode -> theta[j];

            /*Set the pulse's mass/width equal to the proposed value*/
            subnode -> theta[j] = exp(pRE[j]);

            /*Save the mean_contrib for that pulse*/
            for(i=0;i<N;i++)
                old_contrib[i]=subnode->mean_contrib[i];

            /*Recalculate that pulse's mean_contrib assuming proposed mass/width */
            mean_contribution(subnode,ts,parms,N,subject->basehalf[1]);

            /*Calculate likelihood assuming proposed mass/width */
            plikelihood = likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);

            like_ratio = plikelihood - current_like;

            /*Compute the log of the ratio between the two likelihoods*/

            /*Calculate log rho; set alpha equal to min(0,log rho) */
            alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0:temp;

            /*If log U < log rho, accept the proposed value, increase acceptance*/
            /*counter */
            if (log(kiss(seed)) < alpha) {
                tmp1[j]++;
            }

            /*Otherwise, reject the proposed value, set pulse's mass/width back to */
            /*saved current value, and set pulse's mean_contrib back to saved value */
            else {
                subnode->theta[j] = old_val;
                for (i=0;i<N;i++)
                    subnode->mean_contrib[i] = old_contrib[i];
                }
         } /*end of loop through mass & width */

    /*Advance to next pulse*/
    subnode = subnode->succ;
    k++;
     /*free memory*/
    free(pRE);
    free(old_contrib);

      } /*end of loop through pulses*/

    /*advance to next subject*/
    subject=subject->succ;
    parms->subindex++;

} /*end of loop through subjects*/
  /*Set counter equal to temporary vector values*/
  arem = tmp1[0];
  arew = tmp1[1];
  free(tmp1);

}

/*********************************************************************/
              /*START OF error_squared SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*error_squared: This subroutine calculates the sum of squared error. It uses
                the list of pulses and the common parameters to evaluate the
                likelihood and calculate the sum of the squared error between
                observed concentration and expected under the current parameters
    ARGUMENTS: double **ts; this is the matrix of observed data (a column of
                  times and S columns of log(concentration)
               Subject_type *sublist; this is the current list of subjects that exist
               Common_parms *parms; the current values of the common parameters
               int N; the number of observations in **ts;
    RETURNS: double ssq; this value is needed for one of the parameters
             of the posterior distribution of model error variance  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j: generic counters
 *mean: mean concentration under current parameters and pulses
 ssq: sum of squared differences between observed log(concentration) and expected
    values
 *subject: current list of subjects and their qualities

SUBROUTINES USED
 mean_concentration: found in birthdeath.c; sums across mean_contrib vectors for
    each pulse; returns this vector  */
 /**************************************************************************/

double error_squared(double **ts,Subject_type *sublist,Common_parms *parms,int N) {
 int i,j;
 double *mean,ssq;
 Subject_type *subject;
 double *mean_concentration(Node_type *,Common_parms *,int,Node_type *,double);

 subject=sublist->succ;
 ssq = 0;
 for(j=0;j<parms->numsub;j++){

    mean = mean_concentration(subject->list,parms,N,subject->list,subject->basehalf[0]);
    for (i=0;i<N;i++)
        ssq += (ts[i][j+1] - mean[i])*(ts[i][j+1] - mean[i]);

    subject=subject->succ;
    free(mean);
 }
 
 return ssq;
}

/*********************************************************************/
              /*START OF adjust_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust_acceptance: this adjusts the proposal variances based on the inputted
                    acceptance rate and proposal variance. If the acceptance rate is
                    too high or too low, proposal variance will be adjusted.
    ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
                  acceptance counter divided by the attempt counter
               double *X; the current proposal variance
    RETURNS: None; update to the proposal variance is made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new proposal variance based on inputs

SUBROUTINES USED
 None  */
 /**************************************************************************/

void adjust_acceptance(double x,double *X)
{
  double y;
  
  y = 1. + 1000.*(x-.35)*(x-.35)*(x-.35);
  if (y < .9)
    y = .9;
  if (y > 1.1)
    y = 1.1;

  *X *= y;
}

/*********************************************************************/
              /*START OF adjust2_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust2_acceptance: this adjusts the proposal variance-covriance matrix for
                     baseline and halflife based on the inputted acceptance rate
                     and proposal matrix. If the acceptance rate is too high or
                     too low, proposal variance will be adjusted.
    ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
                  acceptance counter divided by the attempt counter
               double **X; the current proposal variance-covariance matrix
               double corr; the correlation between the two proposal variances
    RETURNS: None; update to the proposal variance is made internally*/
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new diagonal elements of proposal variance-covariance matrix based on inputs

SUBROUTINES USED
 None  */
 /**************************************************************************/

void adjust2_acceptance(double x,double **X,double corr)
{
  double y;
  
  y = 1. + 1000.*(x-.25)*(x-.25)*(x-.25);
  if (y < .90) {
    y = .90;
       X[0][0] *= y;
    X[1][1] *= y;
    X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
    }
  if (y > 1.1) {
    y = 1.1;
       X[0][0] *= y;
    X[1][1] *= y;
    X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
    }
}

