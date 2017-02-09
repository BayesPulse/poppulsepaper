/*******************************************************************/
/*************************linklistv2.c *****************************/
/*******************************************************************/

#include "deconvolution_main.h"

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
 arew: counter that increments when an individual pulse width estimate is accetped in the MH draw
 nrew: counter that increments for each loop of the individual pulse width MH draw
 afem: counter that increments when an overall mean pulse mass estimate is accetped in the MH draw
 nfem: counter that increments for each loop of the overall mean pulse mass MH draw
 afew: counter that increments when an overall mean pulse width estimate is accetped in the MH draw
 nfew: counter that increments for each loop of the overall mean pulse width MH draw
 arevm: counter that increments when an overall pulse mass standard deviation estimate is accetped in the MH draw
 nrevm: counter that increments for each loop of the overall pulse mass standard deviation MH draw
 arevw: counter that increments when an overall pulse width standard deviation estimate is accetped in the MH draw
 nrevw: counter that increments for each loop of the overall pulse width standard deviation MH draw
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
long arevm = 0;
long nrevm = 0;
long arevw = 0;
long nrevw = 0;
extern double fitstart;
extern double fitend;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 mcmc 
 mh_time 
 mh_mu_delta
 draw_fixed_effects
 draw_re_precision
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
void mcmc(Node_type *list,Common_parms *parms,double **ts,long iter,int N,
          Priors *priors,unsigned long *seed,char *file1,char *file2,double **pmd_var)
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               Common_parms *parms; the current values of the common parameters;
               double **ts; this is the matrix of observed data (a column of 
                  times and a column of log(concentration);
               long iter; the number of iterations to run;
               int N; the number of observations in **ts;
               Priors *priors; the parameters of the prior distributions;
               unsigned long *seed; seed values needed for the randon number generator;
               char *file1; the output file name for common parameters
               char *file2; the output file name for pulse specific parameters
               double **pmd_var; proposal variance-covariance matrix for baseline and halflife
               Hyper_priors *hyper; hyperpriors; not used
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k,l: generic counters
 num_node: counter of num_node for loops
 num_node2: number of pulses
 NN=50: Output every NNth iteration to files
 NNN=5000: Output every NNNth interation to screen or output file
 vrem: proposal variance for individual masses
 vrew: proposal variance for individual widths
 vm: proposal variance for overall mean mass
 vw: proposal variance for overall mean width
 vmv: proposal variance for st. dev. of masses
 vwv: proposal variance for st. dev. of widths
 vt: proposal variance for individual pulse locations
 ssq: sum of squared differences between log(concentration) and expected value
 *likeli: value of the likelihood
 **pmd_vch: cholesky decomposed proposal variance-covariance matrix for
    baseline and halflife
 *ss:
 *new_node: index used to count number of nodes following birthdeath procedure
 *common: name given to inputted file name for common parameters
 *parm: name given to inputted file name for pulse parameters
 *hyper: name given to inputted file name for hyperprior parameters (not used)

 SUBROUTINES USED
 print_list: found in hash.c; traverses the linked list and prints the contents
             of each node
 inverse_gamma: found in randgen.c; draws from the inverse gamma distribution
 error_squared: found in this file; calculates sum of squared error
 mh_mu_delta: found in this file; draws baseline and halflife
 mh_time: found in this file; draws individual pulse locations
 birth_death: found in birthdeath.c; runs birthdeath process
 draw_fixed_effects: found in this file; draws overall mean mass and width
 draw_re_precision: found in this file; draws st. dev. of mass and width
 draw_random_effects: found in this file; draws individual masses and widths
 adjust_acceptance: found in this file; adjusts proposal variance
 adjust_acceptance2: found in this file; adjusts 2-d proposal variance*/
 /***********************************************************************/

void mcmc(Node_type *list,Common_parms *parms,double **ts,long iter,int N,
          Priors *priors,unsigned long *seed,char *file1,char *file2,double propvar[])
{
 int i,j,k,l,num_node,num_node2,NN=50,NNN=5000;
 double vrem,vrew,vmv,vwv,vt,ssq,*likeli,**pmd_var,**pmd_vch;
 char *ss;
 Node_type *new_node;
 FILE *common,*parm,*hyperp;

 void print_list(Node_type *);
 double inverse_gamma(double alpha,double beta,unsigned long *seed);
 double error_squared(double **ts,Node_type *list,Common_parms *parms,int N);
 void mh_mu_delta(Node_type *,Common_parms *,Priors *,double **,
           double *,int,int,unsigned long *,double **);
 void mh_time(Node_type *,Common_parms *,double **,
                  double *,int,unsigned long *,double);
 void birth_death(Node_type *list,double **ts,Common_parms *parms,int N,
                  double *,unsigned long *seed,int);
 void draw_fixed_effects(Node_type *,Priors *,Common_parms *,unsigned long *);
 void draw_re_precision(Node_type *,Priors *,Common_parms *,double,double,unsigned long *);
 void draw_random_effects(double **ts,Node_type *,Common_parms *,int,double *,double,double,unsigned long *);
    int cholesky_decomp(double **, int);
 void adjust_acceptance(double,double *);
 void adjust2_acceptance(double,double **,double);
 
 ss = (char *)calloc(7,sizeof(char));
 common = fopen(file1,"w");
 parm = fopen(file2,"w");
 /* hyperp = fopen(file3,"a");*/

 /*proposal variance for pulse locations*/
 vt = propvar[6];
 
 /*proposal variances for individual masses and individual widths*/
 vrem = propvar[4];
 vrew = propvar[5];
printf("vrew=%lf\n",vrew);

 /*proposal variances for overall stdev of mass & width*/
 vmv = propvar[2];
 vwv = propvar[3];

 pmd_var = (double **)calloc(2,sizeof(double *));
 for (i=0;i<2;i++)
    pmd_var[i] = (double *)calloc(2,sizeof(double));

 pmd_var[0][0] = propvar[0];
 pmd_var[1][1] = propvar[1];
 pmd_var[0][1] = pmd_var[1][0] = -0.90*sqrt(pmd_var[0][0])*sqrt(pmd_var[1][1]);

 free(ss);

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
       
 /*Allocate memory for likelihood */
 likeli = (double *)calloc(1,sizeof(double));
 for (i=0;i<iter;i++) {
     /*        printf("interation %d\n",i);*/
/* run the birth-death algorithm to select pulses */
   
 
   birth_death(list,ts,parms,N,likeli,seed,i);



   num_node2= 0;
   new_node = list->succ;
   while (new_node != NULL) {
     new_node = new_node->succ;
     num_node2++;
   }

/**************************************************/
/* draw the fixed effects*/
            
 /*1*/   draw_fixed_effects(list,priors,parms,seed);

/* draw sigma_1 and sigma_2 --- the random effects precision */

  /*2*/ draw_re_precision(list,priors,parms,vmv,vwv,seed);

/* draw the random effects */
 /*3*/  draw_random_effects(ts,list,parms,N,likeli,vrem,vrew,seed);

/* draw the times */
 /*4*/  mh_time(list,parms,ts,likeli,N,seed,vt);

/* Metropolis-Hastings steps for baseline and halflife */

 /*5*/ mh_mu_delta(list,parms,priors,ts,likeli,N,num_node2,seed,pmd_vch);

/* draw the model variance from the inverse Gamma distribution */
   
   ssq = error_squared(ts,list,parms,N);
   parms->sigma = inverse_gamma(priors->alpha + N/2,priors->beta + 0.5*ssq,seed);
   parms->lsigma = log(parms->sigma);

/*   printf("a = %lf\n",priors->alpha + N/2);
   printf("ssq = %lf\n",ssq);
   fflush(stdout);*/
   /***************************************************************/

   if (!(i%NN)) {
     num_node = 0;
     new_node = list->succ;
     while (new_node != NULL) {
       fprintf(parm,"%d %d %d %lf %lf %lf\n",i/NN,num_node2,num_node,new_node->time,new_node->theta[0],new_node->theta[1]);
       num_node++;
       new_node = new_node->succ;
     }
     fprintf(common,"%d %lf %lf %lf %lf %lf %lf %lf \n",num_node2,parms->md[0],parms->theta[0],parms->theta[1],parms->md[1],parms->sigma,parms->re_precision[0],parms->re_precision[1]);
     fflush(parm);
     fflush(common);
   }
  
   if (!(i%NNN)) {
     printf("\n\n");
     printf("iter = %d likelihood = %lf\n",i,*likeli);
     printf("mu %.2lf A %.2lf s %.2lf d %.4lf  v %.4le\n",
            parms->md[0],parms->theta[0],parms->theta[1],parms->md[1],parms->sigma);
     print_list(list);
     printf("pct rem = %.2lf pct rew = %.2lf pct time = %.2lf pct md = %.2lf revm = %.2lf revw = %.2lf",
            (double)arem/(double)nrem,(double)arew/(double)nrew,(double)atime/(double)ntime,(double)adelta/(double)ndelta,(double)arevm/(double)nrevm,(double)arevw/(double)nrevw);
     fflush(stdout);
   }

    if (!(i%500) && i<25000 && i >0) {
      adjust2_acceptance((double)adelta/(double)ndelta,pmd_var,-0.90);
      adjust_acceptance((double)atime/(double)ntime,&vt);
      adjust_acceptance((double)arem/(double)nrem,&vrem);
      adjust_acceptance((double)arew/(double)nrew,&vrew);
      adjust_acceptance((double)arevm/(double)nrevm,&vmv);
      adjust_acceptance((double)arevw/(double)nrevw,&vwv);
      
      adelta = ndelta = 0;
      atime = ntime = 0;
      arem = nrem = 0;
      arew = nrew = 0;
      arevm = nrevm = 0;
      arevw = nrevw = 0;
	
    
      /***cholesky decompose the three proposal var-covar matrices***/      
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
 /* fclose(hyperp);*/
 free(likeli);
} /*End of MCMC */


/*********************************************************************/
                    /*START OF mh_time SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mh_time: this runs the M-H draws for each individual pulse location;
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               Common_parms *parms; the current values of the common parameters;
               double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               double *like; the current value of the likelihood;
               int N; the number of observations in **ts;
               unsigned long *seed; seed values needed for the randon number generator;
               double v; the proposal variance for pulse location
    RETURNS: None; all updates are made internally
 *********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *node: current list of pulses and their qualities
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
 time_diff1: part of evaluation of prior_ratio portion of log(rho)
 time_diff1_new: same as time_diff1, but with proposed pulse location
 time_diff2: part of evaluation of prior_ratio portion of log(rho)
 time_diff2_new: same as time_diff1, but with proposed pulse location

 SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
    pulse with specific parameters
 likelihood: found in birthdeath.c; calculates and returns likelihood under
    current parameters
 print_list: found in hash.c; traverses the linked list and prints the contents
             of each node*/
/**************************************************************************/

void mh_time(Node_type *list,Common_parms *parms,double **ts,
             double *like,int N,unsigned long *seed,double v) {
 int i;
 Node_type *node;
 double alpha,plikelihood,ptime,like_ratio,prior_ratio,tmp,current_time,*tmp_mean_contrib;
 double time_diff1, time_diff1_new, time_diff2, time_diff2_new;
 double rnorm(double,double,unsigned long *);
 double kiss(unsigned long *);
 void mean_contribution(Node_type *,double **,Common_parms *,int);
 double likelihood(Node_type *,double **,Common_parms *,int,Node_type *);
 void print_list(Node_type *list);

 /*Allocate Memory*/
 tmp_mean_contrib = (double *)calloc(N,sizeof(double));


 node = list->succ;

 while (node!=NULL) {

   /*Increase denominator of acceptance rate for time*/
   ntime++;

   /*Compute proposal time*/
   ptime = rnorm(node->time,v,seed);

   /*Only proceed if our proposed time is reasonable */
   if (ptime <=fitend && ptime > fitstart) {

        /*If we're at the first pulse, use these as the first numerator/
         denominator of the prior ratio*/
        if (node->pred != NULL) {
            time_diff1_new = ptime - node->pred->time;
            time_diff1 = node->time - node->pred->time;
        }

        /*Otherwise, use these as the first numerator/denominator of
          the prior ratio*/
        else {
            time_diff1_new = ptime - fitstart;
            time_diff1 = node->time - fitstart;
        }

        /*If we're at the last pulse, use these as the second numerator/
         denominator of the prior ratio*/
        if (node->succ != NULL) {
            time_diff2_new = node->succ->time - ptime;
            time_diff2 = node->succ->time - node->time;
        }

        /*Otherwise, use these as the second numerator/denominator*/
        else {
            time_diff2_new = fitend - ptime;
            time_diff2_new = fitend - node->time;
        }
            
        /*Combine it all for the prior ratio*/
	prior_ratio = (mmm-1)*log((time_diff1_new*time_diff2_new)/(time_diff1 * time_diff2));

        /*Save current time*/
        current_time = node->time;

        /*Set pulse's time to the proposal value*/
        node->time = ptime;

        /*Save mean_contrib of this pulse*/
        for (i=0;i<N;i++)
          tmp_mean_contrib[i] = node->mean_contrib[i];

        /*Recalculate mean_contrib of this pulse assuming proposed value*/
        mean_contribution(node,ts,parms,N);

        /*Calculate the likelihood under the proposed value*/
        plikelihood = likelihood(list,ts,parms,N,list);

        /*Calculate the likelihood ratio*/
        like_ratio = plikelihood - *like;

        /*Calculate log rho; set alpha equal to min(0,log rho)*/
        alpha = (0 < (tmp = (prior_ratio + like_ratio))) ? 0:tmp;

        /*If log U < log rho, we accept proposed value. Increase acceptance
         count by one and set likelihood equal to likelihood under proposal*/
        if (log(kiss(seed)) < alpha) {
          atime++;
          *like = plikelihood;
        }

        /*Otherwise, we reject, and we have to revert back to current values*/
        else {

            /*Set pulse's time back to current time*/
            node->time = current_time;

            /*Set mean_contrib of this pulse back to current value*/
            for (i=0;i<N;i++)
               node->mean_contrib[i] = tmp_mean_contrib[i];
        }
   } /*End of if time is feasible statement*/

   /*Advance one pulse*/
   node = node->succ;

 }/*End of loop through pulses*/

 /*Free Memory*/
 free(tmp_mean_contrib);

}

/*********************************************************************/
                    /*START OF mh_mu_delta SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mh_mu_delta: this runs the M-H draw for baseline and halflife (drawn together);
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               Common_parms *parms; the current values of the common parameters;
               Priors *priors; the parameters of the prior distributions;
               double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               double *like; the current value of the likelihood;
               int N; the number of observations in **ts;
               int num_node; the current number of pulses;
               unsigned long *seed; seed values needed for the randon number generator;
               double **var; the proposal variance-covariance matrix for
                  baseline and halflife
    RETURNS: None; all updates are made internally
 *********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 *node: current list of pulses and their qualities
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

SUBROUTINES USED
 rmvnorm: found in randgen.c; draws from the multivariate normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
    pulse with specific parameters
 likelihood: found in birthdeath.c; calculates and returns likelihood under
    current parameters*/
/**************************************************************************/

void mh_mu_delta(Node_type *list,Common_parms *parms,Priors *priors,double **ts,
           double *like,int N,int num_node,unsigned long *seed,double **var) 
{
  int i,j,k;
  Node_type *node;
  double alpha,plikelihood,*pmd,like_ratio,priorb_old,priorb_new,priorh_old,priorh_new,prior_ratio,currentmd[2],**currentmc,temp,current_decay;
  int rmvnorm(double *,double **,int,double *,unsigned long *,int);
  double kiss(unsigned long *);
  void mean_contribution(Node_type *,double **,Common_parms *,int);
  double likelihood(Node_type *,double **,Common_parms *,int,Node_type *);

  /*Allocate memory*/
  currentmc = (double **)calloc(num_node,sizeof(double *));
  for (i=0;i<num_node;i++)
    currentmc[i] = (double *)calloc(N,sizeof(double));
     
  pmd = (double *)calloc(2,sizeof(double));

  /*Increase denominator of acceptance rate for b and hl*/
  ndelta++;

  /*Draw proposal values for b and hl*/
  rmvnorm(pmd,var,2,parms->md,seed,1);

  /*Only proceed if we draw reasonable values*/
  if (pmd[0] > 0 && pmd[1] > 3) {

    /*Compute ratio of prior densities*/
    priorb_old = parms->md[0] - priors->meanbh[0];
    priorb_old *= 0.5*priorb_old/priors->varbh[0];
    priorb_new = pmd[0]-priors->meanbh[0];
    priorb_new *= 0.5*priorb_new/priors->varbh[0];
    priorh_old = parms->md[1] - priors->meanbh[1];
    priorh_old *= 0.5*priorh_old/priors->varbh[1];
    priorh_new = pmd[1]-priors->meanbh[1];
    priorh_new *= 0.5*priorh_new/priors->varbh[1];
    
    prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

    /*Save current values of b and hl and set new current values equal to 
      proposal values */
    for (k=0;k<2;k++) {
      currentmd[k] = parms->md[k];
      parms->md[k] = pmd[k];
    }

    /*Save current decay rate; calculate new current decay rate based on
      proposal value of hl */
    current_decay = parms->decay;
    parms->decay = log(2)/parms->md[1];

    /*Save current mean_contrib for each pulse; calculate new current
      mean_contrib for each pulse based on proposed values of b and hl*/
    i = 0;
    node = list->succ;
    while (node != NULL) {
      for (j=0;j<N;j++)
	currentmc[i][j] = node->mean_contrib[j];
        mean_contribution(node,ts,parms,N);
      i++;
      node = node->succ;
    }

    /*Calculate proposed likelihood and then calculate likelihood ratio*/
    plikelihood = likelihood(list,ts,parms,N,list);
    like_ratio = plikelihood - *like;

    /*Calculate log rho; set alpha equal to min(0,log rho) */
    alpha = (0 < (temp = (prior_ratio+like_ratio))) ? 0:temp;

    /*If log U < log rho, increase acceptance rate by 1 and set the
      likelihood equal to the likelihood under proposed values */
    if (log(kiss(seed)) < alpha) {
      adelta++;
      *like = plikelihood;
    }

    /*Otherwise, we need to revert values back to their current state */
    else {

      /*Set b and hl back equal to current values*/  
      parms->md[0] = currentmd[0];
      parms->md[1] = currentmd[1];

      /*Set mean_contrib back equal to current values for each pulse*/
      i = 0;
      node = list->succ;
      while (node != NULL) {
	for (j=0;j<N;j++) {
	  node->mean_contrib[j] = currentmc[i][j];
        }
	i++;
	node = node->succ;
      }

      /*Set decay rate back equal to current value*/
      parms->decay = current_decay;

    } /*End of if else statement*/

  }

  /*Free memory */
  for (i=0;i<num_node;i++)
    free(currentmc[i]);
  free(currentmc);
  free(pmd);
 
}

/*********************************************************************/
               /*START OF draw_fixed_effects SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_fixed_effects: this runs the M-H draw for overall mean mass and width;
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               Priors *priors; the parameters of the prior distributions;
               Common_parms *parms; the current values of the common parameters;
               unsigned long *seed; seed values needed for the randon number generator;
               double v1; the proposal variance for overall mean mass;
               double v2; the proposal variance for overall mean width;
    RETURNS: None; all updates are made internally */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
 new_sum: part of evaluation of RE_ratio portion of log(rho)
 old_sum: same as new_sum but with current (not proposed) value
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 prior_ratio: comparison between prior distributions under current and proposed values
 old_comp: part of evaluation of RE_ratio portion of log(rho)
 new_comp: same as old_comp but with proposed value
 RE_ratio: comparison between random effects distributions under current and
    proposed values
 temp: value used in evaluation of log(rho)
 alpha: minimum of 0 and log(rho); used to determine acceptance in M-H
 *new_mean: proposed values of overall mean mass and mean width
 *tmp1: current acceptance counters are saved to this vector so we can increment
    inside a loop
 *node: current list of pulses and their qualities

SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution
 kiss: found in randgen.c; draws from U(0,1) distribution
 print_list: found in hash.c; traverses the linked list and prints the contents
             of each node*/
/**************************************************************************/

void draw_fixed_effects(Node_type *list,Priors *priors,Common_parms *parms,unsigned long *seed)
{
/* declare variables */
   int j, numnode;
   double gmean,gvar,resum;
   Node_type *node;
   
/*declare functions */
   double rnorm(double,double,unsigned long *);
   double kiss(unsigned long *);
   void print_list(Node_type *);

    /*Draw new values of mu_a and mu_w via Gibbs Sampler*/
    for (j=0;j<2;j++){

        resum=0.0;
        numnode=0;
        node=list->succ;

        /*This while loop counts the pulses and gets the sum of r.e.*/
        while(node != NULL){
            numnode++;
            resum += log(node->theta[j]);
            node=node->succ;
            } /*end of loop through pulses*/

/*I added this on 11 Sep 2012; numnode is 0 to start with and this causes
 gvar to be very large, sometimes leading to values of mua and muw that are 
 too large or too small; this in turn leads to later extreme values of A-ik 
 and sigma-p-ik which causes precision issues in mean-contrib*/
        if(numnode == 0){
            numnode = 1;
        }
        gmean = (priors->fe_precision[j]*resum);
        gmean += (parms->re_precision[j]*parms->re_precision[j]*priors->fe_mean[j]);
        gmean /= (priors->fe_precision[j]*numnode) + (parms->re_precision[j]*parms->re_precision[j]);

        gvar = priors->fe_precision[j]*parms->re_precision[j]*parms->re_precision[j];
        gvar /= (priors->fe_precision[j]*numnode) + (parms->re_precision[j]*parms->re_precision[j]);

        parms->theta[j] = rnorm(gmean,sqrt(gvar),seed);
    }

}

/*********************************************************************/
               /*START OF draw_re_precision SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_re_precision: this runs the M-H draw for overall standard deviations of
                    mass and width;
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               Priors *priors; the parameters of the prior distributions;
               Common_parms *parms; the current values of the common parameters;
               double v1; the proposal variance for overall st-dev of mass;
               double v2; the proposal variance for overall st-dev of width;
               unsigned long *seed; seed values needed for the randon number generator;
    RETURNS: None; all updates are made internally */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 j: generic counter
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
 *node: current list of pulses and their qualities

SUBROUTINES USED
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution */
/**************************************************************************/

void draw_re_precision(Node_type *list,Priors *priors,Common_parms *parms,double v1,double v2,unsigned long *seed)
{
  int j;
  double *tmp1,*new_var,prop_new,prop_old,psum,pcomp,prior_old,prior_new;
  double prior_ratio,prop_ratio,temp,alpha;
  Node_type *node;
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);

  /*Allocate Memory*/
  new_var = (double *)calloc(2,sizeof(double));
  tmp1 = (double *)calloc(2,sizeof(double));

  /*Add 1 to the counters for acceptance rates of sigma_1 and sigma_2*/
  nrevm++;
  nrevw++;

  /*Assign current acceptance counts to temporary vector*/
  tmp1[0] = arevm;
  tmp1[1] = arevw;

  /*Draw proposed values for sigma_1 and sigma_2*/
  new_var[0]=rnorm(parms->re_precision[0],v1,seed);
  new_var[1]=rnorm(parms->re_precision[1],v2,seed);
  
  /*Accept or Reject sigma_1, then accept or reject for sigma_2*/
  for (j=0; j<2; j++){
      /*Calculate ratio of Cauchy priors*/
      /*prior_new = new_var[j]/priors->re_var[j];
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
      node = list->succ;
      prior_ratio = 1;
      while (node != NULL){
        prior_ratio *= parms->re_precision[j]/new_var[j];

        pcomp = log(node->theta[j]) - parms->theta[j];
        pcomp *= pcomp;
        psum += pcomp;
        node = node->succ;
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
      if(new_var[j]>0 && new_var[j]<priors->re_var[j]){

          /*Compute log rho, and set alpha equal to min(log rho,0)*/
 
          alpha = (0<(temp=(prop_ratio+prior_ratio)))?0:temp;

          /*If log(U) < log rho, accept the proposed value*/
          /*Increase acceptance count by 1*/
          if(log(kiss(seed))<alpha){
              tmp1[j]++;
              parms->re_precision[j]=new_var[j];
          }
      }
  }

  /*Set acceptance count equal to temp vector components*/
  arevm = tmp1[0];
  arevw = tmp1[1];

  free(new_var);
  free(tmp1);

}

/*********************************************************************/
              /*START OF draw_random_effects SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*draw_random_effects: this runs the M-H draw for individual pulse masses and
                      widths;
    ARGUMENTS: double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               Node_type *list; this is the current list of pulses that exist;
               Common_parms *parms; the current values of the common parameters;
               int N; the number of observations in **ts;
               double *like; the current value of the likelihood;
               double v1; the proposal variance for individual pulse masses;
               double v2; the proposal variance for individual pulse widths;
               unsigned long *seed; seed values needed for the randon number generator;
    RETURNS: None; all updates are made internally */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j: generic counters
 temp: value used in evaluation of log(rho)
 prior_old: part of evaluation of prior_ratio portion of log(rho)
 prior_new: same as prior_old, but with proposed value
 old_var: current value of the random effect is saved to this variable so we
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
 *node: current list of pulses and their qualities

SUBROUTINES USED
 mean_contribution: found in birthdeath.c; evaluates mean_contrib for a
    pulse with specific parameters
 likelihood: found in birthdeath.c; calculates and returns likelihood under
    current parameters
 kiss: found in randgen.c; draws from U(0,1) distribution
 rnorm: found in randgen.c; draws from the normal distribution */
/**************************************************************************/

void draw_random_effects(double **ts,Node_type *list,Common_parms *parms,int N,double *like,double v1, double v2,unsigned long *seed)
{
  int i,j;
  double temp,prior_old,prior_new,old_val,*tmp1,*pRE,prior_ratio,like_ratio,alpha,plikelihood,*old_contrib;
  Node_type *node;
  void mean_contribution(Node_type *,double **,Common_parms *,int);
  double likelihood(Node_type *,double **,Common_parms *,int,Node_type *);
  double kiss(unsigned long *);
  double rnorm(double,double,unsigned long *);

  /*Allocate Memory*/
  tmp1 = (double *)calloc(2,sizeof(double));
  old_contrib = (double *)calloc(N,sizeof(double));
  pRE = (double *)calloc(2,sizeof(double));

  /*Set acceptance counts equal to temporary vector*/
  tmp1[0] = arem;
  tmp1[1] = arew;

  /*Go to start of node list*/
  node = list->succ;

  /*Go through each existing pulse*/
  while (node != NULL) {
    /*Increase the denominators of the acceptance rates*/
    nrem++;
    nrew++;

    /*Draw proposed values of current pulse's mass and width*/
    pRE[0]=rnorm(log(node->theta[0]),v1,seed);
    pRE[1]=rnorm(log(node->theta[1]),v2,seed);

    /*Determine if we accept or reject proposed pulse mass then determine*/
    /*if we accept or reject proposed pulse width*/
    for(j=0;j<2;j++){
      /*We can only accept if we draw a non-negative mass/width */

      /*Compute the log of the ratio of the priors*/
      prior_old = log(node->theta[j]) - parms->theta[j];
      prior_old *= 0.5*prior_old;
      prior_new = pRE[j] - parms->theta[j];
      prior_new *= 0.5*prior_new;
      prior_ratio = prior_old - prior_new;
      prior_ratio /= parms -> re_precision[j];
      prior_ratio /= parms -> re_precision[j];

      /*Save the current value of mass/width*/
      old_val = node -> theta[j];

      /*Set the pulse's mass/width equal to the proposed value*/
      node -> theta[j] = exp(pRE[j]);

      /*Save the mean_contrib for that pulse*/
      for(i=0;i<N;i++)
          old_contrib[i]=node->mean_contrib[i];

      /*Recalculate that pulse's mean_contrib assuming proposed mass/width */
      mean_contribution(node,ts,parms,N);

      /*Calculate likelihood assuming proposed mass/width */
      plikelihood = likelihood(list,ts,parms,N,list);

      /*Compute the log of the ratio between the two likelihoods*/
      like_ratio = plikelihood - *like;

      /*Calculate log rho; set alpha equal to min(0,log rho) */
      alpha = (0 < (temp = (prior_ratio + like_ratio))) ? 0:temp;

      /*If log U < log rho, accept the proposed value, increase acceptance*/
      /*counter and set current likelihood equal to likelihood under proposed*/
      /*mass/width */
      if (log(kiss(seed)) < alpha) {
	tmp1[j]++;
	*like = plikelihood;
      }

      /*Otherwise, reject the proposed value, set pulse's mass/width back to */
      /*saved current value, and set pulse's mean_contrib back to saved value */
      else {
	node->theta[j] = old_val;
	for (i=0;i<N;i++)
	  node->mean_contrib[i] = old_contrib[i];
      }
    
    } /*end of loop through mass & width */

    /*Advance to next pulse*/
    node = node->succ;
    
  } /*end of loop through pulses*/

  /*Set counter equal to temporary vector values*/
  arem = tmp1[0];
  arew = tmp1[1];

  /*free memory*/
  free(tmp1);
  free(pRE);
  free(old_contrib);
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
                  times and a column of log(concentration)
               Node_type *list; this is the current list of pulses that exist
               Common_parms *parms; the current values of the common parameters
               int N; the number of observations in **ts;
    RETURNS: double ssq; this value is needed for one of the parameters
             of the posterior distribution of model error variance */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *mean: mean concentration under current parameters and pulses
 ssq: sum of squared differences between observed log(concentration) and expected
    values

SUBROUTINES USED
 mean_concentration: found in birthdeath.c; sums across mean_contrib vectors for
    each pulse; returns this vector  */
 /**************************************************************************/

double error_squared(double **ts,Node_type *list,Common_parms *parms,int N) {
 int i;
 double *mean,ssq;
 double *mean_concentration(Node_type *list,Common_parms *parms,int N,Node_type *node_out,double **ts);

 mean = mean_concentration(list,parms,N,list,ts);
ssq = 0;
 for (i=0;i<N;i++)
   ssq += (ts[i][1] - mean[i])*(ts[i][1] - mean[i]);

 free(mean);
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
    RETURNS: None; update to the proposal variance is made internally */
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
    RETURNS: None; update to the proposal variance is made internally  */
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
