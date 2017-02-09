/*******************************************************************/
/*********************MAIN BDMCMC PROGRAM***************************/
/*******************************************************************/

#include "deconvolution_main.h"

#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT 64
#else
#define ENVIRONMENT 32
#endif
#endif

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT 64
#else
#define ENVIRONMENT 32
#endif
#endif
/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 * M: Sets the environment precision and used for all random number generation
    * Needed specfically for the KISS subroutine
 * mmm: Order statistic used for distribution of pulse locations.
    * This is inputted by the user and is typically 3
 * fitstart: The first time in hours that a pulse may occur
 * fitend: The last time in hours that a pulse may occur

*********************************************************************/

double M;
int mmm;
double fitstart;
double fitend;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 start_position: This subroutine initializes all the parameters in the study.
    Initial baseline and halflife are user-inputted; all others are set in
    the subroutine. This subroutine also creates our list of nodes and inserts
    one new node with parameters equal to 0.
    ARGUMENTS: Common_parms *parms; The values in this data structure are
               modified in this routine;
               double pmean; starting value for baseline inputted by user in
               main; must be greater than 0
               double pdelta; starting value for halflife inputted by user in
               main; must be greater than 0
    RETURNS: Node_type *list; returns the initialized node list
 **********************************************************************/

/**********************************************************************
 MAIN PROGRAM STARTS HERE

 VARIABLE DEFINITIONS
    i: generic counter
    *N: number of observations in inputted file; determined by read_data_file
        subroutine
    iter: number of iterations to run the MCMC subroutine called in main;
          this variable is inputted by the user
    *seed: 3 element vector with seed values for the random number generator;
           read in from an external file
    **ts: matrix containing column of time points and column of hormone
          concentrations
    **pmd_var: Variance-Covariance matrix of the proposal distribution of
               baseline and halflife; matrix is defined in main
    pmean: Starting value for baseline; inputted by the user
    pdelta: Starting value for halflife; inputted by the user
    *fseed: The pointer to the file that contains the random numbers for the
            random number generator
    *parms: Contains the parameters that are common throughout the model
    *priors: Contains all parameters of the prior distributions used
    *hyper: Contains all the hyperpriors (though this is not used)
    *list: Contains the list of nodes and their characteristics

 SUBROUTINES USED
    **read_data_file: Found in format_data.c; scans inputted data file and
                      returns the matrix of data (time and concentration)
    rnorm: Found in randgen.c; draws from the normal distribution
    rgamma: Found in randgen.c; draws from the gamma distribution
    destroy_list: Found in hash.c; frees all memory used by the nodes
    mcmc: Found in linklistv2.c; This runs birth-death mcmc routine
    *start_position: Found in this file; discussed above 
 ************************************************************************/

int main(int argc,char *argv[])
{
    
/* Declare variables*/
  int i,*N,iter,a;

  char datafile[90];
  char common1[100];
  char parm1[100];
  unsigned long *seed;
  double **ts,**pmd_var,propvar[7],priormu1,priorvar1;
  double priormu2,priorvar2,priormub,priorvarb,priormuh,priorvarh;
  double prioralpha, priorbeta, priora1,priora2, priorr;
  double svmu1, svmu2, svbase, svhalf, svevar, svsig1, svsig2;
  double vrem,vrew,vm,vw,vmv,vwv,vt;
  double time[9],mass[9],width[9];

  FILE *fseed, *finput;
  Common_parms *parms;
  Priors *priors;
  Node_type *list,*new_node;

/* Declare subroutines */
  double **read_data_file(char *,int *);
  double rnorm(double,double,unsigned long *);
  double rgamma(double,double,unsigned long *);
  void destroy_list(Node_type *);
  void mcmc(Node_type *,Common_parms *,double **,long,int,Priors *,unsigned long *,
          char *,char *,double []);
  Node_type *start_position(Common_parms *,double,double);
  Node_type *initialize_node(void);
  void insert_node(Node_type *,Node_type *);
  void mean_contribution(Node_type *,double **,Common_parms *,int);

  if (argc != 2 ) exit(0);

    
/*************************************************/
/* M is used in the KISS random number generator */
    
    M = exp(-ENVIRONMENT*log(2.0));
/*************************************************/

    
/* READ IN THE SEED FILE THAT MUST BE IN THE DIRECTORY WHERE THE PROGRAM IS BEING RUN */
    seed = (unsigned long *)calloc(3,sizeof(unsigned long));
    fseed = fopen("seed.dat","r");
    fscanf(fseed,"%lu %lu %lu\n",&seed[0],&seed[1],&seed[2]);
    fclose(fseed);

/* Open the input file and read in the file names, prior parameters and starting values */
    finput = fopen(argv[1],"r");

    fscanf(finput,"%s \n", datafile);   /* read in the data file */
    fscanf(finput,"%s %s \n", common1, parm1);  /*read in the base files names for the two major output files (common parameters and pulse specific parameters) */
    fscanf(finput,"%d \n", &iter);  /* The total run time for the MCMC */


/*read in the hormonal time series*/
    N = (int *)calloc(1,sizeof(int));
    ts = read_data_file(datafile,N);
    
    mmm = 3;  /*specifies the 3rd order statistics for the pulse location model*/


    
/* Start reading in the parameters that define the prior distributions*/
    
    parms = (Common_parms *)calloc(1,sizeof(Common_parms));

/*Create the boundaries of pulse locations just slightly before and after data collection*/
    fitend = ts[*N-1][0]+ ts[0][0] * 2;  /*search 2 units farther in time*/
    fitstart = -ts[0][0] * 4;  /*search 4 units in the past*/

    priors = (Priors *)calloc(1,sizeof(Priors));
    priors->re_var = (double *)calloc(2,sizeof(double));
    priors->fe_precision = (double *)calloc(2,sizeof(double));
    
    fscanf(finput,"%lf %lf\n", &priormu1, &priorvar1); /*Prior mean and variance for the mean mass on the log scale*/
    fscanf(finput,"%lf %lf\n", &priormu2, &priorvar2); /*Prior mean and variance for the mean pulse width on the log scale*/
    fscanf(finput,"%lf %lf\n", &priormub, &priorvarb); /*Prior mean and variance for the baseline parameter*/
    fscanf(finput,"%lf %lf\n", &priormuh, &priorvarh); /*Prior mean and variance for the halflife parameter*/
    fscanf(finput,"%lf %lf\n", &prioralpha, &priorbeta); /*Parameters in the inverse gamma for the model error variance*/
    fscanf(finput,"%lf %lf\n", &priora1, &priora2); /* Maximums in the Uniform priors for the pulse mass and width pulse-to-pulse standard deviations*/
    fscanf(finput,"%lf\n", &priorr); /*Prior pulse rate in the Poisson prior on the pulse number*/

/*Set the values in the data structure of the priors*/
    priors->meanbh[0] = priormub;
    priors->meanbh[1] = priormuh;
    priors->varbh[0] = priorvarb;
    priors->varbh[1] = priorvarh;

    priors->fe_mean = (double *)calloc(2,sizeof(double));

    priors->fe_mean[0] = priormu1;
    priors->fe_mean[1] = priormu2;
    priors->fe_precision[0] = priorvar1;
    priors->fe_precision[1] = priorvar2;
    priors->re_var[0] = priora1;  /*not really variances, SD stored--not used as a matrix*/
    priors->re_var[1] = priora2;  /*not really variances, SD stored*/
    
    priors->alpha = prioralpha;
    priors->beta = priorbeta;

    parms->nprior = priorr;

/* Read in the starting values for all the parameters */
    fscanf(finput,"%lf %lf\n", &svmu1, &svmu2); /*Starting value mean pulse mass, mean pulse width */
    fscanf(finput,"%lf %lf\n", &svbase, &svhalf); /*Starting value baseline and half-life */
    fscanf(finput,"%lf\n", &svevar); /*Starting value error variance*/
    fscanf(finput,"%lf %lf\n", &svsig1, &svsig2); /*Starting value standard deviation of pulse-to-pulse SD of mass and widths--in random effects distribution*/

/*Read in the proposal variances for the MH parts of the MCMC algorithm*/
    fscanf(finput,"%lf %lf\n", &propvar[0], &propvar[1]);  /* Variance for baseline and half-life */
    fscanf(finput,"%lf %lf\n", &propvar[2], &propvar[3]);  /* Variance for pulse-to-pulse SD for log pulse mass and width */
    fscanf(finput,"%lf %lf %lf\n", &propvar[4], &propvar[5], &propvar[6]);  /* ?? */

    fclose(finput);  /*close the input file*/

    /*Set the starting values for  mean pulse mass (mu1) and mean pulse width (mu2)--both on log scale */
    parms->theta[0] = svmu1;
    parms->theta[1] = svmu2;

    /*Set the starting values for baseline and half life*/
    parms->md[0] = svbase;
    parms->md[1] = svhalf;
    parms->decay = log(2)/parms->md[1];

    /*Set the starting value for the error variance***/
    parms->sigma = svevar;
    parms->lsigma = log(parms->sigma);

    parms->re_precision = (double *)calloc(2,sizeof(double));
    
    /*Initialize the starting values for the pulse-to-pulse variation parameters for mass (1) and width (2)*/
    parms->re_precision[0] = svsig1;
    parms->re_precision[1] = svsig2;
    
    
    /*Initialize the pulse link list so an initial pulse can be added in the first MCMC step*/
    list = initialize_node();

    
    /*enter the BDMCMC loop--most of the work is in this function*/

    mcmc(list,parms,ts,iter,*N,priors,seed,common1,parm1,propvar);
  
    /* MCMC loop is complete and do clean up now*/
    
    destroy_list(list);
    /**************************/
    
    /* save the current random number as the seed for the next simulation */
    fseed = fopen("seed.dat","w");
    fprintf(fseed,"%lu %lu %lu\n",seed[0],seed[1],seed[2]);
    fclose(fseed); 
    /**********************************************************************/

    /* deallocate resources */

    free(seed);
    for (i=0;i<*N;i++)
      free(ts[i]);
    free(ts);
    free(N);

    free(priors->fe_mean);

    free(priors->fe_precision);

    free(parms->re_precision);
    free(priors->re_var);

    
    free(priors);
    free(parms);
/************************/
  return 0;
}

