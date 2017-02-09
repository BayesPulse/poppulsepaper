#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

typedef struct node_tag {
  struct node_tag *succ;
  struct node_tag *pred;
  double time;  /*time of the pulse in minutes*/
  double theta[2];  /*mass of pulse, width of pulse(variance)*/
/*  double width;*/
  double *mean_contrib;  /*mean conc. for this pulse for each time point*/
} Node_type;


/*mean of the mean and variance random effects of each pulse*/
/*inverse of the variance of the mean and variance pulse random effects*/

typedef struct {
  double sigma;  /*model error variance for the individual*/
  double lsigma; /*log of hte model error variance for the individual*/
  double theta[2];  /*individual level mean log pulse mass and mean log pulse width*/
  double *re_precision; /*individual level pulse-to-pulse SD:  NOT PRECISION!*/
  int numsub; /* number of subjects total*/
  int subindex; /*counter*/
  int iter; /*iteration counter*/
  double nprior; /*the average number of pulses--parameter in the prior on the number of pulses*/
} Common_parms;

typedef struct {
  double meanbh[2]; /*prior mean on baseline and halflife */
  double varbh[2]; /*variance on prior of baseline and halflife */
  double *fe_mean;
  double *fe_precision;
  double *re_var;
  double alpha;  /* prior parameters for model error */                     
  double beta;
 
} Priors;

typedef struct {
    double hmean[2];
    double hvar[2];
    double prec[2];
    double meanmeanbh[2];
    double meanvarbh[2];
    double varbh[2];

} Hyper_priors;

typedef struct subject_tag {
    struct subject_tag *succ;
    struct subject_tag *pred;
    Node_type *list;
    double theta[2];
    double basehalf[2];
    double decay;
    char *common;
    char *pulse;
    FILE *csub;
    FILE *psub;
} Subject_type;
