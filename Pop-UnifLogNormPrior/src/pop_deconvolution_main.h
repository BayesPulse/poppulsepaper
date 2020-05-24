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
  double kappa[2];  /*scale in the t-distribution of the pulse mass and width*/
  double *mean_contrib;  /*mean conc. for this pulse for each time point*/
} Node_type;


/*mean of the mean and variance random effects of each pulse*/
/*inverse of the variance of the mean and variance pulse random effects*/

typedef struct {
  double sigma;  /*model error variance for the individual*/
  double lsigma; /*log of hte model error variance for the individual*/
  double theta[2];  /*individual level mean pulse mass and log pulse width:see defined below in subj structure.  Careful with use.  Double check what is going on*/
  double *re_precision; /*individual level pulse-to-pulse SD:  NOT PRECISION!*/
  int numsub; /* number of subjects total*/
  int subindex; /*counter*/
  int iter; /*iteration counter*/
  double nprior; /*the average number of pulses--parameter in the prior on the number of pulses*/
} Common_parms;

typedef struct {
  double meanbh[2]; /*prior mean on baseline and halflife */
  double varbh[2]; /*variance on prior of baseline and halflife */
  double *fe_mean; /*population mean pulse mass and width*/
  double *fe_precision; /*sub-to-sub SD (check) */
  double *re_var; /* Maxium values in the Unif prior defining the pulse-to-pulse SD of pulse mass and width*/
  double alpha;  /* prior parameters for model error */                     
  double beta;
 
} Priors;

typedef struct {
    double hmean[2]; /* means in the prior for population mean pulse mass and width*/
    double hvar[2]; /*Variances in the prior for the population mean pulse mass and width*/
    double prec[2]; /* Maximum value of the pulse-to-pulse SD of the individual pulse mass and width random effects*/
    double meanmeanbh[2]; /*means in the priors for the population mean baseline and halflife*/
    double meanvarbh[2]; /*Variances in the priors for the population mean baseline and half-life*/
    double varbh[2]; /*upper bound on the uniform prior for the sd of the sub-to-sub baseline and half-lives*/

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
