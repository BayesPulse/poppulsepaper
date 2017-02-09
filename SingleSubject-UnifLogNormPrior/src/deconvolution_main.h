#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

typedef struct node_tag {
  struct node_tag *succ;
  struct node_tag *pred;
  double time;
  double theta[2];
  double *mean_contrib;
} Node_type;


/*mean of the mean and variance random effects of each pulse*/
/*inverse of the variance of the mean and variance pulse random effects*/

typedef struct {
  double md[2];
  double sigma;
  double lsigma;
  double theta[2];  
  double *re_precision;
  double decay;
  double nprior;
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
