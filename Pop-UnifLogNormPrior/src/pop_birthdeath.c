/*******************************************************************/
/************************birthdeath.c ******************************/
/*******************************************************************/

#include "pop_deconvolution_main.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 fitstart: The first time in hours that a pulse may occur
 fitend: The last time in hours that a pulse may occur
 mmm: Order statistic used for distribution of pulse locations.
      This is inputted by the user and is typically 3

*********************************************************************/

extern double fitstart;
extern double fitend;
extern int mmm;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 birth_death
 mean_contribution
 *mean_concentration
 likelihood
 *calc_death_rate

**********************************************************************/

/*********************************************************************/
                 /*START OF birth_death SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*birth_death: this runs the birth-death step of the BDMCMC process
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               Common_parms *parms; the current values of the common parameters;
               int N; the number of observations in **ts;
               double *likeli; the current value of the likelihood;
               unsigned long *seed; seed values needed for the randon number generator;
               int iter; which iteration are we on;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j,k: generic counters
 remove: if a death occurs, remove is a draw from the multinomial distribution
    and represents which pulse we will remove
 num_node: number of pulses
 flag: if a birth occurs, we must draw parameters for this new pulse; flag will
    be set to 0 until a positive mass and a positive width are drawn
 aaa: this counter ensures we do not run the birthdeath loop too many times
 max_num_node=60: this variable is set to 60 and ensures a death if we exceed
    60 pulses
 S: the sum of exponential draws; the birthdeath loop stops when S exceeds T
 Birth_rate: a constant birth rate needed for the process
 T=1: stop birthdeath loop when S exceeds T=1
 full_likelihood: current likelihood before a birth or death occurs
 max: part of the calculation for Death_rate
 Death_rate: total death rate
 *death_rate: vector of death rates for each pulse
 position: if a birth occurs, the pulse's location is position, drawn from a
    uniform distribution
 *partial_likelihood: a vector of likelihoods where the ith element is the
    likelihood if pulse i were removed
 new_theta: used when drawing new pulse's mass and width
 new_tsd : the std in the t-prior for drawing new pulse masses and widths
 new_kappa: used when drawing the new pulsess massess and widths
 *node: counter through pulses
 *new_node: if a birth occurs, we give the new pulse this name

 SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution
 rmultinomial: found in randgen.c; draws from the multinomial distribution
 kiss: found in randgen.c; draws from U(0,1) distribution
 rexp: found in randgen.c; draws from the exponential distribution
 runif_atob: found in randgen.c; draws from U(a,b) distribution
 mean_contribution: found in this file; evaluates mean_contrib for a
    pulse with specific parameters
 likelihood: found in this file; calculates and returns likelihood under
    current parameters Node_type *initialize_node(void);
 *initialize_node: found in hash.c; allocates memory for a new pulse and creates it;
 insert_node: found in hash.c; inserts a newly created node in the appropriate
    spot in the linked list
 delete_node: found in hash.c; eliminates a node and links the neighbors
 *calc_death_rate: found in this file; creates a vector where the ith element
    represents the death rate for pulse i   */
 /***********************************************************************/

void birth_death(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                 unsigned long *seed,int iter)
{
  int i,j,k,remove,num_node,flag,aaa,max_num_node=60;
  double S,Birth_rate,T=1.,full_likelihood,full_likelihood2,max,new_kappa, new_tsd;
  double Death_rate,*death_rate,position,*partial_likelihood,new_theta;
  Node_type *node,*new_node,*list;
  Subject_type *subject;
  double rnorm(double,double,unsigned long *);void print_list(Node_type *);
  long rmultinomial(double *,long,unsigned long *);
  double kiss(unsigned long *);
  double rexp(double, unsigned long *);
  double runif_atob(unsigned long *,double,double);
    double rgamma(double,double,unsigned long *);
  void mean_contribution(Node_type *,double **,Common_parms *,int,double);
  double likelihood(Subject_type *list,double **ts,Common_parms *parms,int N,
                    Node_type *node_out);
  double likelihood2(Node_type *list,double **ts,Common_parms *parms,int N,
                  Node_type *node_out, double baseline);
  Node_type *initialize_node(void);
  void insert_node(Node_type *new_node,Node_type *list);
  void delete_node(Node_type *node,Node_type *list);
  double *calc_death_rate(Node_type *,int,double *,double,double,double);


  /*For first 100 iterations, birth rate is 10, then it becomes 1;
   This is to force more births early on.*/
 /* if (iter < 100) {
    Birth_rate = 10;
  }
  else {
    Birth_rate = 1;
  }
*/
    
  Birth_rate = parms->nprior;
    
  subject=sublist->succ;

  parms->subindex = 0;

  while (subject != NULL){

     list = subject->list;
      aaa = 0;

      S = 0.0;

  /*Save Likelihood*/
  /*full_likelihood = *likeli;
*/
  /*Go until the loop is broken*/

      while (1) {

          /*This counter keeps this loop from running too many times*/
          aaa++;
  
          /*Count number of pulses*/
          num_node = 0;
          node = list->succ;
          while (node != NULL) {
              num_node++;
              node = node->succ;
          }

    
    /*Allocate memory for partial likelihood vector*/
    partial_likelihood = (double *)calloc(num_node,sizeof(double));

    /*Calculate the likelihood if pulse i is removed*/
    i = 0;
    node = list->succ;

  /*  printf("num_node = %d \n",num_node);*/

    while (node != NULL) {
        partial_likelihood[i] = likelihood2(subject->list,ts,parms,N,node,subject->basehalf[0]);
      /*printf("partial likelihood = %lf \n",partial_likelihood[i]);*/
      i++;
      node = node->succ;
    }

    full_likelihood2=likelihood2(subject->list,ts,parms,N,subject->list,subject->basehalf[0]);

    /* CALCULATE DEATH RATE FOR EACH COMPONENT */
    death_rate = NULL;

    death_rate = calc_death_rate(list,num_node,partial_likelihood,full_likelihood2,Birth_rate,parms->nprior);

/*    node = list->succ;
    i=0;
    while (node != NULL) {
      printf("death_rate = %lf \n",death_rate[i]);
      i++;
      node = node->succ;
    }
*/
    /*The next portion computes D = sum(d_i); This is a little more complicated
     than summing them up because of precision issues.*/
    if (death_rate != NULL) {
      Death_rate = death_rate[0];
      for (i=1;i<num_node;i++) {
          max = (Death_rate > death_rate[i]) ? Death_rate:death_rate[i];
          Death_rate = log(exp(Death_rate-max) + exp(death_rate[i]-max)) + max;
      }

      for (i=0;i<num_node;i++) {
          death_rate[i] -= Death_rate;
          death_rate[i] = exp(death_rate[i]);
      }

      for (i=1;i<num_node;i++)
          death_rate[i] += death_rate[i-1];

        Death_rate = exp(Death_rate);
    }
    else 
      Death_rate = 0;

  
    free(partial_likelihood);


    if (num_node <= 1) Death_rate = 0;

    /* Draw from exp(B+D) and add to current S */
    S += rexp(Birth_rate + Death_rate,seed);

    /*If S exceeds T or if we've run this too many times, break*/
    if (S > T)  break;
    if (aaa > 5000) break;

    /* SIMULATE JUMP TYPE (BIRTH OR DEATH) */
    /* If we only have 0 or 1 pulses, we definitely have a birth*/
    if (num_node <= 1)
      max = 1.1;

    /* If we have too many pulses (60), we definitely have a death*/
    else if (num_node >= max_num_node)
      max = -0.1;
 
    /* Otherwise, set max = B/(B+D) */
    else
      max = Birth_rate/(Birth_rate + Death_rate);

/*    printf("Birth_rate = %lf \n",Birth_rate);
   printf("Death_rate = %lf \n",Death_rate);
*/
    /*If U < B/(B+D), a birth occurs */
    if (kiss(seed) < max) {  

  
        /*simulate a pulse position uniformly */
        position = runif_atob(seed,fitstart,fitend);
	
        node = list->succ;

        /* initialize a new node */
        new_node = initialize_node();
        new_node->time = position;
        
      /* simulate random effects from the prior */
        for (j=0;j<2;j++) {
            new_theta = -1;
            new_kappa = rgamma(2,2,seed);
            new_tsd = sqrt((parms->re_precision[0] * parms->re_precision[0]) / new_kappa);
            while (new_theta < 0 ) {
                new_theta = rnorm(subject->theta[j],new_tsd,seed);
            }
            new_node->theta[j] = new_theta;
            new_node->kappa[j] = new_kappa;
        }
        
        new_node->mean_contrib = (double *)calloc(N,sizeof(double));
        mean_contribution(new_node,ts,parms,N,subject->basehalf[1]);
          /* insert the new node; this routine also makes sure it's
           located properly */
        insert_node(new_node,list);
    	
      }

    /*Otherwise, a death occurs */
    else { 
      /* Pick a node to remove based on the multinomial */
      
      remove = (int)rmultinomial(death_rate,(long)num_node,seed) + 1;
      node = list;
      for (i=0;i<remove;i++)
        node = node->succ;
      delete_node(node,list);
    }
    free(death_rate);
  /*  fflush(stdout);*/
  } /*End of While Loop*/

  /*Update likelihood before leaving Birthdeath.c*/
  /**likeli = full_likelihood;
*/
  
  free(tmp);
/*free(death_rate);*/

  subject=subject->succ;
  parms->subindex++;

  }

  }

/*********************************************************************/
               /*START OF mean_contribution SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mean_contribution: this updates a pulse's mean_contrib vector based on current
                     values of parameters
    ARGUMENTS: Node_type *node; what pulse are we updating;
               double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               Common_parms *parms; the current values of the common parameters;
               int N; the number of observations in **ts;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 x,y,z,tmp: these are all variables that are part of the arithmetic needed in
    calculating the mean contribution

 SUBROUTINES USED
 erf: this is a subroutine used to help in integration  */
 /***********************************************************************/

void mean_contribution(Node_type *node,double **ts,Common_parms *parms,int N,double halflife)
{
/* This function updates a pulse's mean_contrib vector based on inputted parms*/
  int i;
  double x,y,z,tmp;
  double erf(double);

/*  printf("theta1 = %lf \n",node->theta[1]);
  printf("halflife = %lf \n",halflife);
  printf("time = %lf \n",node->time);
  printf("theta0 = %lf \n",node->theta[0]);*/
  
  z = node->theta[1]*0.6931472/halflife;
  y = 0.6931472*(0.5*z/halflife + node->time/halflife);
  z += node->time;
  tmp = sqrt(2.*node->theta[1]);
  for (i=0;i<N;i++) {
    x = (ts[i][0] - z)/tmp;
    x = 0.5*(1.+erf(x));
    if(x == 0) node->mean_contrib[i] = 0;
    else node->mean_contrib[i] = node->theta[0]*x*exp(y - ts[i][0]*0.6931472/halflife);
   
  }
 /* printf("mc55 = %lf \n",node->mean_contrib[55]);

fflush(stdout);
*/
}

/*********************************************************************/
               /*START OF mean_concentration SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*mean_concentration: this takes each pulse's mean_contrib vector and sums
                      across them
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               Common_parms *parms; the current values of the common parameters;
               int N; the number of observations in **ts;
               Node_type *node_out; if we want, we can ignore a pulse;
               double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
    RETURNS: x, the vector of sums  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *x: the vector of sums
 *node: the counter of pulses

 SUBROUTINES USED
 rnorm: found in randgen.c; draws from the normal distribution  */
/***********************************************************************/

double *mean_concentration(Node_type *list,Common_parms *parms,int N,Node_type *node_out,double baseline)
{
 /* This function sums the mean_contrib vector across pulses */
 /* The result is a vector of sums */
 int i;
 double *x;
 Node_type *node;
 double rnorm(double,double,unsigned long *);

 x = (double *)calloc(N,sizeof(double));

 /* add the contribution to the mean from each pulse */
 node = list->succ;
 while (node != NULL) {
   if (node != node_out) {
     for (i=0;i<N;i++)
       x[i] += node->mean_contrib[i];
   }
/*printf("mc15 = %lf \n",node->mean_contrib[15]);
printf("x45 = %lf \n",x[45]);
*/
   node = node->succ;
 }

/* add the baseline contribution */
 for (i=0;i<N;i++) {
   x[i] += baseline;
   x[i] = log(x[i]);
 }

 return x;
}

/*********************************************************************/
                   /*START OF likelihood SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*likelihood: computes the current likelihood using the observed log-concentrations
              and mean concentration
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               double **ts; this is the matrix of observed data (a column of
                  times and a column of log(concentration);
               Common_parms *parms; the current values of the common parameters;
               int N; the number of observations in **ts;
               Node_type *node_out; if we want, we can ignore a pulse;

    RETURNS: x, the likelihood computed based on the inputs  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 x=0: likelihood to be calculated; initialized at zero
 *mean: vector of sums of mean_contribs calculated using mean_concentration

 SUBROUTINES USED
 mean_concentration: found in this file; computes a vector of sums of mean_contribs  */
 /***********************************************************************/

double likelihood(Subject_type *sublist,double **ts,Common_parms *parms,int N,
                  Node_type *node_out)
{
  /* This function computes the likelihood under inputted parameters */
  /* The output is a scalar */
  int i,j;
  double x=0,*mean;
  Subject_type *subject;

  double *mean_concentration(Node_type *,Common_parms *,int,Node_type *,double );

  /* Sum across mean_contribs */
  subject=sublist->succ;

  for (j=0;j<parms->numsub;j++){

      
      /*13*/
   /* if(j==parms->subindex)
        mean = mean_concentration(subject->list,parms,N,node_out,subject->basehalf[0]);
      printf("mean eq = %lf \n",mean);
*/
  /*  if(j!=parms->subindex)*/
      mean = mean_concentration(subject->list,parms,N,subject->list,subject->basehalf[0]);


    for (i=0;i<N;i++)
        x += (ts[i][j+1]-mean[i])*(ts[i][j+1]-mean[i]);


    subject=subject->succ;

    free(mean);

  }

  x /= (-2.0*parms->sigma);
  x += -0.5*N*parms->numsub*(1.8378771+parms->lsigma);

  
  return x;
}

double likelihood2(Node_type *list,double **ts,Common_parms *parms,int N,
                  Node_type *node_out,double baseline)
{
  /* This function computes the likelihood under inputted parameters */
  /* The output is a scalar */
  int i,j;
  double x=0,*mean;

  double *mean_concentration(Node_type *,Common_parms *,int,Node_type *,double );

  j=parms->subindex;

  /* Sum across mean_contribs */
  mean = mean_concentration(list,parms,N,node_out,baseline);
  for (i=0;i<N;i++)
    x += (ts[i][j+1]-mean[i])*(ts[i][j+1]-mean[i]);

  x /= (-2.0*parms->sigma);
  x += -0.5*N*(1.8378771+parms->lsigma);
  free(mean);
  return x;
}



/*********************************************************************/
                /*START OF calc_death_rate SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*calc_death_rate: calculates a vector of death rates, one for each existing pulse
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               int num_node; current number of pulses
               double *partial_likelihood; vector of likelihoods, where the ith
                  element represents the likelihood with the ith pulse removed
               double full_likelihood; value of the full likelihood
               double Birth_rate; value of the birth rate

    RETURNS: death_rate; a vector where the ith element represents the death
             rate of the ith pulse  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j: generic counters
 x: variable for death rate
 *death_rate: vector of death rates
 coef_denom, coef_num: part of the calculation of death rate
 *node: counter for pulses

 SUBROUTINES USED
  None  */
 /***********************************************************************/
double *calc_death_rate(Node_type *list,int num_node,double *partial_likelihood,double full_likelihood,double Birth_rate,double r)
{
  /* This function calculates the death rate vector */
  int i,j;
  double x,*death_rate;
  double coef_denom,coef_num;
  Node_type *node;

  if (num_node > 1) {
    death_rate = (double *)calloc(num_node,sizeof(double));
    node = list->succ;
    i = 0;

    /* Calculate the coefficient of the distribution of the taus conditional
     on number of pulses. In this portion, I have an extra num_node in the
     numerator */
    coef_num = 1;
    for (j=1;j<mmm;j++)
      coef_num *= j;
    coef_denom = mmm*num_node;
    for (j=1;j<mmm;j++)
      coef_denom *= (mmm*num_node+j);
    while (node != NULL) {
      if (mmm > 0)

          /*This computes the portion of death rate due to the Poisson prior on
         N, the birth rate, and the likelihood */
       	x = log(num_node*Birth_rate/r) + partial_likelihood[i] - full_likelihood;
     else
	x = log(Birth_rate/num_node) + partial_likelihood[i] - full_likelihood;

      /*Now, compute the portion of the death rate due to distribution of the
       taus conditional on number of pulses. */
      if (mmm > 1) {
	if (i==0)

            /*If we are on the first pulse*/
            x += (mmm-1)*log(fitend-fitstart) +
	    (mmm-1)*log((node->succ->time - fitstart)/((node->succ->time-node->time)*(node->time-fitstart)))+log(coef_num/coef_denom);

        else if (i>0) {

            /*If we are not on the first pulse */
            if (node->succ)
	    x += (mmm-1)*log(fitend-fitstart) + 
	      (mmm-1)*log((node->succ->time - node->pred->time)/((node->succ->time-node->time)*(node->time-node->pred->time)))+log(coef_num/coef_denom);

            /*If we are not on the first or last pulse */
            else
	    x += (mmm-1)*log(fitend-fitstart) + 
	      (mmm-1)*log((fitend - node->pred->time)/((fitend-node->time)*(node->time-node->pred->time))) + log(coef_num/coef_denom);
        }
      }

      /*Save to death rate vector */
      death_rate[i] = x;



      /*Advance to next pulse*/
      node = node->succ;

      /*Advance counter 1*/
      i++;
    }
  }

  else {

      /*if we have 0 or 1 pulses*/
      if (num_node==0) /*if we don't have any pulses, no death rate*/
      return NULL;

      else {

          /*If we have 1 pulse, don't kill it*/
          death_rate = (double *)calloc(num_node,sizeof(double));
          death_rate[0] = -1e300;
    }
  }
  return death_rate;  
}


