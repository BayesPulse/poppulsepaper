/*******************************************************************/
/************************format_data.c *****************************/
/*******************************************************************/

#include "deconvolution_main.h"

/*********************************************************************/
                /*START OF read_data_file SUBROUTINE*/
             /*This is the only subroutine in this file*/
/*********************************************************************/
/*********************************************************************/
/*read_data_file: this scans the inputted data file and saves the second two
                  columns; they contain the time and the concentration. It
                  log-transforms the concentration before outputting the data
    ARGUMENTS: char *datafile; the name of inputted data file
               int *N; the variable where we will store the number of observations
    RETURNS: data; the Nx2 matrix of times and log-concentrations    */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
   int i,j;
  char *labels;
  double x,y,**data;
  FILE *inf;
 i: generic counter
 j: where we store the observation number
 *labels: this is where we store the column names
 x: we store the data in column 2 here while we count the number of rows
 y: we store the data in column 3 here while we count the number of rows
 **data: we read the entries in columns 2 and 3 (minus the header) into the matrix
    data
 *inf: infile variable

 SUBROUTINES USED
    None
**************************************************************************/

double **read_data_file(char *datafile,int *N)
{
  int i,j;
  char *labels;
  double x,y,**data;
  FILE *inf;

/* scans the data file and count the number of records */
/* the first column is the observation number          */
/* the second column is the observation time           */
/* the third column is the observation itself          */

  inf = fopen(datafile,"r");
  if (inf == NULL) {
    printf("file %s does not exist\n",datafile);
    free(N);
    exit(0);
  }

  *N=0;
 /* not interested in the column names, read them in and forget about them */
  labels = (char *)calloc(50,sizeof(char *));
  fscanf(inf,"%s %s %s\n",labels,labels,labels);
 /**************************************************************************/
  while (fscanf(inf,"%d %lf %lf\n",&j,&x,&y) != EOF) 
    (*N)++;

  rewind(inf);
  fscanf(inf,"%s %s %s\n",labels,labels,labels);
  free(labels);
/*******************************************************/

/* allocate memory to hold the hormonal time series.                          */
/* data is an Nx2 matrix, each row contains a recording of the time series.   */
/* The first column is the time the measurement was taken, the second column  */
/* records the log hormonal concentration (note: we are not using the obs. #) */

  data = (double **)calloc(*N,sizeof(double *));
  for (i=0;i<*N;i++)
    data[i] = (double *)calloc(2,sizeof(double));
/******************************************************************************/

/*************************************************************/
/* read the time series into memory                          */
/*************************************************************/
     for (i=0;i<*N;i++) {
    	fscanf(inf,"%d %lf %lf\n",&j,&data[i][0],&data[i][1]);
     data[i][1] = log(data[i][1]);
  }

  fclose(inf);
/************************************/

  return data;
}


