/*******************************************************************/
/****************************hash.c ********************************/
/******************************************************************************/
/* this module implements the linked list.  in this implementation, the first */
/* node in the linked list is a dummy node.  one must initialize the list with*/
/* the dummy node as the first.                                               */
/******************************************************************************/
/*******************************************************************/

#include "pop_deconvolution_main.h"

/*******************************************************************
*******************GLOBAL VARIABLE DEFINITIONS**********************

 fitstart: The first time in hours that a pulse may occur

*********************************************************************/

extern double fitstart;

/********************************************************************/
/*SUBROUTINES THAT EXIST IN THIS PROGRAM

 *initialize_node
 *sequential_search
 insert_node
 delete_node
 print_list
 destroy_list

**********************************************************************/

/*********************************************************************/
               /*START OF initialize_node SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*initialize_node: this allocates memory for a new pulse and gives it some
                   initial parameters
    ARGUMENTS: None
    RETURNS: p, the created pulse  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *p: initialized pulse

 SUBROUTINES USED
  None  */
/***********************************************************************/

Subject_type *initialize_subject(void)
{
  int i;
  Subject_type *s;
  Node_type *initialize_node(void);

/* initialize memory for a subject */
  if ((s = (Subject_type *)malloc(sizeof(Subject_type))) == NULL) {
    printf("initialize_subject: out of memory, can't allocate node\n");
    exit(1);
  }
       
/* all nodes initialized with a pointer to the NULL pointer     */
/* and time = 0, the parameters all initialized to zero        */
/* it is up to the programmer to insure that the correct values */
/* are inserted into the node at the appropriate time           */

  s->succ = s->pred = NULL;
  s->list = initialize_node();
  s->decay = 0;
  for (i=0;i<2;i++)
    s->theta[i] = 1;
  for (i=0;i<2;i++)
    s->basehalf[i] = 0;

  return s;
}

Node_type *initialize_node(void)
{
  int i;
  Node_type *p;

/* initialize memory for a node */
    if ((p = (Node_type *)malloc(sizeof(Node_type))) == NULL) {
    printf("initialize_node: out of memory, can't allocate node\n");
    exit(1);
  }

/* all nodes initialized with a pointer to the NULL pointer     */
/* and time = 0, the parameters all initialized to zero        */
/* it is up to the programmer to insure that the correct values */
/* are inserted into the node at the appropriate time           */

  p->succ = p->pred = NULL;
  p->mean_contrib = NULL;
  p->time = fitstart;
    for (i=0;i<2;i++) {
        p->theta[i] = 1;
        p->kappa[i] = 1;
    }

  return p;
}

/*********************************************************************/
               /*START OF sequential_search SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*sequential_search: finds where a pulse's location fits within an established
                     linked list
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
               double time; this is the location of a new pulse; note that we
                  do not actually input a new pulse, we just input its location
    RETURNS: *loc, the newly created pulse's prececessor  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 *loc: the newly created pulse's predecessor
 *locm1: the newly created pulse's successor

 SUBROUTINES USED
  None  */
/***********************************************************************/

Node_type *sequential_search(Node_type *list,double time)
{
  Node_type *loc,*locm1;

/* searches through the linked list "list" for the first occassion of */
/* "list->time" that is less than "time".                             */
/* the new node will be inserted between "locm1" and "loc".  "locm1"  */
/* preceeds "loc" in "list" by one position                           */

  locm1 = list;
  for (loc=list;loc;loc=loc->succ){
    if (loc->time >= time)
      return locm1;
    else
      locm1 = loc;
  }
  return locm1;
}

/*********************************************************************/
               /*START OF insert_node SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*insert_node: integrates a newly created pulse into the linked list
    ARGUMENTS: Node_type *new_node; the newly created pulse;
               Node_type *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 *node: new_node's predecessor

 SUBROUTINES USED
  sequential_search: found in this file; identifies inputted pulse's predecessor  */
/***********************************************************************/

void insert_node(Node_type *new_node,Node_type *list)
{
  Node_type *node;

/* insert the node "new_node" in the linked list "list" */
/* find the position in which to insert "new_node"      */

  node = sequential_search(list,new_node->time);
/********************************************************/

/* reassign pointers so that everyone is pointing to the correct node */

  new_node->pred = node;
  new_node->succ = node->succ;
  if (node->succ != NULL)
    (node->succ)->pred = new_node;
  node->succ = new_node;
}


void insert_subject(Subject_type *new_subj,Subject_type *sublist)
{

/* insert the node "new_node" in the linked list "list" */
/* find the position in which to insert "new_node"      */

/********************************************************/

/* reassign pointers so that everyone is pointing to the correct node */
  if (sublist->succ != NULL)
    (sublist->succ)->pred = new_subj;
  new_subj->pred = sublist;
  new_subj->succ = sublist->succ;
  sublist->succ = new_subj;

}

/*********************************************************************/
               /*START OF delete_node SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*delete_node: frees memory associated with inputted pulse and reassigns
               pointers of remaining pulses
    ARGUMENTS: Node_type *node; the pulse to be deleted;
               Node_type *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
  None

 SUBROUTINES USED
  None  */
/***********************************************************************/

void delete_node(Node_type *node,Node_type *list)
{
  free(node->mean_contrib);
  if (node->succ == NULL)
    (node->pred)->succ = node->succ;
  else {
    (node->pred)->succ = node->succ;
    (node->succ)->pred = node->pred;
  }
  free(node);
}

/*********************************************************************/
               /*START OF print_list SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*print_list: prints information about all pulses in the linked list
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i: generic counter
 *node: counter through pulses

 SUBROUTINES USED
  None  */
/***********************************************************************/

void print_list(Node_type *list)
{
  int i;
  Node_type *node;

/* traverses the linked list and prints the contents of each node */
  i = 1;
  node = list->succ;
  while (node != NULL) {
    printf("%2d %8.4lf %8.4lf %8.4lf \n",i,node->time,node->theta[0],node->theta[1]);
    node = node->succ;
    i++;
  }
}

/*********************************************************************/
               /*START OF destroy_list SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*destroy_list: frees memory throughout the linked list
    ARGUMENTS: Node_type *list; this is the current list of pulses that exist;
    RETURNS: None; all updates are made internally  */
/*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 *loc: counter through pulses

 SUBROUTINES USED
  None  */
/***********************************************************************/

void destroy_list(Node_type *list)
{
  Node_type *loc;

  while (list->succ != NULL) {
    loc = list->succ;
    free(list->mean_contrib);
    free(list);
    list = loc;
  }
  free(list->mean_contrib);
  free(list);
}

void destroy_sublist(Subject_type *sublist)
{
  Subject_type *loc;
  void destroy_list(Node_type *);

  while (sublist->succ != NULL) {
    loc = sublist->succ;
    destroy_list(sublist->list);
    free(sublist);
    sublist = loc;
  }
  destroy_list(sublist->list);
  free(sublist);
}
