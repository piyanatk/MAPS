#include <stdio.h>
#include <stdlib.h>
#include "list_proc.h"

/*
 * Create list closed into ring with values from val[n] array.
 */

struct list_int *list_create(/* nelem,  */
			     int n, int val[n]) 
{
  struct list_int *list, *ptr;
  int i;
  if (n == 0) return NULL;
  list = (struct list_int *) malloc(sizeof(struct list_int));
  if (list == NULL) {printf("list_create: malloc error\n"); exit(1);}
  list->num = val[0];
  ptr = list;
  for (i = 1; i < n; i++) {
    ptr->next = (struct list_int *) malloc(sizeof(struct list_int));
    if (ptr->next == NULL) {printf("list_create: malloc error\n"); exit(1);}
    ptr = ptr->next;
    ptr->num = val[i];
  }
  ptr->next = list;
  return list;
}

/*
 * Create one-atom list closed into ring with the value from val.
 */

struct list_int *list_create_atom(/* num  */
			     int val) 
{
  struct list_int *list;
  list = (struct list_int *) malloc(sizeof(struct list_int));
  if (list == NULL) {printf("list_create_atom: malloc error\n"); exit(1);}
  list->num = val;
  list->next = list;
  return list;
}


/*
 * List duplication.
 * The returned list is closed into ring, i.e.
 * its last element points at its first one.
 */
struct list_int *list_copy(/* *orig */
			   struct list_int *orig) 
{
  /* int i = 0; */
  struct list_int *copy, *porig, *pcopy;
  
  if (orig == NULL) return NULL;

  copy = (struct list_int *) malloc(sizeof(struct list_int));
  if (copy == NULL) {printf("list_copy: malloc error\n"); exit(1);}
  porig = orig;
  pcopy = copy;

  while (porig->next != orig) {  /* Use this to only copy CLOSED lists! */
    pcopy->num = porig->num;
    porig = porig->next;
    pcopy->next = (struct list_int *) malloc(sizeof(struct list_int));
    if (pcopy->next == NULL) {printf("list_copy: malloc error\n"); exit(1);}
    pcopy = pcopy->next;
    /*printf("i = %d, porig = %p, pcopy = %p\n", i, porig, pcopy); */
    /* i++; */
  } while (porig->next != orig);  /* Use this to only copy CLOSED lists! */

  pcopy->num = porig->num;
  pcopy->next = copy;  /* Close the list into ring */
  return copy;
}

/*
 * Insert atom with integer value num after ptr
 */
void list_insert(struct list_int *ptr, int num) {
  struct list_int *atom;
  //next = ptr->next;
  atom = (struct list_int *) malloc(sizeof(struct list_int));
  if (atom == NULL) {printf("list_insert: malloc error\n"); exit(1);}
  atom->num = num;
  atom->next = ptr->next;
  ptr->next = atom;
}

/*
 * Delete list and free the memory.
 * The list MUST be closed into ring!
 */
void list_del(/* *list */
			   struct list_int *list) 
{
  struct list_int *plst, *next;
  plst = list;
  next = NULL;

  while (next != list) {
    next = plst->next;
    free(plst);
    plst = next;
  }
}

/*
 * Print list after str.
 */
void list_print(/* str, *list */
		char *str,        
		struct list_int *list) 
{
  struct list_int *ptr;
  ptr = list;
  printf("%s", str);
  if (list == NULL) {
    printf("nil\n");
    return;
  }
  do {
    printf("%d ", ptr->num);
    ptr = ptr->next;
  } while (ptr != list);
  printf("\n");
}
