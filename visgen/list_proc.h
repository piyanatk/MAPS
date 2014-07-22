/*************************************************************************
 * list_proc.h                                                           *
 *                                                                       *
 * Function prototypes and data structures for list_proc.c, the set of   *
 * list processing primitives.                                           *
 *                                                                       *
 * Created 15 December 2010 by L. Benkevitch                             *
 *                                                                       *
 *************************************************************************/
#define TRUE 1
#define FALSE 0


struct list_int {
  int num;                  /* Integer value of an atom */
  struct list_int *next;    /* Pointer to next atom in list */
};

struct list_int *list_create(int n, int val[n]);
struct list_int *list_create_atom(int num);
struct list_int *list_copy(struct list_int *orig);
void list_insert(struct list_int *ptr, int num);
void list_del(struct list_int *list);
void list_print(char *str, struct list_int *list);
