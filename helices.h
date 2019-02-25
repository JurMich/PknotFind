#include "rna.h"
#include "energies.h"

#ifndef HELICES_H_
#define HELICES_H_

int MAX_SKEW;
int MAX_ALPHA;
int MAX_BETA;

/* Four-dimensional table containing value Z for helices at i,j,k, delta 
 * */
TYPE **** HELIX_PART_FCI;

void compute_helix_partition(plain_sequence * rna);


/* Table containing poset for every helix */

typedef struct poset_node poset_node;
typedef struct poset_children poset_children;


/* nodes of a poset */
struct poset_node{
					int k;
					int l;
					int status;  /* points out whether we already visited the node */
					poset_children * children;
					poset_children * last_child;
				};


/* wrapper that links nodes with children  */
struct poset_children{
					poset_children * next_child;
					poset_node * child;
				};

poset_node *** POSET_TAB;


#endif