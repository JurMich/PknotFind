#include "rna.h"
#include "energies.h"

#ifndef HELICES_H_
#define HELICES_H_

int MAX_SKEW;
int MAX_ALPHA;
int MAX_BETA;
int MAX_LIST; /* maximum length of a list */

void compute_helix_partition(plain_sequence * rna);

/* nodes of the helix list */

typedef struct helices_list helices_list;
typedef struct helix helix;

/* table of linked lists of all helices */
helices_list *** HELIX_PART_FCI;

struct helices_list{
					int count;
					TYPE total_pf;
					helix *first_hx;
};

struct helix{
			int k;
			int l;
			TYPE helix_pf;
			TYPE goodness; /* estimates quality of the helix */
			helix *next_hx; /* next helix in linked list */
};
#endif