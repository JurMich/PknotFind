#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"
#include "helices.h"
#include "energies.h"

/* ###############################
 * #   Helix table constructor   #
 * ###############################
 */


/* initializes partition fci table */
void initialize_helix_partition(plain_sequence * rna){
	
	int 		i, j, k, delta;
	printf("%d \n ", rna->size+2);
	HELIX_PART_FCI = (TYPE****)malloc((rna->size+2)*sizeof(TYPE***));
	for(i = 0; i <= rna->size; i++){
		HELIX_PART_FCI[i] = (TYPE***)malloc((rna->size+2)*sizeof(TYPE**));
		for(j = 0; j <= rna->size; j++){
			HELIX_PART_FCI[i][j] = (TYPE**)malloc((rna->size+2)*sizeof(TYPE*));
			for(k = 0; k <= rna->size; k++){
				HELIX_PART_FCI[i][j][k] = (TYPE*)malloc((2*MAX_SKEW+1)*sizeof(TYPE));
				for(delta = 0; delta<=2*MAX_SKEW; delta++){
					HELIX_PART_FCI[i][j][k][delta] = 0.;
				}	
			}
		}	
	}
}


/* computes partition fci table for every quadruplet i,j,k,delta */
void compute_helix_partition(plain_sequence * rna){
	
	int 	i, j, k, l, dist, alpha, beta, delta;
	int 	beta_min, beta_max, delta_min, delta_max; /* pre-calculate limits of beta and delta respectively */
	TYPE 	part_fci;
	int 	theta = 3;
	
	initialize_helix_partition(rna);
	for(dist = theta; dist < rna->size; dist++){
		for(i = 1; i <= rna->size - dist; i++){
			j = i + dist;
			if(pairable(i, j, rna)){ /* don't continue if i and j cannot be paired */
				HELIX_PART_FCI[i][j][i][MAX_SKEW] = 1.;
				for(k = i+1; k <= j-1; k++){
					delta_min = (i - k < -MAX_SKEW ? -MAX_SKEW : (i - k));
					delta_max = (i + j - 2*k > MAX_SKEW ? MAX_SKEW : (i + j - 2*k));
					
					for(delta = delta_min; delta <= delta_max; delta++){
						l = i + j - k - delta;	
						if(pairable(k, l, rna)){
							part_fci = 0;
							for(alpha = 0; alpha <= MAX_SKEW; alpha++){
								/* setting limits for beta, which depend on choices of alpha and delta */
								beta_min = ((delta + alpha - MAX_SKEW >= 0) ? (delta + alpha - MAX_SKEW) : 0);
								beta_max = ((delta + alpha + MAX_SKEW <= MAX_SKEW) ? (delta + alpha + MAX_SKEW) : MAX_SKEW);
								for(beta = beta_min; beta <= beta_max; beta++){
									if((alpha == 0) && (beta == 0)){
										part_fci += energy_to_bfactor(stacking_energy(i,j)) * HELIX_PART_FCI[i+1][j-1][k][delta+MAX_SKEW];	
									}
									else{
										part_fci += energy_to_bfactor(internal_loop_energy(i, i+alpha+1,j-beta-1,j))
										 * HELIX_PART_FCI[i+alpha+1][j-beta-1][k][delta - beta + alpha + MAX_SKEW];
									}
								}	
							}
							HELIX_PART_FCI[i][j][k][delta+MAX_SKEW] += part_fci;
						}
					}
				}	
			}
		}
	}	
}


/* ###############################
 * #   Helix poset constructor   #
 * ###############################
 */

 
/* creates a poset node */ 
poset_node * construct_poset_node(int k, int l){
	
	poset_node * new_node = (poset_node*)malloc(sizeof(poset_node));
	new_node->k = k;
	new_node->l = l;
	new_node->status = 1;
	new_node->children = NULL;
	new_node->last_child = NULL;
	return new_node;
}


/* builds an initial table with clique posets */
void initiate_poset_table(plain_sequence * rna){
	
	int 			i, j, k, delta;
	poset_node 		*node_tmp;
	
	POSET_TAB = (poset_node***)malloc((rna->size+2)*sizeof(poset_node**));
	for(i = 1; i <= rna->size; i++){
		POSET_TAB[i] = (poset_node**)malloc((rna->size+2)*sizeof(poset_node*));
		for(j = 1+1; j <= rna->size; j++){
			if(HELIX_PART_FCI[i][j][i][0]==1.0){
				POSET_TAB[i][j] = construct_poset_node(i,j);
				for(k = 0; k < rna->size; k++){
					for	(delta = -MAX_SKEW; delta <= MAX_SKEW; delta++){
							
					}
				}
			}	
		}
	}
}
