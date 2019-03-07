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
	HELIX_PART_FCI = (helices_list***)malloc((rna->size+2)*sizeof(helices_list**));
	for(i = 0; i <= rna->size+1; i++){
		HELIX_PART_FCI[i] = (helices_list**)malloc((rna->size+2)*sizeof(helices_list*));
		for(j = 0; j <= rna->size+1; j++){
			HELIX_PART_FCI[i][j] = (helices_list*)malloc(sizeof(helices_list));
			HELIX_PART_FCI[i][j]->count = 0;
			HELIX_PART_FCI[i][j]->total_pf = 0;
			HELIX_PART_FCI[i][j]->first_hx = NULL;
		}	
	}
}

TYPE compute_goodness(int i, int j, int k, int l, TYPE pf){
	return 	log(pf)/(j-l+k-i+2);
}

/* builds a helix object */
helix* new_helix(int i, int j, int k, int l, TYPE helix_pf){
	helix *new_hx = (helix*)malloc(sizeof(helix));
	new_hx->k = k;
	new_hx->l = l;
	new_hx->helix_pf = helix_pf;
	new_hx->goodness = compute_goodness(i,j,k,l,helix_pf); /* adjust later */
	new_hx->next_hx = NULL;
	return new_hx;
}

/* insert the helix in a manner so they are ordered in decreasing goodness */
void insert_helix(int i, int j, helix *insert_hx, helix *last_inserted){
	helix	*current_hx, *tmp_hx;
	/* if last goodness is better, then continue from it, otherwise start from beginning */
	if(!last_inserted){ /* first insert */
		HELIX_PART_FCI[i][j]->first_hx = insert_hx;
		return;
		
	}else if(insert_hx->goodness > last_inserted->goodness){
		current_hx = HELIX_PART_FCI[i][j]->first_hx;
		
	}else{
		current_hx = last_inserted;	
		
	}
	if(HELIX_PART_FCI[i][j]->first_hx->goodness < insert_hx->goodness){
		tmp_hx = HELIX_PART_FCI[i][j]->first_hx ;
		HELIX_PART_FCI[i][j]->first_hx  = insert_hx;
		insert_hx->next_hx = tmp_hx;
		return;	
	}
	
	while(current_hx->next_hx && (current_hx->next_hx->goodness > insert_hx->goodness)){
		current_hx = current_hx->next_hx;
	}
	tmp_hx = current_hx->next_hx;
	current_hx->next_hx = insert_hx;
	insert_hx->next_hx = tmp_hx;
	HELIX_PART_FCI[i][j]->count++;
	return;
}


void reinsert_helix(int i, int j, helix* reinsert_hx, int position){
	helix	*prev_node1, *prev_node2, *tmp_hx;
	helix 	*curr_helix = HELIX_PART_FCI[i][j]->first_hx;
	int count = 0;
	
	/* get first preceeding node */
	if(curr_helix->goodness < reinsert_hx->goodness){
		count = 1;
		prev_node1 = NULL;	
	}
	else{
		while(curr_helix){
			count++;
			if(curr_helix->next_hx->goodness <  reinsert_hx->goodness){
				prev_node1 = curr_helix;
				tmp_hx = curr_helix->next_hx;
				break;	
			}
			curr_helix = curr_helix->next_hx;
		}
	}
	
	while(count < position){
		if(count == (position-1)){
			prev_node2 = curr_helix;
			break;
		}
		curr_helix = curr_helix->next_hx;
		count++;
	}
	
	prev_node2->next_hx = reinsert_hx->next_hx;
	if(prev_node1){
		tmp_hx = prev_node1->next_hx;
		prev_node1->next_hx = reinsert_hx;
		reinsert_hx->next_hx = tmp_hx;
	}else{
		tmp_hx = HELIX_PART_FCI[i][j]->first_hx;
		HELIX_PART_FCI[i][j]->first_hx = reinsert_hx;
		reinsert_hx->next_hx = tmp_hx;
	}	
}


/* adds partition function to an existing helix, otherwise it creates a new node in list with given helix */
helix* add_or_insert(int i, int j, int k, int l, TYPE helix_pf, helix *last_inserted){
	helix 		*new_hx = NULL;
	helix 		*curr_helix = HELIX_PART_FCI[i][j]->first_hx;
	int 		count = 0;
	TYPE		last_goodness = 0;
	while(curr_helix){
		count++;										  /* used later for reorder purposes */
		if((curr_helix->k == k) && (curr_helix->l == l)){ /* helix in list exists we continue */
			curr_helix->helix_pf += helix_pf;
			curr_helix->goodness = compute_goodness(i, j, k, l, curr_helix->helix_pf);
			if((curr_helix->goodness > last_goodness) && (count > 1)){
				reinsert_helix(i, j, curr_helix, count);
			}
			return curr_helix;
		}
		last_goodness = curr_helix->goodness;
		curr_helix = curr_helix->next_hx;
	}
	new_hx = new_helix(i, j, k, l, helix_pf);
	insert_helix(i, j, new_hx, last_inserted);
	return new_hx;	
}


/* removes extra helices */
void remove_overlimit_hx(int i, int j){
	if(HELIX_PART_FCI[i][j]->count > MAX_LIST){
		HELIX_PART_FCI[i][j]->count = MAX_LIST;
		int 	count = 0;
		helix 	*current_hx = HELIX_PART_FCI[i][j]->first_hx;
		helix 	*tmp_hx;
		while(count< MAX_LIST){ /* pass through the part of the list that stays */
			count++;
			current_hx = current_hx->next_hx;
		}
		tmp_hx = current_hx; /* reset the pointer of the last list */
		current_hx = current_hx->next_hx;
		tmp_hx->next_hx = NULL;
			
		while(current_hx){ /* erase the rest of the list */
			tmp_hx = current_hx;
			current_hx = current_hx->next_hx;
			free(tmp_hx);
		}
	}
}


/* computes the total partition function of HELIX_PART_FCI[i][j]*/
void get_total_pf(int i, int j){
	helix	*curr_helix = HELIX_PART_FCI[i][j]->first_hx;
	while(curr_helix){
		HELIX_PART_FCI[i][j]->total_pf += curr_helix->helix_pf;
		curr_helix = curr_helix->next_hx;	
	}
}


/* computes partition fci table for every quadruplet i,j,k,delta */
void compute_helix_partition(plain_sequence * rna){
	
	int 	i, j, k, l, dist, alpha, beta, delta;
	int 	beta_min, beta_max, delta_min, delta_max; /* pre-calculate limits of beta and delta respectively */
	TYPE 	part_fci;
	int 	theta = 3;
	helix	*new_hx, *curr_helix, *last_hx;
	
	initialize_helix_partition(rna);
	for(dist = theta; dist < rna->size; dist++){
		for(i = 1; i <= rna->size - dist; i++){
			j = i + dist;						
			if(pairable(i, j, rna)){ /* don't continue if i and j cannot be paired */
				last_hx = NULL;
				new_hx = new_helix(i, j, i, j, 1.);
				insert_helix(i,j, new_hx, last_hx);
				last_hx = new_hx;
				
				part_fci = 0;
				for(alpha = 0; alpha <= MAX_SKEW; alpha++){
					/* setting limits for beta, which depend on choices of alpha and delta */
					for(beta = 0; beta <= MAX_SKEW; beta++){
						if((j-beta-1 > i+alpha+1) && (j-beta-1 > 0) && (i+alpha+1 <= rna->size)){
							if((alpha == 0) && (beta == 0)){ /* stack */
								curr_helix = HELIX_PART_FCI[i+1][j-1]->first_hx;
								while(curr_helix){
									new_hx = add_or_insert(i, j, curr_helix->k, curr_helix->l,
									 energy_to_bfactor(stacking_energy(i,j)) * curr_helix->helix_pf, last_hx);
									if(new_hx){
										last_hx = new_hx;
									}			
									curr_helix = curr_helix->next_hx;
								}
							}
							else{ /* bulge or internal loop */
								curr_helix = HELIX_PART_FCI[i+alpha+1][j-beta-1]->first_hx;
								while(curr_helix){
									if(abs(j-curr_helix->l-(curr_helix->k-i))<=MAX_SKEW && 
									 abs(j-beta-curr_helix->l-(curr_helix->k-alpha-i))<=MAX_SKEW){ /* check if skew of i-k l-j is within limits */
										new_hx = add_or_insert(i, j, curr_helix->k, curr_helix->l,
										 energy_to_bfactor(internal_loop_energy(i, i+alpha+1,j-beta-1,j)) * curr_helix->helix_pf, last_hx);
										if(new_hx){
											last_hx = new_hx;
										}
									}			
									curr_helix = curr_helix->next_hx;
								}
							}
						}
					}
				}	
			}
			remove_overlimit_hx(i, j);
			get_total_pf(i, j);
		}
	}
}

