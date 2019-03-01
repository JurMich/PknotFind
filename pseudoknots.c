#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/utils/basic.h>

#include "rna.h"
#include "helices.h"
#include "energies.h"
#include "pseudoknots.h"
#include "boltzmann_sampling_pk.h"


/* #####################################################################
 * ###                                                               ###
 * ###             PARTITION FUNCTION COMPUTATIONS                   ###
 * ###                                                               ###
 * #####################################################################
 */

/* computes pseudoknot partition function for interval i,j */
void exp_E_pseudoknot (	vrna_fold_compound_t *fcom,
						int i,
						int j ){
	
	int				*my_iindx;						
	int 			i_prime, j_prime, k, l, k2, l2;
	int				delta1, delta1_min, delta1_max; 
	int				delta2, delta2_min, delta2_max; 
	int				h11, h12, h21, h22; /* paired regions of helices */
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk1, qtemp;
	FLT_OR_DBL 		*q_m, *scale, *expMLbase;

	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q_m 		= matrices->qm;
	scale 		= matrices->scale;
	expMLbase 	= matrices->expMLbase;
		
//	printf("%f %d \n", expMLbase[0], my_iindx[1]);
	
	qpk1 = 0.;
	
	/* get all points of helix H1 */
	
	
	for(j_prime = i+1; j_prime < j; j_prime++){ /* H1 right our */
		for(k = i; k < j_prime - 1; k++){ /* H1 left in */
			delta1_min = (i - k < -MAX_SKEW ? -MAX_SKEW : (i - k));
			delta1_max = (i + j_prime - 2*k > MAX_SKEW ? MAX_SKEW : (i + j_prime - 2*k));
			for(delta1 = delta1_min; delta1 <= delta1_max; delta1++){ /* H1 right in */
				l = i + j_prime - k - delta1;
				
				/* continue if something interesting exists for helix H1 */
				if(HELIX_PART_FCI[i][j_prime][k][delta1+MAX_SKEW]>0){
					
					qtemp =  HELIX_PART_FCI[i][j_prime][k][delta1];				
					
					/* now pick H2 so that it doesn't cover with H1 */
					for(i_prime = k + 1; i_prime < l; i_prime++){ /* H2 left out */
						for(k2 = i_prime; k2 < l; k2++){ /* H2 left in */
							delta2_min = (i_prime - k2 < -MAX_SKEW ? -MAX_SKEW : (i_prime - k2));
							delta2_max = (i_prime + j - j_prime - k2 - 1 > MAX_SKEW ? MAX_SKEW : (i_prime + j - k2 - j_prime - 1));
							for(delta2 = delta2_min; delta2 <= delta2_max; delta2++){ /* H1 right in */
								l2 = i_prime + j - k2 - delta2;
								
								/* continue if something interesting exists for helix H2 too */
								if(HELIX_PART_FCI[i_prime][j][k2][delta2+MAX_SKEW] > 0){
									//printf("(%d,%d,%d,%d)<->(%d,%d,%d,%d) %d,%d \n", i,j_prime,k,delta1,i_prime,j,k2,delta2, delta1+MAX_SKEW, delta2+MAX_SKEW);
									//printf("Ls (%d,%d)\n", l,l2);
									
									/* H1 x H2 */
									qtemp = 	HELIX_PART_FCI[i][j_prime][k][delta1+MAX_SKEW] * 
												HELIX_PART_FCI[i_prime][j][k2][delta2+MAX_SKEW];
												
									//printf("deltas %d %d (%d, %d) \n", delta1_min, delta1_max,  delta2_min, delta2_max);
									//~ printf("QTEMP: %f %f\ %f \n", HELIX_PART_FCI[i][j_prime][k][delta1+MAX_SKEW] * 
												//~ HELIX_PART_FCI[i_prime][j][k2][delta2+MAX_SKEW], 
												//~ HELIX_PART_FCI[i][j_prime][k][delta1+MAX_SKEW], 
												//~ HELIX_PART_FCI[i_prime][j][k2][delta2+MAX_SKEW]);
									
									/* all limits have been found, now include motifs for all
									 * three intermediate sectors (which can be unpaired) */
									
									//printf("maximum beforu %f, \n", qtemp);
									
									if(i_prime - k - 1 > 0){
										/* either there is some motif or only unpaired bases */
										qtemp *= (q_m[my_iindx[k + 1] - (i_prime - 1)] + expMLbase[i_prime - k - 1]); 
									}									
										
									if(l - k2 - 1 > 0){
										/* same as above */
										qtemp *= (q_m[my_iindx[k2 + 1] - (l - 1)] + expMLbase[l -  k2 - 1]); 
									}
									
									if(l2 - j_prime - 1 > 0){
										/* same as two before */
										qtemp *= (q_m[my_iindx[j_prime + 1] - (l2 - 1)] + expMLbase[l2 - j_prime - 1]); 

									}	
									
									//printf("maximum %f, \n", qtemp);
									/* scale values the segments 'consumed' by helices */
									h11 = k - i + 1;
									h12 = k2 - i_prime + 1;
									h21 = j_prime - l + 1;
									h22 = j - l2 + 1;		
									
									//printf("helices length %d \n", h11 + h12 + h21 + h22);				
								 	
								 	qtemp *= scale[h11 + h12 + h21 + h22];
								 	
								 	/* sum up to other possibilities */
								 	qpk1 += qtemp;
								 	qpk1 *= exp(-PK_ENERGY/RT); /* penalty of H-type pseudoknot */
								}
							} 	
						}
					}
				}
			}			
		}
	}
	
	/* contribution for intervals included in [i,j] (computed beforehand) and qpk1,
	 * this also includes limit cases 
	 * 
	 * NOTE: THIS IS INCORRECT and leads to ambiguity. Leaving it here as warning for
	 * future generations
	if(j - i < 3){ 
		qpki1j = 0.; 
		qpkij1 = 0.;
		qpki1j1 = 0.;
	}else{
		qpki1j = qpk[my_iindx[i+1] - j]; 
		qpkij1 = qpk[my_iindx[i] - (j-1)];
		qpki1j1 = qpk[my_iindx[i+1] - (j-1)];
	}
		
	//printf("qpk(%d, %d) %f %f %f %f \n ", i, j, qpki1j, qpkij1, qpki1j1, qpk1);
	qpk[my_iindx[i] - j] = qpki1j + qpkij1 - qpki1j1 + qpk1; */
	qpk[my_iindx[i] - j] = qpk1; 
}


FLT_OR_DBL exp_E_pseudoknot_new (	vrna_fold_compound_t *fcom,
								int i,
								int j){
									
	return qpk[fcom->iindx[i] - j];
}


/* computes partition function for a pseudoknot appearing in an external loop (f) */
FLT_OR_DBL exp_E_pseudoknot_f	(	vrna_fold_compound_t *fcom,
									int i,
									int j ){
	
	int				*my_iindx;						
	int 			k;
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk_extloop;
	FLT_OR_DBL 		*q;						
	
	qpk_extloop = 0.;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q 			= matrices->q;
	
	
	/* start computing - pseudoknot occupies entire region */
	//qpk_extloop += qpk[my_iindx[i] - j]; 
	
	for(k = i+1; k<=j; k++){
		//printf("pseudoknot F check: (%d,%d,%d) %f %f\n", i,j,k, qpk[my_iindx[i] - k], q[my_iindx[k] - j]);
		qpk_extloop += qpk[my_iindx[i] - (k - 1)] * q[my_iindx[k] - j];
	}						
	//printf("pseudoknot total: %f \n", qpk_extloop);
	
	return qpk_extloop;
}


/* computes a partition fci for a pseudoknot appearing between paired bases *
 * Also computes corresponding partition function for a given interval [i,j] */

FLT_OR_DBL exp_E_pseudoknot_c (	vrna_fold_compound_t *fcom,
								int i,
								int j){
	/* i,j is paired as external pair */
	
	
	int				*my_iindx;						
	int 			k, l, unp1, unp2;
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk_intloop, qtemp;
	FLT_OR_DBL 		*expMLbase;

	
	qpk_intloop = 0.;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
	
	matrices 	= fcom->exp_matrices;
	expMLbase 	= matrices->expMLbase;
	
	
	/* compute pseudoknots appearing hested in base pair (i,j) */
	for(k = i + 1; k < j - 1; k++){
		for(l = k + MIN_TURN; l < j; l++){
			/* we condsider pseudoknots as multiloops so this comes with penalties for unpaired bases */
			unp1 = k - i - 1;
			unp2 = j - l - 1;
			qtemp = expMLbase[unp1] * expMLbase[unp2]; 
			qtemp *= qpk[my_iindx[k] - l];
			
			//printf("Pseudoknot C check (%d,%d,%d,%d): %f, unp: %f \n", i,j,k,l, qpk[my_iindx[k] - l], expMLbase[unp1]);
			
			qpk_intloop += qtemp;
		} 
	}
	//printf("C: (%d,%d) %f\n", i,j, qpk_intloop);
	
	return qpk_intloop;						
}


/* computes partition function for a pseudoknot appearing in multilooped region (1+ helices) */
FLT_OR_DBL exp_E_pseudoknot_m (	vrna_fold_compound_t *fcom,
								int i,
								int j){
									
	int				*my_iindx;						
	int 			k, unp1;
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk_multiloop, qtemp;
	FLT_OR_DBL 		*q_m, *expMLbase;
	
	qpk_multiloop = 0.;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
	
	matrices 	= fcom->exp_matrices;
	expMLbase 	= matrices->expMLbase;
	q_m 		= matrices->qm;
	
	//qpk_multiloop = qpk[my_iindx[i] - j];
	
	for(k = i+1; k < j; k++){
		/* case A where i - k-1 is unpaired */
		unp1 = k - i;
		qtemp = expMLbase[unp1] * qpk[my_iindx[k] - j];
		
		/* case B where i - k-1 is another instance of qm */
		qtemp += q_m[my_iindx[i] - (k - 1)] * qpk[my_iindx[k] - j];
		
		qpk_multiloop += qtemp;
		

	}
	
	return qpk_multiloop;									
}


/* computes partition function for a pseudoknot appearing in multilooped region (1 helix) */					
FLT_OR_DBL exp_E_pseudoknot_m1 (	vrna_fold_compound_t *fcom,
									int i,
									int j){
	int				*my_iindx;						
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
	
	return qpk[my_iindx[i] - j];								
}


/* initiates auxilliary grammar and everything related.
 * Notably creates DP matrix for pseudoknots and returns it.*/
void init_aux_grammar (	plain_sequence * rna,
						vrna_fold_compound_t * fcom
					 ){

	/* init pseudoknot matrix */
	int 				size, n, i;
	vrna_gr_aux_t 		*aux_grammar; /* auxiliary pseudoknot grammar */

	n = rna->size;
	qpk = NULL;
	size = ((n + 1) * (n + 2))/2;
	
	qpk = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);	
	for(i = 0; i < size; i++)
		qpk[i] = 0.;
		
	/* create aux grammar and add rules to it */
	fcom->aux_grammar = (vrna_gr_aux_t *)vrna_alloc(sizeof(vrna_gr_aux_t));
	//fcom->aux_grammar->cb_aux_exp_f = exp_E_pseudoknot_f;
	fcom->aux_grammar->cb_aux_exp_c = exp_E_pseudoknot_c;
	//fcom->aux_grammar->cb_aux_exp_m = exp_E_pseudoknot_m;
	//fcom->aux_grammar->cb_aux_exp_m1 = exp_E_pseudoknot_m1;
			
	//fcom->aux_grammar->cb_aux_exp_new = exp_E_pseudoknot_new;
}


/* #####################################################################
 * ###                                                               ###
 * ###                  TRACEBACK COMPUTATIONS                       ###
 * ###                                                               ###
 * #####################################################################
 */

/* draws pair to structure according to choices */
void draw_pk_pair(	int 	i,
					int		j,
					char*	pstruc,
					int		type_closure){
	if(type_closure == 1){
		pstruc[i - 1] = '[';
		pstruc[j - 1] = ']';	
	}
	else if(type_closure == 2){
		pstruc[i - 1] = '{';
		pstruc[j - 1] = '}';
	}
	else if(type_closure == 3){
		pstruc[i - 1] = '{';
		pstruc[j - 1] = '}';
	}
	else if(type_closure == 4){
		pstruc[i - 1] = '<';
		pstruc[j - 1] = '>';
	}
						
}


/* function that recursively traces the arcs for helix that:
 * - has (i,j) as opening arc
 * - has (k,l) as closing arc 
 * - also, int_closure indicates symbol to indicate pseudoknots */
void traceback_helix_rec(	int i,
							int j,
							int k,
							int l,
							char* pstruc,
							int type_closure){
	
	draw_pk_pair(i, j, pstruc, type_closure);
	if(i == k){		/* if i == k then obligatorily j == l since there can't be base triplet! */
		//printf("Arc founded (%d,%d)\n", i,j);
		return;	
	}
	
	int 			alpha, beta, delta;
	int 			beta_min, beta_max; 
	plain_sequence *rna;
	FLT_OR_DBL 		r, qpk1;

	delta = j - l - k + i;
	qpk1 = 0.;
	r = ((double)rand()) / ((double)RAND_MAX) * HELIX_PART_FCI[i][j][k][delta + MAX_SKEW];
	
	for(alpha = 0; alpha <= MAX_SKEW; alpha++){
		/* setting limits for beta, which depend on choices of alpha and delta */
		beta_min = ((delta + alpha - MAX_SKEW >= 0) ? (delta + alpha - MAX_SKEW) : 0);
		beta_max = ((delta + alpha + MAX_SKEW <= MAX_SKEW) ? (delta + alpha + MAX_SKEW) : MAX_SKEW);
		for(beta = beta_min; beta <= beta_max; beta++){
			if((alpha == 0) && (beta == 0)){
				qpk1 += energy_to_bfactor(stacking_energy(i,j)) * HELIX_PART_FCI[i+1][j-1][k][delta+MAX_SKEW];	
			}
			else if((i + alpha + 1) < (j - beta - 1)){
				qpk1 += energy_to_bfactor(internal_loop_energy(i, i+alpha+1,j-beta-1,j))
				 * HELIX_PART_FCI[i+alpha+1][j-beta-1][k][delta - beta + alpha + MAX_SKEW];
			}
			if(qpk1 > r){
				goto ARC_FOUND;	
			}
		}	
	}	
	ARC_FOUND:
	
	//printf("Arc (%d,%d)\n", i,j);	
	
	traceback_helix_rec(i + alpha + 1, j - beta - 1, k, l, pstruc, type_closure);	
	
}


/* backtracks 'empty' parts inside pseudoknots */
void traceback_interspace(	int 		i,
							int 		j,
							char* 		pstruc,
							vrna_fold_compound_t  *fcom	
						){
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		r, p_unp;
	FLT_OR_DBL 		*q_m, *expMLbase;
	int 			*my_iindx;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q_m 		= matrices->qm;
	expMLbase 	= matrices->expMLbase;
	
	p_unp = expMLbase[j - i - 1];
	
	//printf("total %f unp %f vs. something %f \n", q_m[my_iindx[i+1] - (j-1)] + p_unp, p_unp, q_m[my_iindx[i+1] - (j-1)]);	
	
	if(j - i - 1 <= 0){
		if(j - i - 1 < 0)
			//printf("Negative space in pseudoknots: this shouldn't happen\n");
		return;
	}
	
	r = (((double)rand()) / ((double)RAND_MAX)) * (q_m[my_iindx[i+1] - (j-1)] + p_unp);
	if(r > p_unp){  /* do not leave the space empty */
		backtrack_qm(i + 1 , j - 1, pstruc, fcom);	
		//printf("backtracked \n");
	}						
}
							

/* main traceback function, which picks up helices then arcs for them */
void traceback_pks(	int	i,
					int j,
					char* pstruc,
					vrna_fold_compound_t  *fcom
					){
	
	int 		i1, j1, k1, l1, i2, j2, k2, l2;
	int			delta1, delta1_min, delta1_max; 
	int			delta2, delta2_min, delta2_max; 
	int			h11, h12, h21, h22;
	int 		*my_iindx;
	
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk1, qtemp, r;
	FLT_OR_DBL 		*q_m, *scale, *expMLbase;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q_m 		= matrices->qm;
	scale 		= matrices->scale;
	expMLbase 	= matrices->expMLbase;
	
	i1 = i; /* i is the start of H1 */
	j2 = j; /* j is the end of H2 */
	qpk1 = 0.;
	
	/*
	 * ### STEP A: CHOOSE HELICES H1 AND H2 ###
	 */
	 
	/* randomly pick two helices we vant */
	r = (((double)rand()) / ((double)RAND_MAX)) * qpk[my_iindx[i] - j];	
	//printf("random %f %f\n", r, qpk[my_iindx[i] - j]);
	
	for(j1 = i1+1; j1 < j2; j1++){ /* H1 right out */
		for(k1 = i1; k1 < j1 - 1; k1++){ /* H1 left in */
			delta1_min = (i1 - k1 < -MAX_SKEW ? -MAX_SKEW : (i1 - k1));
			delta1_max = (i1 + j1 - 2*k1 > MAX_SKEW ? MAX_SKEW : (i1 + j1 - 2*k1));
			for(delta1 = delta1_min; delta1 <= delta1_max; delta1++){ /* H1 right in */
				l1 = i1 + j1 - k1 - delta1;
				
				/* continue if something interesting exists for helix H1 */
				if(HELIX_PART_FCI[i1][j1][k1][delta1+MAX_SKEW]>0){
					
					qtemp =  HELIX_PART_FCI[i1][j1][k1][delta1];				
					
					/* now pick H2 so that it doesn't cover with H1 */
					for(i2 = k1 + 1; i2 < l1; i2++){ /* H2 left out */
						for(k2 = i2; k2 < l1; k2++){ /* H2 left in */
							delta2_min = (i2 - k2 < -MAX_SKEW ? -MAX_SKEW : (i2 - k2));
							delta2_max = (i2 + j2 - j1 - k2 - 1 > MAX_SKEW ? MAX_SKEW : (i2 + j2 - k2 - j1 - 1));
							for(delta2 = delta2_min; delta2 <= delta2_max; delta2++){ /* H1 right in */
								l2 = i2 + j2 - k2 - delta2;
								
								/* continue if something interesting exists for helix H2 too */
								if(HELIX_PART_FCI[i2][j2][k2][delta2+MAX_SKEW] > 0){
									
									
									/* H1 x H2 */
									qtemp = 	HELIX_PART_FCI[i1][j1][k1][delta1+MAX_SKEW] * 
												HELIX_PART_FCI[i2][j2][k2][delta2+MAX_SKEW];
												
									/* all limits have been found, now include motifs for all
									 * three intermediate sectors (which can be unpaired) */
									
									
									if(i2 - k1 - 1 > 0){
										/* either there is some motif or only unpaired bases */
										qtemp *= (q_m[my_iindx[k1 + 1] - (i2 - 1)] + expMLbase[i2 - k1 - 1]); 
									}									
										
									if(l1 - k2 - 1 > 0){
										/* same as above */
										qtemp *= (q_m[my_iindx[k2 + 1] - (l1 - 1)] + expMLbase[l1 -  k2 - 1]); 
									}
									
									if(l2 - j1 - 1 > 0){
										/* same as two before */
										qtemp *= (q_m[my_iindx[j1 + 1] - (l2 - 1)] + expMLbase[l2 - j1 - 1]); 

									}	
									
									/* scale values the segments 'consumed' by helices */
									h11 = k1 - i1 + 1;
									h12 = k2 - i2 + 1;
									h21 = j1 - l1 + 1;
									h22 = j2 - l2 + 1;		
									
									//printf("helices length %d \n", h11 + h12 + h21 + h22);				
								 	
								 	qtemp *= scale[h11 + h12 + h21 + h22];
								 	
								 	/* sum up to other possibilities */
								 	qpk1 += qtemp;
								 	qpk1 *= PK_ENERGY;
								 	
								 	//printf("H1:((%d,%d) - (%d,%d)),  H2:((%d,%d) - (%d,%d)), r: %f, total: %f, this segment: %f, \n", i1, j1, k1, l1, i2, j2, k2, l2, r, qpk1, qtemp);
								 	
								 	/* helices chosen, get out of loops*/
								 	if(qpk1 >= r){
										goto HELICES_CHOSEN;	
									}
								}
							} 	
						}
					}
				}
			}			
		}
	}
	
	HELICES_CHOSEN: 
	//printf("test ok\n");	
	
	//printf("chosen helices are H1:((%d,%d) - (%d,%d) - delta: %d),  H2:((%d,%d) - (%d,%d) - delta: %d)\n", i1, j1, k1, l1, delta1, i2, j2, k2, l2, delta2);		
	
	/*
	 * ### STEP B: CHOOSE ARCS IN HELCIES H1 AND H2 ###
	 */		
	 
	/* B1: Choose arcs for H1 */
	traceback_helix_rec(i1, j1, k1, l1, pstruc, 1);	

	/* B2: Choose arcs for H2 */
	traceback_helix_rec(i2, j2, k2, l2, pstruc, 2);	
	
	/* traceback everything else */
	traceback_interspace(k1, i2, pstruc, fcom);
	traceback_interspace(k2, l1, pstruc, fcom);
	traceback_interspace(j1, l2, pstruc, fcom);
}
