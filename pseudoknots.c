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
void exp_E_pseudoknot(	vrna_fold_compound_t *fcom,
						int i,
						int j ){
	
	int				*my_iindx;						
	int 			i_prime, j_prime;
	int				h11, h12, h21, h22; /* paired regions of helices */
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk1, qtemp;
	FLT_OR_DBL 		*q_m, *scale, *expMLbase;
	helix			*hx_H1, *hx_H2;

	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q_m 		= matrices->qm;
	scale 		= matrices->scale;
	expMLbase 	= matrices->expMLbase;
	
	qpk1 = 0.;
	
	/* get j' for H1 */
	
	if((j-i+1) > MIN_PK){
		for(j_prime = i+1; j_prime < j; j_prime++){ /* H1 - j' */
			hx_H1 = HELIX_PART_FCI[i][j_prime]->first_hx;
			while(hx_H1){
				/* get i' for H2 */
				for(i_prime = hx_H1->k+1; i_prime < hx_H1->l; i_prime++){ /* H2 - i' */
					hx_H2 = HELIX_PART_FCI[i_prime][j]->first_hx;
					while(hx_H2){
						//printf("H2: (%d, %d,%d,%d)", i_prime,hx_H2->k,hx_H2->l,i_prime);
						if((hx_H1->k < i_prime) && (hx_H2->k < hx_H1->l)&&(hx_H2->l > j_prime)){ /* H1 and H2 are compatible, proceed with pseudoknots */
							//printf("Hs (%d,%d,%d,%d) - (%d,%d, %d,%d)\n"); 
							qtemp = hx_H1->helix_pf*hx_H2->helix_pf;
							if(i_prime - hx_H1->k - 1 > 0){
								/* either there is some motif or only unpaired bases */
								qtemp *= (q_m[my_iindx[hx_H1->k + 1] - (i_prime - 1)] + expMLbase[i_prime - hx_H1->k - 1]); 
							}									
															
							if(hx_H1->l - hx_H2->k - 1 > 0){
								/* same as above */
								qtemp *= (q_m[my_iindx[hx_H2->k + 1] - (hx_H1->l - 1)] + expMLbase[hx_H2->k  -  hx_H1->l - 1]); 
							}
										
							if(hx_H2->l - j_prime - 1 > 0){
								/* same as two before */
								qtemp *= (q_m[my_iindx[j_prime + 1] - (hx_H2->l - 1)] + expMLbase[hx_H2->l - j_prime - 1]); 

							}
							/* scale values the segments 'consumed' by helices */
							h11 = hx_H1->k - i + 1;
							h12 = hx_H2->k - i_prime + 1;
							h21 = j_prime - hx_H1->l + 1;
							h22 = j - hx_H2->l + 1;		
										
							//printf("helices length %d \n", h11 + h12 + h21 + h22);				
										
							//printf("scale %f %d\n", scale[h11 + h12 + h21 + h22], h11 + h12 + h21 + h22);		 	
							qtemp *= scale[h11 + h12 + h21 + h22];
							qtemp *= exp(-PK_ENERGY/RT); /* penalty of H-type pseudoknot */
										
							/* sum up to other possibilities */
							qpk1 += qtemp;
						}
						hx_H2 = hx_H2->next_hx;
					}			
				}
				hx_H1 = hx_H1->next_hx;
			}
		}
	}
	qpk[my_iindx[i] - j] = qpk1; 
	//printf("qpk(%d,%d)=%f\n", i,j, qpk[my_iindx[i] - j]);
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
	printf("pseudoknot total: %f \n", qpk_extloop);
	
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
		pstruc[i - 1] = '<';
		pstruc[j - 1] = '>';
	}
	else if(type_closure == 4){
		pstruc[i - 1] = '<';
		pstruc[j - 1] = '>';
	}
						
}


/* traces back helices */
void traceback_helix(	int 	i, 
						int 	j, 
						helix* 	hx,
						char*	pstruc, 
						int 	type_closure,
						double 	*pk_energy
					){
											
	FLT_OR_DBL		r, qpk1, qtemp;
	int				i_prime, j_prime, alpha, beta;
	helix			*curr_helix;
	
	draw_pk_pair(i, j, pstruc, type_closure);
	r =  (((double)rand()) / ((double)RAND_MAX)) * hx->helix_pf;
	qpk1 = 0;
	
	if((i==hx->k) && (j==hx->l)){
		return;	
	}
	else{
		for(alpha = 0; alpha <= MAX_SKEW; alpha++){
			/* setting limits for beta, which depend on choices of alpha and delta */
			for(beta = 0; beta <= MAX_SKEW; beta++){
				if((alpha == 0) && (beta == 0)){ /* stack */
					curr_helix = HELIX_PART_FCI[i+1][j-1]->first_hx;
					while(curr_helix){
						if((curr_helix->k==hx->k) && (curr_helix->l==hx->l)){
							qpk1 += energy_to_bfactor(stacking_energy(i,j)) * curr_helix->helix_pf;
							i_prime = i+1;
							j_prime = j-1;
							if(qpk1 > r){
								printf("IL: %f\n", *pk_energy);
								*pk_energy *= energy_to_bfactor(stacking_energy(i,j));
								printf("stack: %f -> %f \n", *pk_energy, energy_to_bfactor(stacking_energy(i,j)));
								traceback_helix(i_prime, j_prime, curr_helix, pstruc, type_closure, pk_energy);
								return;	
							}
							
						}
						curr_helix = curr_helix->next_hx;			
					}
				}
				else{ /* bulge or internal loop */
					curr_helix = HELIX_PART_FCI[i+alpha+1][j-beta-1]->first_hx;
					while(curr_helix){
						if((curr_helix->k==hx->k) && (curr_helix->l==hx->l)){
							qpk1 += energy_to_bfactor(internal_loop_energy(i, i+alpha+1,j-beta-1,j)) * curr_helix->helix_pf;
							i_prime = i+alpha+1;
							j_prime = j-beta-1;
							if(qpk1 > r){
								printf("IL: %f\n", *pk_energy);
								*pk_energy *= energy_to_bfactor(internal_loop_energy(i, i+alpha+1,j-beta-1,j));
								printf("IL: %f -> %f \n", *pk_energy, energy_to_bfactor(internal_loop_energy(i, i+alpha+1,j-beta-1,j)));
								traceback_helix(i_prime, j_prime, curr_helix, pstruc, type_closure, pk_energy);	
								return;
							}
						}
						curr_helix = curr_helix->next_hx;
					}
				}
			}
		}		
	}
}


/* backtracks 'empty' parts inside pseudoknots */
void traceback_interspace(	int 		i,
							int 		j,
							char* 		pstruc,
							vrna_fold_compound_t  *fcom,
							double 		*pk_energy	
						){
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		r;
	FLT_OR_DBL 		*q_m, *expMLbase;
	int 			*my_iindx;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q_m 		= matrices->qm;
	expMLbase 	= matrices->expMLbase;
	
	if(j - i - 1 <= 0){
		if(j - i - 1 < 0)
			//printf("Negative space in pseudoknots: this shouldn't happen\n");
		return;
	}
	
	r = (((double)rand()) / ((double)RAND_MAX)) * (q_m[my_iindx[i+1] - (j-1)] + expMLbase[j - i - 1]);
	if(r > expMLbase[j - i - 1]){  /* do not leave the space empty */
		backtrack_qm(i + 1 , j - 1, pstruc, fcom, pk_energy);	
		//printf("backtracked \n");
	}else{
		*pk_energy *= expMLbase[j - i - 1];	
	}						
}
						

/* main traceback function, which picks up helices then arcs for them */
void traceback_pks(	int	i,
					int j,
					char* pstruc,
					vrna_fold_compound_t  *fcom,
					double 				*pk_energy
					){
	
	int 		i_prime, j_prime;
	int			h11, h12, h21, h22;
	int 		*my_iindx;
	helix		*hx_H1, *hx_H2;
	
	vrna_mx_pf_t 	* matrices;
	FLT_OR_DBL		qpk1, qtemp, r, pf_H1, pf_H2;
	FLT_OR_DBL 		*q_m, *scale, *expMLbase;
	
	/* extract data from fold compound */	
	my_iindx 	= fcom->iindx;
		
	matrices 	= fcom->exp_matrices;
	q_m 		= matrices->qm;
	scale 		= matrices->scale;
	expMLbase 	= matrices->expMLbase;
	
	qpk1 = 0.;
	 
	/* randomly pick two helices we vant */
	r = (((double)rand()) / ((double)RAND_MAX)) * qpk[my_iindx[i] - j];	
	
	if((j-i+1) > MIN_PK){
		for(j_prime = i+1; j_prime < j; j_prime++){ /* H1 - j' */
			hx_H1 = HELIX_PART_FCI[i][j_prime]->first_hx;
			while(hx_H1){
		
				/* get i' for H2 */
				
				for(i_prime = hx_H1->k+1; i_prime < hx_H1->l; i_prime++){ /* H2 - i' */
					hx_H2 = HELIX_PART_FCI[i_prime][j]->first_hx;
					while(hx_H2){
						
						if((hx_H1->k < i_prime) && (hx_H2->k < hx_H1->l)&&(hx_H2->l > j_prime)){ /* H1 and H2 are compatible, proceed with pseudoknots */
							pf_H1 = hx_H1->helix_pf;
							pf_H2 = hx_H2->helix_pf;
							qtemp = pf_H1 * pf_H2;
							
							//printf("%d, %d, %f, %f XX\n", i,j, hx_H1->helix_pf, hx_H2->helix_pf);
							if(i_prime - hx_H1->k - 1 > 0){
								/* either there is some motif or only unpaired bases */
								qtemp *= (q_m[my_iindx[hx_H1->k + 1] - (i_prime - 1)] + expMLbase[i_prime - hx_H1->k - 1]); 
							}									
											
							if(hx_H1->l - hx_H2->k - 1 > 0){
								/* same as above */
								qtemp *= (q_m[my_iindx[hx_H2->k + 1] - (hx_H1->l - 1)] + expMLbase[hx_H2->k  -  hx_H1->l - 1]); 
							}
										
							if(hx_H2->l - j_prime - 1 > 0){
								/* same as two before */
								qtemp *= (q_m[my_iindx[j_prime + 1] - (hx_H2->l - 1)] + expMLbase[hx_H2->l - j_prime - 1]); 

							}
							/* scale values the segments 'consumed' by helices */
							h11 = hx_H1->k - i + 1;
							h12 = hx_H2->k - i_prime + 1;
							h21 = j_prime - hx_H1->l + 1;
							h22 = j - hx_H2->l + 1;		
										
							//printf("helices length %d \n", h11 + h12 + h21 + h22);				
													
							qtemp *= scale[h11 + h12 + h21 + h22];
							qtemp *= exp(-PK_ENERGY/RT); /* penalty of H-type pseudoknot */
							
							//printf("tamper %f \n", qtemp);
										
							/* sum up to other possibilities */
							qpk1 += qtemp;
							if(qpk1 >= r){
								printf("pkfound %f \n", *pk_energy);
								*pk_energy *= exp(-PK_ENERGY/RT) * scale[h11 + h12 + h21 + h22];
								printf("pkfound -> cont %.9f %.09f %f %d %f\n", *pk_energy, exp(-PK_ENERGY/RT), scale[h11 + h12 + h21 + h22], h11 + h12 + h21 + h22, scale[1]);
								goto HELICES_CHOSEN;	/* move to selected helices */
							}
						}
						hx_H2 = hx_H2->next_hx;
					}			
				}
				hx_H1 = hx_H1->next_hx;
			}
		}
	}
	 
	HELICES_CHOSEN:
	
	/* traceback helices */
	
	//printf("(%d,%d,%d,%d) - (%d,%d,%d,%d)\n", i, j_prime, hx_H1->k, hx_H1->l, i_prime, j, hx_H2->k, hx_H2->l);
	traceback_helix(i,j_prime, hx_H1, pstruc, 1, pk_energy);
	traceback_helix(i_prime, j, hx_H2, pstruc, 2, pk_energy);
	
	/* traceback intermediate sectors */
	traceback_interspace(hx_H1->k, i_prime, pstruc, fcom, pk_energy);
	traceback_interspace(hx_H2->k, hx_H1->l, pstruc, fcom, pk_energy);
	traceback_interspace(j_prime, hx_H2->l, pstruc, fcom, pk_energy);
}

