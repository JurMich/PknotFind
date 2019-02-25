#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/loop_energies.h>

#include "rna.h"
#include "energies.h" 


/* sets ViennaRNA energy model */
void set_E_fold_model(plain_sequence * rna){	
	double				mfe;
	vrna_md_t			md;
	
	/* transforms sequence rna into type accepted by ViennaRNA */
	char *rna_seq = (char*) malloc(sizeof(char)*(rna->size+1));
	strncpy(rna_seq, rna->label+1, rna->size);
	rna_seq[rna->size] = '\0';
	
	/* creates model */
	set_model_details(&md);
	
	md.uniq_ML = 1;
	md.temperature = TEMP;  
	MIN_TURN = md.min_loop_size;
	
	E_fold_cp = vrna_fold_compound(rna_seq, &md, VRNA_OPTION_PF);
	
	mfe = vrna_mfe(E_fold_cp, NULL);
	vrna_exp_params_rescale(E_fold_cp, &mfe); /* rescale scale params */
}


/* this converts base into a number compatible with ViennaRNA's
 * vrna_md_t->pair/ptype definition*/
int base_to_int(char base){
  switch(base){
	case 'a': return 1;
	case 'A': return 1;
	case 'c': return 2;
	case 'C': return 2;
	case 'g': return 3;
	case 'G': return 3;
	case 'u': return 4;
	case 'U': return 4;
	default : return 0;  
  }	
}


/* obtains a type as defined by ViennaRNA of a pair entered in parameters.*/
int get_type(char base_5, char base_3){ /* base_5 - from 5' end, base_3 - from 3' end */
	int b_5, b_3, type;
	b_5 = base_to_int(base_5);
	b_3 = base_to_int(base_3);
	type = E_fold_cp->params->model_details.pair[b_5][b_3];
	return type;
}


/* defines whether base pairs can be paired or no */
int pairable(	int base_5, 
				int base_3, 
				plain_sequence * rna){
					
					
	/* not distant enough */
	if((base_3 - base_5 - 1) < TURN) 
		return 0;
		
	if(!get_type(rna->label[base_5], rna->label[base_3])){
		return 0;	
	}
	return 1;
}


TYPE energy_to_bfactor(double energy){
	return (exp(-energy/RT));
}


/* Energy of (i+1,j-1) stacking over (i,j), i+1<j-1 */
double stacking_energy(int i, int j){
	double stack_energy;
	stack_energy = vrna_E_stack(E_fold_cp, i, j);  /*ViennaRNA*/
	return stack_energy;
}


/* Energy of internal/bulge loop defined by two base pairs (i,l) and (j,k), i<k<l<j */
double internal_loop_energy(int i, int k, int l, int j){
	double it_energy;
	int 		ptype1, ptype2, size1, size2;
	short		*S, *S2;
	vrna_md_t	*md;
	
	md          = &(E_fold_cp->params->model_details);
	S			= E_fold_cp->sequence_encoding;
	S2			= E_fold_cp->sequence_encoding2;
  
	size1 = k-i-1;
	size2 = j-l-1;
	//ptype1 = get_type(rna->label[i], rna->label[j]);
	//ptype2 = get_type(rna->label[l], rna->label[k]);
	
	ptype1 = vrna_get_ptype_md(S2[i], S2[j], md);
    ptype2 = vrna_get_ptype_md(S2[l], S2[k], md);
	it_energy = E_IntLoop(size1, size2, ptype1, ptype2, S[i], S[j], S[k], S[l], E_fold_cp->params);
	return it_energy;
}