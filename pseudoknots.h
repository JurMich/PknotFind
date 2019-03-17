#include <ViennaRNA/data_structures.h>

#include "energies.h"

#ifndef PSEUDOKNOTS_H_
#define PSEUDOKNOTS_H_

/* this matrix */
FLT_OR_DBL * qpk; /* matrix containing the values for pseudoknots */
FLT_OR_DBL PK_ENERGY; /* penalty for creating a pseudoknot */
int MIN_PK; /* minimum pseudoknot size */


/**
 * @brief computes enerky of all pseudoknots that can be created within
 * the interval [i,j] (it is exactly delinited by positions i an j)
 * and all possible objects that can contain (such as hairpins and other)
 * and stores them into DP matrix directly (so no return). You MUST have
 * computed the same for [i+1,j], [i,j-1] and [i+1,j-1] beforehand (DP
 * hopefully ensures that).
 * 
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * @param i: start of interval, included
 * @param j: send of interval, included 
 * 
 * **/

void exp_E_pseudoknot (	vrna_fold_compound_t *fcom,
						int i,
						int j);


/**
 * @brief an extension for external loop case of ViennaRNA library, it 
 * computes partition function over those cases that have a pseudoknot s
 * starting somewhere in external loop for interval [i,j]. 
 * 
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * @param i: start of interval, included
 * @param j: send of interval, included 
 * 
 * @return : a partition function of all pseudoknots (including 
 * containing elements) starting at external loop in [i,j].
 * 
 * **/ 
						
FLT_OR_DBL exp_E_pseudoknot_f (	vrna_fold_compound_t *fcom,
								int i,
								int j);


/**
 * @brief an extension for cases contained within a base pair, it returns
 * all pseudoknots that can be created within the base interval [i,j]. 
 * This function shouldn be called if base pair (i,j) can't be created.
 * 
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * @param i: start of interval, included
 * @param j: send of interval, included 
 * 
 * @return : a partition function of all pseudoknots (including 
 * containing elements) within a base pair (i,j).
 * 
 * **/ 						
						
FLT_OR_DBL exp_E_pseudoknot_c (	vrna_fold_compound_t *fcom,
								int i,
								int j);


/**
 * @brief an extension for cases contained within a multiloop (multiple
 * branches part), it returns a partition function of all cases where 
 * there is a pseudoknot in all but the last position in a multiloop.
 * 
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * @param i: start of interval, included
 * @param j: send of interval, included 
 * 
 * @return : a partition function of all pseudoknots (including 
 * containing elements) within all but the last branches of a multiloop.
 * 
 * **/ 	

FLT_OR_DBL exp_E_pseudoknot_m (	vrna_fold_compound_t *fcom,
								int i,
								int j);


/**
 * @brief an extension for cases contained within a multiloop (single
 * branch part), it returns a partition function of all cases where 
 * there is a pseudoknot in the last position in a multiloop.
 * 
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * @param i: start of interval, included
 * @param j: send of interval, included 
 * 
 * @return : a partition function of all pseudoknots (including 
 * containing elements) within the last branch of a multiloop.
 * 
 * **/
						
FLT_OR_DBL exp_E_pseudoknot_m1 (	vrna_fold_compound_t *fcom,
									int i,
									int j);


/**
 * @brief returns boltzmann factor of pseudoknots spanning over [i,j]
 * 
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * @param i: start of interval, included
 * @param j: send of interval, included 
 * 
 * @return : a partition function of all pseudoknots (including 
 * containing elements) over [i,j].
 * 
 * **/


FLT_OR_DBL exp_E_pseudoknot_new (	vrna_fold_compound_t *fcom,
								int i,
								int j); 

/**
 * @brief this function prepares DP pseudoknot matrix and also auxilliary 
 * grammar object (notably associating previous functions to their 
 * callbacks) for computations.
 * 
 * @param rna: an rna sequence in plain_sequence format read from fasta.
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * 
 * **/			
 																

void init_aux_grammar(	plain_sequence * rna,
						vrna_fold_compound_t * fc);
						

/**
 * @brief traces back a pseudoknot for an interval [i,j] (starting and
 * and ending at i and j respectively).
 * 
 * @param i: start of pseudoknot's first helix
 * @param j: end of pseudoknot's second helix
 * @param pstruc: secondary structure string where pseudoknot will be
 * drawn
 * @param fcom: fold compound object with all necessary matrices. Note 
 * that pk matrix is a global variable
 * 
 * **/

void traceback_pks(	int	i,
				int j,
				char* pstruc,
				vrna_fold_compound_t  *fcom,
				double *pk_energy);

#endif