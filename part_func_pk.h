#ifndef VIENNA_RNA_PACKAGE_PART_FUNC_H
#define VIENNA_RNA_PACKAGE_PART_FUNC_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/centroid.h>
#include <ViennaRNA/equilibrium_probs.h>
#include <ViennaRNA/boltzmann_sampling.h>

/**
 *  @name Basic global partition function interface
 *  @{
 */

/**
 *  @brief Compute the partition function @f$Q@f$ for a given RNA sequence, or sequence alignment
 *
 *  If @a structure is not a NULL pointer on input, it contains on
 *  return a string consisting of the letters " . , | { } ( ) " denoting
 *  bases that are essentially unpaired, weakly paired, strongly paired without
 *  preference, weakly upstream (downstream) paired, or strongly up-
 *  (down-)stream paired bases, respectively.
 *  If the model's compute_bpp is set to 0 base pairing probabilities will not
 *  be computed (saving CPU time), otherwise after calculations took place #pr will
 *  contain the probability that bases @a i and @a j pair.
 * 
 *  @note This function is polymorphic. It accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE.
 *
 *  @note This function may return #INF / 100. in case of contradicting constraints
 *        or numerical over-/underflow. In the latter case, a corresponding warning
 *        will be issued to @p stdout.
 *
 *  @see #vrna_fold_compound_t, vrna_fold_compound(), vrna_pf_fold(), vrna_pf_circfold(),
 *        vrna_fold_compound_comparative(), vrna_pf_alifold(), vrna_pf_circalifold(),
 *        vrna_db_from_probs(), vrna_exp_params(), vrna_aln_pinfo()
 *
 *  @param[in,out]  vc              The fold compound data structure
 *  @param[in,out]  structure       A pointer to the character array where position-wise pairing propensity
 *                                  will be stored. (Maybe NULL)
 *  @return         The Gibbs free energy of the ensemble (@f$G = -RT \cdot \log(Q) @f$) in kcal/mol
 */
float vrna_pf_pk(vrna_fold_compound_t *vc, char *structure);

#endif
