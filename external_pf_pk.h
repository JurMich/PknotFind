#ifndef VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_H
#define VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/params/basic.h>

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

/**
 *  @brief  Evaluate a stem branching off the exterior loop (Boltzmann factor version)
 *
 *  Given a base pair @f$(i,j)@f$ encoded by @em type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *  If either of the adjacent nucleotides @f$(i - 1)@f$ and @f$(j+1)@f$ must not
 *  contribute stacking energy, the corresponding encoding must be @f$-1@f$.
 *
 *  @see vrna_E_ext_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters (Boltzmann factor version)
 *  @return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
FLT_OR_DBL
vrna_exp_E_ext_stem_pk(unsigned int      type,
                    int               n5d,
                    int               n3d,
                    vrna_exp_param_t  *p);


struct vrna_mx_pf_aux_el_s *
vrna_exp_E_ext_fast_init_pk(vrna_fold_compound_t *fc);


void
vrna_exp_E_ext_fast_rotate_pk(struct vrna_mx_pf_aux_el_s *aux_mx);


void
vrna_exp_E_ext_fast_free_pk(struct vrna_mx_pf_aux_el_s *aux_mx);


FLT_OR_DBL
vrna_exp_E_ext_fast_pk(vrna_fold_compound_t        *fc,
                    int                         i,
                    int                         j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx);


void
vrna_exp_E_ext_fast_update_pk(vrna_fold_compound_t       *fc,
                           int                        j,
                           struct vrna_mx_pf_aux_el_s *aux_mx);

#endif
