#ifndef VIENNA_RNA_PACKAGE_LOOPS_MULTIBRANCH_H
#define VIENNA_RNA_PACKAGE_LOOPS_MULTIBRANCH_H

#include <ViennaRNA/utils/basic.h>
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

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/* End basic interface */
/**@}*/


/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */


/**
 *  @brief  Auxiliary helper arrays for fast exterior loop computations
 *
 *  @see vrna_exp_E_ml_fast_init(), vrna_exp_E_ml_fast_rotate(),
 *  vrna_exp_E_ml_fast_free(), vrna_exp_E_ml_fast()
 */
typedef struct vrna_mx_pf_aux_ml_s *vrna_mx_pf_aux_ml_t;


FLT_OR_DBL
vrna_exp_E_mb_loop_fast(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j,
                        vrna_mx_pf_aux_ml_t   aux_mx);


vrna_mx_pf_aux_ml_t
vrna_exp_E_ml_fast_init_pk(vrna_fold_compound_t *fc);


void
vrna_exp_E_ml_fast_rotate_pk(vrna_mx_pf_aux_ml_t aux_mx);


void
vrna_exp_E_ml_fast_free_pk(vrna_mx_pf_aux_ml_t aux_mx);


const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm_pk(struct vrna_mx_pf_aux_ml_s *aux_mx);


const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm1_pk(struct vrna_mx_pf_aux_ml_s *aux_mx);


FLT_OR_DBL
vrna_exp_E_ml_fast_pk(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j,
                   vrna_mx_pf_aux_ml_t  aux_mx);


/* End partition function interface */



/**
 * @}
 */

#endif
