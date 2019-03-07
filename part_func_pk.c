/*
 *   				THIS FILE IS TAKEN FROM ViennaRNA AND
 * 					IT IS DISTRIBUTED UNDER ITS LICENSE!
 * 			
 *                partiton function for RNA secondary structures
 * 
 *                CREDITED TO:
 * 
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 * 
 *				BE SURE TO VISIT: https://www.tbi.univie.ac.at/RNA/
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* #defines FLT_MAX ... */
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe.h"
#include "external_pf_pk.h"
#include "multibranch_pf_pk.h"
#include "part_func_pk.h"
#include "pseudoknots.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
 
PRIVATE int
fill_arrays_pk(vrna_fold_compound_t *fc);

/*
 #################################
 # OTHER DECLARATIONS            #
 #################################
 */

const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm_pk(struct vrna_mx_pf_aux_ml_s *aux_mx);

FLT_OR_DBL
vrna_exp_E_ext_fast_pk(vrna_fold_compound_t        *fc,
                    int                         i,
                    int                         j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx);

FLT_OR_DBL
vrna_exp_E_ml_fast_pk(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j,
                   vrna_mx_pf_aux_ml_t  aux_mx);
                   
void
vrna_exp_E_ml_fast_rotate_pk(vrna_mx_pf_aux_ml_t aux_mx);


void
vrna_exp_E_ml_fast_free_pk(vrna_mx_pf_aux_ml_t aux_mx);


void
vrna_exp_E_ext_fast_rotate_pk(struct vrna_mx_pf_aux_el_s *aux_mx);


void
vrna_exp_E_ext_fast_free_pk(struct vrna_mx_pf_aux_el_s *aux_mx);



/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_pf_pk(vrna_fold_compound_t  *fc,
        char                  *structure)
{
  int               n;
  FLT_OR_DBL        Q;
  double            free_energy;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  free_energy = (float)(INF / 100.);

  if (fc) {
    /* make sure, everything is set up properly to start partition function computations */
    if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_PF)) {
      vrna_message_warning("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
      return free_energy;
    }

    n         = fc->length;
    params    = fc->exp_params;
    matrices  = fc->exp_matrices;
    md        = &(params->model_details);

#ifdef _OPENMP
    /* Explicitly turn off dynamic threads */
    omp_set_dynamic(0);
#endif

#ifdef SUN4
    nonstandard_arithmetic();
#elif defined(HP9)
    fpsetfastmode(1);
#endif

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_PF_PRE, fc->auxdata);

    /* call user-defined grammar pre-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_PRE, fc->aux_grammar->data);

    if (!fill_arrays_pk(fc)) {
#ifdef SUN4
      standard_arithmetic();
#elif defined(HP9)
      fpsetfastmode(0);
#endif
      return (float)(INF / 100.);
    }

    /* calculate base pairing probability matrix (bppm)  */
    if (md->compute_bpp) {
      vrna_pairing_probs(fc, structure);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

      /*
       *  Backward compatibility:
       *  This block may be removed if deprecated functions
       *  relying on the global variable "pr" vanish from within the package!
       */
      pr = matrices->probs;

#endif
    }

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_PF_POST, fc->auxdata);

    /* call user-defined grammar post-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_POST, fc->aux_grammar->data);

    switch (md->backtrack_type) {
      case 'C':
        Q = matrices->qb[fc->iindx[1] - n];
        break;

      case 'M':
        Q = matrices->qm[fc->iindx[1] - n];
        break;

      default:
        Q = (md->circ) ? matrices->qo : matrices->q[fc->iindx[1] - n];
        break;
    }

    /* ensemble free energy in Kcal/mol              */
    if (Q <= FLT_MIN)
      vrna_message_warning("pf_scale too large");

    free_energy = (-log(Q) - n * log(params->pf_scale)) *
                  params->kT /
                  1000.0;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      free_energy /= fc->n_seq;

#ifdef SUN4
    standard_arithmetic();
#elif defined(HP9)
    fpsetfastmode(0);
#endif
  }

  return free_energy;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE int
fill_arrays_pk(vrna_fold_compound_t *fc)
{
  unsigned char       *hard_constraints;
  int                 n, i, j, k, ij, d, *my_iindx, *jindx, with_gquad, turn,
                      with_ud, hc_decompose, *pscore;
  FLT_OR_DBL          temp, temp2, Qmax, qbt1, *q, *qb, *qm, *qm1, *q1k, *qln, *scale;
  double              kTn, max_real;
  vrna_ud_t           *domains_up;
  vrna_md_t           *md;
  vrna_hc_t           *hc;
  vrna_mx_pf_t        *matrices;
  vrna_mx_pf_aux_el_t aux_mx_el;
  vrna_mx_pf_aux_ml_t aux_mx_ml;
  vrna_exp_param_t    *pf_params;

  n                 = fc->length;
  my_iindx          = fc->iindx;
  jindx             = fc->jindx;
  pscore            = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->pscore : NULL;
  matrices          = fc->exp_matrices;
  pf_params         = fc->exp_params;
  kTn               = pf_params->kT / 10.;  /* kT in cal/mol */
  hc                = fc->hc;
  domains_up        = fc->domains_up;
  q                 = matrices->q;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  q1k               = matrices->q1k;
  qln               = matrices->qln;
  scale 			= matrices->scale;
  md                = &(pf_params->model_details);
  with_gquad        = md->gquad;
  turn              = md->min_loop_size;
  hard_constraints  = hc->matrix;

  with_ud = (domains_up && domains_up->exp_energy_cb && (!(fc->type == VRNA_FC_TYPE_COMPARATIVE)));
  Qmax    = 0;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (with_ud && domains_up->exp_prod_cb)
    domains_up->exp_prod_cb(fc, domains_up->data);

  /* no G-Quadruplexes for comparative partition function (yet) */
  if (with_gquad && (!(fc->type == VRNA_FC_TYPE_COMPARATIVE))) {
    free(fc->exp_matrices->G);
    fc->exp_matrices->G = get_gquad_pf_matrix(fc->sequence_encoding2,
                                              fc->exp_matrices->scale,
                                              fc->exp_params);
  }

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(fc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(fc);

  /*array initialization ; qb,qm,q
   * qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  for (d = 0; d <= turn; d++)
    for (i = 1; i <= n - d; i++) {
      j       = i + d;
      ij      = my_iindx[i] - j;
      qb[ij]  = 0.0;
    }

  for (j = turn + 2; j <= n; j++) {
    for (i = j - turn - 1; i >= 1; i--) {
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[jindx[j] + i];
      qbt1          = 0;
      
      //exp_E_pseudoknot(fc, i ,j);
      exp_E_pseudoknot(fc, i ,j);
      //printf("qpk1 <-> %f (%d,%d) (%d,%d)\n", qpk[my_iindx[1] - n], 1, n, i, j);

      if (hc_decompose) {
        /* process hairpin loop(s) */
        qbt1 += vrna_exp_E_hp_loop(fc, i, j);
        /* process interior loop(s) */
        qbt1 += vrna_exp_E_int_loop(fc, i, j);
        /* process multibranch loop(s) */
        qbt1 += vrna_exp_E_mb_loop_fast(fc, i, j, aux_mx_ml);

        if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_c))
          qbt1 += fc->aux_grammar->cb_aux_exp_c(fc, i, j, fc->aux_grammar->data);

        if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
          qbt1 *= exp(pscore[jindx[j] + i] / kTn);
      }

      qb[ij] = qbt1;

      /* Multibranch loop */
      temp = vrna_exp_E_ml_fast_pk(fc, i, j, aux_mx_ml);

      /* apply auxiliary grammar rule for multibranch loop case */
      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m))
        temp += fc->aux_grammar->cb_aux_exp_m(fc, i, j, fc->aux_grammar->data);

      qm[ij] = temp;

      if (qm1) {
        temp = vrna_exp_E_ml_fast_qqm_pk(aux_mx_ml)[i]; /* for stochastic backtracking and circfold */

        /* apply auxiliary grammar rule for multibranch loop (M1) case */
        if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m1))
          temp += fc->aux_grammar->cb_aux_exp_m1(fc, i, j, fc->aux_grammar->data);

        qm1[jindx[j] + i] = temp;
      }

      /* Exterior loop */
      temp = vrna_exp_E_ext_fast_pk(fc, i, j, aux_mx_el);
      
      /* apply auxiliary grammar rule for exterior loop case */
      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_f)){
        temp += fc->aux_grammar->cb_aux_exp_f(fc, i, j, fc->aux_grammar->data);
	}

      q[ij] = temp;	

      if (temp > Qmax) {
        Qmax = temp;
        if (Qmax > max_real / 10.)
          vrna_message_warning("Q close to overflow: %d %d %g", i, j, temp);
      }

      if (temp >= max_real) {
        vrna_message_warning("overflow while computing partition function for segment q[%d,%d]\n"
                             "use larger pf_scale", i, j);

        vrna_exp_E_ml_fast_free_pk(aux_mx_ml);
        vrna_exp_E_ext_fast_free_pk(aux_mx_el);

        return 0; /* failure */
      }
    }

    /* rotate auxiliary arrays */
    vrna_exp_E_ext_fast_rotate_pk(aux_mx_el);
    vrna_exp_E_ml_fast_rotate_pk(aux_mx_ml);
  }

  /* prefill linear qln, q1k arrays */
  if (q1k && qln) {
    for (k = 1; k <= n; k++) {
      q1k[k]  = q[my_iindx[1] - k];
      qln[k]  = q[my_iindx[k] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;
  }

  /* free memory occupied by auxiliary arrays for fast exterior/multibranch loops */
  vrna_exp_E_ml_fast_free(aux_mx_ml);
  vrna_exp_E_ext_fast_free(aux_mx_el);

  return 1;
}
