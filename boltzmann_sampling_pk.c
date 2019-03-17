/*
 *   				THIS FILE IS TAKEN FROM ViennaRNA AND
 * 					IT IS DISTRIBUTED UNDER ITS LICENSE!
 * 
 *                 CREDITED TO:
 * 			
 *                partiton function for RNA secondary structures
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/alphabet.h"
#include "boltzmann_sampling_pk.h"
#include "pseudoknots.h"
#include "external_pf_pk.h"
#include "multibranch_pf_pk.h"

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
PRIVATE char *info_set_uniq_ml =
  "Activate unique multiloop decomposition by setting the"
  " uniq_ML field of the model details structure to a non-zero"
  " value before running vrna_pf()!";


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void  backtrack(int                   i,
                        int                   j,
                        char                  *pstruc,
                        vrna_fold_compound_t  *vc,
                        double			  *pk_energy);


void  backtrack_qm(int                  i,
                           int                  j,
                           char                 *pstruc,
                           vrna_fold_compound_t *vc,
                           double			  *pk_energy
                           );


PRIVATE void  backtrack_qm1(int                   i,
                            int                   j,
                            char                  *pstruc,
                            vrna_fold_compound_t  *vc,
                            double			  *pk_energy);


PRIVATE void  backtrack_qm2(int                   u,
                            int                   n,
                            char                  *pstruc,
                            vrna_fold_compound_t  *vc,
                            double			  *pk_energy);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/*
 * stochastic backtracking in pf_fold arrays
 * returns random structure S with Boltzman probabilty
 * p(S) = exp(-E(S)/kT)/Z
 */
PUBLIC char *
vrna_pbacktrack_pk(vrna_fold_compound_t *vc)
{
  char    *structure  = NULL;
  double  prob        = 1.;

  if (vc) {
    if (!vc->exp_params) {
      vrna_message_warning("vrna_pbacktrack: DP matrices are missing! Call vrna_pf() first!");
      return NULL;
    } else if (!vc->exp_params->model_details.uniq_ML) {
      vrna_message_warning("vrna_pbacktrack: Unique multiloop decomposition is unset!");
      vrna_message_info(stderr, info_set_uniq_ml);
      return NULL;
    }

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return vrna_pbacktrack5(vc, vc->length);

        break;

      default:
        vrna_message_warning("unrecognized fold compound type");
        return structure;
        break;
    }
  }

  return structure;
}


PUBLIC char *
vrna_pbacktrack5_pk(vrna_fold_compound_t *vc,
                 int                  length,
                 double				  *pk_energy)
{
  FLT_OR_DBL        r, qt, q_temp, qkl, qpkl;
  int               i, j, ij, n, k, u, type, start, paired;
  char              *pstruc;
  char              *ptype;
  int               *my_iindx, *jindx, hc_decompose, *hc_up_ext;
  FLT_OR_DBL        *q, *qb, *q1k, *qln, *scale;
  unsigned char     *hard_constraints;
  short             *S1, *S2;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  vrna_exp_param_t  *pf_params;

  n = vc->length;

  pf_params = vc->exp_params;
  md        = &(vc->exp_params->model_details);
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;
  matrices  = vc->exp_matrices;
  
  ptype = vc->ptype;

  hc  = vc->hc;
  sc  = vc->sc;
  S1  = vc->sequence_encoding;
  S2  = vc->sequence_encoding2;

  hard_constraints  = hc->matrix;
  hc_up_ext         = hc->up_ext;
  paired 			= 0; /* needed to take pseudoknots into account */

  if (length > n) {
    vrna_message_warning("vrna_pbacktrack5: 3'-end exceeds sequence length");
    return NULL;
  } else if (length < 1) {
    vrna_message_warning("vrna_pbacktrack5: 3'-end too small");
    return NULL;
  } else if ((!matrices) || (!matrices->q) || (!matrices->qb) || (!matrices->qm) || (!pf_params)) {
    vrna_message_warning("vrna_pbacktrack5: DP matrices are missing! Call vrna_pf() first!");
    return NULL;
  } else if ((!vc->exp_params->model_details.uniq_ML) || (!matrices->qm1)) {
    vrna_message_warning("vrna_pbacktrack5: Unique multiloop decomposition is unset!");
    vrna_message_info(stderr, info_set_uniq_ml);
    return NULL;
  }

  q     = matrices->q;
  qb    = matrices->qb;
  q1k   = matrices->q1k;
  qln   = matrices->qln;
  scale = matrices->scale;

  pstruc = vrna_alloc((length + 1) * sizeof(char));

  for (i = 0; i < length; i++)
    pstruc[i] = '.';

  if (!(q1k && qln)) {
    matrices->q1k = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    matrices->qln = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    q1k           = matrices->q1k;
    qln           = matrices->qln;
    for (k = 1; k <= n; k++) {
      q1k[k]  = q[my_iindx[1] - k];
      qln[k]  = q[my_iindx[k] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;
  }

#ifdef VRNA_WITH_BOUSTROPHEDON
  j = length;
  while (j > 1) {
    /* find j position of first pair */
    for (; j > 1; j--) {
      if (hc_up_ext[j]) {
        r       = vrna_urn() * q1k[j];
        q_temp  = q[my_iindx[1] - j + 1] * scale[1];

        if (sc) {
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[j][1];

          if (sc->exp_f)
            q_temp *= sc->exp_f(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (r > q_temp)
          break;                /* j is paired */
        *pk_energy *= q_temp/q[my_iindx[1] - j + 1];
      }
    }
    if (j <= md->min_loop_size + 1)
      break;         /* no more pairs */

    /* now find the pairing partner i */
    r = vrna_urn() * (q1k[j] - q_temp);
    u = j - 1;

    for (qt = 0, k = 1; k < j; k++) {
      /* apply alternating boustrophedon scheme to variable i */
      i = (int)(1 + (u - 1) * ((k - 1) % 2)) +
          (int)((1 - (2 * ((k - 1) % 2))) * ((k - 1) / 2));
      ij            = my_iindx[i] - j;
      hc_decompose  = hard_constraints[jindx[j] + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        type  = vrna_get_ptype_md(S2[i], S2[j], md);
        qkl   = qb[ij] * exp_E_ExtLoop(type,
                                       (i > 1) ? S1[i - 1] : -1,
                                       (j < n) ? S1[j + 1] : -1,
                                       pf_params);

        if (i > 1) {
          qkl *= q1k[i - 1];
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(1, j, i - 1, i, VRNA_DECOMP_EXT_EXT_STEM, sc->data);
        } else {
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        qt += qkl;
        if (qt > r){
		  *pk_energy *= qkl/qb[ij];
		  if (i > 1) {
			*pk_energy /= q1k[i - 1];
		  }
          paired = 1;
          break;           /* j is paired */
		}
		qpkl = qpk[ij];
		if(j < length)
			qpkl *= qln[j + 1];
		
		qt += qpkl;
		
		
		if (qt > r)		   /* j is involved in pseudoknot */
          break;				
      }
    }
    if (k == j)
      vrna_message_error("backtracking failed in ext loop");

	if(paired){ /* basic pair */
		backtrack(i, j, pstruc, vc, pk_energy);
	}else{ /* pseudoknot */
		traceback_pks(i, j, pstruc, vc, pk_energy);
	}
    j = i - 1;
  }
#else
  start = 1;
  while (start < length) {
    /* find i position of first pair */
    for (i = start; i < length; i++) {
      if (hc_up_ext[i]) {
        r       = vrna_urn() * qln[i];
        q_temp  = qln[i + 1] * scale[1];

        if (sc) {
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][1];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, length, i + 1, length, VRNA_DECOMP_EXT_EXT, sc->data);
        }

        if (r > q_temp)
          break;                /* i is paired */
        *pk_energy *= q_temp/qln[i + 1];
      }
    }
    if (i >= length){
	  *pk_energy *= scale[1];
      break;              /* no more pairs */
	}
	
    /* now find the pairing partner j */
    r = vrna_urn() * (qln[i] - q_temp);
    for (qt = 0, j = i + 1; j <= length; j++) {
      ij            = my_iindx[i] - j;
      type          = vrna_get_ptype(jindx[j] + i, ptype);
      hc_decompose  = hard_constraints[jindx[j] + i];
      if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        qkl = qb[ij] * exp_E_ExtLoop(type,
                                     (i > 1) ? S1[i - 1] : -1,
                                     (j < n) ? S1[j + 1] : -1,
                                     pf_params);

        if (j < length) {
          qkl *= qln[j + 1];
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(i, length, j, j + 1, VRNA_DECOMP_EXT_STEM_EXT, sc->data);
        } else {
          if (sc)
            if (sc->exp_f)
              qkl *= sc->exp_f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
        }

        qt += qkl;
        if (qt > r){
		  *pk_energy *= qkl/qb[ij];
		  if (j < length) {
			*pk_energy /= qln[j + 1];
		  }
		  paired = 1;	
          break;           /* j is paired */
		}			
      }
      qpkl = qpk[ij];
      if(j < length)
        qpkl *= qln[j + 1];
		
      qt += qpkl;
			
      if (qt > r)		   /* j is involved in pseudoknot */
        break;	
    }
    if (j == length + 1)
      vrna_message_error("backtracking failed in ext loop");

    start = j + 1;
    if(paired){
		backtrack(i, j, pstruc, vc, pk_energy);
		paired = 0;		  /* must reset  the value after the backtrak is performed */
	}else{ /* pseudoknot */
		traceback_pks(i, j, pstruc, vc, pk_energy);
	}
	if(start == length){
		*pk_energy *= scale[1];
	}
  }
#endif
  return pstruc;
}



void backtrack_qm(int                  i,
             int                  j,
             char                 *pstruc,
             vrna_fold_compound_t *vc,
             double			  *pk_energy)
{
	
  /* divide multiloop into qm and qm1  */
  FLT_OR_DBL    qmt, r, q_temp;
  int           k, u, cnt, span, turn;
  FLT_OR_DBL    *qm, *qm1, *expMLbase;
  int           *my_iindx, *jindx, *hc_up_ml;
  vrna_sc_t     *sc;
  vrna_hc_t     *hc;

  vrna_mx_pf_t  *matrices = vc->exp_matrices;

  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  hc        = vc->hc;
  sc        = vc->sc;
  hc_up_ml  = hc->up_ml;

  qm        = matrices->qm;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;

  turn = vc->exp_params->model_details.min_loop_size;
  


	//printf("qmtot (%d,%d) %f\n",i,j, qm[my_iindx[i] - j] );

  while (j > i) {
    /* now backtrack  [i ... j] in qm[] */
    r   = vrna_urn() * qm[my_iindx[i] - j];
    
    qmt = qm1[jindx[j] + i];
    k   = cnt = i;
    if (qmt < r) {
      for (span = j - i, cnt = i + 1; cnt <= j; cnt++) {
#ifdef VRNA_WITH_BOUSTROPHEDON
        k = (int)(i + 1 + span * ((cnt - i - 1) % 2)) +
            (int)((1 - (2 * ((cnt - i - 1) % 2))) * ((cnt - i) / 2));
#else
        k = cnt;
#endif
        q_temp  = 0.;
        u       = k - i;
        /* [i...k] is unpaired */
        if (hc_up_ml[i] >= u) {
          q_temp += expMLbase[u] * qm1[jindx[j] + k];

          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][u];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);
          }

          qmt += q_temp;
        }

        /* split between k-1, k */
       
        q_temp = qm[my_iindx[i] - (k - 1)] * qm1[jindx[j] + k];
		
		if (sc)
          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

        qmt += q_temp;
        

        if (qmt >= r)
          break;

      }
    }

    if (cnt > j)
      vrna_message_error("backtrack failed in qm");

	
	backtrack_qm1(k, j, pstruc, vc, pk_energy);

    if (k < i + turn){
	  *pk_energy *= expMLbase[k-i];
      break;            /* no more pairs */
	}

    u = k - i;
    /* check whether we make the decision to leave [i..k-1] unpaired */
    if (hc_up_ml[i] >= u) {
      q_temp = expMLbase[u];

      if (sc) {
        if (sc->exp_energy_up)
          q_temp *= sc->exp_energy_up[i][u];

        if (sc->exp_f)
          q_temp *= sc->exp_f(i, k - 1, i, k - 1, VRNA_DECOMP_ML_UP, sc->data);
      }

      r = vrna_urn() * (qm[my_iindx[i] - (k - 1)] + q_temp);
      if (q_temp >= r){
		*pk_energy *= q_temp;  
		break;
		
	  }
    }

    j = k - 1;
  }
}


PRIVATE void
backtrack_qm1(int                   i,
              int                   j,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              double			  *pk_energy)
{
  /* i is paired to l, i<l<j; backtrack in qm1 to find l */
  int               ii, l, il, type, turn, paired;
  FLT_OR_DBL        qt, r, q_temp;
  FLT_OR_DBL        *qm1, *qb, *expMLbase;
  vrna_mx_pf_t      *matrices;
  int               u, *my_iindx, *jindx, *hc_up_ml;
  char              *ptype;
  unsigned char     *hard_constraints;
  short             *S1;
  vrna_sc_t         *sc;
  vrna_hc_t         *hc;
  vrna_exp_param_t  *pf_params;


  pf_params = vc->exp_params;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  ptype = vc->ptype;

  sc                = vc->sc;
  hc                = vc->hc;
  hc_up_ml          = hc->up_ml;
  hard_constraints  = hc->matrix;

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm1       = matrices->qm1;
  expMLbase = matrices->expMLbase;
  S1        = vc->sequence_encoding;

  turn = pf_params->model_details.min_loop_size;

  paired 	= 0;

  r   = vrna_urn() * qm1[jindx[j] + i];
  
  ii  = my_iindx[i];
  for (qt = 0., l = j; l > i + turn; l--) {
    il = jindx[l] + i;
    if (hard_constraints[il] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
      u = j - l;
      if (hc_up_ml[l + 1] >= u) {
        type    = vrna_get_ptype(il, ptype);
        q_temp  = qb[ii - l]
                  * exp_E_MLstem(type, S1[i - 1], S1[l + 1], pf_params)
                  * expMLbase[j - l];

        if (sc) {
          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[l + 1][j - l];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, i, l, VRNA_DECOMP_ML_STEM, sc->data);
        }
		
        qt += q_temp;
        if (qt >= r){
		  paired = 1;
		  *pk_energy *= q_temp/qb[ii - l];
          break;
        
		}			
      } else {
        l = i + turn;
        break;
      }
    }
    
    qt += qpk[my_iindx[i] - l] * expMLbase[j - l];
	if (qt >= r){
	  *pk_energy *= expMLbase[j - l];
      break;
	}
  }
  if (l < i + turn + 1)
    vrna_message_error("backtrack failed in qm1");

  if(paired){	
	backtrack(i, l, pstruc, vc, pk_energy);
  }else{
	traceback_pks(i, l, pstruc, vc, pk_energy);  
  }	  
}


PRIVATE void
backtrack_qm2(int                   k,
              int                   n,
              char                  *pstruc,
              vrna_fold_compound_t  *vc,
              double			  *pk_energy)
{
  FLT_OR_DBL  qom2t, r;
  int         u, turn;
  FLT_OR_DBL  *qm1, *qm2;
  int         *jindx;
  vrna_sc_t   *sc;

  jindx = vc->jindx;
  qm1   = vc->exp_matrices->qm1;
  qm2   = vc->exp_matrices->qm2;
  turn  = vc->exp_params->model_details.min_loop_size;
  sc    = vc->sc;

  r = vrna_urn() * qm2[k];
  /* we have to search for our barrier u between qm1 and qm1  */
  if ((sc) && (sc->exp_f)) {
    for (qom2t = 0., u = k + turn + 1; u < n - turn - 1; u++) {
      qom2t +=  qm1[jindx[u] + k] *
                qm1[jindx[n] + (u + 1)] *
                sc->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

      if (qom2t > r)
        break;
    }
  } else {
    for (qom2t = 0., u = k + turn + 1; u < n - turn - 1; u++) {
      qom2t += qm1[jindx[u] + k] * qm1[jindx[n] + (u + 1)];
      if (qom2t > r)
        break;
    }
  }

  if (u == n - turn)
    vrna_message_error("backtrack failed in qm2");

  backtrack_qm1(k, u, pstruc, vc, pk_energy);
  backtrack_qm1(u + 1, n, pstruc, vc, pk_energy);
}


PRIVATE void
backtrack(int                   i,
          int                   j,
          char                  *pstruc,
          vrna_fold_compound_t  *vc,
          double			  *pk_energy)
{
  char              *ptype, *sequence;
  unsigned char     *hard_constraints, hc_decompose;
  vrna_exp_param_t  *pf_params;
  FLT_OR_DBL        *qb, *qm, *qm1, *scale, *expMLbase;	
  FLT_OR_DBL        r, qbt1, qt, q_temp;
  vrna_mx_pf_t      *matrices;
  int               *my_iindx, *jindx, *hc_up_int;
  vrna_sc_t         *sc;
  vrna_hc_t         *hc;
  short             *S1;

  sequence  = vc->sequence;
  pf_params = vc->exp_params;
  ptype     = vc->ptype;
  S1        = vc->sequence_encoding;
  my_iindx  = vc->iindx;
  jindx     = vc->jindx;

  sc                = vc->sc;
  hc                = vc->hc;
  hc_up_int         = hc->up_int;
  hard_constraints  = hc->matrix;

  matrices  = vc->exp_matrices;
  qb        = matrices->qb;
  qm        = matrices->qm;
  qm1       = matrices->qm1;
  scale     = matrices->scale;
  expMLbase = matrices->expMLbase;

  int noGUclosure = pf_params->model_details.noGUclosure;
  int turn        = pf_params->model_details.min_loop_size;
  int *rtype      = &(pf_params->model_details.rtype[0]);

  qbt1 = 0.;

  do {
    int           k, l, kl, u, u1, u2, max_k, min_l;
    unsigned char type;
    k = i;
    l = j;

    pstruc[i - 1] = '(';
    pstruc[j - 1] = ')';

    r             = vrna_urn() * qb[my_iindx[i] - j];
    type          = vrna_get_ptype(jindx[j] + i, ptype);
    hc_decompose  = hard_constraints[jindx[j] + i];

    /* hairpin contribution */
    qbt1 = vrna_exp_E_hp_loop(vc, i, j);
    
    //printf("hairpino %f (%d,%d)\n", qbt1, i,j);

    if (qbt1 >= r){
	  *pk_energy *= qbt1;
      return;            /* found the hairpin we're done */
	}

    if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
      /* interior loop contributions */
      max_k = i + MAXLOOP + 1;
      max_k = MIN2(max_k, j - turn - 2);
      max_k = MIN2(max_k, i + 1 + hc_up_int[i + 1]);
      for (k = i + 1; k <= max_k; k++) {
        u1    = k - i - 1;
        min_l = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
        kl    = my_iindx[k] - j + 1;
        for (u2 = 0, l = j - 1; l >= min_l; l--, kl++, u2++) {
          if (hc_up_int[l + 1] < u2)
            break;

          if (hard_constraints[jindx[l] + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
            unsigned int type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

            /* add *scale[u1+u2+2] */
            q_temp = qb[kl]
                     * scale[u1 + u2 + 2]
                     * exp_E_IntLoop(u1,
                                     u2,
                                     type,
                                     type_2,
                                     S1[i + 1],
                                     S1[j - 1],
                                     S1[k - 1],
                                     S1[l + 1],
                                     pf_params);

            if (sc) {
              if (sc->exp_energy_up)
                q_temp *= sc->exp_energy_up[i + 1][u1]
                          * sc->exp_energy_up[l + 1][u2];

              if (sc->exp_energy_bp)
                q_temp *= sc->exp_energy_bp[jindx[j] + i];

              if (sc->exp_energy_stack) {
                if ((i + 1 == k) && (j - 1 == l)) {
                  q_temp *= sc->exp_energy_stack[i]
                            * sc->exp_energy_stack[k]
                            * sc->exp_energy_stack[l]
                            * sc->exp_energy_stack[j];
                }
              }

              if (sc->exp_f)
                q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
            }

            qbt1 += q_temp;
            if (qbt1 >= r){
			  *pk_energy *= q_temp/qb[kl];
              break;
            }      
          }
          qbt1 += expMLbase[k - i - 1] * expMLbase[j - l - 1] * qpk[kl];
          if (qbt1 >= r){ /* found a pseudoknot, we're done (after we check it out) */
			*pk_energy *= expMLbase[k - i - 1] * expMLbase[j - l - 1];
			traceback_pks(k, l, pstruc, vc, pk_energy);
            return;      
              
          }
        }
        if (qbt1 >= r)
          break;
      }
      if (k <= max_k) {
        i = k;
        j = l;
      } else {
        /* interior loop contributions did not exceed threshold, so we break */
        break;
      }
    } else {
      /* must not be interior loop, so we break out */
      break;
    }
  } while (1);

  /* backtrack in multi-loop */
  {
    int         k, ii, jj, tt;
    FLT_OR_DBL  closingPair;
    tt          = rtype[vrna_get_ptype(jindx[j] + i, ptype)];
    closingPair = pf_params->expMLclosing
                  * exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params)
                  * scale[2];
    if (sc) {
      if (sc->exp_energy_bp)
        q_temp *= sc->exp_energy_bp[jindx[j] + i];

      if (sc->exp_f)
        closingPair *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    i++;
    j--;
    /* find the first split index */
    ii  = my_iindx[i];  /* ii-j=[i,j] */
    jj  = jindx[j];     /* jj+i=[j,i] */

    if ((sc) && (sc->exp_f)) {
      for (qt = qbt1, k = i + 1; k < j; k++) {
        q_temp =  qm[ii - (k - 1)] *
                  qm1[jj + k] *
                  closingPair *
                  sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);

        qt    += q_temp;
        qbt1  += q_temp;

        if (qt >= r){
		  *pk_energy *= closingPair * sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
          break;
	    }
      }
    } else {
      for (qt = qbt1, k = i + 1; k < j; k++) {
        q_temp =  qm[ii - (k - 1)] *
                  qm1[jj + k] *
                  closingPair;

        qt    += q_temp;
        qbt1  += q_temp;

        if (qt >= r){
		  *pk_energy *= closingPair;
          break;
        }
      }
    }
    if (k >= j)
      vrna_message_error("backtrack failed, can't find split index ");

    backtrack_qm1(k, j, pstruc, vc, pk_energy);

    j = k - 1;
    backtrack_qm(i, j, pstruc, vc, pk_energy);
  }
}
