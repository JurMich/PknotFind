#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include "energies.h"
#include "rna.h"
#include "helices.h"
#include "pseudoknots.h"
#include "part_func_pk.h"

char * nameFileInFasta;
double TEMPSCALE;
int N_SAMPLES;

void read_parameters(int argc, char **argv){
	char char_read;
	while((char_read = getopt(argc, argv, "i:j:p:k:l:")) != EOF){
		switch(char_read){
			case 'i' : /* input sequence in Fasta format*/
				nameFileInFasta=(char*) malloc ((strlen(optarg)+1) * sizeof(char));
				strcpy (nameFileInFasta, optarg);
				break;
			case 'j' : 	
				MAX_SKEW=atoi(optarg);
				if (MAX_SKEW<0){
					printf("\nBad parameter value (-j): This value should be positive.\n\n");  
					exit(EXIT_SUCCESS);
				}
				break;
			case 'p' :
				N_SAMPLES=atoi(optarg);
				if (MAX_SKEW<0){
					printf("\nBad parameter value (-p): This value should be positive.\n\n");  
					exit(EXIT_SUCCESS);
				}
				break;
			case 'k' :
				PK_ENERGY=atoi(optarg);
				if (PK_ENERGY<0){
					printf("\nBad parameter value (-k): This value should be positive.\n\n");  
					exit(EXIT_SUCCESS);
				}
				break;
			case 'l' :
				MAX_LIST=atoi(optarg);
				if (MAX_LIST<0){
					printf("\nBad parameter value (-l): This value should be positive.\n\n");  
					exit(EXIT_SUCCESS);
				}
				break;
			case 'm' :
				MAX_LIST=atoi(optarg);
				if (MAX_LIST<0){
					printf("\nBad parameter value (-m): This value should be positive.\n\n");  
					exit(EXIT_SUCCESS);
				}
				break;
		}		
	}	
}

int main(int argc, char **argv){
	srand ( time(NULL) );
	char RNAname[200];
	nameFileInFasta = " ";
	
	/* default parameters */
	TEMP = 37;
	TEMPSCALE = 1;
	RT = TEMPSCALE*0.0019872370936902486 * (273.15 + TEMP) * 100;
	MAX_SKEW = 2;
	MAX_LIST = 10;
	MIN_PK = 6;
	PK_ENERGY = 900.;
	N_SAMPLES = 100;
	
	read_parameters(argc, argv);
	
	/* create sequence then model */
	plain_sequence  * rna_seq;
	rna_seq= (plain_sequence *) get_plain_sequence(nameFileInFasta, RNAname);
	set_E_fold_model(rna_seq);
	
	printf("%d\n",MAX_LIST);
	compute_helix_partition(rna_seq);
	int i,j,k, delta, d;
	
	//~ for(i = 1; i<rna_seq->size; i++){
		//~ for(j = i; j<=rna_seq->size; j++){
			//~ helix *HX = HELIX_PART_FCI[i][j]->first_hx;
			//~ if(HX)
				//~ printf("Start of Magneto %p (%d,%d)\n", HELIX_PART_FCI[i][j]->first_hx, i,j);
			//~ while(HX){
				//~ printf("H: %f (%f) (%d,%d)-(%d,%d)\n", HX->goodness, HX->helix_pf, i,j, HX->k, HX->l);
				//~ HX = HX->next_hx;
			//~ }
		//~ }
	//~ }
	
	init_aux_grammar(rna_seq, E_fold_cp);
	float floatzel;
	floatzel = vrna_pf_pk(E_fold_cp, NULL);
	
	vrna_mx_pf_t 	* matrix;
	FLT_OR_DBL		* all, *qb, *qm, *qm1;
	matrix = E_fold_cp->exp_matrices;
	all = matrix->q;
	qb = matrix->qb;
	qm = matrix->qm;
	qm1 = matrix->qm1;
	
	//~ for(i = 1; i<rna_seq->size; i++){
		//~ for(j = i; j<=rna_seq->size; j++){
			//~ printf("q(%d,%d) = %f, qm(%d,%d) = %f\n",i,j, all[E_fold_cp->iindx[i] - j], qm[E_fold_cp->iindx[i] - j]);

		//~ }
	//~ }

	char * pstruc;
	pstruc = (char*)malloc(sizeof(char) * (rna_seq->size + 2));
	pstruc[rna_seq->size+1] = '\0';
	
	char * pstruc2;

	display_plain_sequence(rna_seq, NULL);

	for(int a = 1; a <=N_SAMPLES; a++){
		pstruc2 = vrna_pbacktrack5(E_fold_cp, rna_seq->size);
		printf("%s \n", pstruc2);
	}
	
	printf("Total partition function: q = %f \n", all[E_fold_cp->iindx[1] - rna_seq->size]*pow(E_fold_cp->exp_params->pf_scale, rna_seq->size));
	
}

