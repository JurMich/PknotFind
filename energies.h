#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include "rna.h"

#ifndef ENERGIES_H_
#define ENERGIES_H_

#define TYPE double

double TEMP; /*temperature */
double TEMPSCALE; /*temperature scaling factor*/
double RT; /*Boltzmann's cte*/

int MIN_TURN;

// parameters related with Vienna package
vrna_fold_compound_t *E_fold_cp;

/* converts int to a base */
int base2int(char base);

/* motif energies functions */
void set_E_fold_model(plain_sequence * rna);
int get_type(char base_5, char base_3);
int pairable(int base_5, int base_3, plain_sequence * rna);
TYPE energy_to_bfactor(double energy);
double stacking_energy(int i, int j);
double internal_loop_energy(int i, int k, int l, int j);

#endif