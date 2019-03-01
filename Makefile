CC = gcc
#CFLAGS = -O3 -Wall
#-ansi removed from nextt line
override CFLAGS +=  -IViennaRNA -std=gnu99 -Wall -pedantic -ggdb -g -O3 -fstack-protector-all 
LIBS =   -lRNA -lm -lmpfr -lgmp -fopenmp 

PKnotFind: main.o  rna.o energies.o helices.o pseudoknots.o external_pf_pk.o multibranch_pf_pk.o part_func_pk.o boltzmann_sampling_pk.o
	$(CC) $(CFLAGS) -o PKnotFind main.o rna.o energies.o helices.o pseudoknots.o external_pf_pk.o multibranch_pf_pk.o part_func_pk.o boltzmann_sampling_pk.o $(LIBS)

main.o: main.c  rna.h  helices.h pseudoknots.h external_pf_pk.h multibranch_pf_pk.h part_func_pk.h boltzmann_sampling_pk.h 

boltzmann_sampling_pk.o: pseudoknots.h external_pf_pk.h multibranch_pf_pk.h boltzmann_sampling_pk.h boltzmann_sampling_pk.c

part_func_pk.o: pseudoknots.h external_pf_pk.h multibranch_pf_pk.h part_func_pk.h part_func_pk.c 

multibranch_pf_pk.o: multibranch_pf_pk.h pseudoknots.h multibranch_pf_pk.c

external_pf_pk.o: external_pf_pk.h pseudoknots.h external_pf_pk.c

pseudoknots.o: pseudoknots.c pseudoknots.h

helices.o: helices.c helices.h

energies.o: energies.c energies.h

rna.o: rna.c rna.h

clean:
	rm -f *.o ; rm RNANR