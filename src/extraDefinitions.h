//
// Created by Ethan Kang on 12/4/19.
//

#ifndef MYSIMULATOR_EXTRADEFINITIONS_H
#define MYSIMULATOR_EXTRADEFINITIONS_H

#define NUM_CAN_ML_GENOTYPES 10
#define FiniteMu 3
#define FiniteMuDel 4
#define FiniteMuInDel 5

#define INVALID (-1)
#define DELETED (-2)
#define ADOExtd (-3)

#define SINGLE_ADO 1
#define LOCUS_DROPOUT 2

#define LOW_VAF 2

typedef struct {
  int code;
  int numExistingAlleles;

  /*
   * At most three alleles can exist on two chromosomes.
   * The third allele is the inserted allele if there is insertions, otherwise
   * it is marked with INVALID.
   * Any allele being deleted is marked with DELETED.
   * On the tips, an ADO on an allele is marked with ADOExtd.
   * The order of the alleles follows the convention of the genotypes.
   * */
  int *alleles;

  /*
   * including ancestors
   * mutations, back mutations, potential back mutations, deletions, insertions
   */
  unsigned int *pastEvolutionaryEvents;
} DataStr;

typedef struct {
  long id;
  int start;
  int end; // does not include this site
  int numSubSites;
  int numSubSNVs;

  DataStr **data;
  int *hasParallelMutations;    // [numSubSites]
  unsigned **parallelMutations; // [numSubSites][numCells]
  int **summary; // [numSubSites][numMus, numBackMus, numDels, numIns]
  int *isSNV;    // [numSubSites]
  int numMutatedSubSites[3];
} ThreadStr;

typedef struct {
  int numReads;
  int *readCount;
  double shapeController;
  int trueCode, observedCode, MLCode;
  int trueNumExistingAlleles, observedNumExistingAlleles, MLNumExistingAlleles;
  int *trueAlleles, *observedAlleles, *MLAlleles;
  unsigned int *pastEvolutionaryEvents; // mutations, back mutations, potential
                                        // back mutations, deletions, insertions
  int hasADO, numADO;
  int hasAmplificationError, fixedAmplificationError;
  double amplificationError;
  double *log10GenotypeLikelihood, *scaledLog10GenotypeLikelihood;
  int uniqueGenotypeLikelihoodNr;
  double *uniqueLog10GenotypeLikelihood, *uniqueScaledLog10GenotypeLikelihood;
} CellSiteStrExtd;

typedef struct {
  CellSiteStrExtd *site;
  double sizeFactor;
} CellStrExtd;

#endif // MYSIMULATOR_EXTRADEFINITIONS_H
