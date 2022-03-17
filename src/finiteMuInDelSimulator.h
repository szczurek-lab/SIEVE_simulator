//
// finiteMuInDelSimulator.h
// CellCoal
//
// Created by Senbai Kang on 11/20/19.
//

#ifndef MYSIMULATOR_FINITEMUINDELSIMULATOR_H
#define MYSIMULATOR_FINITEMUINDELSIMULATOR_H

#include "definitions.h"
#include "hashMap.h"
#include <limits.h>
#include <omp.h>
#include <pthread.h>
#include <string.h>

typedef struct {
  long id;
  TreeNode *root;
  int numGenotypes;
  long int *seed;
} TreeRootRandomizer;

typedef struct {
  long id;
  CellStrExtd *cellStr;
} ThreadedCellStr;

typedef struct {
  FILE *fp;
  char fileName[MAX_NAME];
} RawDataFile;

typedef struct {
  long id;
  unsigned int numCells;
  long *seed;

  unsigned int *cellIndices;
  RawDataFile *files;
} ThreadedRawDataPrinter;

int GetMaxInt(const int *arr, size_t size);
int GetShapeController(const int *alleles, int numExistingAlleles);
void Normalize(double *arr, size_t size);
void PrepareGlobalFilesExtd(int argc, char **argv);
void PrepareSeparateFoldersAndFiles(int replicate);
void PrepareSeparateFilesExtd(int replicate);
void PrepareRawReadSeqDataFiles(int replicate);
void InitializeGenomesExtd(TreeNode *p, long int *seed);
extern int EigenREVMD(double RootExtd[], double CijkExtd[]);
extern int EigenREVMDI(double RootExtd[], double CijkExtd[]);
void EvolveSitesOnTreeExtd(TreeNode *treeRoot, long int *seed);
void SimulateFiniteMuInDel(TreeNode *treeRoot, long int *seed);
void ThreadedSimulateFiniteMuInDel(TreeNode *treeRoot, int nrOfGenotypes,
                                   long int *seed);
void *ThreadedSimulateFiniteMuInDelForSite(void *arg);
void SimulateFiniteMuInDelForSite(TreeNode *treeRoot, ThreadStr *threadStr,
                                  int site, int nrOfGenotypes, double **EPij,
                                  long int *seed);
// void ThreadedLikelihoodRatioTest(CellStrExtd *cellStr);
// void *ThreadedLikelihoodRatioTestForSite(void *arg);
int ChooseUniformStateExtd(const double *prob, long int *seed);
int CountTrueVariantsExtd();
void GenerateReadCountsExtd(int replicate, long int *seed);
void SiteReadCountsExtd(CellSiteStrExtd *c, double sizeFactor, int i, int j,
                        double *probs, double **ngsEij, double **ampEij,
                        long int *seed);
void RandomGaussian(int size, double mean, double variance,
                    const long int *seed, double *out);
int GenerateSequencingCoverage(int numExistingAlleles, double alleleCoverage,
                               double sizeFactor, double alleleRawVariance,
                               long int *seed);
int RandomBetaBinomial(double mean, double overdispersion, int numTrials,
                       long int *seed);
void RandomDirichletMultinomial(double *probs, double overdispersion,
                                int numTrials, int pivot, size_t size,
                                long int *seed, int *readCounts);
void GenotypeLikelihoodsExtd(CellSiteStrExtd *c, int i, int j, double **ngsEij,
                             double **ampEij);
int CountAllelesInMLGenotypesExtd(long int *seed);
int CountTrueSNVs();
int ValidateFilterAndDetermineVariantNucleotide(int numAltAllelesInCalls,
                                                int site, long int *seed);
void ProcessBackgroundInformation(int replicate);
void ProcessTotalBackgroundInformation(int replicate);
void ProcessCoverageProb(int replicate);
void PrintTumCellNames(FILE *fp, char separator, int includeNormalCell);
void PrintSNVSitesNames(FILE *fp);
void PrintSNVAllelicInformation(FILE *fp);
void PrintAdoStates(FILE *fp);
void PrintParallelMut(FILE *fp);
void PrintSizeFactors(FILE *fp);
void PrintTrueSNV(FILE *fpAl, FILE *fpNr, FILE *fpMt);
void PrintEvolutionaryEventsSummary(FILE *fp);
void PrintRawData(RawDataFile *files);
void PrintRawDataFASTQ();
void *ThreadedPrintRawDataFASTQ(void *arg);
char GetANucleotide(int *readCounts, double emptyProb, long *seed);
void PrintTrueSNVCov(FILE *fp);
void PrintMPileUp(FILE *fpNM, FILE *fpCN);
void PrintForMonoVar(int replicate, FILE *fpNMTrueVal, FILE *fpCNTrueVal,
                     FILE *fpNMDFTVal, FILE *fpCNDFTVal);
void PrintForSciPhi(int replicate, FILE *fpNM, FILE *fpCN);
void PrintSNVGenotypesExtd(FILE *fp, int toStderrFlag);
int GenerateAnAlternativeNucleotide(int referenceAllele, long int *seed);
void PrintFullGenotypesExtd(FILE *fp, int toStderrFlag);
void PrintVCFHead(FILE *fp, int simpler);
void PrintVCFExtd(FILE *fp, char missing, int phased, int simpler);
void PrintReads(FILE *fp);
void PrintForSiFit(FILE *fpCmdTrueVal, FILE *fpCmdDFTVal, FILE *fpData,
                   char *dataName);
void PrintForCellPhy(FILE *fpFASTACmd, FILE *fpFASTAData, char *FASTADataName,
                     FILE *fpGLVCFCmd, FILE *fpGLVCFData, char *GLVCFDataName,
                     FILE *fpEPVCFCmd, FILE *fpEPVCFData, char *EPVCFDataName);
void PrintSiteInfoExtd(FILE *fp, int i);
void PrintRunInformationExtd(FILE *fp);
void PrintSummaryToLog(FILE *fp, int replicate, int toStderrFlag);

unsigned int workingMode;

/**
 * 0: read length of both ends
 * 1: read length of internal fragments
 */
static int readLengths[2] = {0, 150};
int numReadBlocks;
int *startSiteBlocks;

/**
 * for two OBSERVED alleles and MuInDel model
 * 0: no OBSERVED number of existing alleles larger than 2 in this read block
 * 1: otherwise
 */
const double noReadProb[] = {0.0, 0.05};

/**
 * cells; read blocks
 */
double **noReadProbBlocks;

pthread_t *threads;
int numThreads;
int *sitesPointsForThreads;
ThreadStr *threadStr;
int fixedTree;

ThreadedRawDataPrinter *rawDataPrinters;

// always print out
char cellNamesFile[MAX_NAME], SNVSitesNamesFile[MAX_NAME],
    SNVSitesNamesDir[MAX_NAME], SNVAllelicSeqInfoFile[MAX_NAME],
    SNVAllelicSeqInfoDir[MAX_NAME];
FILE *fpSNVSitesNames, *fpSNVAllelicSeqInfo;

// only under working mode 1
char trueSNVGenotypesAlFile[MAX_NAME], trueSNVGenotypesNrFile[MAX_NAME],
    trueSNVTernaryMatrixFile[MAX_NAME], trueSNVDir[MAX_NAME],
    evolutionaryEventsFile[MAX_NAME], rawSeqDataRefFile[MAX_NAME],
    rawSeqDataDir[MAX_NAME], MPileUpDataNMFile[MAX_NAME],
    MPileUpDataCNFile[MAX_NAME], MPileUpDataDir[MAX_NAME],
    parallelMutDir[MAX_NAME], parallelMutFile[MAX_NAME], adoStatesDir[MAX_NAME],
    adoStatesFile[MAX_NAME], trueSNVCovFile[MAX_NAME], bamFilesDir[MAX_NAME],
    MonoVarDir[MAX_NAME], bamFileNamesIncNMMonoVar[MAX_NAME],
    bamFileNamesExcNMMonoVar[MAX_NAME], MonoVarCMDIncNMTrueVal[MAX_NAME],
    MonoVarCMDExcNMTrueVal[MAX_NAME], MonoVarCMDIncNMDFTVal[MAX_NAME],
    MonoVarCMDExcNMDFTVal[MAX_NAME], SciPhiDir[MAX_NAME],
    bamFileNamesIncNMSciPhi[MAX_NAME], bamFileNamesExcNMSciPhi[MAX_NAME],
    SciPhiCMDIncNM[MAX_NAME], SciPhiCMDExcNM[MAX_NAME];
FILE *fpTrueSNVGenotypesAl, *fpTrueSNVGenotypesNr, *fpTrueSNVTernaryMatrix,
    *fpMPileUpDataNM, *fpMPileUpDataCN, *fpParallelMut, *fpAdoStates,
    *fpTrueSNVCov, *fpEvolutionaryEvents, *fpMonoVarCMDIncNMTrueVal,
    *fpMonoVarCMDExcNMTrueVal, *fpMonoVarCMDIncNMDFTVal,
    *fpMonoVarCMDExcNMDFTVal, *fpBamFileNamesIncNMMonoVar,
    *fpBamFileNamesExcNMMonoVar, *fpSciPhiCMDIncNM, *fpSciPhiCMDExcNM,
    *fpBamFileNamesIncNMSciPhi, *fpBamFileNamesExcNMSciPhi;

/* size: numCells + 2
 * 0: reference sequence;
 * 1 to numCells + 1: cancerous cells
 * numCells + 2: healthy cell
 */
RawDataFile *rawDataFiles;

/**
 * Tree structure represented as a matrix with non-negative elements.
 * Each element represents the relative parent-child relationship between its
 * row and column. 0: no parent-child relationship. 1: the row node itself. 2:
 * the column node is a direct child of the row node. 3: the column node is a
 * grandchild of the row node. 4: ... and so on.
 */
unsigned **treeStructure;

int numMutatedSites[4]; /* 0: the number of mutated and recorded sites;
                         * 1: the number of sites wrongly classified as mutated;
                         * 2: the number of mutated sites with low vaf
                         * (does not belong either to SNV sites or to
                         * background sites);
                         * 3: the number of mutated but not recorded sites.
                         */
char **cellNamesExtd;
int ADOtype;
static int filterSNVs = NO, filterLowVAF = YES, cellGroupMode = 0;
int filterCellThreshold, filterAltNucleotideReadsThreshold;
static double meanLikelihoodRatioTestFilter = 0.25;
double MijMD[15][15];  // mutation, deletion
double MijMDI[65][65]; // mutation, deletion, insertion
double allelicCovMean, allelicCovVar, allelicRawVarMean, allelicRawVarVar;
double *allelicCov, *allelicRawVar; // for all sites
double *sizeFactors;
double meanSizeFactor, varSizeFactor;
double shapeController[2] = {0.0, 0.0};
static int cumNucleotides[4] = {1, 2, 3, 4};
static int fixedAmplificationError = 1;
static double numericalStabilizer = 1e-6;
static int doPrintHealthyTip = NO;
static int doCoalescenceTree = YES;
static int doPrintFullNucleotideReads =
    NO; // save reads in a form of: 0, (A, C, G, T), or: 1, (c, m)

// save the data set in a form of being processed by SiFit
static int doPrintForSiFit = NO, doPrintTreeWithoutTrunk = NO,
           doPrintDataType = 0; // 1 for binary, 2 for ternary

// save the data set in a form of being processed by CellPhy
static int doPrintForCellPhy = NO, doPrintTreeWithHealthyTip = NO,
           doPrintForCellPhyMode = 0; // 1: FASTA, 2: vcf, 3: both

FILE *fpTreeWithoutTrunk, *fpTreeWithHealthyTip, *fpCellNames, *fpSizeFactors,
    *fpReadsTrueSNVs, *fpReadsCalledSNVs, *fpSiFitCellNames,
    *fpDataSiFitTrueSNVs, *fpDataSiFitCalledSNVs, *fpCMDSiFitTrueSNVsTrueVal,
    *fpCMDSiFitCalledSNVsTrueVal, *fpCMDSiFitTrueSNVsDFTVal,
    *fpCMDSiFitCalledSNVsDFTVal, *fpFASTADataCellPhyTrueSNVs,
    *fpFASTACMDCellPhyTrueSNVs, *fpFASTADataCellPhyCalledSNVs,
    *fpFASTACMDCellPhyCalledSNVs, *fpGLVCFDataCellPhyTrueSNVs,
    *fpGLVCFCMDCellPhyTrueSNVs, *fpGLVCFDataCellPhyCalledSNVs,
    *fpGLVCFCMDCellPhyCalledSNVs, *fpEPVCFDataCellPhyTrueSNVs,
    *fpEPVCFCMDCellPhyTrueSNVs, *fpEPVCFDataCellPhyCalledSNVs,
    *fpEPVCFCMDCellPhyCalledSNVs;
char treeWithTrunkDir[MAX_NAME], treeWithoutTrunkDir[MAX_NAME],
    treeWithoutTrunkFile[MAX_NAME], treeWithoutTrunkFileCopy[MAX_NAME],
    treeWithHealthyTipDir[MAX_NAME], treeWithHealthyTipFile[MAX_NAME],
    dataFromTrueSNVsDir[MAX_NAME], dataFromCalledSNVsDir[MAX_NAME],
    sizeFactorsDir[MAX_NAME], sizeFactorsFile[MAX_NAME], readsDir[MAX_NAME],
    readsTrueSNVsFile[MAX_NAME], readsCalledSNVsFile[MAX_NAME],
    SiFitDir[MAX_NAME], SiFitCellNames[MAX_NAME],
    SiFitDataTrueSNVsFile[MAX_NAME], SiFitDataCalledSNVsFile[MAX_NAME],
    SiFitCMDTrueSNVsFileTrueVal[MAX_NAME],
    SiFitCMDCalledSNVsFileTrueVal[MAX_NAME],
    SiFitCMDTrueSNVsFileDFTVal[MAX_NAME],
    SiFitCMDCalledSNVsFileDFTVal[MAX_NAME], CellPhyDir[MAX_NAME],
    CellPhyDataTrueSNVsFASTAFile[MAX_NAME],
    CellPhyDataTrueSNVsFASTAFileCopy[MAX_NAME],
    CellPhyDataCalledSNVsFASTAFile[MAX_NAME],
    CellPhyDataCalledSNVsFASTAFileCopy[MAX_NAME],
    CellPhyCMDTrueSNVsFASTAFile[MAX_NAME],
    CellPhyCMDCalledSNVsFASTAFile[MAX_NAME],
    // corresponding to GL mode of CELLPHY
    CellPhyDataTrueSNVsVCFGLFile[MAX_NAME],
    CellPhyDataTrueSNVsVCFGLFileCopy[MAX_NAME],
    CellPhyDataCalledSNVsVCFGLFile[MAX_NAME],
    CellPhyDataCalledSNVsVCFGLFileCopy[MAX_NAME],
    CellPhyCMDTrueSNVsVCFGLFile[MAX_NAME],
    CellPhyCMDCalledSNVsVCFGLFile[MAX_NAME],
    // corresponding to EP17 mode of CELLPHY
    CellPhyDataTrueSNVsVCFEPFile[MAX_NAME],
    CellPhyDataTrueSNVsVCFEPFileCopy[MAX_NAME],
    CellPhyDataCalledSNVsVCFEPFile[MAX_NAME],
    CellPhyDataCalledSNVsVCFEPFileCopy[MAX_NAME],
    CellPhyCMDTrueSNVsVCFEPFile[MAX_NAME],
    CellPhyCMDCalledSNVsVCFEPFile[MAX_NAME];

double CijkMD[15 * 15 * 15];
double CijkMDI[65 * 65 * 65];
double RootMD[15];
double RootMDI[65];
extern double QijMD[15 * 15];
extern double QijMDI[65 * 65];
static const int candidateMLGenotypesMD[NUM_CAN_ML_GENOTYPES] = {0, 1, 2, 3, 4,
                                                                 5, 6, 7, 8, 9};
static const int candidateMLGenotypesMDI[NUM_CAN_ML_GENOTYPES] = {
    40, 42, 43, 44, 48, 49, 50, 51, 52, 53};
static const char genotypesInCellPhy[NUM_CAN_ML_GENOTYPES] = {
    'A', 'M', 'R', 'W', 'C', 'S', 'Y', 'G', 'K', 'T'};
static double insertionRate = 0.0;
struct HashMapIntToIntArr *genotypeHashMapIntToIntArr;
struct HashMapIntArrToInt *genotypeHashMapIntArrToInt;
int genotypeArrMD[15][3];
int genotypeArrMDI[65][3];
int numINS;
double cumNumINS, meanNumINS, varNumINS, cumNumINSSq;
int *mutatedSites, *insertedSites, *deletedSites;
int *numPureVariants; /* nVariants, nMutations, nDeletions, nInsertions*/
CellStrExtd *cellExtd;
int *backgroundSize; /* coverage (c), reads for MLVariantAllele (m), c - m */
int **backgroundInfo;
int *maxBackgroundSize;
long **backgroundInfoAll; /* store total background information for more than 1
                             data sets of a fixed tree */
long completeDeletions, completeADOs, noCoverage, noReadsPositiveCoverage,
    zeroCoverage, zeroProportionVAF, completeDeletionsAll, completeADOsAll,
    noCoverageAll, noReadsPositiveCoverageAll, zeroCoverageAll,
    zeroProportionVAFAll;
/* for SNV sites
 * 1. the index of each row is the number of existing alleles
 * 2. the fist column is the total number of zero coverage occurrences of a
 * specific number of existing alleles for the current data set (for different
 * unfixed trees) or for all data sets (for a fixed tree)
 * 3. the second column has the same meaning as the first column does, only
 * those sites counted with no reads (representing non-coverage) */
int allelicZeroCoverageProportion[4][2], allelicZeroCoverageProportionAll[4][2];
/* the probability of a successful trial for SNV sites
 * 1. the first dimension is the number of existing alleles (0, 1, 2, 3).
 * 2. the second dimension contains the minimum (0) and the maximum (1) number
 * of probability or number of successful trials, and the corresponding mean (2)
 * as well as the variance (3). */
double probNBInfo[4][4], probNBInfoAll[4][4];
/* details of the probability of a successful trial for SNV sites
 * 1. the first dimension is the number of existing alleles (0, 1, 2, 3).
 * 2. the second dimension is the concrete value of probability for each site in
 * each cell. The space will be dynamically allocated at run time. */
double *probNBDetails[4], *probNBDetailsAll[4];

/*********************
 *      Methods      *
 *********************/

void AllocateTreeStructure(int numNodes) {
  treeStructure = (unsigned **)calloc(numNodes + 1, sizeof(unsigned *));
  for (int i = 0; i < numNodes; i++) {
    treeStructure[i] = (unsigned *)calloc(numNodes + 1, sizeof(unsigned));
  }
} // AllocateTreeStructure

void InitTreeStructure(int numNodes) {
  for (int i = 0; i < numNodes; i++) {
    for (int j = 0; j < numNodes; j++) {
      if (i == j) {
        treeStructure[i][i] = 1;
      } else {
        treeStructure[i][j] = 0;
      }
    }
  }
} // InitTreeStructure

void SetTreeStructure(const TreeNode *node, int numNodes) {
  if (node == NULL)
    return;

  SetTreeStructure(node->left, numNodes);
  SetTreeStructure(node->right, numNodes);

  if (node->left != NULL) {
    for (int i = 0; i < numNodes; i++) {
      if (treeStructure[node->left->label][i] > 0) {
        treeStructure[node->label][i] +=
            treeStructure[node->left->label][i] + 1;
      }
    }
  }

  if (node->right != NULL) {
    for (int i = 0; i < numNodes; i++) {
      if (treeStructure[node->right->label][i] > 0) {
        treeStructure[node->label][i] +=
            treeStructure[node->right->label][i] + 1;
      }
    }
  }
} // SetTreeStructure

void PrintTreeStructure(int numNodes) {
  fprintf(stdout, "\n |");
  for (int i = 0; i < numNodes; i++) {
    fprintf(stdout, "%d ", i);
  }
  fprintf(stdout, "\n");

  char sep[2 * numNodes + 2];
  memset(sep, '_', sizeof(sep) + 1);
  sep[2 * numNodes + 1] = '\0';
  fprintf(stdout, "%s\n", sep);

  for (int i = 0; i < numNodes; i++) {
    fprintf(stdout, "%d|", i);
    for (int j = 0; j < numNodes; j++) {
      fprintf(stdout, "%d ", treeStructure[i][j]);
    }
    fprintf(stdout, "\n");
  }
} // PrintTreeStructure

void FreeTreeStructure(int numNodes) {
  for (int i = 0; i < numNodes; i++) {
    free(treeStructure[i]);
    treeStructure[i] = NULL;
  }
  free(treeStructure);
  treeStructure = NULL;
} // AllocateTreeStructure

void SetReadLengths(int cellNr, int siteNr) {
  if (siteNr % readLengths[1] == 0) {
    readLengths[0] = readLengths[1];
    numReadBlocks = siteNr / readLengths[0];
  } else {
    int readsNr = (int)(((double)siteNr) / readLengths[1] - 1);
    readLengths[0] = (siteNr - readsNr * readLengths[1]) / 2;
    numReadBlocks = readsNr + 2;
  }

  startSiteBlocks = (int *)malloc((numReadBlocks + 1) * sizeof(int));
  startSiteBlocks[0] = 0;
  for (int i = 0; i < numReadBlocks; i++) {
    startSiteBlocks[i + 1] =
        startSiteBlocks[i] +
        (i == 0 || i == numReadBlocks - 1 ? readLengths[0] : readLengths[1]);
  }

  noReadProbBlocks = (double **)malloc((cellNr + 2) * sizeof(double *));
  for (int i = 0; i < cellNr + 1; i++) {
    noReadProbBlocks[i] =
        (double *)malloc((numReadBlocks + 1) * sizeof(double));
    memset(noReadProbBlocks[i], noReadProb[0], numReadBlocks);
  }
}

/*
 * Initialize the hash map store the information of mapping codes to genotypes.
 * For example, 0 -> A/A
 * */
void InitGenotypeHashMapMD() {
  genotypeHashMapIntToIntArr = CreateHashMapIntToIntArr(3);
  genotypeHashMapIntArrToInt = CreateHashMapIntArrToInt(3);

  for (int i = 0; i < 15; ++i) {
    InsertKeyValIntToIntArr(genotypeHashMapIntToIntArr, i, genotypeArrMD[i]);
    InsertKeyValIntArrToInt(genotypeHashMapIntArrToInt, genotypeArrMD[i], i);
  }
}

/*
 * mutation, deletion
 * define the substitution rate matrix internally.
 * users can configure the mutation rate, deletion rate and as
 * they wish. mutationRate is defined as the rate of a reference nucleotide
 * mutating to any one of the three alternative nucleotides.
 */
void InitRateMatrixMD(double mutationRate, double deletionRate) {
  // normalize the deletionRate and insertionRate regarding to the mutationRate
  deletionRate /= (mutationRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < 15; ++i)
    for (int j = 0; j < 15; ++j)
      MijMD[i][j] = 0.0;

  // 0    1    2    3    4    5    6    7    8    9    10   11   12   13   14
  // A/A, A/C, A/G, A/T, C/C, C/G, C/T, G/G, G/T, T/T, A/-, C/-, G/-, T/-, -

  // 0, A/A
  MijMD[0][1] = MijMD[0][2] = MijMD[0][3] = 1.0 / 3;
  MijMD[0][10] = deletionRate;

  genotypeArrMD[0][0] = 0;
  genotypeArrMD[0][1] = 0;
  genotypeArrMD[0][2] = INVALID;

  // 1, A/C
  MijMD[1][0] = MijMD[1][2] = MijMD[1][3] = MijMD[1][4] = MijMD[1][5] =
      MijMD[1][6] = 1.0 / 6;
  MijMD[1][10] = MijMD[1][11] = deletionRate / 2;

  genotypeArrMD[1][0] = 0;
  genotypeArrMD[1][1] = 1;
  genotypeArrMD[1][2] = INVALID;

  // 2, A/G
  MijMD[2][0] = MijMD[2][1] = MijMD[2][3] = MijMD[2][5] = MijMD[2][7] =
      MijMD[2][8] = 1.0 / 6;
  MijMD[2][10] = MijMD[2][12] = deletionRate / 2;

  genotypeArrMD[2][0] = 0;
  genotypeArrMD[2][1] = 2;
  genotypeArrMD[2][2] = INVALID;

  // 3, A/T
  MijMD[3][0] = MijMD[3][1] = MijMD[3][2] = MijMD[3][6] = MijMD[3][8] =
      MijMD[3][9] = 1.0 / 6;
  MijMD[3][10] = MijMD[3][13] = deletionRate / 2;

  genotypeArrMD[3][0] = 0;
  genotypeArrMD[3][1] = 3;
  genotypeArrMD[3][2] = INVALID;

  // 4, C/C
  MijMD[4][1] = MijMD[4][5] = MijMD[4][6] = 1.0 / 3;
  MijMD[4][11] = deletionRate;

  genotypeArrMD[4][0] = 1;
  genotypeArrMD[4][1] = 1;
  genotypeArrMD[4][2] = INVALID;

  // 5, C/G
  MijMD[5][1] = MijMD[5][2] = MijMD[5][4] = MijMD[5][6] = MijMD[5][7] =
      MijMD[5][8] = 1.0 / 6;
  MijMD[5][11] = MijMD[5][12] = deletionRate / 2;

  genotypeArrMD[5][0] = 1;
  genotypeArrMD[5][1] = 2;
  genotypeArrMD[5][2] = INVALID;

  // 6, C/T
  MijMD[6][1] = MijMD[6][3] = MijMD[6][4] = MijMD[6][5] = MijMD[6][8] =
      MijMD[6][9] = 1.0 / 6;
  MijMD[6][11] = MijMD[6][13] = deletionRate / 2;

  genotypeArrMD[6][0] = 1;
  genotypeArrMD[6][1] = 3;
  genotypeArrMD[6][2] = INVALID;

  // 7, G/G
  MijMD[7][2] = MijMD[7][5] = MijMD[7][8] = 1.0 / 3;
  MijMD[7][12] = deletionRate;

  genotypeArrMD[7][0] = 2;
  genotypeArrMD[7][1] = 2;
  genotypeArrMD[7][2] = INVALID;

  // 8, G/T
  MijMD[8][2] = MijMD[8][3] = MijMD[8][5] = MijMD[8][6] = MijMD[8][7] =
      MijMD[8][9] = 1.0 / 6;
  MijMD[8][12] = MijMD[8][13] = deletionRate / 2;

  genotypeArrMD[8][0] = 2;
  genotypeArrMD[8][1] = 3;
  genotypeArrMD[8][2] = INVALID;

  // 9, T/T
  MijMD[9][3] = MijMD[9][6] = MijMD[9][8] = 1.0 / 3;
  MijMD[9][13] = deletionRate;

  genotypeArrMD[9][0] = 3;
  genotypeArrMD[9][1] = 3;
  genotypeArrMD[9][2] = INVALID;

  // 10, A/-
  MijMD[10][11] = MijMD[10][12] = MijMD[10][13] = 1.0 / 6;
  MijMD[10][14] = deletionRate / 2;

  genotypeArrMD[10][0] = 0;
  genotypeArrMD[10][1] = DELETED;
  genotypeArrMD[10][2] = INVALID;

  // 11, C/-
  MijMD[11][10] = MijMD[11][12] = MijMD[11][13] = 1.0 / 6;
  MijMD[11][14] = deletionRate / 2;

  genotypeArrMD[11][0] = 1;
  genotypeArrMD[11][1] = DELETED;
  genotypeArrMD[11][2] = INVALID;

  // 12, G/-
  MijMD[12][10] = MijMD[12][11] = MijMD[12][13] = 1.0 / 6;
  MijMD[12][14] = deletionRate / 2;

  genotypeArrMD[12][0] = 2;
  genotypeArrMD[12][1] = DELETED;
  genotypeArrMD[12][2] = INVALID;

  // 13, T/-
  MijMD[13][10] = MijMD[13][11] = MijMD[13][12] = 1.0 / 6;
  MijMD[13][14] = deletionRate / 2;

  genotypeArrMD[13][0] = 3;
  genotypeArrMD[13][1] = DELETED;
  genotypeArrMD[13][2] = INVALID;

  // 14, -
  genotypeArrMD[14][0] = genotypeArrMD[14][1] = DELETED;
  genotypeArrMD[14][2] = INVALID;

  InitGenotypeHashMapMD();
}

/*
 * Initialize the hash map store the information of mapping codes to genotypes.
 * For example, 0 -> A/AA
 * */
void InitGenotypeHashMapMDI() {
  genotypeHashMapIntToIntArr = CreateHashMapIntToIntArr(13);
  genotypeHashMapIntArrToInt = CreateHashMapIntArrToInt(13);

  for (int i = 0; i < 65; ++i) {
    InsertKeyValIntToIntArr(genotypeHashMapIntToIntArr, i, genotypeArrMDI[i]);
    InsertKeyValIntArrToInt(genotypeHashMapIntArrToInt, genotypeArrMDI[i], i);
  }
}

/*
 * mutation, deletion, insertion
 * define the substitution rate matrix internally.
 * users can configure the mutation rate, deletion rate and insertion rate as
 * they wish. mutationRate is defined as the rate of a reference nucleotide
 * mutating to any one of the three alternative nucleotides.
 */
/*
void initRateMatrix() {
  // diagonal terms
  MijMDI[0][0] = MijMDI[1][1] = MijMDI[2][2] = MijMDI[3][3] =
MijMDI[4][4] = MijMDI[5][5] = MijMDI[6][6] = MijMDI[7][7] =
MijMDI[8][8] = MijMDI[9][9] = MijMDI[10][10] = MijMDI[11][11] =
MijMDI[12][12] = MijMDI[13][13] = MijMDI[14][14] = 0.0;
  // 0/00
  MijMDI[0][1] = 2/3.0; MijMDI[0][2] = 1/3.0; MijMDI[0][3] =
MijMDI[0][4] = MijMDI[0][5] = 0.0; MijMDI[0][6] = 2 * deletionRate
/ 3.0; MijMDI[0][7] = deletionRate / 3.0; MijMDI[0][8] =
MijMDI[0][9] = MijMDI[0][10] = MijMDI[0][11] = MijMDI[0][12] =
MijMDI[0][13] = MijMDI[0][14] =0.0;
  // 0/01
  MijMDI[1][0] = 1/9.0; MijMDI[1][2] = 0.0; MijMDI[1][3] =
MijMDI[1][4] = 1/3.0; MijMDI[1][5] = 0.0; MijMDI[1][6] =
deletionRate / 3.0; MijMDI[1][7] = 0.0; MijMDI[1][8] = MijMDI[1][8]
= deletionRate / 3.0; MijMDI[1][10] = MijMDI[1][11] = MijMDI[1][12]
= MijMDI[1][13] = MijMDI[1][14] = 0.0;
  // 1/00
  MijMDI[2][0] = 1/9.0; MijMDI[2][1] = MijMDI[2][3] = 0.0;
MijMDI[2][4] = 2/3.0; MijMDI[2][5] = MijMDI[2][6] = 0.0;
MijMDI[2][7] = deletionRate / 3.0; MijMDI[2][8] = 2 * deletionRate
/ 3.0; MijMDI[2][9] = MijMDI[2][10] = MijMDI[2][11] =
MijMDI[2][12] = MijMDI[2][13] = MijMDI[2][14] = 0.0;
  // 0/11
  MijMDI[3][0] = 0.0; MijMDI[3][1] = 2/9.0; MijMDI[3][2] = 0.0;
MijMDI[3][4] = 0.0; MijMDI[3][5] = 1/3.0; MijMDI[3][6] =
MijMDI[3][7] = 0.0; MijMDI[3][8] = 2 * deletionRate / 3.0;
MijMDI[3][9] =  MijMDI[3][10] = 0.0; MijMDI[3][11] = deletionRate
/ 3.0; MijMDI[3][12] = MijMDI[3][13] = MijMDI[3][14] = 0.0;
  // 1/01
  MijMDI[4][0] = 0.0; MijMDI[4][1] = MijMDI[4][2] = 1/9.0;
MijMDI[4][3] = 0.0; MijMDI[4][5] = 1/3.0; MijMDI[4][6] =
MijMDI[4][7] = 0.0; MijMDI[4][8] = MijMDI[4][9] = MijMDI[4][10]
= deletionRate / 3.0; MijMDI[4][11] = MijMDI[4][12] = MijMDI[4][13]
= MijMDI[4][14] = 0.0;
  // 1/11
  MijMDI[5][0] = MijMDI[5][1] = MijMDI[5][2] = 0.0; MijMDI[5][3]
= 1/9.0; MijMDI[5][4] = 2/9.0; MijMDI[5][6] = MijMDI[5][7] =
MijMDI[5][8] = MijMDI[5][9] = 0.0; MijMDI[5][10] = 2 * deletionRate
/ 3.0; MijMDI[5][11] = deletionRate / 3.0; MijMDI[5][12] =
MijMDI[5][13] = MijMDI[5][14] = 0.0;
  // 0/0
  MijMDI[6][0] = 2 * insertionRate / 3.0; MijMDI[6][1] =
MijMDI[6][2] = MijMDI[6][3] = MijMDI[6][4] = MijMDI[6][5] =
MijMDI[6][7] = 0.0; MijMDI[6][8] = 2/3.0; MijMDI[6][9] =
MijMDI[6][10] = MijMDI[6][11] = 0.0; MijMDI[6][12] = 2 *
deletionRate / 3.0; MijMDI[6][13] = MijMDI[6][14] = 0.0;
  // -/00
  MijMDI[7][0] = MijMDI[7][1] = MijMDI[7][2] = MijMDI[7][3] =
MijMDI[7][4] = MijMDI[7][5] = MijMDI[7][6] = MijMDI[7][8] = 0.0;
MijMDI[7][9] = 2/3.0; MijMDI[7][10] = MijMDI[7][11] = 0.0;
MijMDI[7][12] = 2 * deletionRate / 3.0; MijMDI[7][13] =
MijMDI[7][14] = 0.0;
  // 0/1
  MijMDI[8][0] = MijMDI[8][1] = 0.0; MijMDI[8][2] = MijMDI[8][3]
= insertionRate / 3.0; MijMDI[8][4] = MijMDI[8][5] = 0.0;
MijMDI[8][6] = 1/9.0; MijMDI[8][7] = MijMDI[8][9] = 0.0;
MijMDI[8][10] = 1/3.0; MijMDI[8][11] = 0.0; MijMDI[8][12] =
MijMDI[8][13] = deletionRate / 3.0; MijMDI[8][14] = 0.0;
  // -/01
  MijMDI[9][0] = MijMDI[9][1] = MijMDI[9][2] = MijMDI[9][3] =
MijMDI[9][4] = MijMDI[9][5] = MijMDI[9][6] = 0.0; MijMDI[9][7] =
1/9.0; MijMDI[9][8] = MijMDI[9][10] = 0.0; MijMDI[9][11] =1/3.0;
MijMDI[9][12] = MijMDI[9][13] = deletionRate / 3.0; MijMDI[9][14] =
0.0;
  // 1/1
  MijMDI[10][0] = MijMDI[10][1] = MijMDI[10][2] = MijMDI[10][3]
= MijMDI[10][4] = 0.0; MijMDI[10][5] = 2 * insertionRate / 3.0;
MijMDI[10][6] = MijMDI[10][7] = 0.0; MijMDI[10][8] = 2/9.0;
MijMDI[10][9] = MijMDI[10][11] = MijMDI[10][12] = 0.0;
MijMDI[10][13] = 2 * deletionRate / 3.0; MijMDI[10][14] = 0.0;
  // -/11
  MijMDI[11][0] = MijMDI[11][1] = MijMDI[11][2] = MijMDI[11][3]
= MijMDI[11][4] = MijMDI[11][5] = MijMDI[11][6] = MijMDI[11][7]
= MijMDI[11][8] = 0.0; MijMDI[11][9] = 2/9.0; MijMDI[11][10] =
MijMDI[11][12] = 0.0; MijMDI[11][13] = 2 * deletionRate / 3.0;
MijMDI[11][14] = 0.0;
  // 0/-
  MijMDI[12][0] = MijMDI[12][1] = MijMDI[12][2] = MijMDI[12][3]
= MijMDI[12][4] = MijMDI[12][5] = MijMDI[12][6] = 0.0;
MijMDI[12][7] = insertionRate / 3.0; MijMDI[12][8] = MijMDI[12][9] =
MijMDI[12][10] = MijMDI[12][11] = 0.0; MijMDI[12][13] = 1/3.0;
MijMDI[12][14] = deletionRate / 3.0;
  // 1/-
  MijMDI[13][0] = MijMDI[13][1] = MijMDI[13][2] = MijMDI[13][3]
= MijMDI[13][4] = MijMDI[13][5] = MijMDI[13][6] = MijMDI[13][7]
= MijMDI[13][8] = MijMDI[13][9] = MijMDI[13][10] = 0.0;
MijMDI[13][11] = insertionRate / 3.0; MijMDI[13][12] = 1/9.0;
MijMDI[13][14] = deletionRate / 3.0;
  // -
  MijMDI[14][0] = MijMDI[14][1] = MijMDI[14][2] = MijMDI[14][3]
= MijMDI[14][4] = MijMDI[14][5] = MijMDI[14][6] = MijMDI[14][7]
= MijMDI[14][8] = MijMDI[14][9] = MijMDI[14][10] =
MijMDI[14][11] = MijMDI[14][12] = MijMDI[14][13] = 0.0;
}
*/
void InitRateMatrixMDI(double mutationRate, double deletionRate,
                       double insertionRate) {
  // normalize the deletionRate and insertionRate regarding to the mutationRate
  deletionRate /= (mutationRate + 1e-15);
  insertionRate /= (mutationRate + 1e-15);

  // initialize all entries to 0.0
  for (int i = 0; i < 65; ++i)
    for (int j = 0; j < 65; ++j)
      MijMDI[i][j] = 0.0;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 0, A/AA
  MijMDI[0][1] = MijMDI[0][2] = MijMDI[0][3] = 2 / 9.0;
  MijMDI[0][4] = MijMDI[0][5] = MijMDI[0][6] = 1 / 9.0;
  MijMDI[0][40] = 2 * deletionRate / 3.0;
  MijMDI[0][41] = deletionRate / 3.0;

  genotypeArrMDI[0][0] = 0;
  genotypeArrMDI[0][1] = 0;
  genotypeArrMDI[0][2] = 0;

  // 1, A/AC
  MijMDI[1][0] = MijMDI[1][2] = MijMDI[1][3] = MijMDI[1][7] = MijMDI[1][8] =
      MijMDI[1][9] = MijMDI[1][13] = MijMDI[1][16] = MijMDI[1][19] = 1 / 9.0;
  MijMDI[1][40] = MijMDI[1][42] = MijMDI[1][45] = deletionRate / 3.0;

  genotypeArrMDI[1][0] = 0;
  genotypeArrMDI[1][1] = 0;
  genotypeArrMDI[1][2] = 1;

  // 2, A/AG
  MijMDI[2][0] = MijMDI[2][1] = MijMDI[2][3] = MijMDI[2][8] = MijMDI[2][10] =
      MijMDI[2][11] = MijMDI[2][14] = MijMDI[2][17] = MijMDI[2][20] = 1 / 9.0;
  MijMDI[2][40] = MijMDI[2][43] = MijMDI[2][46] = deletionRate / 3.0;

  genotypeArrMDI[2][0] = 0;
  genotypeArrMDI[2][1] = 0;
  genotypeArrMDI[2][2] = 2;

  // 3, A/AT
  MijMDI[3][0] = MijMDI[3][1] = MijMDI[3][2] = MijMDI[3][9] = MijMDI[3][11] =
      MijMDI[3][12] = MijMDI[3][15] = MijMDI[3][18] = MijMDI[3][21] = 1 / 9.0;
  MijMDI[3][40] = MijMDI[3][44] = MijMDI[3][47] = deletionRate / 3.0;

  genotypeArrMDI[3][0] = 0;
  genotypeArrMDI[3][1] = 0;
  genotypeArrMDI[3][2] = 3;

  // 4, C/AA
  MijMDI[4][0] = MijMDI[4][5] = MijMDI[4][6] = 1 / 9.0;
  MijMDI[4][13] = MijMDI[4][14] = MijMDI[4][15] = 2 / 9.0;
  MijMDI[4][41] = deletionRate / 3.0;
  MijMDI[4][42] = 2 * deletionRate / 3.0;

  genotypeArrMDI[4][0] = 1;
  genotypeArrMDI[4][1] = 0;
  genotypeArrMDI[4][2] = 0;

  // 5, G/AA
  MijMDI[5][0] = MijMDI[5][4] = MijMDI[5][6] = 1 / 9.0;
  MijMDI[5][16] = MijMDI[5][17] = MijMDI[5][18] = 2 / 9.0;
  MijMDI[5][41] = deletionRate / 3.0;
  MijMDI[5][43] = 2 * deletionRate / 3.0;

  genotypeArrMDI[5][0] = 2;
  genotypeArrMDI[5][1] = 0;
  genotypeArrMDI[5][2] = 0;

  // 6, T/AA
  MijMDI[6][0] = MijMDI[6][4] = MijMDI[6][5] = 1 / 9.0;
  MijMDI[6][19] = MijMDI[6][20] = MijMDI[6][21] = 2 / 9.0;
  MijMDI[6][41] = deletionRate / 3.0;
  MijMDI[6][44] = 2 * deletionRate / 3.0;

  genotypeArrMDI[6][0] = 3;
  genotypeArrMDI[6][1] = 0;
  genotypeArrMDI[6][2] = 0;

  // 7, A/CC
  MijMDI[7][1] = MijMDI[7][8] = MijMDI[7][9] = 2 / 9.0;
  MijMDI[7][22] = MijMDI[7][28] = MijMDI[7][34] = 1 / 9.0;
  MijMDI[7][42] = 2 * deletionRate / 3.0;
  MijMDI[7][54] = deletionRate / 3.0;

  genotypeArrMDI[7][0] = 0;
  genotypeArrMDI[7][1] = 1;
  genotypeArrMDI[7][2] = 1;

  // 8, A/CG
  MijMDI[8][1] = MijMDI[8][2] = MijMDI[8][7] = MijMDI[8][9] = MijMDI[8][10] =
      MijMDI[8][11] = MijMDI[8][23] = MijMDI[8][29] = MijMDI[8][35] = 1 / 9.0;
  MijMDI[8][42] = MijMDI[8][43] = MijMDI[8][55] = deletionRate / 3.0;

  genotypeArrMDI[8][0] = 0;
  genotypeArrMDI[8][1] = 1;
  genotypeArrMDI[8][2] = 2;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 9, A/CT
  MijMDI[9][1] = MijMDI[9][3] = MijMDI[9][7] = MijMDI[9][8] = MijMDI[9][11] =
      MijMDI[9][12] = MijMDI[9][24] = MijMDI[9][30] = MijMDI[9][36] = 1 / 9.0;
  MijMDI[9][42] = MijMDI[9][44] = MijMDI[9][56] = deletionRate / 3.0;

  genotypeArrMDI[9][0] = 0;
  genotypeArrMDI[9][1] = 1;
  genotypeArrMDI[9][2] = 3;

  // 10, A/GG
  MijMDI[10][2] = MijMDI[10][8] = MijMDI[10][11] = 2 / 9.0;
  MijMDI[10][25] = MijMDI[10][31] = MijMDI[10][37] = 1 / 9.0;
  MijMDI[10][43] = 2 * deletionRate / 3.0;
  MijMDI[10][57] = deletionRate / 3.0;

  genotypeArrMDI[10][0] = 0;
  genotypeArrMDI[10][1] = 2;
  genotypeArrMDI[10][2] = 2;

  // 11, A/GT
  MijMDI[11][2] = MijMDI[11][3] = MijMDI[11][8] = MijMDI[11][9] =
      MijMDI[11][10] = MijMDI[11][12] = MijMDI[11][26] = MijMDI[11][32] =
          MijMDI[11][38] = 1 / 9.0;
  MijMDI[11][43] = MijMDI[11][44] = MijMDI[11][58] = deletionRate / 3.0;

  genotypeArrMDI[11][0] = 0;
  genotypeArrMDI[11][1] = 2;
  genotypeArrMDI[11][2] = 3;

  // 12, A/TT
  MijMDI[12][3] = MijMDI[12][9] = MijMDI[12][11] = 2 / 9.0;
  MijMDI[12][27] = MijMDI[12][33] = MijMDI[12][39] = 1 / 9.0;
  MijMDI[12][44] = 2 * deletionRate / 3.0;
  MijMDI[12][59] = deletionRate / 3.0;

  genotypeArrMDI[12][0] = 0;
  genotypeArrMDI[12][1] = 3;
  genotypeArrMDI[12][2] = 3;

  // 13, C/AC
  MijMDI[13][1] = MijMDI[13][4] = MijMDI[13][14] = MijMDI[13][15] =
      MijMDI[13][16] = MijMDI[13][19] = MijMDI[13][22] = MijMDI[13][23] =
          MijMDI[13][24] = 1 / 9.0;
  MijMDI[13][42] = MijMDI[13][45] = MijMDI[13][48] = deletionRate / 3.0;

  genotypeArrMDI[13][0] = 1;
  genotypeArrMDI[13][1] = 0;
  genotypeArrMDI[13][2] = 1;

  // 14, C/AG
  MijMDI[14][2] = MijMDI[14][4] = MijMDI[14][13] = MijMDI[14][15] =
      MijMDI[14][17] = MijMDI[14][20] = MijMDI[14][23] = MijMDI[14][25] =
          MijMDI[14][26] = 1 / 9.0;
  MijMDI[14][42] = MijMDI[14][46] = MijMDI[14][49] = deletionRate / 3.0;

  genotypeArrMDI[14][0] = 1;
  genotypeArrMDI[14][1] = 0;
  genotypeArrMDI[14][2] = 2;

  // 15, C/AT
  MijMDI[15][3] = MijMDI[15][4] = MijMDI[15][13] = MijMDI[15][14] =
      MijMDI[15][18] = MijMDI[15][21] = MijMDI[15][24] = MijMDI[15][26] =
          MijMDI[15][27] = 1 / 9.0;
  MijMDI[15][42] = MijMDI[15][47] = MijMDI[15][50] = deletionRate / 3.0;

  genotypeArrMDI[15][0] = 1;
  genotypeArrMDI[15][1] = 0;
  genotypeArrMDI[15][2] = 3;

  // 16, G/AC
  MijMDI[16][1] = MijMDI[16][5] = MijMDI[16][13] = MijMDI[16][17] =
      MijMDI[16][18] = MijMDI[16][19] = MijMDI[16][28] = MijMDI[16][29] =
          MijMDI[16][30] = 1 / 9.0;
  MijMDI[16][43] = MijMDI[16][45] = MijMDI[16][49] = deletionRate / 3.0;

  genotypeArrMDI[16][0] = 2;
  genotypeArrMDI[16][1] = 0;
  genotypeArrMDI[16][2] = 1;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 17, G/AG
  MijMDI[17][2] = MijMDI[17][5] = MijMDI[17][14] = MijMDI[17][16] =
      MijMDI[17][18] = MijMDI[17][20] = MijMDI[17][29] = MijMDI[17][31] =
          MijMDI[17][32] = 1 / 9.0;
  MijMDI[17][43] = MijMDI[17][46] = MijMDI[17][51] = deletionRate / 3.0;

  genotypeArrMDI[17][0] = 2;
  genotypeArrMDI[17][1] = 0;
  genotypeArrMDI[17][2] = 2;

  // 18, G/AT
  MijMDI[18][3] = MijMDI[18][5] = MijMDI[18][15] = MijMDI[18][16] =
      MijMDI[18][17] = MijMDI[18][21] = MijMDI[18][30] = MijMDI[18][32] =
          MijMDI[18][33] = 1 / 9.0;
  MijMDI[18][43] = MijMDI[18][47] = MijMDI[18][52] = deletionRate / 3.0;

  genotypeArrMDI[18][0] = 2;
  genotypeArrMDI[18][1] = 0;
  genotypeArrMDI[18][2] = 3;

  // 19, T/AC
  MijMDI[19][1] = MijMDI[19][6] = MijMDI[19][13] = MijMDI[19][16] =
      MijMDI[19][20] = MijMDI[19][21] = MijMDI[19][34] = MijMDI[19][35] =
          MijMDI[19][36] = 1 / 9.0;
  MijMDI[19][44] = MijMDI[19][45] = MijMDI[19][50] = deletionRate / 3.0;

  genotypeArrMDI[19][0] = 3;
  genotypeArrMDI[19][1] = 0;
  genotypeArrMDI[19][2] = 1;

  // 20, T/AG
  MijMDI[20][2] = MijMDI[20][6] = MijMDI[20][14] = MijMDI[20][17] =
      MijMDI[20][19] = MijMDI[20][21] = MijMDI[20][35] = MijMDI[20][37] =
          MijMDI[20][38] = 1 / 9.0;
  MijMDI[20][44] = MijMDI[20][46] = MijMDI[20][52] = deletionRate / 3.0;

  genotypeArrMDI[20][0] = 3;
  genotypeArrMDI[20][1] = 0;
  genotypeArrMDI[20][2] = 2;

  // 21, T/AT
  MijMDI[21][3] = MijMDI[21][6] = MijMDI[21][15] = MijMDI[21][18] =
      MijMDI[21][19] = MijMDI[21][20] = MijMDI[21][36] = MijMDI[21][38] =
          MijMDI[21][39] = 1 / 9.0;
  MijMDI[21][44] = MijMDI[21][47] = MijMDI[21][53] = deletionRate / 3.0;

  genotypeArrMDI[21][0] = 3;
  genotypeArrMDI[21][1] = 0;
  genotypeArrMDI[21][2] = 3;

  // 22, C/CC
  MijMDI[22][7] = MijMDI[22][28] = MijMDI[22][34] = 1 / 9.0;
  MijMDI[22][13] = MijMDI[22][23] = MijMDI[22][24] = 2 / 9.0;
  MijMDI[22][48] = 2 * deletionRate / 3.0;
  MijMDI[22][54] = deletionRate / 3.0;

  genotypeArrMDI[22][0] = 1;
  genotypeArrMDI[22][1] = 1;
  genotypeArrMDI[22][2] = 1;

  // 23, C/CG
  MijMDI[23][8] = MijMDI[23][13] = MijMDI[23][14] = MijMDI[23][22] =
      MijMDI[23][24] = MijMDI[23][25] = MijMDI[23][26] = MijMDI[23][29] =
          MijMDI[23][35] = 1 / 9.0;
  MijMDI[23][48] = MijMDI[23][49] = MijMDI[23][55] = deletionRate / 3.0;

  genotypeArrMDI[23][0] = 1;
  genotypeArrMDI[23][1] = 1;
  genotypeArrMDI[23][2] = 2;

  // 24, C/CT
  MijMDI[24][9] = MijMDI[24][13] = MijMDI[24][15] = MijMDI[24][22] =
      MijMDI[24][23] = MijMDI[24][26] = MijMDI[24][27] = MijMDI[24][30] =
          MijMDI[24][36] = 1 / 9.0;
  MijMDI[24][48] = MijMDI[24][50] = MijMDI[24][56] = deletionRate / 3.0;

  genotypeArrMDI[24][0] = 1;
  genotypeArrMDI[24][1] = 1;
  genotypeArrMDI[24][2] = 3;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 25, C/GG
  MijMDI[25][10] = MijMDI[25][31] = MijMDI[25][37] = 1 / 9.0;
  MijMDI[25][14] = MijMDI[25][23] = MijMDI[25][26] = 2 / 9.0;
  MijMDI[25][49] = 2 * deletionRate / 3.0;
  MijMDI[25][57] = deletionRate / 3.0;

  genotypeArrMDI[25][0] = 1;
  genotypeArrMDI[25][1] = 2;
  genotypeArrMDI[25][2] = 2;

  // 26, C/GT
  MijMDI[26][11] = MijMDI[26][14] = MijMDI[26][15] = MijMDI[26][23] =
      MijMDI[26][24] = MijMDI[26][25] = MijMDI[26][27] = MijMDI[26][32] =
          MijMDI[26][38] = 1 / 9.0;
  MijMDI[26][49] = MijMDI[26][50] = MijMDI[26][58] = deletionRate / 3.0;

  genotypeArrMDI[26][0] = 1;
  genotypeArrMDI[26][1] = 2;
  genotypeArrMDI[26][2] = 3;

  // 27, C/TT
  MijMDI[27][12] = MijMDI[27][33] = MijMDI[27][39] = 1 / 9.0;
  MijMDI[27][15] = MijMDI[27][24] = MijMDI[27][26] = 2 / 9.0;
  MijMDI[27][50] = 2 * deletionRate / 3.0;
  MijMDI[27][59] = deletionRate / 3.0;

  genotypeArrMDI[27][0] = 1;
  genotypeArrMDI[27][1] = 3;
  genotypeArrMDI[27][2] = 3;

  // 28, G/CC
  MijMDI[28][7] = MijMDI[28][22] = MijMDI[28][34] = 1 / 9.0;
  MijMDI[28][16] = MijMDI[28][29] = MijMDI[28][30] = 2 / 9.0;
  MijMDI[28][49] = 2 * deletionRate / 3.0;
  MijMDI[28][54] = deletionRate / 3.0;

  genotypeArrMDI[28][0] = 2;
  genotypeArrMDI[28][1] = 1;
  genotypeArrMDI[28][2] = 1;

  // 29, G/CG
  MijMDI[29][8] = MijMDI[29][16] = MijMDI[29][17] = MijMDI[29][23] =
      MijMDI[29][28] = MijMDI[29][30] = MijMDI[29][31] = MijMDI[29][32] =
          MijMDI[29][35] = 1 / 9.0;
  MijMDI[29][49] = MijMDI[29][51] = MijMDI[29][55] = deletionRate / 3.0;

  genotypeArrMDI[29][0] = 2;
  genotypeArrMDI[29][1] = 1;
  genotypeArrMDI[29][2] = 2;

  // 30, G/CT
  MijMDI[30][9] = MijMDI[30][16] = MijMDI[30][18] = MijMDI[30][24] =
      MijMDI[30][28] = MijMDI[30][29] = MijMDI[30][32] = MijMDI[30][33] =
          MijMDI[30][36] = 1 / 9.0;
  MijMDI[30][49] = MijMDI[30][52] = MijMDI[30][56] = deletionRate / 3.0;

  genotypeArrMDI[30][0] = 2;
  genotypeArrMDI[30][1] = 1;
  genotypeArrMDI[30][2] = 3;

  // 31, G/GG
  MijMDI[31][10] = MijMDI[31][25] = MijMDI[31][37] = 1 / 9.0;
  MijMDI[31][17] = MijMDI[31][29] = MijMDI[31][32] = 2 / 9.0;
  MijMDI[31][51] = 2 * deletionRate / 3.0;
  MijMDI[31][57] = deletionRate / 3.0;

  genotypeArrMDI[31][0] = 2;
  genotypeArrMDI[31][1] = 2;
  genotypeArrMDI[31][2] = 2;

  // 32, G/GT
  MijMDI[32][11] = MijMDI[32][17] = MijMDI[32][18] = MijMDI[32][26] =
      MijMDI[32][29] = MijMDI[32][30] = MijMDI[32][31] = MijMDI[32][33] =
          MijMDI[32][38] = 1 / 9.0;
  MijMDI[32][51] = MijMDI[32][52] = MijMDI[32][58] = deletionRate / 3.0;

  genotypeArrMDI[32][0] = 2;
  genotypeArrMDI[32][1] = 2;
  genotypeArrMDI[32][2] = 3;

  // 33, G/TT
  MijMDI[33][12] = MijMDI[33][27] = MijMDI[33][39] = 1 / 9.0;
  MijMDI[33][18] = MijMDI[33][30] = MijMDI[33][32] = 2 / 9.0;
  MijMDI[33][52] = 2 * deletionRate / 3.0;
  MijMDI[33][59] = deletionRate / 3.0;

  genotypeArrMDI[33][0] = 2;
  genotypeArrMDI[33][1] = 3;
  genotypeArrMDI[33][2] = 3;

  // 34, T/CC
  MijMDI[34][7] = MijMDI[34][22] = MijMDI[34][28] = 1 / 9.0;
  MijMDI[34][19] = MijMDI[34][35] = MijMDI[34][36] = 2 / 9.0;
  MijMDI[34][50] = 2 * deletionRate / 3.0;
  MijMDI[34][54] = deletionRate / 3.0;

  genotypeArrMDI[34][0] = 3;
  genotypeArrMDI[34][1] = 1;
  genotypeArrMDI[34][2] = 1;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 35, T/CG
  MijMDI[35][8] = MijMDI[35][19] = MijMDI[35][20] = MijMDI[35][23] =
      MijMDI[35][29] = MijMDI[35][34] = MijMDI[35][36] = MijMDI[35][37] =
          MijMDI[35][38] = 1 / 9.0;
  MijMDI[35][50] = MijMDI[35][52] = MijMDI[35][55] = deletionRate / 3.0;

  genotypeArrMDI[35][0] = 3;
  genotypeArrMDI[35][1] = 1;
  genotypeArrMDI[35][2] = 2;

  // 36, T/CT
  MijMDI[36][9] = MijMDI[36][19] = MijMDI[36][21] = MijMDI[36][24] =
      MijMDI[36][30] = MijMDI[36][34] = MijMDI[36][35] = MijMDI[36][38] =
          MijMDI[36][39] = 1 / 9.0;
  MijMDI[36][50] = MijMDI[36][53] = MijMDI[36][56] = deletionRate / 3.0;

  genotypeArrMDI[36][0] = 3;
  genotypeArrMDI[36][1] = 1;
  genotypeArrMDI[36][2] = 3;

  // 37, T/GG
  MijMDI[37][10] = MijMDI[37][25] = MijMDI[37][31] = 1 / 9.0;
  MijMDI[37][20] = MijMDI[37][35] = MijMDI[37][38] = 2 / 9.0;
  MijMDI[37][52] = 2 * deletionRate / 3.0;
  MijMDI[37][57] = deletionRate / 3.0;

  genotypeArrMDI[37][0] = 3;
  genotypeArrMDI[37][1] = 2;
  genotypeArrMDI[37][2] = 2;

  // 38, T/GT
  MijMDI[38][11] = MijMDI[38][20] = MijMDI[38][21] = MijMDI[38][26] =
      MijMDI[38][32] = MijMDI[38][35] = MijMDI[38][36] = MijMDI[38][37] =
          MijMDI[38][39] = 1 / 9.0;
  MijMDI[38][52] = MijMDI[38][53] = MijMDI[38][58] = deletionRate / 3.0;

  genotypeArrMDI[38][0] = 3;
  genotypeArrMDI[38][1] = 2;
  genotypeArrMDI[38][2] = 3;

  // 39, T/TT
  MijMDI[39][12] = MijMDI[39][27] = MijMDI[39][33] = 1 / 9.0;
  MijMDI[39][21] = MijMDI[39][36] = MijMDI[39][38] = 2 / 9.0;
  MijMDI[39][53] = 2 * deletionRate / 3.0;
  MijMDI[39][59] = deletionRate / 3.0;

  genotypeArrMDI[39][0] = 3;
  genotypeArrMDI[39][1] = 3;
  genotypeArrMDI[39][2] = 3;

  // 40, A/A
  MijMDI[40][0] = 2 * insertionRate / 3.0;
  MijMDI[40][42] = MijMDI[40][43] = MijMDI[40][44] = 2 / 9.0;
  MijMDI[40][60] = 2 * deletionRate / 3.0;

  genotypeArrMDI[40][0] = 0;
  genotypeArrMDI[40][1] = 0;
  genotypeArrMDI[40][2] = INVALID;

  // 41, -/AA
  MijMDI[41][45] = MijMDI[41][46] = MijMDI[41][47] = 2 / 9.0;
  MijMDI[41][60] = 2 * deletionRate / 3.0;

  genotypeArrMDI[41][0] = DELETED;
  genotypeArrMDI[41][1] = 0;
  genotypeArrMDI[41][2] = 0;

  // 42, A/C
  MijMDI[42][4] = MijMDI[42][7] = insertionRate / 3.0;
  MijMDI[42][40] = MijMDI[42][43] = MijMDI[42][44] = MijMDI[42][48] =
      MijMDI[42][49] = MijMDI[42][50] = 1 / 9.0;
  MijMDI[42][60] = MijMDI[42][61] = deletionRate / 3.0;

  genotypeArrMDI[42][0] = 0;
  genotypeArrMDI[42][1] = 1;
  genotypeArrMDI[42][2] = INVALID;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 43, A/G
  MijMDI[43][5] = MijMDI[43][10] = insertionRate / 3.0;
  MijMDI[43][40] = MijMDI[43][42] = MijMDI[43][44] = MijMDI[43][49] =
      MijMDI[43][51] = MijMDI[43][52] = 1 / 9.0;
  MijMDI[43][60] = MijMDI[43][62] = deletionRate / 3.0;

  genotypeArrMDI[43][0] = 0;
  genotypeArrMDI[43][1] = 2;
  genotypeArrMDI[43][2] = INVALID;

  // 44, A/T
  MijMDI[44][6] = MijMDI[44][12] = insertionRate / 3.0;
  MijMDI[44][40] = MijMDI[44][42] = MijMDI[44][43] = MijMDI[44][50] =
      MijMDI[44][52] = MijMDI[44][53] = 1 / 9.0;
  MijMDI[44][60] = MijMDI[44][63] = deletionRate / 3.0;

  genotypeArrMDI[44][0] = 0;
  genotypeArrMDI[44][1] = 3;
  genotypeArrMDI[44][2] = INVALID;

  // 45, -/AC
  MijMDI[45][41] = MijMDI[45][46] = MijMDI[45][47] = MijMDI[45][54] =
      MijMDI[45][55] = MijMDI[45][56] = 1 / 9.0;
  MijMDI[45][60] = MijMDI[45][61] = deletionRate / 3.0;

  genotypeArrMDI[45][0] = DELETED;
  genotypeArrMDI[45][1] = 0;
  genotypeArrMDI[45][2] = 1;

  // 46, -/AG
  MijMDI[46][41] = MijMDI[46][45] = MijMDI[46][47] = MijMDI[46][55] =
      MijMDI[46][57] = MijMDI[46][58] = 1 / 9.0;
  MijMDI[46][60] = MijMDI[46][62] = deletionRate / 3.0;

  genotypeArrMDI[46][0] = DELETED;
  genotypeArrMDI[46][1] = 0;
  genotypeArrMDI[46][2] = 2;

  // 47, -/AT
  MijMDI[47][41] = MijMDI[47][45] = MijMDI[47][46] = MijMDI[47][56] =
      MijMDI[47][58] = MijMDI[47][59] = 1 / 9.0;
  MijMDI[47][60] = MijMDI[47][63] = deletionRate / 3.0;

  genotypeArrMDI[47][0] = DELETED;
  genotypeArrMDI[47][1] = 0;
  genotypeArrMDI[47][2] = 3;

  // 48, C/C
  MijMDI[48][22] = 2 * insertionRate / 3.0;
  MijMDI[48][42] = MijMDI[48][49] = MijMDI[48][50] = 2 / 9.0;
  MijMDI[48][61] = 2 * deletionRate / 3.0;

  genotypeArrMDI[48][0] = 1;
  genotypeArrMDI[48][1] = 1;
  genotypeArrMDI[48][2] = INVALID;

  // 49, C/G
  MijMDI[49][25] = MijMDI[49][28] = insertionRate / 3.0;
  MijMDI[49][42] = MijMDI[49][43] = MijMDI[49][48] = MijMDI[49][50] =
      MijMDI[49][51] = MijMDI[49][52] = 1 / 9.0;
  MijMDI[49][61] = MijMDI[49][62] = deletionRate / 3.0;

  genotypeArrMDI[49][0] = 1;
  genotypeArrMDI[49][1] = 2;
  genotypeArrMDI[49][2] = INVALID;

  // 50, C/T
  MijMDI[50][27] = MijMDI[50][34] = insertionRate / 3.0;
  MijMDI[50][42] = MijMDI[50][44] = MijMDI[50][48] = MijMDI[50][49] =
      MijMDI[50][52] = MijMDI[50][53] = 1 / 9.0;
  MijMDI[50][61] = MijMDI[50][63] = deletionRate / 3.0;

  genotypeArrMDI[50][0] = 1;
  genotypeArrMDI[50][1] = 3;
  genotypeArrMDI[50][2] = INVALID;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 51, G/G
  MijMDI[51][31] = 2 * insertionRate / 3.0;
  MijMDI[51][43] = MijMDI[51][49] = MijMDI[51][52] = 2 / 9.0;
  MijMDI[51][62] = deletionRate / 3.0;

  genotypeArrMDI[51][0] = 2;
  genotypeArrMDI[51][1] = 2;
  genotypeArrMDI[51][2] = INVALID;

  // 52, G/T
  MijMDI[52][33] = MijMDI[52][37] = insertionRate / 3.0;
  MijMDI[52][43] = MijMDI[52][44] = MijMDI[52][49] = MijMDI[52][50] =
      MijMDI[52][51] = MijMDI[52][53] = 1 / 9.0;
  MijMDI[52][62] = MijMDI[52][63] = deletionRate / 3.0;

  genotypeArrMDI[52][0] = 2;
  genotypeArrMDI[52][1] = 3;
  genotypeArrMDI[52][2] = INVALID;

  // 53, T/T
  MijMDI[53][39] = 2 * insertionRate / 3.0;
  MijMDI[53][44] = MijMDI[53][50] = MijMDI[53][52] = 2 / 9.0;
  MijMDI[53][63] = 2 * deletionRate / 3.0;

  genotypeArrMDI[53][0] = 3;
  genotypeArrMDI[53][1] = 3;
  genotypeArrMDI[53][2] = INVALID;

  // 54, -/CC
  MijMDI[54][45] = MijMDI[54][55] = MijMDI[54][56] = 2 / 9.0;
  MijMDI[54][61] = 2 * deletionRate / 3.0;

  genotypeArrMDI[54][0] = DELETED;
  genotypeArrMDI[54][1] = 1;
  genotypeArrMDI[54][2] = 1;

  // 55, -/CG
  MijMDI[55][45] = MijMDI[55][46] = MijMDI[55][54] = MijMDI[55][56] =
      MijMDI[55][57] = MijMDI[55][58] = 1 / 9.0;
  MijMDI[55][61] = MijMDI[55][62] = deletionRate / 3.0;

  genotypeArrMDI[55][0] = DELETED;
  genotypeArrMDI[55][1] = 1;
  genotypeArrMDI[55][2] = 2;

  // 56, -/CT
  MijMDI[56][45] = MijMDI[56][47] = MijMDI[56][54] = MijMDI[56][55] =
      MijMDI[56][58] = MijMDI[56][59] = 1 / 9.0;
  MijMDI[56][61] = MijMDI[56][63] = deletionRate / 3.0;

  genotypeArrMDI[56][0] = DELETED;
  genotypeArrMDI[56][1] = 1;
  genotypeArrMDI[56][2] = 3;

  // 57, -/GG
  MijMDI[57][46] = MijMDI[57][55] = MijMDI[57][58] = 2 / 9.0;
  MijMDI[57][62] = 2 * deletionRate / 3.0;

  genotypeArrMDI[57][0] = DELETED;
  genotypeArrMDI[57][1] = 2;
  genotypeArrMDI[57][2] = 2;

  // 58, -/GT
  MijMDI[58][46] = MijMDI[58][47] = MijMDI[58][55] = MijMDI[58][56] =
      MijMDI[58][57] = MijMDI[58][59] = 1 / 9.0;
  MijMDI[58][62] = MijMDI[58][63] = deletionRate / 3.0;

  genotypeArrMDI[58][0] = DELETED;
  genotypeArrMDI[58][1] = 2;
  genotypeArrMDI[58][2] = 3;

  // 0     1     2     3     4     5     6     7     8     9     10    11
  // A/AA, A/AC, A/AG, A/AT, C/AA, G/AA, T/AA, A/CC, A/CG, A/CT, A/GG, A/GT,

  // 12    13    14    15    16    17    18    19    20    21    22    23
  // A/TT, C/AC, C/AG, C/AT, G/AC, G/AG, G/AT, T/AC, T/AG, T/AT, C/CC, C/CG,

  // 24    25    26    27    28    29    30    31    32    33    34    35
  // C/CT, C/GG, C/GT, C/TT, G/CC, G/CG, G/CT, G/GG, G/GT, G/TT, T/CC, T/CG,

  // 36    37    38    39    40   41    42   43   44   45    46    47    48
  // T/CT, T/GG, T/GT, T/TT, A/A, -/AA, A/C, A/G, A/T, -/AC, -/AG, -/AT, C/C,

  // 49   50   51   52   53   54    55    56    57    58    59    60   61
  // C/G, C/T, G/G, G/T, T/T, -/CC, -/CG, -/CT, -/GG, -/GT, -/TT, A/-, C/-,

  // 62   63   64
  // G/-, T/-, -

  // 59, -/TT
  MijMDI[59][47] = MijMDI[59][56] = MijMDI[59][58] = 2 / 9.0;
  MijMDI[59][63] = 2 * deletionRate / 3.0;

  genotypeArrMDI[59][0] = DELETED;
  genotypeArrMDI[59][1] = 3;
  genotypeArrMDI[59][2] = 3;

  // 60, A/-
  MijMDI[60][41] = insertionRate / 3.0;
  MijMDI[60][61] = MijMDI[60][62] = MijMDI[60][63] = 1 / 9.0;
  MijMDI[60][64] = deletionRate / 3.0;

  genotypeArrMDI[60][0] = 0;
  genotypeArrMDI[60][1] = DELETED;
  genotypeArrMDI[60][2] = INVALID;

  // 61, C/-
  MijMDI[61][54] = insertionRate / 3.0;
  MijMDI[61][60] = MijMDI[61][62] = MijMDI[61][63] = 1 / 9.0;
  MijMDI[61][64] = deletionRate / 3.0;

  genotypeArrMDI[61][0] = 1;
  genotypeArrMDI[61][1] = DELETED;
  genotypeArrMDI[61][2] = INVALID;

  // 62, G/-
  MijMDI[62][57] = insertionRate / 3.0;
  MijMDI[62][60] = MijMDI[62][61] = MijMDI[62][63] = 1 / 9.0;
  MijMDI[62][64] = deletionRate / 3.0;

  genotypeArrMDI[62][0] = 2;
  genotypeArrMDI[62][1] = DELETED;
  genotypeArrMDI[62][2] = INVALID;

  // 63, T/-
  MijMDI[63][59] = insertionRate / 3.0;
  MijMDI[63][60] = MijMDI[63][61] = MijMDI[63][62] = 1 / 9.0;
  MijMDI[63][64] = deletionRate / 3.0;

  genotypeArrMDI[63][0] = 3;
  genotypeArrMDI[63][1] = DELETED;
  genotypeArrMDI[63][2] = INVALID;

  // 64, -
  genotypeArrMDI[64][0] = genotypeArrMDI[64][1] = DELETED;
  genotypeArrMDI[64][2] = INVALID;

  InitGenotypeHashMapMDI();
}

/**
 * Get the base name of a folder.
 * E.g., .../path/to/me/, where 'me' is a folder, then 'me' will be returned.
 *
 * @param path a path
 * @return base name of a folder
 */
char *GetFolderBaseName(char *path) {
  int separatorIndex = -1, i = 0;
  while (1) {
    if (i > MAX_NAME) {
      fprintf(stderr, "Error! The path doesn't end with '%c'!", '\0');
      exit(-1);
    } else {
      if (path[i] == '\0')
        break;
      else if (path[i] == '/')
        separatorIndex = i;
    }
    i++;
  }

  if (separatorIndex == -1)
    return path;
  else
    return path + separatorIndex + 1;
}

int GetMaxInt(const int *arr, size_t size) {
  int max = arr[0];

  for (size_t i = 1; i < size; i++)
    if (max < arr[i])
      max = arr[i];

  return max;
}

void Normalize(double *arr, size_t size) {
  double sum = 0.0;

  for (size_t i = 0; i < size; i++)
    sum += arr[i];

  if (sum == 1.0)
    return;

  for (size_t i = 0; i < size; i++)
    arr[i] = arr[i] / sum;
}

/* determine whether to make a coalescence tree or not */
int DoCoalescenceTree(int numDataSets, int replicate, int fixedTree,
                      int noisy) {
  if (numDataSets == 1) {
    return 1;
  } else {
    if (replicate == 0) {
      return 1;
    } else {
      if (fixedTree == 1) {
        if (noisy > 2) {
          fprintf(stderr,
                  "\n Use the fixed tree which has been coalesced before.");
        }

        return 0;
      } else {
        return 1;
      }
    }
  }
}

int GetExistingAllelesNrMD(int code) {
  if (code >= 0 && code < 10) {
    return 2;
  } else if (code < 14 && code >= 10) {
    return 1;
  } else if (code == 14) {
    return 0;
  } else {
    return -1;
  }
}

int GetExistingAllelesNrMDI(int code) {
  if (code >= 0 && code < 40) {
    return 3;
  } else if (code < 60 && code >= 40) {
    return 2;
  } else if (code < 64 && code >= 60) {
    return 1;
  } else if (code == 64) {
    return 0;
  } else {
    return -1;
  }
}

int WhichRefGenotype(int nucleotide, int MuInDelModel) {
  if (nucleotide == A) {
    if (MuInDelModel == YES) {
      return 40;
    } else {
      return 0;
    }
  } else if (nucleotide == C) {
    if (MuInDelModel == YES) {
      return 48;
    } else {
      return 4;
    }
  } else if (nucleotide == G) {
    if (MuInDelModel == YES) {
      return 51;
    } else {
      return 7;
    }
  } else if (nucleotide == T) {
    if (MuInDelModel == YES) {
      return 53;
    } else {
      return 9;
    }
  } else {
    fprintf(stderr, "\nERROR in WhichRefGenotype: nucleotide = %c\n",
            nucleotide);
    exit(-1);
  }
}

int NumOverlap(const int *paGenotype, const int *chGenotype, int maxCopies,
               int paStart, int chStart, int *overlap) {
  int num = 0, flag_1, flag_2;
  int localOverlap[maxCopies];
  int *tmp;

  if (overlap == NULL) {
    tmp = localOverlap;
  } else {
    tmp = overlap;
  }

  for (int i = 0; i < maxCopies; i++) {
    tmp[i] = -1;
  }

  for (int i = chStart; i < maxCopies; ++i) {
    flag_1 = 0;
    if (chGenotype[i] >= 0) {
      for (int j = paStart; j < maxCopies; ++j) {
        flag_2 = 0;
        if (paGenotype[j] >= 0) {
          if (chGenotype[i] == paGenotype[j]) {
            for (int k = 0; k < maxCopies; ++k) {
              if (tmp[k] == -1) {
                tmp[k] = j;
                num++;
                flag_1 = 1;
                break;
              } else {
                if (tmp[k] == j) {
                  flag_2 = 1;
                  break;
                } else {
                  continue;
                }
              }
            }
            if (flag_1 == 1)
              break;
            if (flag_2 == 1)
              continue;
          }
        }
      }
    }
  }

  return num;
}

void AllocateThreads(int threadsNr) {
  threads = (pthread_t *)malloc(threadsNr * sizeof(pthread_t));
  if (!threads) {
    fprintf(stderr, "Could not allocate the threads structure\n");
    exit(-1);
  }
}

void InitThreadsStructure(int cellNr, int siteNr, int threadsNr) {
  threadStr = (ThreadStr *)malloc((threadsNr + 1) * sizeof(ThreadStr));
  if (!threadStr) {
    fprintf(stderr, "Could not allocate the threadStr structure\n");
    exit(-1);
  }
  sitesPointsForThreads = (int *)calloc(threadsNr + 2, sizeof(int));
  if (!sitesPointsForThreads) {
    fprintf(stderr, "Could not allocate the sitesPointsForThreads structure\n");
    exit(-1);
  }

  int nrSitesPerThread = siteNr / threadsNr;
  int start = 0, end = nrSitesPerThread;
  for (long i = 0; i < threadsNr; i++) {
    threadStr[i].id = i;
    threadStr[i].start = start;
    threadStr[i].end = i == threadsNr - 1 ? siteNr : end;
    threadStr[i].numSubSites = threadStr[i].end - threadStr[i].start;

    /* allocate genotype data to be stored in data[cell][site] */
    /* add 1 more size when allocating to make sure free() correctly works */
    /* only allocate space for the first dataset, and release it after all
     * dataset processed to save space and time */
    threadStr[i].data = (DataStr **)calloc(2 * cellNr + 2, sizeof(DataStr *));
    if (!threadStr[i].data) {
      fprintf(stderr, "Could not allocate the data structure in thread %ld\n",
              i);
      exit(-1);
    }

    for (int j = 0; j < 2 * cellNr + 1; j++) {
      threadStr[i].data[j] =
          (DataStr *)calloc(threadStr[i].numSubSites + 1, sizeof(DataStr));
      if (!threadStr[i].data[j]) {
        fprintf(stderr,
                "Could not allocate the data[] structure in thread %ld for "
                "cell %d\n",
                i, j);
        exit(-1);
      }

      for (int k = 0; k < threadStr[i].numSubSites; k++) {
        threadStr[i].data[j][k].pastEvolutionaryEvents =
            (unsigned int *)calloc(6, sizeof(unsigned int));
        if (!threadStr[i].data[j][k].pastEvolutionaryEvents) {
          fprintf(stderr,
                  "Could not allocate the data[][].pastEvolutionaryEvents"
                  " structure in thread %ld for cell %d and site %d\n",
                  i, j, k);
          exit(-1);
        }
      }
    }

    threadStr[i].hasParallelMutations =
        (int *)calloc(threadStr[i].numSubSites + 1, sizeof(int));
    threadStr[i].parallelMutations =
        (unsigned **)calloc(threadStr[i].numSubSites + 1, sizeof(unsigned *));
    threadStr[i].summary =
        (int **)calloc(threadStr[i].numSubSites + 1, sizeof(int *));
    for (int j = 0; j < threadStr[i].numSubSites; j++) {
      threadStr[i].summary[j] = (int *)calloc(5, sizeof(int));
      threadStr[i].parallelMutations[j] =
          (unsigned *)calloc(cellNr + 2, sizeof(unsigned));
    }

    threadStr[i].isSNV =
        (int *)calloc(threadStr[i].numSubSites + 1, sizeof(int));

    sitesPointsForThreads[i + 1] = threadStr[i].end;

    start = end;
    end += nrSitesPerThread;
  }
}

void InitRawDataPrinters(int cellNr, int threadsNr, long *seed) {
  rawDataPrinters = (ThreadedRawDataPrinter *)malloc(
      (threadsNr + 1) * sizeof(ThreadedRawDataPrinter));

  int numCellsPerThread = 0;
  if (cellNr + 1 > threadsNr)
    numCellsPerThread = (cellNr + 1) / threadsNr;
  else
    numCellsPerThread = 1;

  int numCellsLeft = cellNr + 1;
  int cellIndex = 0;
  for (long i = 0; i < threadsNr; i++) {
    rawDataPrinters[i].id = i;
    rawDataPrinters[i].seed = seed;

    if (i < threadsNr - 1) {
      if (numCellsLeft > 0 && numCellsLeft >= numCellsPerThread)
        rawDataPrinters[i].numCells = numCellsPerThread;
      else if (numCellsLeft > 0 && numCellsLeft < numCellsPerThread)
        rawDataPrinters[i].numCells = numCellsLeft;
      else
        rawDataPrinters[i].numCells = 0;
    } else if (i == threadsNr - 1) {
      if (numCellsLeft > 0)
        rawDataPrinters[i].numCells = numCellsLeft;
      else
        rawDataPrinters[i].numCells = 0;
    }

    if (rawDataPrinters[i].numCells > 0) {

      rawDataPrinters[i].cellIndices = (unsigned int *)malloc(
          (rawDataPrinters[i].numCells + 1) * sizeof(unsigned int));
      rawDataPrinters[i].files = (RawDataFile *)malloc(
          (rawDataPrinters[i].numCells + 1) * sizeof(RawDataFile));

      for (int j = 0; j < rawDataPrinters[i].numCells; j++) {
        rawDataPrinters[i].cellIndices[j] = cellIndex + j;
      }

    } else {
      rawDataPrinters[i].cellIndices = NULL;
      rawDataPrinters[i].files = NULL;
    }

    numCellsLeft -= rawDataPrinters[i].numCells;
    if (numCellsLeft < 0)
      numCellsLeft = 0;

    cellIndex += rawDataPrinters[i].numCells;
  }
}

void PrepareThreadedRawReadSeqDataFiles(int threadNr) {
  for (long i = 0; i < threadNr; i++) {
    for (int j = 0; j < rawDataPrinters[i].numCells; j++) {
      RawDataFile tmp = rawDataFiles[rawDataPrinters[i].cellIndices[j] + 1];
      rawDataPrinters[i].files[j].fp = tmp.fp;
      strcpy(rawDataPrinters[i].files[j].fileName, tmp.fileName);
    }
  }
}

int GetThreadIndex(int siteIndex) {
  for (int i = 0; i < numThreads; i++) {
    if (sitesPointsForThreads[i] <= siteIndex &&
        siteIndex < sitesPointsForThreads[i + 1]) {
      return i;
    }
  }

  fprintf(stderr,
          "Could not find the index of thread (from 0 - %d) according "
          "to the given site %d \n",
          sitesPointsForThreads[numThreads], siteIndex);
  exit(-1);
}

void AllocateCellNamesExtd(int cellNr) {
  cellNamesExtd = (char **)calloc(cellNr + 2, sizeof(char *));
  if (!cellNamesExtd) {
    fprintf(stderr, "Could not allocate tipNames (%ld)\n",
            (cellNr + 2) * sizeof(char *));
    exit(-1);
  }
  for (int i = 0; i < cellNr + 1; i++) {
    cellNamesExtd[i] = (char *)calloc(MAX_NAME, sizeof(char));
    if (!cellNamesExtd[i]) {
      fprintf(stderr, "Could not allocate tipNames[%d] (%ld)\n", i,
              (MAX_NAME) * sizeof(char));
      exit(-1);
    }
  }
}

void AllocateCellStructureExtd(int cellNr, int siteNr) {
  cellExtd = (CellStrExtd *)calloc(cellNr + 2, sizeof(CellStrExtd));
  if (!cellExtd) {
    fprintf(stderr, "Could not allocate cellExtd (%ld)\n",
            (cellNr + 1) * sizeof(CellStrExtd));
    exit(-1);
  }

  for (int i = 0; i < cellNr + 1; ++i) {
    cellExtd[i].site =
        (CellSiteStrExtd *)calloc(siteNr + 1, sizeof(CellSiteStrExtd));
    if (!cellExtd[i].site) {
      fprintf(stderr, "Could not allocate the cellExtd[%d].site structure\n",
              i);
      exit(-1);
    }

    for (int j = 0; j < siteNr; ++j) {
      cellExtd[i].site[j].hasADO = NO;
      cellExtd[i].site[j].numADO = 0;

      cellExtd[i].site[j].readCount = (int *)calloc(5, sizeof(int));
      if (!cellExtd[i].site[j].readCount) {
        fprintf(stderr,
                "Could not allocate the cellExtd[%d].site[%d].readCount "
                "structure\n",
                i, j);
        exit(-1);
      }

      cellExtd[i].site[j].observedAlleles = (int *)calloc(5, sizeof(int));
      if (!cellExtd[i].site[j].observedAlleles) {
        fprintf(stderr,
                "Could not allocate the cellExtd[%d].site[%d].observedAlleles "
                "structure\n",
                i, j);
        exit(-1);
      }

      cellExtd[i].site[j].log10GenotypeLikelihood =
          (double *)calloc(NUM_CAN_ML_GENOTYPES + 1, sizeof(double));
      if (!cellExtd[i].site[j].log10GenotypeLikelihood) {
        fprintf(
            stderr,
            "Could not allocate the cellExtd[%d].site[%d].genotypeLikelihood "
            "structure\n",
            i, j);
        exit(-1);
      }

      cellExtd[i].site[j].scaledLog10GenotypeLikelihood =
          (double *)calloc(NUM_CAN_ML_GENOTYPES + 1, sizeof(double));
      if (!cellExtd[i].site[j].scaledLog10GenotypeLikelihood) {
        fprintf(stderr,
                "Could not allocate the "
                "cellExtd[%d].site[%d].scaledGenotypeLikelihood "
                "structure\n",
                i, j);
        exit(-1);
      }
    }
  }
}

int GetReadBlockIndex(int site) {
  for (int i = 0; i < numReadBlocks; i++) {
    if (site < startSiteBlocks[i + 1])
      return i;
  }
}

int ValueInArray(int value, const int *array, int size) {
  for (int i = 0; i < size; i++)
    if (array[i] == value)
      return 1;

  return 0;
}

int AddValueToArray(int value, int *array, int size, int marker) {
  for (int i = 0; i < size; i++)
    if (array[i] == marker) {
      array[i] = value;
      return 0;
    }

  return 1;
}

void CountMutationNum(const int *chAlleles, int startPos, int maxCopies,
                      int referenceNuc, const int *paAlleles,
                      const int *overlap, unsigned int *mutationCounter) {
  int exOverlapPos[maxCopies], compareFlag;

  for (int i = 0; i < maxCopies; i++)
    exOverlapPos[i] = -1;

  for (int i = startPos; i < maxCopies; i++) {
    compareFlag = 1;

    if (chAlleles[i] < 0)
      continue;

    for (int j = 0; j < maxCopies && overlap[j] >= 0; j++) {
      if (chAlleles[i] == paAlleles[overlap[j]]) {
        if (ValueInArray(j, exOverlapPos, maxCopies) == 0) {
          if (AddValueToArray(j, exOverlapPos, maxCopies, -1) == 0) {
            compareFlag = 0;
            break;
          } else {
            fprintf(stderr,
                    "ERROR! Cannot add %d in an array of length "
                    "%d as it is already full: [",
                    j, maxCopies);
            for (int k = 0; k < maxCopies; k++) {
              fprintf(stderr, "%d", exOverlapPos[k]);

              if (k < maxCopies - 1)
                fprintf(stderr, ",");
              else
                fprintf(stderr, "]\n");
            }
            exit(1);
          }
        } else
          continue;
      }
    }

    if (compareFlag == 1) {
      mutationCounter[0]++;

      if (chAlleles[i] == referenceNuc)
        mutationCounter[1]++;
    }
  }
}

void GetMuDelInsNr(ThreadStr *thStr, int anccell, int cell, int site,
                   int maxCopies, int referenceNuc) {
  int overlapNr, overlap[maxCopies];

  if (thStr->data[cell][site].numExistingAlleles !=
      thStr->data[anccell][site].numExistingAlleles) {
    if (thStr->data[anccell][site].numExistingAlleles == 3) {
      // allSites[site].numDeletions++;
      // numDEL++;
      thStr->summary[site][3]++;
      thStr->data[cell][site].pastEvolutionaryEvents[3]++;

      if (thStr->data[cell][site].numExistingAlleles == 2) {
        if (thStr->data[cell][site].alleles[0] < 0) {
          // e.g., A/AA -> -/AA
          overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                                 thStr->data[cell][site].alleles, maxCopies, 1,
                                 1, overlap);
          // allSites[site].numMutations +=
          // (thStr->data[cell][site].numExistingAlleles - overlapNr); numMU +=
          // (thStr->data[cell][site].numExistingAlleles - overlapNr);
          thStr->summary[site][0] +=
              (thStr->data[cell][site].numExistingAlleles - overlapNr);
          CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                           referenceNuc, thStr->data[anccell][site].alleles,
                           overlap,
                           thStr->data[cell][site].pastEvolutionaryEvents);
        } else {
          overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                                 thStr->data[cell][site].alleles, maxCopies, 0,
                                 0, overlap);
          if (overlapNr == 2) {
            if (thStr->data[anccell][site].alleles[0] !=
                    thStr->data[cell][site].alleles[0] &&
                thStr->data[anccell][site].alleles[0] !=
                    thStr->data[cell][site].alleles[1]) {
              // e.g., A/CG -> C/G
              // Ambiguity will raise.
              // allSites[site].numMutations++;
              // numMU++;
              thStr->summary[site][0]++;
              thStr->data[cell][site].pastEvolutionaryEvents[0]++;
              if (ValueInArray(referenceNuc, thStr->data[cell][site].alleles,
                               2) == 1)
                thStr->data[cell][site].pastEvolutionaryEvents[2]++;
            }
          } else {
            // allSites[site].numMutations +=
            // (thStr->data[cell][site].numExistingAlleles - overlapNr); numMU
            // += (thStr->data[cell][site].numExistingAlleles - overlapNr);
            thStr->summary[site][0] +=
                (thStr->data[cell][site].numExistingAlleles - overlapNr);
            CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        }
      } else if (thStr->data[cell][site].numExistingAlleles == 1) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        if (overlapNr == 0) {
          // allSites[site].numMutations++;
          // numMU++;
          thStr->summary[site][0]++;
          CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                           referenceNuc, thStr->data[anccell][site].alleles,
                           overlap,
                           thStr->data[cell][site].pastEvolutionaryEvents);
        }
      } else {
        // allSites[site].numDeletions += 2;
        // numDEL += 2;
        thStr->summary[site][3] += 2;
        thStr->data[cell][site].pastEvolutionaryEvents[3] += 2;
      }
    } else if (thStr->data[anccell][site].numExistingAlleles == 2) {
      if (thStr->data[cell][site].numExistingAlleles == 3) {
        if (thStr->data[anccell][site].alleles[0] >= 0) {
          // allSites[site].numInsertions++;
          // numINS++;
          thStr->summary[site][4]++;
          thStr->data[cell][site].pastEvolutionaryEvents[4]++;

          overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                                 thStr->data[cell][site].alleles, maxCopies, 0,
                                 0, overlap);
          if (thStr->data[anccell][site].alleles[0] ==
              thStr->data[anccell][site].alleles[1]) {
            // e.g., G/G
            if (overlapNr == 2) {
              if (thStr->data[cell][site].alleles[0] !=
                      thStr->data[cell][site].alleles[1] ||
                  thStr->data[cell][site].alleles[0] !=
                      thStr->data[cell][site].alleles[2] ||
                  thStr->data[cell][site].alleles[1] !=
                      thStr->data[cell][site].alleles[2]) {
                // e.g., G/G -> C/GG or G/CG
                // allSites[site].numMutations++;
                // numMU++;
                thStr->summary[site][0]++;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            } else if (overlapNr == 1) {
              if (thStr->data[cell][site].alleles[1] ==
                  thStr->data[cell][site].alleles[2]) {
                // e.g., G/G -> G/AA
                // allSites[site].numMutations++;
                // numMU++;
                thStr->summary[site][0]++;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., G/G -> A/CG
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            } else {
              if (thStr->data[cell][site].alleles[1] ==
                  thStr->data[cell][site].alleles[2]) {
                // e.g., G/G -> A/CC
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., G/G -> A/CT
                // allSites[site].numMutations += 3;
                // numMU += 3;
                thStr->summary[site][0] += 3;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            }
          } else {
            // e.g., C/T
            if (overlapNr == 2) {
              if (thStr->data[cell][site].alleles[1] !=
                  thStr->data[cell][site].alleles[2]) {
                if ((thStr->data[cell][site].alleles[0] !=
                     thStr->data[anccell][site].alleles[0]) &&
                    (thStr->data[cell][site].alleles[0] !=
                     thStr->data[anccell][site].alleles[1])) {
                  // e.g., C/T -> A/CT
                  // Ambiguity could raise.
                  // allSites[site].numMutations += 2;
                  // numMU += 2;
                  thStr->summary[site][0] += 2;
                  thStr->data[cell][site].pastEvolutionaryEvents[0] += 2;
                  if (referenceNuc == thStr->data[cell][site].alleles[0])
                    thStr->data[cell][site].pastEvolutionaryEvents[1] += 1;
                  if (referenceNuc == thStr->data[cell][site].alleles[1] ||
                      referenceNuc == thStr->data[cell][site].alleles[2])
                    thStr->data[cell][site].pastEvolutionaryEvents[2] += 1;
                } else {
                  // e.g., C/T -> C/AT
                  // allSites[site].numMutations++;
                  // numMU++;
                  thStr->summary[site][0]++;
                  CountMutationNum(
                      thStr->data[cell][site].alleles, 0, maxCopies,
                      referenceNuc, thStr->data[anccell][site].alleles, overlap,
                      thStr->data[cell][site].pastEvolutionaryEvents);
                }
              }
            } else if (overlapNr == 1) {
              if (thStr->data[cell][site].alleles[1] !=
                  thStr->data[cell][site].alleles[2]) {
                // e.g., C/T -> C/AG
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., C/T -> A/TT
                // allSites[site].numMutations++;
                // numMU++;
                thStr->summary[site][0]++;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            } else {
              if (thStr->data[cell][site].alleles[1] !=
                  thStr->data[cell][site].alleles[2]) {
                // e.g., C/T -> G/AG
                // allSites[site].numMutations += 3;
                // numMU += 3;
                thStr->summary[site][0] += 3;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, maxCopies, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              } else {
                // e.g., C/T -> G/AA
                // allSites[site].numMutations += 2;
                // numMU += 2;
                thStr->summary[site][0] += 2;
                CountMutationNum(
                    thStr->data[cell][site].alleles, 0, 2, referenceNuc,
                    thStr->data[anccell][site].alleles, overlap,
                    thStr->data[cell][site].pastEvolutionaryEvents);
              }
            }
          }
        }
      } else if (thStr->data[cell][site].numExistingAlleles == 1) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        if (overlapNr == 0) {
          // e.g., C/T -> A/-
          // allSites[site].numMutations++;
          // numMU++;
          thStr->summary[site][0]++;
          CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                           referenceNuc, thStr->data[anccell][site].alleles,
                           overlap,
                           thStr->data[cell][site].pastEvolutionaryEvents);
        }
      } else {
        // allSites[site].numDeletions += 2;
        // numDEL += 2;
        thStr->summary[site][3] += 2;
        thStr->data[cell][site].pastEvolutionaryEvents[3] += 2;
      }
    } else if (thStr->data[anccell][site].numExistingAlleles == 1) {
      // e.g., C/-
      if (thStr->data[cell][site].numExistingAlleles == 2 &&
          thStr->data[cell][site].alleles[0] == DELETED) {
        // allSites[site].numInsertions++;
        // numINS++;
        thStr->summary[site][4]++;
        thStr->data[cell][site].pastEvolutionaryEvents[4]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        if (overlapNr == 0) {
          if (thStr->data[cell][site].alleles[1] ==
              thStr->data[cell][site].alleles[2]) {
            // e.g., C/- -> -/AA
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, 2,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          } else {
            // e.g., C/- -> -/GT
            // allSites[site].numMutations += 2;
            // numMU += 2;
            thStr->summary[site][0] += 2;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        } else {
          // e.g., C/- -> -/CT
          if (thStr->data[cell][site].alleles[1] !=
              thStr->data[cell][site].alleles[2]) {
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        }
      } else if (thStr->data[cell][site].numExistingAlleles == 0) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;
      }
    }
  } else {
    if (thStr->data[anccell][site].numExistingAlleles == 3) {
      if (thStr->data[anccell][site].alleles[0] !=
          thStr->data[cell][site].alleles[0]) {
        // e.g., A/CC -> G/??
        // allSites[site].numMutations++;
        // numMU++;
        thStr->summary[site][0]++;
        thStr->data[cell][site].pastEvolutionaryEvents[0]++;
        if (referenceNuc == thStr->data[cell][site].alleles[0])
          thStr->data[cell][site].pastEvolutionaryEvents[1]++;
      }

      overlapNr =
          NumOverlap(thStr->data[anccell][site].alleles,
                     thStr->data[cell][site].alleles, maxCopies, 1, 1, overlap);
      // allSites[site].numMutations +=
      // (thStr->data[anccell][site].numExistingAlleles - 1 - overlapNr); numMU
      // += (thStr->data[anccell][site].numExistingAlleles - 1 - overlapNr);
      thStr->summary[site][0] +=
          (thStr->data[anccell][site].numExistingAlleles - 1 - overlapNr);
      CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                       referenceNuc, thStr->data[anccell][site].alleles,
                       overlap, thStr->data[cell][site].pastEvolutionaryEvents);
    } else if (thStr->data[anccell][site].numExistingAlleles == 1) {
      // allSites[site].numMutations++;
      // numMU++;
      thStr->summary[site][0]++;
      CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                       referenceNuc, thStr->data[anccell][site].alleles,
                       overlap, thStr->data[cell][site].pastEvolutionaryEvents);
    } else {
      if (thStr->data[anccell][site].alleles[0] != DELETED &&
          thStr->data[cell][site].alleles[0] == DELETED) {
        // allSites[site].numDeletions++;
        // numDEL++;
        thStr->summary[site][3]++;
        thStr->data[cell][site].pastEvolutionaryEvents[3]++;
        // allSites[site].numInsertions++;
        // numINS++;
        thStr->summary[site][4]++;
        thStr->data[cell][site].pastEvolutionaryEvents[4]++;

        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 1,
                               overlap);
        if (overlapNr == 2) {
          if (thStr->data[anccell][site].alleles[0] !=
              thStr->data[anccell][site].alleles[1]) {
            // e.g., G/T -> -/GT
            // Ambiguity might raise.
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            thStr->data[cell][site].pastEvolutionaryEvents[0]++;
            if (ValueInArray(referenceNuc, thStr->data[cell][site].alleles,
                             3) == 1)
              thStr->data[cell][site].pastEvolutionaryEvents[2]++;
          }
        } else if (overlapNr == 1) {
          if (thStr->data[cell][site].alleles[1] !=
              thStr->data[cell][site].alleles[2]) {
            // e.g., G/T -> -/AT
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        } else {
          if (thStr->data[cell][site].alleles[1] ==
              thStr->data[cell][site].alleles[2]) {
            // e.g., G/T -> -/AA
            // allSites[site].numMutations++;
            // numMU++;
            thStr->summary[site][0]++;
            CountMutationNum(thStr->data[cell][site].alleles, 1, 2,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          } else {
            // e.g., G/T -> -/AC
            // allSites[site].numMutations += 2;
            // numMU += 2;
            thStr->summary[site][0] += 2;
            CountMutationNum(thStr->data[cell][site].alleles, 1, maxCopies,
                             referenceNuc, thStr->data[anccell][site].alleles,
                             overlap,
                             thStr->data[cell][site].pastEvolutionaryEvents);
          }
        }
      } else if (!(thStr->data[anccell][site].alleles[0] == DELETED &&
                   thStr->data[cell][site].alleles[0] != DELETED)) {
        overlapNr = NumOverlap(thStr->data[anccell][site].alleles,
                               thStr->data[cell][site].alleles, maxCopies, 0, 0,
                               overlap);
        // allSites[site].numMutations +=
        // (thStr->data[anccell][site].numExistingAlleles - overlapNr); numMU +=
        // (thStr->data[anccell][site].numExistingAlleles - overlapNr);
        thStr->summary[site][0] +=
            (thStr->data[anccell][site].numExistingAlleles - overlapNr);
        CountMutationNum(thStr->data[cell][site].alleles, 0, maxCopies,
                         referenceNuc, thStr->data[anccell][site].alleles,
                         overlap,
                         thStr->data[cell][site].pastEvolutionaryEvents);
      }
    }
  }
}

int GetAllele(double randomNr, int num) {
  int alleleNr = -1, intArray[num];

  for (int i = 0; i < num; ++i) {
    intArray[i] = i + 1;
  }

  for (int i = 0; i < num; ++i) {
    if ((randomNr * num) < intArray[i]) {
      alleleNr = intArray[i] - 1;
      break;
    }
  }

  return alleleNr;
}

void GetGenotypeStructure(const int *alleles, int *out) {
  for (int i = 0; i < 4; ++i)
    out[i] = 0;

  for (int i = 0; i < 3; ++i)
    if (alleles[i] >= 0)
      out[alleles[i]]++;
}

/* ******************** GetShapeController ******************** */
int GetShapeController(const int *alleles, int numExistingAlleles) {
  if (numExistingAlleles == 0) {
    fprintf(stderr,
            "Could not get shape controller value for 0 existing alleles.\n");
    exit(-1);
  }

  int result = -1, maxAlleles, genotypeStructure[4];

  GetGenotypeStructure(alleles, genotypeStructure);
  maxAlleles = GetMaxInt(genotypeStructure, 4);

  if (maxAlleles == numExistingAlleles) {
    result = 0;
  } else if (maxAlleles * 2 == numExistingAlleles) {
    result = 1;
  }

  if (result == -1) {
    fprintf(stderr, "Error!");
    exit(-1);
  }

  return result;
}

/******************** ProcessZeroCoverageAndZeroProportion *******************/
void ProcessZeroCoverageAndZeroProportion(int cell, int site,
                                          int MLVariantNucleotide,
                                          int numDataSets) {
  int coverageFlag;

  /* analyze the reason for zero coverage in negative binomial distribution */
  if (cellExtd[cell].site[site].numReads == 0) /* if there is zero coverage */
  {
    coverageFlag = 0;
    zeroCoverage++;
    if (fixedTree == YES && numDataSets > 1) {
      zeroCoverageAll++;
    }

    if (cellExtd[cell].site[site].trueNumExistingAlleles ==
        0) /* if there is complete deletions */
    {
      coverageFlag = 1;
      completeDeletions++;
      if (fixedTree == YES && numDataSets > 1) {
        completeDeletionsAll++;
      }
    } else /* if there are alleles after evolving along the tree */
    {
      if (cellExtd[cell].site[site].observedNumExistingAlleles ==
          0) /* if there is complete ADO */
      {
        coverageFlag = 1;
        completeADOs++;
        if (fixedTree == YES && numDataSets > 1) {
          completeADOsAll++;
        }
      } else /* if there are alleles after evolution and ADO */
      {
        coverageFlag = 1;
        noCoverage++;
        allelicZeroCoverageProportion
            [cellExtd[cell].site[site].observedNumExistingAlleles][1]++;
        if (fixedTree == YES && numDataSets > 1) {
          noCoverageAll++;
        }
      }
    }

    if (coverageFlag == 0) {
      fprintf(
          stderr,
          "WARNING! Unknown reason for zero coverage at site %d in cell %d!",
          site, cell);
    }
  }

  /* analyze the reason for zero proportion in beta binomial distribution */
  if (cellExtd[cell].site[site].readCount[MLVariantNucleotide] ==
      0) /* if there is zero reads for the site-wise variant nucleotide */
  {
    zeroProportionVAF++;
    if (fixedTree == YES && numDataSets > 1)
      zeroProportionVAFAll++;

    if (cellExtd[cell].site[site].numReads >
        0) /* if there is positive coverage */
    {
      noReadsPositiveCoverage++;
      if (fixedTree == YES && numDataSets > 1)
        noReadsPositiveCoverageAll++;
    }
  }
}

char WhichNucExtd(int nucleotide) {
  if (nucleotide >= 0) {
    if (nucleotide == A)
      return ('A');
    else if (nucleotide == C)
      return ('C');
    else if (nucleotide == G)
      return ('G');
    else if (nucleotide == T)
      return ('T');
    else
      return ('N');
  } else {
    return ('-');
  }
}

/******************** ComputeMeanVarianceDoubleArr *******************/
void ComputeMeanVarianceDoubleArr(const double *arr, int size, double *meanAddr,
                                  double *varAddr, int unbiasedVar) {
  double sum, distSquared, mean, var;

  sum = distSquared = 0.0;

  for (int i = 0; i < size; ++i) {
    sum += arr[i];
  }

  mean = sum / size;

  for (int i = 0; i < size; ++i) {
    distSquared += pow(arr[i] - mean, 2);
  }

  if (unbiasedVar) {
    var = distSquared / (size - 1);
  } else {
    var = distSquared / size;
  }

  *meanAddr = mean;
  *varAddr = var;
}

void PrintGenotypeInLettersToFile(FILE *fp, int numExistingAlleles,
                                  const int *alleles, int withSpaceFlag) {
  if ((numExistingAlleles == 3) ||
      (numExistingAlleles == 2 && alleles[0] < 0)) {
    if (fp) {
      if (withSpaceFlag)
        fprintf(fp, " ");

      fprintf(fp, "%c|%c%c", WhichNucExtd(alleles[0]), WhichNucExtd(alleles[1]),
              WhichNucExtd(alleles[2]));
    }
  } else if ((numExistingAlleles == 2 && alleles[0] >= 0) ||
             numExistingAlleles == 1) {
    if (fp) {
      if (withSpaceFlag)
        fprintf(fp, " ");

      fprintf(fp, "%c|%c", WhichNucExtd(alleles[0]), WhichNucExtd(alleles[1]));
    }
  } else if (numExistingAlleles == 0) {
    if (fp) {
      if (withSpaceFlag)
        fprintf(fp, " ");

      fprintf(fp, "./.");
    }
  }
}

void GetGenotypeInNumbersToFileGivenExistingAlleles(int referenceAllele,
                                                    int stopPosition,
                                                    const int *alleles,
                                                    char missing,
                                                    char *results) {
  /* sanity check */
  if (stopPosition < 2 || stopPosition > 3) {
    fprintf(stderr,
            "ARGUMENT ERROR! Stop position should be either 2 or 3! "
            "(%d is given)\n",
            stopPosition);
    exit(-1);
  }

  for (int i = 0; i < stopPosition; ++i) {
    if (alleles[i] < 0) {
      results[i] = missing;
    } else if (alleles[i] >= 0 && alleles[i] == referenceAllele) {
      results[i] = '0';
    } else {
      results[i] = '1';
    }
  }
}

void changeToUnphasedInNumbers(int startPos, char missing, char *results) {
  // sanity check
  if (results[startPos] == missing || results[startPos + 1] == missing)
    return;

  if (results[startPos] == '1' && results[startPos + 1] == '0') {
    results[startPos] = '0';
    results[startPos + 1] = '1';
  }
}

void PrintGenotypeInNumbers(FILE *fp, int stopPosition, char separator,
                            const char *results) {
  for (int i = 0; i < stopPosition; i++) {
    if (i == 1)
      fprintf(fp, "%c", separator);

    fprintf(fp, "%c", results[i]);
  }
}

void PrintGenotypeInTernaryMatrix(FILE *fp, int startPos, int endPos,
                                  char missing, const char *results) {
  int counts[] = {0, 0, 0}; // for '0', '1', and "missing" character

  for (int i = startPos; i < endPos; i++) {
    if (results[i] == '0') {
      counts[0]++;
    } else if (results[i] == '1') {
      counts[1]++;
    } else if (results[i] == missing) {
      counts[2]++;
    }
  }

  if (endPos - startPos != counts[0] + counts[1] + counts[2]) {
    fprintf(stderr, "Unrecognized genotype");
    return;
  }

  const int comparedTo = endPos - startPos - counts[2];
  int code;

  if (counts[0] == comparedTo)
    code = 0;
  else if (counts[1] == comparedTo)
    code = 2;
  else
    code = 1;

  fprintf(fp, "%d", code);
}

void PrintGenotypeInNumbersToFile(FILE *fpNr, FILE *fpMt, int referenceAllele,
                                  int numExistingAlleles, const int *alleles,
                                  int phased, char missing,
                                  int withWhiteSpaceFlag, int useWhiteSpace) {
  if (withWhiteSpaceFlag) {
    if (fpNr)
      fprintf(fpNr, "%c", useWhiteSpace == 1 ? ' ' : '\t');

    if (fpMt)
      fprintf(fpMt, "%c", useWhiteSpace == 1 ? ' ' : '\t');
  }

  char results[] = {'0', '0', '0', '\0'};

  if ((numExistingAlleles == 3) ||
      (numExistingAlleles == 2 && alleles[0] < 0)) {

    GetGenotypeInNumbersToFileGivenExistingAlleles(referenceAllele, 3, alleles,
                                                   missing, results);

    if (!phased)
      changeToUnphasedInNumbers(1, missing, results);

    if (fpNr)
      PrintGenotypeInNumbers(fpNr, 3, phased > 0 ? '|' : '/', results);

    if (fpMt)
      PrintGenotypeInTernaryMatrix(fpMt, 0, 3, missing, results);

  } else if ((numExistingAlleles == 2 && alleles[0] >= 0) ||
             numExistingAlleles == 1) {

    GetGenotypeInNumbersToFileGivenExistingAlleles(referenceAllele, 2, alleles,
                                                   missing, results);

    if (!phased)
      changeToUnphasedInNumbers(0, missing, results);

    if (fpNr)
      PrintGenotypeInNumbers(fpNr, 2, phased > 0 ? '|' : '/', results);

    if (fpMt)
      PrintGenotypeInTernaryMatrix(fpMt, 0, 2, missing, results);

  } else if (numExistingAlleles == 0) {

    if (fpNr)
      fprintf(fpNr, "%c%c%c", missing, phased > 0 ? '|' : '/', missing);

    if (fpMt)
      fprintf(fpMt, "3");
  }
}

void FreeFiniteMuInDelSimulator(int size) {
  for (int i = 0; i < size; ++i) {
    struct NodeIntToIntArr *tmp = NULL;

    while (genotypeHashMapIntToIntArr->list[i]) {
      tmp = genotypeHashMapIntToIntArr->list[i];
      genotypeHashMapIntToIntArr->list[i] =
          genotypeHashMapIntToIntArr->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(genotypeHashMapIntToIntArr->list);
  genotypeHashMapIntToIntArr->list = NULL;
  free(genotypeHashMapIntToIntArr);
  genotypeHashMapIntToIntArr = NULL;

  for (int i = 0; i < size; ++i) {
    struct NodeIntArrToInt *tmp = NULL;

    while (genotypeHashMapIntArrToInt->list[i]) {
      tmp = genotypeHashMapIntArrToInt->list[i];
      genotypeHashMapIntArrToInt->list[i] =
          genotypeHashMapIntArrToInt->list[i]->next;
      free(tmp);
      tmp = NULL;
    }
  }
  free(genotypeHashMapIntArrToInt->list);
  genotypeHashMapIntArrToInt->list = NULL;
  free(genotypeHashMapIntArrToInt);
  genotypeHashMapIntArrToInt = NULL;
}

void ProbMtxSanityCheck(double **Pij, int nrOfGenotypes) {
  double rowSum;
  int flag = 0;

  for (int i = 0; i < nrOfGenotypes; i++) {
    rowSum = 0.0;
    for (int j = 0; j < nrOfGenotypes; j++) {
      if (Pij[i][j] < 0.0 && fabs(Pij[i][j]) < 1e-8) {
        Pij[i][j] = 0.0;
      }
      rowSum += Pij[i][j];
    }
    if (fabs(rowSum - 1.0) > 1e-8) {
      fprintf(stderr,
              "\n!!! PROBABILITY MATRIX ERROR: "
              "The sum of row %d is %.8f, not 1. \n\n",
              i, rowSum);
      flag = 1;
    }
    if (flag == YES) {
      exit(-1);
    }
  }
}

void BackgroundInfoSanityCheck(int numOfSites, int numOfCells, int numOfSNVs) {
  const int expected = (numOfSites - numOfSNVs) * numOfCells;
  int cnt[3] = {0, 0, 0}; // coverage, variant, normal

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < backgroundSize[i] + 1; j++) {
      cnt[i] += backgroundInfo[i][j];
    }

    if (cnt[i] != expected) {
      fprintf(stderr,
              "Unmatched background information size (for %s); %d expected, %d "
              "observed.",
              (i == 0 ? "coverage" : (i == 1 ? "variant" : "normal")), expected,
              cnt[i]);
    }
  }
}

void BubbleSortIntArrMinToMax(int *arr, int startPosition, int endPosition) {
  while (startPosition < endPosition) {
    for (int i = startPosition; i < endPosition - 1; ++i) {
      if (arr[i + 1] < arr[i]) {
        int tmp;
        tmp = arr[i + 1];
        arr[i + 1] = arr[i];
        arr[i] = tmp;
      }
    }
    startPosition++;
  }
}

void BubbleSortIntArrMaxToMin(int *arr, int startPosition, int endPosition) {
  while (startPosition < endPosition) {
    for (int i = startPosition; i < endPosition - 1; ++i) {
      if (arr[i + 1] > arr[i]) {
        int tmp;
        tmp = arr[i + 1];
        arr[i + 1] = arr[i];
        arr[i] = tmp;
      }
    }
    endPosition--;
  }
}

int BinarySearch(const int *arr, int start, int end, int val) {
  if (start == end)
    return -1;

  int midPoint = start + (end - start) / 2;

  if (arr[midPoint] == val)
    return midPoint;
  else if (arr[midPoint] < val)
    return BinarySearch(arr, midPoint + 1, end, val);
  else
    return BinarySearch(arr, start, midPoint, val);
}

#endif // MYSIMULATOR_FINITEMUINDELSIMULATOR_H
