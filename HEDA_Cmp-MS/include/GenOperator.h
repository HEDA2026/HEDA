
#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H

void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void GnrTskSchLst_HGA(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& chromosomes);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_LWSGA(chromosome& chrom);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
double IFBSI(chromosome& chrom);
void LBCRI(chromosome& chrom);
int SltChr(vector<double>& A);
void Crs_TS(chromosome& ch1, chromosome& ch2, bool Flag , int& stg);
bool Crs_IL(chromosome& ch1, chromosome& ch2,int& stg);
void Mutation(vector<chromosome>& newPopulations,int& stg);
void MtnSS_TS(chromosome& a);
void MtnSS_IL(chromosome& ch);
void SltCrs(vector<chromosome>& Pop,vector<chromosome>& newPopulation,vector<double>& A,int& stg);

#endif //CSTCHANGE_CROSSOVER_H
