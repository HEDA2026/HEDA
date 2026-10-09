
#ifndef CSTCHANGE_GENERATEACHROM_H
#include "common.h"
#define CSTCHANGE_GENERATEACHROM_H

void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
vector<int> GnrSS_Lvl();
vector<double> GnrDecimalsByAscend();
double DcdEvl(chromosome& chrom, bool isForward);
double NrmDcd(chromosome& chrom, bool isForward);
double HrsDcd_CTP(chromosome& chrom);
double HrsDcd_EFT_ADBRKGA(chromosome& ch);
double GnrMS_Evl(chromosome& chrom);
double IFBDI(chromosome& ch);
void LBCAI(chromosome& ch);
chromosome GnrChr_HEFT_b(vector<double> Rank_b);
chromosome GnrPrtByRank_Rnd(vector<double>& Rank);
chromosome GnrPrtByRank_EFT(vector<double>& Rank);
void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskId, int& RscId, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime);
void RepairMapAndGnrRscAlcLst(chromosome& ch);
int FindNearestRscId(int TaskId, double value);
void RepairPriorityAndGnrSchOrd(chromosome& chrom);
void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime);
void InitProModelOfResAlc(vector<vector<double> >& PMR);
void InitProModelOfTskSch(vector<vector<double> >& PMS, vector<int>& NumOfAncestors,vector<int>& NumOfNonDescendants, vector<double>& Rank_b);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double> >& PMR);
chromosome GnrTskLstOfChr_prp(vector<vector<double> >& PMS, vector<double>& eta_TSO);
void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom);
void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom);
//new
void W_Cal_Average_S(vector<double>& w);
void Calculate_Rank_b_S(vector<double>& RankList, vector<double>& ExeTime);
void Calculate_OCT_S(vector<vector<double>>& OCT);
chromosome GnrChr_PEFT_S(vector<vector<double>> OCT);
chromosome GnrChr_HEFT_b_S(vector<double> Rank_b);
void InitProModel(vector<vector<double>>& PMS, chromosome BstChr);
void GeneratePermutation(chromosome& Chrom,vector<vector<double>>& PMS);
void HybridDecodingMechanism(chromosome& Chrom, vector<vector<double>>& OCT);
void LocalIntensification(chromosome& BestChrom, vector<chromosome>& ElitePop, vector<vector<double>>& OCT);
void UpdateProModel(vector<vector<double>>& PMS, vector<chromosome>& ElitePop);
void SeletRsc_OHEFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscId, double& FinalStartTime, double& FinalEndTime, vector<vector<double>>& OCT);
double GnrMS_Evl_S(chromosome& ch);
double GnrMS_EFT_S(chromosome& ch);
double GnrMS_HEFT_S(chromosome& ch);
double GnrMS_OEFT_S(chromosome& ch, vector<vector<double>>& OCT);
double GnrMS_OHEFT_S(chromosome& ch, vector<vector<double>>& OCT);
void HEFT2EFT(chromosome& Chrom, chromosome& Chrom_HEFT);
chromosome PathRelinking(int Y, chromosome& BestChrom, chromosome& Chrom_Elite, vector<vector<double>>& OCT);
void SeletRsc_EFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime);
void SeletRsc_EFT_1_S(chromosome& ch, vector<double>& avlTime, int& TaskIndex, int& RscId, double& FinalStartTime, double& FinalEndTime);
void SeletRsc_OEFT_S(chromosome& ch, vector<double>& avlTime, int& TaskIndex, int& RscId, double& FinalStartTime, double& FinalEndTime, vector<vector<double>>& OCT);

void IntPop(vector<vector<chromosome>>& populations, int& stg);
//void HEFT2EFT(chromosome& Chrom, chromosome& Chrom_HEFT);
void GnrInsertSetAndFollowTask(vector<int>& InsertSet, vector<int>& FollowTaskSet ,chromosome& ch);
bool transfer(vector<int>& InsertSet,vector<int >& FollowTaskSet, chromosome& Chrom_HEFT);

void IntChr(chromosome& chrom);
vector<int> GnrSS_TS();
vector<int> GnrSS_Lvl();
double GnrMS_Evl(chromosome& chrom);
chromosome GnrChr_HEFT(vector<double> rnk);
chromosome GnrChr_Lvl_EFT();
chromosome GnrChr_TS_EFT();
chromosome GnrChr_TS_Rnd();
void IntPop(vector<vector<chromosome>>& populations, int& stg);
void ChrExc(vector<vector<chromosome>>& populations);

bool produce(chromosome& chrom_HEFT, chromosome& chrom);
bool isValid(vector<int >& TskList);
vector<chromosome> PRInsert(chromosome& Chrom_R, chromosome& Chrom_S);
vector<chromosome> PRSwap(chromosome& Chrom_R, chromosome& Chrom_S);
void Decode_R2GA(chromosome_R2GA& chrom) ;
chromosome_R2GA GnrChr_HEFT_To_R2GA(vector<double> Rank_b);
double DcdEvl__R2GA(chromosome& ch, bool IsFrw);
// void Decode_R2GA(chromosome_R2GA& chrom) ;
//void HEFT2EFT(chromosome& Chrom, chromosome& Chrom_HEFT)
#endif //CSTCHANGE_GENERATEACHROM_H
