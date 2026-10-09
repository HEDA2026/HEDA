#include <cstdlib>
#include <GenerateAChrom.h>
#include <unordered_set>
#include "tools.hpp"
#include "common.h"
#include "GenOperator.h"

using namespace std;

//{select two different chromosomes using the tournament method}
//it can only be used in the population where the chromosome have been sorted according fitness from good to bad
void SelectionTournament(int& parent_1, int& parent_2 , int& NumOfChromPerPop) {
    int P1 = rand() % NumOfChromPerPop;
    int P2 = rand() % NumOfChromPerPop;
    while (P1 == P2)
        P2 = rand() % NumOfChromPerPop;
    if (P1 < P2)
        parent_1 = P1;
    else
        parent_1 = P2;

    parent_2 = parent_1;
    while (parent_2 == parent_1) {
        P1 = rand() % NumOfChromPerPop;
        P2 = rand() % NumOfChromPerPop;
        while (P1 == P2)
            P2 = rand() % NumOfChromPerPop;
        if (P1 < P2)
            parent_2 = P1;
        else
            parent_2 = P2;
    }
}

void CrsMS_MP(chromosome& chrom1, chromosome& chrom2) {
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()]; //select a task randomly
        XY_SWAP(chrom1.RscAlcLst[TaskId], chrom2.RscAlcLst[TaskId], int);
    }
}

//{The MS mutation based multiple point}
void MtnMS_MP(chromosome& ch) {
    int gamma = 1 + rand() % (int(comConst.NumOfTsk / 4) + 1);
    while (gamma--) {
        int i = rand() % comConst.NumOfTsk;
        int j = rand() % Tasks[i].ElgRsc.size();
        ch.RscAlcLst[i] = Tasks[i].ElgRsc[j];
    }
}

//{HGA: single point crossover }
void CrsMS_SP(chromosome& ch1, chromosome& ch2) {
    int CrossPoint = rand() % (comConst.NumOfTsk - 1) + 1;
    for (int i = 0; i < CrossPoint; ++i) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

//{HGA:two point crossover, HGA}
void CrsMS_DP(chromosome& ch1, chromosome& ch2) {
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point1 == point2) {
        point2 = rand() % comConst.NumOfTsk;
    }
    if (point1 > point2) {
        XY_SWAP(point1, point2, int);
    }
    for (int i = point1;  i <= point2; ++i ) {
        XY_SWAP(ch1.RscAlcLst[i], ch2.RscAlcLst[i], int);
    }
}

void GnrTskSchLst_HGA(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = ch.RscAlcLst[parent];
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.TskSchLst[i] = TaskIndex;
    }
}

//｛HGA crossover｝
void Crossover_HGA(chromosome& ch1, chromosome& ch2) {
    chromosome x1 = ch1;
    chromosome x2 = ch2;
    chromosome y1 = ch1;
    chromosome y2 = ch2;
    CrsMS_SP(x1, x2);
    GnrTskSchLst_HGA(x1);
    DcdEvl(x1, true);
    GnrTskSchLst_HGA(x2);
    DcdEvl(x2, true);

    CrsMS_DP(y1, y2);
    GnrTskSchLst_HGA(y1);
    DcdEvl(y1, true);
    GnrTskSchLst_HGA(y2);
    DcdEvl(y2, true);
    vector<chromosome> sub;
    sub.push_back(x1);
    sub.push_back(x2);
    sub.push_back(y1);
    sub.push_back(y2);
    sort(sub.begin(), sub.end(), SortPopOnFitValueByAscend);
    ch1 = sub[0];
    ch2 = sub[1];
}

void Mutation_HGA(chromosome& ch) {
    chromosome x = ch;
    chromosome y = ch;
    //{single point mutation on x}
    int point = rand() % comConst.NumOfTsk;
    x.RscAlcLst[point] = Tasks[point].ElgRsc[rand() % Tasks[point].ElgRsc.size()];
    //{double point mutation on y}
    int point1 = rand() % comConst.NumOfTsk;
    int point2 = rand() % comConst.NumOfTsk;
    while (point2 == point1) {
        point2 = rand() % comConst.NumOfTsk;
    }
    y.RscAlcLst[point1] = Tasks[point1].ElgRsc[rand() % Tasks[point1].ElgRsc.size()];
    y.RscAlcLst[point2] = Tasks[point2].ElgRsc[rand() % Tasks[point2].ElgRsc.size()];
    GnrTskSchLst_HGA(x);
    DcdEvl(x, true);
    GnrTskSchLst_HGA(y);
    DcdEvl(y, true);
    if ( y.FitnessValue + PrecisionValue < x.FitnessValue ) {
        ch = y;
    } else {
        ch = x;
    }
}

//load balance improvement for HGA
void RscLoadAdjust_HGA(vector<chromosome>& Pop) {
    vector<double> lb(Pop.size());
#pragma omp parallel for
    for(int n =0 ;n < Pop.size(); ++n){
        //calculate the finish times of all task for each resource and find out the maximum
        vector<double> FT(comConst.NumOfRsc,0);
        for(int j = 0 ;j < Pop[n].RscAlcLst.size(); ++j){
            int IndexRsc = Pop[n].RscAlcLst[j];
            if(FT[IndexRsc] < Pop[n].EndTime[j]){
                FT[IndexRsc] = Pop[n].EndTime[j];
            }
        }
        //{find the minimum in FT}
        double min = InfiniteValue;
        for(int j = 0; j < comConst.NumOfRsc; ++j){
            if(min>FT[j]){
                min = FT[j];
            }
        }
        lb[n] = Pop[n].FitnessValue - min;
    }
    //{sort lb}
    vector<int> IndexLb(Pop.size());
    IndexSortByValueOnAscend(IndexLb,lb);             //sorting chromosome by lb from small to large
    //{According to LB, select the last 50% from small to large to improve}
#pragma omp parallel for
    for (int n = Pop.size() / 2; n < Pop.size(); ++n) {
        chromosome chrom = Pop[IndexLb[n]];
        vector<double> ld (comConst.NumOfRsc,0);
        vector<vector<int>> TSK(comConst.NumOfRsc);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            int RscIndex = chrom.RscAlcLst[i];
            ld[RscIndex] += (1.0 * Tasks[i].length) / Rscs[RscIndex].pc;
            TSK[RscIndex].push_back(i);
        }
        vector<int> Ind(comConst.NumOfRsc);
        IndexSortByValueOnAscend(Ind, ld);                                                           //load sort
        int BigRsc = Ind[Ind.size() - 1];                                                            //obtain the maximum
        int RandTask = TSK[BigRsc][rand() % TSK[BigRsc].size()];                                     //select a task from the Rsc with the largest load
        chrom.RscAlcLst[RandTask] = Tasks[RandTask].ElgRsc[rand() % Tasks[RandTask].ElgRsc.size()];  //reallocation
        GnrTskSchLst_HGA(chrom);
        DcdEvl(chrom, true);
        if ( chrom.FitnessValue + PrecisionValue < Pop[IndexLb[n]].FitnessValue ) {
            Pop[IndexLb[n]] = chrom;
        }
    }
}

//{LWSGA: level swapping (exchange all tasks in the level) }
void Crs_Lvl(chromosome& chrom1, chromosome& chrom2){
    int RandLevel = rand()%TskLstInLvl.size();
    //{Rsc swap}
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        int TaskIndex = TskLstInLvl[RandLevel][i];
        XY_SWAP(chrom1.RscAlcLst[TaskIndex], chrom2.RscAlcLst[TaskIndex], int);
    }
    //{find the start point of level }
    int pos = 0;
    for(int i = 0; i < RandLevel; ++i){
        pos += TskLstInLvl[i].size();
    }
    for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
        XY_SWAP(chrom1.TskSchLst[pos], chrom2.TskSchLst[pos], int);
        ++pos;
    }
}

//{LWSGA: exchange two tasks in each level }
void CrsSS_ExcTskInLvl(chromosome& chrom1, chromosome& chrom2) {
    int p1 = 0, p2 = 0;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        int TaskId = TskLstInLvl[i][rand() % TskLstInLvl[i].size()];
        //{Since the task with small level must be in front of the task with large level, the search can be started from the last recorded position}
        for (int j = p1; j < comConst.NumOfTsk; ++j) {
            if (chrom1.TskSchLst[j] == TaskId) {
                p1 = j;
                break;
            }
        }
        for (int j = p2; j < comConst.NumOfTsk; ++j) {
            if (chrom2.TskSchLst[j] == TaskId) {
                p2 = j;
                break;
            }
        }
        XY_SWAP(chrom1.TskSchLst[p1], chrom1.TskSchLst[p2], int);
        XY_SWAP(chrom2.TskSchLst[p1], chrom2.TskSchLst[p2], int);
    }
}

void Crossover_LWSGA(chromosome& ch1, chromosome& ch2) {
    int method = rand() % 3;
    if (method == 0) {
        Crs_Lvl(ch1, ch2);
    } else if (method == 1) {
        CrsSS_ExcTskInLvl(ch1, ch2);
    } else {
        CrsMS_MP(ch1, ch2);
    }
}

//(mutation: exchange two tasks in level)
void MtnSS_ExcTskInLvl(chromosome& chrom) {
    int RandLevel = rand() % TskLstInLvl.size();
    int t1 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int t2 = TskLstInLvl[RandLevel][rand() % TskLstInLvl[RandLevel].size()];
    int p1 = -1, p2 = -1;
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t1) {
            p1 = i;
            break;
        }
    }
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        if(chrom.TskSchLst[i] == t2) {
            p2 = i;
            break;
        }
    }
    XY_SWAP(chrom.TskSchLst[p1], chrom.TskSchLst[p2], int);
}

//{select a level, rearrange these tasks in this level and reallocate the resources for these tasks in this level}
void Mtn_rebuild_level(chromosome& ch) {
    int SctLvl = rand() % TskLstInLvl.size();
    vector<int> TemTskLst = TskLstInLvl[SctLvl];
    random_shuffle(TemTskLst.begin(), TemTskLst.end()); // rearrange these tasks
    int StartIndex = 0;
    for(int i = 0; i < SctLvl; ++i) {
        StartIndex = StartIndex + TskLstInLvl[i].size();
    }
    for(int i = 0;i < TemTskLst.size(); ++i) {
        int index = StartIndex + i;
        int TaskId = TemTskLst[i];
        ch.TskSchLst[index] = TaskId;
        int RscIndex = Tasks[TaskId].ElgRsc[rand() % Tasks[TaskId].ElgRsc.size()];
        ch.RscAlcLst[TaskId] = RscIndex;
    }
}

void Mutation_LWSGA(chromosome& ch) {
    int method = rand() % 3;
    if (method == 0) {
        MtnSS_ExcTskInLvl(ch);
    } else if (method == 1) {
        MtnMS_MP(ch);
    } else {
        Mtn_rebuild_level(ch);
    }
}

//{calculate the cumulative probabilities for the population whose chromosome have been sorted}
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A , int& NumOfChormPerPop) {
    for (int n = 0; n < NumOfChormPerPop; ++n) {
        A[n] = pow(RtOfSltPrb,NumOfChormPerPop-1-n) * (RtOfSltPrb - 1) / (pow(RtOfSltPrb, NumOfChormPerPop) - 1);
    }
    for (int n = 1; n < NumOfChormPerPop; ++n){
        A[n] = A[n] + A[n - 1];
    }
}

//IFBSI verison2 I/O independent -w
//{Iterative Forward and Backward Scheduling Improvement (IFBSI)}
double IFBSI(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSort(ind, OldChrom.EndTime);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.TskSchLst[comConst.NumOfTsk - 1 - i] = ind[i];
        }
        DcdEvl(NewChrom, IsFrw);
        IsFrw = !IsFrw;
    } while (NewChrom.FitnessValue + PrecisionValue < OldChrom.FitnessValue);
    if (IsFrw) { //the last is backward
        ch = OldChrom;
    } else {
        ch = NewChrom;
    }
    return ch.FitnessValue;
}

//{ Load Balancing with Communication Reduction Improvement (LBCRI)}
void LBCRI(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Id(comConst.NumOfRsc,0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        Id[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSort(ind, Id);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (abs(Id[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    } else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            int TaskIndex = TSK[RscWithMinLd][i];
            ST.insert(Tasks[TaskIndex].children.begin(),Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(),Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==
                    Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){//-w
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Id[ch.RscAlcLst[s]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.RscAlcLst[t[0].first] = RscWithMinLd;
    DcdEvl(ch, true);
    IFBSI(ch);
    if (OldCh.FitnessValue + PrecisionValue < ch.FitnessValue) {
        ch = OldCh;
    }
}

void SltCrs(vector<chromosome>& Pop,vector<chromosome>& newPopulation,vector<double>& A,int& stg) {
    bool flg = true;
    #pragma omp parallel for
    for (int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; n += 2) {
        int ind1 = SltChr(A);
        int ind2 = SltChr(A);
        while (ind1 == ind2) {
            ind2 = SltChr(A);
        }
        chromosome chrom1 = Pop[ind1];
        chromosome chrom2 = Pop[ind2];
        if(stg == 1){                                            //the evolution is in stage 1
            int MOD = rand() % 2;
            if(MOD == 0){
                if (!Crs_IL(chrom1, chrom2, stg)) {   // The improved level crossover -xy4
                    Crs_TS(chrom1, chrom2, flg, stg);
                    flg = !flg;
                }
            }
            if (MOD == 1){
                Crs_TS(chrom1, chrom2, flg, stg);     // The crossover based TS -w  //-xy4
                flg = !flg;
            }
        } else{  //the evolution is in stage 2
            int MOD = rand() % 3;
            if(MOD == 0){
                CrsMS_MP(chrom1, chrom2);                 //The MS crossover based multiple points
            }
            if(MOD == 1){
                if (!Crs_IL(chrom1, chrom2, stg)) {
                    int method = rand() % 2;
                    if(method == 0){
                        Crs_TS(chrom1, chrom2, flg, stg);
                        flg = !flg;
                    } else{
                        CrsMS_MP(chrom1, chrom2);
                    }
                }
            }
            if(MOD == 2){
                Crs_TS(chrom1, chrom2, flg, stg);
                flg = !flg;
            }
        }
        newPopulation[n] = chrom1;
        newPopulation[n + 1] = chrom2;
    }
}

//{select a chromosome using roulette wheel selection scheme}
int SltChr(vector<double>& A) {
    double lambda = RandomDouble(0, 1);
    for (int n = 0; n < A.size(); ++n)  // -xy
        if (lambda <= A[n])
            return n;
}

void Crs_TS(chromosome& ch1, chromosome& ch2, bool Flag , int& stg) {
    int gamma = 1 + rand() % (comConst.NumOfTsk - 1);
    chromosome p1;
    p1 = ch1;
    if (Flag) {//Right crossover;
        int delta = comConst.NumOfTsk-1;
        for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
            bool fd = false;        //mark whether the task at position i in the SS of ch2 is in the left half of the SS of ch1
            for (int j = gamma-1; j >= 0; --j) {
                if (ch2.TskSchLst[i] == ch1.TskSchLst[j]) {
                    fd = true;      //the task at position i in the SS of ch2 is in the left half of the SS of ch1
                    break;
                }
            }
            if(!fd){//flase
                ch1.TskSchLst[delta] = ch2.TskSchLst[i];
                ch1.RscAlcLst[ch2.TskSchLst[i]]=(2-stg)*ch1.RscAlcLst[ch2.TskSchLst[i]]+(stg-1)*ch2.RscAlcLst[ch2.TskSchLst[i]];
                delta--;
                if (delta < gamma){
                    break;
                }
            }
        }
        delta = comConst.NumOfTsk-1;
        for (int i = comConst.NumOfTsk-1; i >= 0 ; --i) {
            bool fd = false;        //mark whether the task at position i in the SS of original ch1 is in the left half of the SS of ch2
            for (int j = gamma-1; j >= 0; --j) {
                if (p1.TskSchLst[i] == ch2.TskSchLst[j]) {
                    fd = true;      //the task at position i in the SS of original ch1 is in the left half of the SS of ch2
                    break;
                }
            }
            if(!fd){
                ch2.TskSchLst[delta] = p1.TskSchLst[i];
                ch2.RscAlcLst[p1.TskSchLst[i]]=(2-stg)*ch2.RscAlcLst[p1.TskSchLst[i]]+(stg-1)*p1.RscAlcLst[p1.TskSchLst[i]];
                delta--;
                if (delta < gamma){
                    break;
                }
            }
        }
    } else {//Left crossover
        int delta = 0;
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            bool fd = false;        //  mark whether the task at position i in the SS of ch2 is in the right half of the SS of ch1
            for (int j = gamma; j < comConst.NumOfTsk; ++j) {
                if(ch2.TskSchLst[i] == ch1.TskSchLst[j]){
                    fd = true;      // the task at position i in the SS of ch2 is in the right half of the SS of ch1
                    break;
                }
            }
            if (!fd ){
                ch1.TskSchLst[delta] = ch2.TskSchLst[i];
                ch1.RscAlcLst[ch2.TskSchLst[i]]=(2-stg)*ch1.RscAlcLst[ch2.TskSchLst[i]]+(stg-1)*ch2.RscAlcLst[ch2.TskSchLst[i]];
                delta++;
                if(delta >= gamma){
                    break;
                }
            }
        }
        delta = 0;
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            bool fd = false;
            for (int j = gamma; j < comConst.NumOfTsk; ++j) {
                if(p1.TskSchLst[i] == ch2.TskSchLst[j]){
                    fd = true;
                    break;
                }
            }
            if (!fd ){
                ch2.TskSchLst[delta] = p1.TskSchLst[i];
                ch2.RscAlcLst[p1.TskSchLst[i]]=(2-stg)*ch2.RscAlcLst[p1.TskSchLst[i]]+(stg-1)*p1.RscAlcLst[p1.TskSchLst[i]];
                delta++;
                if(delta >= gamma){
                    break;
                }
            }
        }
    }
}


bool Crs_IL(chromosome& ch1, chromosome& ch2,int& stg) {
    bool scs = true;
    vector<vector<int>> IsLvl1(TskLstInLvl.size(),vector<int>(2,0));
    vector<vector<int>> IsLvl2(TskLstInLvl.size(),vector<int>(2,0));
    FndLvl(ch1,IsLvl1);
    FndLvl(ch2,IsLvl2);
    vector<int> ComLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        if (IsLvl1[i][0] * IsLvl2[i][0] == 1){
            ComLvl.push_back(i);
        }
    }
    if(ComLvl.empty()){
        scs = false;
    } else{
        int RandLevel = ComLvl[rand() % ComLvl.size()];
        int s1 = IsLvl1[RandLevel][1];
        int s2 = IsLvl2[RandLevel][1];
        for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) { // swap tasks
            XY_SWAP(ch1.TskSchLst[s1], ch2.TskSchLst[s2], int);
            ++s1;
            ++s2;
        }
        if(stg == 2){                                             // swap resources
            for (int i = 0; i < TskLstInLvl[RandLevel].size(); ++i) {
                int TaskIndex = TskLstInLvl[RandLevel][i];
                XY_SWAP(ch1.RscAlcLst[TaskIndex], ch2.RscAlcLst[TaskIndex], int);
            }
        }
    }
    return scs;
}

void Mutation(vector<chromosome>& NewSubPop,int& stg){
    if(stg == 1){
#pragma omp parallel for
        for (int n = 0; n < Parameter_TMGA.NumOfChrInSubPop; ++n) {
            if (RandomDouble(0, 1) < Parameter_TMGA.MutationRate_TMGA) {
                int MOD = rand() % 2;
                if(MOD == 0){
                    MtnSS_IL(NewSubPop[n]);
                }
                if(MOD == 1){
                    MtnSS_TS(NewSubPop[n]);
                }
            }
        }
    } else{
#pragma omp parallel for
        for (int i = 0; i < Parameter_TMGA.NumOfChrInSubPop; ++i) {
            if (RandomDouble(0, 1) < Parameter_TMGA.NumOfChrInSubPop) {
                int MOD = rand() % 3;
                if(MOD == 0){
                    MtnSS_IL(NewSubPop[i]);
                }
                if(MOD == 1){
                    MtnSS_TS(NewSubPop[i]);
                }
                if(MOD == 2){
                    MtnMS_MP(NewSubPop[i]);
                }
            }
        }
    }
}

//{The SS mutation based TS:select a task randomly,then select a different position in its valid range to insert it}
void MtnSS_TS(chromosome& ch) {
    int pos = rand() % comConst.NumOfTsk;
    int TaskID = ch.TskSchLst[pos];
    int str = pos - 1, end = pos + 1;
    while (str > -1 && (find(Tasks[TaskID].parents.begin(), Tasks[TaskID].parents.end(), ch.TskSchLst[str]) ==
                        Tasks[TaskID].parents.end()))
        --str;
    ++str;
    while (end < comConst.NumOfTsk &&
           (find(Tasks[TaskID].children.begin(), Tasks[TaskID].children.end(), ch.TskSchLst[end]) ==
            Tasks[TaskID].children.end()))
        ++end;

    if (end - str <= 1) {
        return;
    }
    int InsertPoint = rand() % (end - str) + str;
    while (InsertPoint == pos) { // select a different position
        InsertPoint = rand() % (end - str) + str;
    }
    if (InsertPoint < pos) {
        for (int i = pos; i > InsertPoint; --i) {
            ch.TskSchLst[i] = ch.TskSchLst[i - 1];
        }
        ch.TskSchLst[InsertPoint] = TaskID;
    } else {
        for (int i = pos; i < InsertPoint; ++i) {
            ch.TskSchLst[i] = ch.TskSchLst[i + 1];
        }
        ch.TskSchLst[InsertPoint] = TaskID;
    }
}

void MtnSS_IL(chromosome& ch) {
    vector<vector<int>> IsLvl(TskLstInLvl.size(),vector<int>(2,0));
    FndLvl(ch,IsLvl);
    vector<int> AvlLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        if (IsLvl[i][0] == 1){
            AvlLvl.push_back(i);
        }
    }
    if (AvlLvl.empty()) return;
    int RandLevel = AvlLvl[rand() % AvlLvl.size()];
    int s = IsLvl[RandLevel][1];
    int gamma1 = rand() % TskLstInLvl[RandLevel].size() + s;
    int gamma2 = rand() % TskLstInLvl[RandLevel].size() + s;
    while (gamma1 == gamma2) {
        gamma2 = rand() % TskLstInLvl[RandLevel].size() + s;
    }
    XY_SWAP(ch.TskSchLst[gamma1], ch.TskSchLst[gamma2], int);
}