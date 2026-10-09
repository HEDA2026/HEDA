#include <cstdlib>
#include "GenerateAChrom.h"
#include "GenOperator.h"
#include "tools.hpp"

//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}

//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] = w[TaskId] + ChildMaxRankc;
        }
    }
}

void Calculate_Rank_t(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w) {
    for(int i =1 ;i < TskLstInLvl.size(); ++i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId = TskLstInLvl[i][j];
            for (int k = 0; k < Tasks[TaskId].parents.size(); ++k) {
                int tem = Tasks[TaskId].parents[k];
                double re = w[tem] + c[tem][TaskId] + RankList[tem];
                if (RankList[TaskId] < re) {
                    RankList[TaskId] = re;
                }
            }
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.RscAlcPart.resize(comConst.NumOfTsk);
    chrom.TskSchPart.resize(comConst.NumOfTsk);
    chrom.VTskSchPart.resize(comConst.NumOfTsk,0.0);
    chrom.VRscAlcPart.resize(comConst.NumOfTsk,0.0);
    chrom.EndTime.resize(comConst.NumOfTsk);
    chrom.StartTime.resize(comConst.NumOfTsk);
    chrom.pc = 0.8;
    chrom.pm = 0.6;
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }
    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}
//{generate a task scheduling order by the levels of tasks from small to large} -xy2
//{Those haveing the same level are ranked arbitrarily among them} -xy2
vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

vector<double> GnrDecimalsByAscend() {
    vector<double> decimals(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        decimals[i] = (i + RandomDouble(0,1)) / comConst.NumOfTsk;
    }
    return decimals;
}

// I/O independent
double DcdEvl(chromosome& ch, bool IsFrw) {
    double makespan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;
        if(IsFrw) {                              //forward-loading
            for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {
                int ParentTask = Tasks[TaskIndex].parents[j];
                int ParentRsc = ch.RscAlcLst[ParentTask];
                double fft = ch.EndTime[ParentTask];
                if(RscIndex != ParentRsc) {
                    fft += ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw)); // -xy
                }
                if (ReadyTime < fft) {
                    ReadyTime = fft;
                }
            }
        } else {                                 //backward-loading
            for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                int ChildTask = Tasks[TaskIndex].children[j];
                int ChildRsc = ch.RscAlcLst[ChildTask];
                double fft = ch.EndTime[ChildTask];
                if(RscIndex != ChildRsc) {
                    fft += ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                }
                if (ReadyTime < fft) {
                    ReadyTime = fft;
                }
            }
        }
        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[RscIndex],ExecutionTime,ReadyTime); //{find an idle time-slot in ITL which can finish the task  at the earliest}
        ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExecutionTime;
        if (makespan < ch.EndTime[TaskIndex]) {
            makespan = ch.EndTime[TaskIndex];
        }
        UpdateITL(ITL[RscIndex],ch.StartTime[TaskIndex],ch.EndTime[TaskIndex]);            //{update ITL}
    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
}

double NrmDcd(chromosome& ch, bool IsFrw) {
    vector<int > upr(comConst.NumOfTsk,-1);
    list<int> RTI;
    if(IsFrw)
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].parents.size();
            if (upr[i]==0)  RTI.push_back(i);
        }
    else
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            upr[i]=Tasks[i].children.size();
            if (upr[i]==0)  RTI.push_back(i);
        }
    //generate resource allocation list and task scheduling order list
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        ch.RscAlcLst[i] = floor(ch.Code_RK[i]);
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            int decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; //小数部分最小的那个优先调度
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        //更新RTI;
        if (IsFrw)
            for (int l = 0; l < Tasks[ch.TskSchLst[i]].children.size(); ++l) {
                int childId = Tasks[ch.TskSchLst[i]].children[l];
                upr[childId] = upr[childId] - 1;
                if (upr[childId]==0)   RTI.push_back(childId);
            }
        else
            for (int l = 0; l < Tasks[ch.TskSchLst[i]].parents.size(); ++l) {
                int parentId = Tasks[ch.TskSchLst[i]].parents[l];
                upr[parentId] = upr[parentId] - 1;
                if (upr[parentId]==0)  RTI.push_back(parentId);
            }
    }
    DcdEvl(ch, IsFrw);
    return ch.FitnessValue;
}

double HrsDcd_CTP(chromosome& ch) {
    vector<double> w(comConst.NumOfTsk, 0);
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<int> ind(comConst.NumOfTsk);
    vector<vector<double>> TransferTime(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    //{calculate the transfer time between tasks when resource(Rsc) allocation has been determined}
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = floor(ch.Code_RK[i]);
        if(Tasks[i].parents.size() !=  0){
            for (int j = 0; j < Tasks[i].parents.size(); ++j) {
                int parent = Tasks[i].parents[j];
                int ParRsc = floor(ch.Code_RK[parent]);
                if(ParRsc != RscIndex){
                    TransferTime[parent][i] = ParChildTranFileSizeSum[parent][i] / VALUE * 8 / XY_MIN(Rscs[RscIndex].bw,Rscs[ParRsc].bw) ;
                }
            }
        }
        w[i] = Tasks[i].length / Rscs[RscIndex].pc;
    }
    Calculate_Rank_b(Rank_b,TransferTime, w);
    IndexSortByValueOnAscend(ind, Rank_b);
    vector<double> Decimals = GnrDecimalsByAscend();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ind[comConst.NumOfTsk - i - 1];
        ch.Code_RK[TaskIndex] = floor(ch.Code_RK[TaskIndex]) + Decimals[i];
    }
    NrmDcd(ch, true);
    return ch.FitnessValue;
}

double IFBDI(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
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

void LBCAI(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Ld(comConst.NumOfRsc,0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        Ld[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSortByValueOnAscend(ind, Ld);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (abs(Ld[RscWithMinLd]) < PrecisionValue) {
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
        if(ST.empty()){
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Ld[ch.RscAlcLst[s]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.RscAlcLst[t[0].first] = RscWithMinLd;
    DcdEvl(ch, true);
    IFBDI(ch);
    if (OldCh.FitnessValue + PrecisionValue < ch.FitnessValue) {
        ch = OldCh;
    }
}

double GnrMS_Evl(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                                   //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskId = ch.TskSchLst[i], RscId = -1;
        double FinalEndTime = InfiniteValue, FinalStartTime = 0;
        SeletRsc_EFT(ch, ITL, TaskId, RscId, FinalStartTime, FinalEndTime);
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId], FinalStartTime, FinalEndTime);     //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

double HrsDcd_EFT_ADBRKGA(chromosome& ch) {
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    vector<int > upr(comConst.NumOfTsk,0.0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1;
        double tmp = 1;
        list<int>::iterator pit;
        for (list<int>::iterator lit = RTI.begin(); lit != RTI.end(); ++lit) {
            double decimal = ch.Code_RK[*lit] - floor(ch.Code_RK[*lit]);
            if (decimal < tmp) {
                tmp = decimal; pit = lit; //小数部分最小的那个优先调度
            }
        }
        ch.TskSchLst[i] = *pit;
        RTI.erase(pit);
        double FinalEndTime = InfiniteValue;
        double FinalStartTime = 0;
        SeletRsc_EFT(ch,ITL,ch.TskSchLst[i],RscId,FinalStartTime,FinalEndTime);
        ch.EndTime[ch.TskSchLst[i]] = FinalEndTime;
        ch.Code_RK[ch.TskSchLst[i]] = RscId + tmp;
        ch.RscAlcLst[ch.TskSchLst[i]] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime); //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
        for (int l = 0; l < Tasks[ch.TskSchLst[i]].children.size(); ++l) {
            int ChildId = Tasks[ch.TskSchLst[i]].children[l];
            upr[ChildId] = upr[ChildId] - 1;
            if (upr[ChildId]==0)   RTI.push_back(ChildId);
        }
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChr_HEFT_b(vector<double> Rank_b) {
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnDescend(TemChrom.TskSchLst, Rank_b);
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

void SeletRsc_EFT(chromosome& ch, vector<set<double>>& ITL, int& TaskId, int& RscId, double& FinalStartTime, double& FinalEndTime) {
    for (int j = 0; j < Tasks[TaskId].ElgRsc.size(); ++j) {
        double ReadyTime = 0;
        int RscIdOfCrnTsk = Tasks[TaskId].ElgRsc[j];
        for (int n = 0; n < Tasks[TaskId].parents.size(); ++n) { //calculate the ready time of the task
            int PrnTskId = Tasks[TaskId].parents[n];
            int RscIdOfPrnTsk = ch.RscAlcLst[PrnTskId];
            double fft = ch.EndTime[PrnTskId];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                double TransferData = ParChildTranFileSizeSum[PrnTskId][TaskId];
                fft += TransferData / VALUE * 8 / (XY_MIN(Rscs[RscIdOfCrnTsk].bw,Rscs[RscIdOfPrnTsk].bw));
            }
            if (ReadyTime + PrecisionValue < fft){
                ReadyTime = fft;
            }
        }
        double ExeTime = Tasks[TaskId].length / Rscs[RscIdOfCrnTsk].pc;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk],ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscId = RscIdOfCrnTsk;
        }
    }
}

double FindIdleTimeSlot(set<double>& ITLofRscId,double& ExeTime,double& ReadyTime){
    set<double>::iterator pre  = ITLofRscId.begin();
    set<double>::iterator post = ITLofRscId.begin();
    ++post;
    while(post != ITLofRscId.end()) {
        if((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post)-ExeTime) {
            return  XY_MAX(*pre, ReadyTime);
        } else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}

void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime){
    if(ITLofRscId.find(StartTime) != ITLofRscId.end()) {
        ITLofRscId.erase(StartTime);
    } else {
        ITLofRscId.insert(StartTime);
    }
    if(ITLofRscId.find(EndTime) != ITLofRscId.end()) {
        ITLofRscId.erase(EndTime);
    } else {
        ITLofRscId.insert(EndTime);
    }
 }

// {in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of the number of child tasks}
chromosome GnrChr_HEFT_Baseline() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++]=TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSortByValueOnAscend(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}


chromosome GnrPrtByRank_Rnd(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        chrom.RscAlcPart[i] =RandomDouble2(0,comConst.NumOfRsc-1); //RandomDouble2(0,comConst.NumOfRsc);//rand() % comConst.NumOfRsc + rand() % 1000 / 1000.0 - 0.5;
    }
    RepairMapAndGnrRscAlcLst(chrom); //GnrRscAlcLst(chrom); //
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    DcdEvl(chrom, true);
    return chrom;
}

chromosome GnrPrtByRank_EFT(vector<double>& Rank) {
    chromosome chrom;
    IntChr(chrom);
    chrom.TskSchPart = Rank;
    RepairPriorityAndGnrSchOrd(chrom);
    GnrMS_Evl(chrom);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.RscAlcPart[i] = chrom.RscAlcLst[i] - 0.5 + (rand() % 10000) / 10000.0;
    }
    return chrom;
}

void RepairMapAndGnrRscAlcLst(chromosome& ch) {
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        int RscId = round(ch.RscAlcPart[i]);
        if(RscId < Tasks[i].ElgRsc[0]) { //超出下限的处理
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            continue;
        }
        if(RscId > Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1]) { //超出上限的处理
            ch.RscAlcPart[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1];
            ch.RscAlcLst[i] = Tasks[i].ElgRsc[Tasks[i].ElgRsc.size()-1];
            continue;
        }
        if(find(Tasks[i].ElgRsc.begin(), Tasks[i].ElgRsc.end(), RscId) == Tasks[i].ElgRsc.end()){ //不存在的处理
            if(Tasks[i].ElgRsc.size() == 1) {
                ch.RscAlcPart[i] = Tasks[i].ElgRsc[0];
                ch.RscAlcLst[i] = Tasks[i].ElgRsc[0];
            } else {
                int TemRscId = FindNearestRscId(i, ch.RscAlcPart[i]);
                ch.RscAlcPart[i] = TemRscId;
                ch.RscAlcLst[i] = TemRscId;
            }
            continue;
        }
        ch.RscAlcLst[i] = RscId;
    }
}

int FindNearestRscId(int TaskId, double value ) {
    for (int j = 0; j < Tasks[TaskId].ElgRsc.size()-1; ++j ){
        if (Tasks[TaskId].ElgRsc[j] < value && value < Tasks[TaskId].ElgRsc[j+1] ) {
            if ( Tasks[TaskId].ElgRsc[j+1] - value < value - Tasks[TaskId].ElgRsc[j] ) {
                return Tasks[TaskId].ElgRsc[j+1];
            } else {
                return Tasks[TaskId].ElgRsc[j];
            }
        }
    }
}

void RepairPriorityAndGnrSchOrd(chromosome& chrom) {
    vector<int> V, Q;
    vector<int> N(comConst.NumOfTsk, -1);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        N[i] = round(chrom.TskSchPart[i]);
    }
    vector<int> upr(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0) {
            Q.push_back(i);
        }
    }
    int MaxV = -1;
    while (V.size() != comConst.NumOfTsk) {
        for (int i = 0; i < Q.size(); ++i) {
            int TaskId = Q[i];
            int MaxP = -1;
            for (int i1 = 0; i1 < Tasks[TaskId].parents.size(); ++i1) {
                if (MaxP < N[Tasks[TaskId].parents[i1]]) {
                    MaxP = N[Tasks[TaskId].parents[i1]];
                }
            }
            if (N[TaskId] <= MaxP) {
                N[TaskId] = MaxP + 1;
            }
            for (int i1 = 0; i1 < V.size(); ++i1) {
                if (N[TaskId] == N[V[i1]])  {
                    N[TaskId] = MaxV + 1;
                    MaxV += 1;
                    break;
                }
            }
            MaxV = XY_MAX(N[TaskId], MaxV);
            V.push_back(TaskId);
        }
        vector<int> TemQ;
        for (int i = 0; i < Q.size(); ++i) {
            int taskId = Q[i];
            for (int i2 = 0; i2 < Tasks[taskId].children.size(); ++i2) {
                int childId = Tasks[taskId].children[i2];
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    TemQ.push_back(childId);
                }
            }
        }
        Q = TemQ;
    }
//    chrom.TskSchPart = N;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        chrom.TskSchPart[i] = N[i];
    }
    IndexSortByValueOnAscend(chrom.TskSchLst, N);
}

void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime){
    Parameter_HPSO.InertiaWeight = 0.1 * (1-(runtime / SchTime)) + 0.9;
    Parameter_HPSO.c1 = 2 * (1-(runtime / SchTime));
    Parameter_HPSO.c2 = 2 * (runtime / SchTime);
    double r1 = RandomDouble(0,1);
    double r2 = RandomDouble(0,1);
    for(int i = 0; i < comConst.NumOfTsk; ++i){
        ch.VTskSchPart[i] = Parameter_HPSO.InertiaWeight * ch.VTskSchPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.TskSchPart[i] - ch.TskSchPart[i])
                               + Parameter_HPSO.c2 * r2 * (Gbest.TskSchPart[i] - ch.TskSchPart[i]);
        ch.TskSchPart[i] += ch.VTskSchPart[i];

        ch.VRscAlcPart[i] = Parameter_HPSO.InertiaWeight * ch.VRscAlcPart[i] + Parameter_HPSO.c1 * r1 * (Pbest.RscAlcPart[i] - ch.RscAlcPart[i])
                               + Parameter_HPSO.c2 * r2 * (Gbest.RscAlcPart[i] - ch.RscAlcPart[i]);
        ch.RscAlcPart[i] += ch.VRscAlcPart[i];
    }
    RepairMapAndGnrRscAlcLst(ch); //GnrRscAlcLst(ch); //
    RepairPriorityAndGnrSchOrd(ch);
}

void InitProModelOfResAlc(vector<vector<double> >& PMR) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j : Tasks[i].ElgRsc) {
            PMR[i][j] =  1.0 / Tasks[i].ElgRsc.size();
        }
    }
}

void InitProModelOfTskSch(vector<vector<double> >& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants, vector<double>& Rank_b) {
    vector<int> STS(comConst.NumOfTsk ,0);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        int left  = NumOfAncestors[i];
        int right = NumOfNonDescendants[i];
        for(int j = left; j < right; ++j) {
            PMS[i][j] = 1;
            ++STS[j];
        }
    }

    for(int j = 0; j < comConst.NumOfTsk; ++j) {
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            PMS[i][j] = PMS[i][j] / STS[j];
        }
    }

//    for(int j = 0; j < comConst.NumOfTsk; ++j) { //initializing PMS based levels
////        int sum = 0;
//        double sum = 0.0;
//        for(int i = 0; i < comConst.NumOfTsk; ++i) {
////            sum = sum + (TskLstInLvl.size() - LevelIdOfTask[i]) * PMS[i][j];
//            sum = sum + Rank_b[i] * PMS[i][j];
//        }
//        for(int i = 0; i < comConst.NumOfTsk; ++i) {
////            PMS[i][j] =  (TskLstInLvl.size() - LevelIdOfTask[i]) * PMS[i][j] / sum;
//            PMS[i][j] = Rank_b[i] * PMS[i][j] / sum;
//        }
//    }

}

void GnrRscLstOfChr(chromosome& chrom, vector<vector<double> >& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double rnd = double(rand()%100) / 100;
        double sum = 0;
        for(int j: Tasks[i].ElgRsc){
            sum += PMR[i][j];
            if(rnd < sum) {
                chrom.RscAlcLst[i] = j;
                break;
            }
        }
    }
}

chromosome GnrTskLstOfChr_prp(vector<vector<double> >& PMS, vector<double>& eta_TSO) {
    chromosome chrom;
    IntChr(chrom);
    vector<int > upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double sum = 0;
        for(int k : RTI){
            sum += PMS[k][i] * eta_TSO[k];
        }
        vector<double> SltProb(comConst.NumOfTsk);
        for (int k : RTI) {
            SltProb[k] = PMS[k][i] * eta_TSO[k] / sum;
        }
        if (RandomDouble(0,1) < Parameter_HEDA.prp) {
            double MaxPrt = -1;
            for (int k : RTI)  {
                if (MaxPrt + PrecisionValue < SltProb[k]) {
                    MaxPrt = SltProb[k];
                    chrom.TskSchLst[i] = k;
                }
            }
        } else {
            double rnd = double(rand()%100) / 100;
            double ProbSum = 0;
            for (int k : RTI) {
                ProbSum += SltProb[k];
                if (rnd + PrecisionValue < ProbSum) {
                    chrom.TskSchLst[i] = k;
                    break;
                }
            }
        }
        RTI.erase(find(RTI.begin(), RTI.end(), chrom.TskSchLst[i]));
        for (int k = 0; k < Tasks[chrom.TskSchLst[i]].children.size(); ++k) {
            upr[Tasks[chrom.TskSchLst[i]].children[k]]--;
            if (upr[Tasks[chrom.TskSchLst[i]].children[k]] == 0){
                RTI.push_back(Tasks[chrom.TskSchLst[i]].children[k]);
            }
        }
    }
    return chrom;
}

void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom){
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            int count = 0;
            if (bstChrom.RscAlcLst[i] == j) {
                count = 1;
            }
            PMR[i][j] = (1 - Parameter_HEDA.theta1) * PMR[i][j] + Parameter_HEDA.theta1 * count;
        }
    }
}

void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom){
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j = 0; j < comConst.NumOfTsk; ++j) {
            int count = 0;
            if(bstChrom.TskSchLst[i] == j) {
                count = 1;
            }
            PMS[j][i] = (1-Parameter_HEDA.theta2) * PMS[j][i] + Parameter_HEDA.theta2 * count;
        }
    }
}
////new
//{calculate the average execution time of tasks}
void W_Cal_Average_S(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0, AllTransferData = Tasks[i].IFileSizeSum + Tasks[i].OFileSizeSum;
        for (int RscId : Tasks[i].ElgRsc)
            sum += Tasks[i].length / Rscs[RscId].pc + AllTransferData / VALUE * 8 / Rscs[RscId].bw;
        w[i] = sum / Tasks[i].ElgRsc.size();
    }
}

//calculate the rank of tasks based on shared database
void Calculate_Rank_b_S(vector<double>& RankList, vector<double>& ExeTime){
    for (int TskId: TskLstInLvl[TskLstInLvl.size()-1]) {
        RankList[TskId] = ExeTime[TskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int TaskId: TskLstInLvl[i]) {
            for (int Child: Tasks[TaskId].children) {
                if(RankList[TaskId] + PrecisionValue < RankList[Child] ){
                    RankList[TaskId] = RankList[Child];
                }
            }
            RankList[TaskId] = RankList[TaskId] + ExeTime[TaskId];
        }
    }
}

void Calculate_OCT_S(vector<vector<double>>& OCT){
    //按层次计算，最后一层的OCT值不用计算，为初始化值0；
    for(int i = TskLstInLvl.size()-2 ;i >= 0 ;--i){
        for (int CurTask: TskLstInLvl[i]) {
            for (int CurRsc: Tasks[CurTask].ElgRsc) {
                double FinalOct = 0;
                for (int CurChild: Tasks[CurTask].children) {
                    double MinOct = 9999999999;
                    double TransferSize = Tasks[CurChild].IFileSizeSum; //the size of input files that need to obtain from share database;
                    for (int ChildRsc: Tasks[CurChild].ElgRsc) {
                        if (CurRsc == ChildRsc) {
                            TransferSize = TransferSize - ParChildTranFileSizeSum[CurTask][CurChild];
                        }
                        double CurOct = OCT[CurChild][ChildRsc] + Tasks[CurChild].length / Rscs[ChildRsc].pc +
                                        (Tasks[CurChild].OFileSizeSum + TransferSize) / VALUE * 8 / Rscs[ChildRsc].bw;
                        if(CurOct + PrecisionValue < MinOct){
                            MinOct = CurOct;
                        }
                    }
                    if(MinOct > FinalOct + PrecisionValue){
                        FinalOct = MinOct;
                    }
                }
                OCT[CurTask][CurRsc] = FinalOct;
            }
        }
    }
}

chromosome GnrChr_PEFT_S(vector<vector<double>> OCT) {
    vector<double> rnk_oct(comConst.NumOfTsk, 0.0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            rnk_oct[i] += OCT[i][j];
        }
        rnk_oct[i] = rnk_oct[i]/Tasks[i].ElgRsc.size();
    }
    chromosome chrom;
    IntChr(chrom);
    chrom.FitnessValue = -1.0;
    vector<int > upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);  a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        list<int>::iterator it = RTI.begin();
        double maxRnk = rnk_oct[*it];
        int tskId = *it;
        ++it;
        for (it; it!=RTI.end(); ++it) {
            if ( maxRnk + PrecisionValue < rnk_oct[*it] ) {
                maxRnk = rnk_oct[*it]; tskId = *it;
            }
        }
        double FinalEndTime = -1, FinalStartTime = -1;
        int  RscId = -1;
        chrom.TskSchLst[i] = tskId;
        SeletRsc_OHEFT_S(chrom, ITL, tskId, RscId, FinalStartTime, FinalEndTime, OCT);
        chrom.EndTime[tskId] = FinalEndTime;
        chrom.StartTime[tskId] = FinalStartTime;
        chrom.RscAlcLst[tskId] = RscId;
        UpdateITL(ITL[RscId], FinalStartTime, FinalEndTime);     //{update ITL}
        chrom.FitnessValue = XY_MAX(chrom.FitnessValue, FinalEndTime);
        RTI.erase(find(RTI.begin(), RTI.end(), tskId));
        for (int childTskId: Tasks[tskId].children) {
            upr[childTskId]--;
            if (upr[childTskId] == 0){
                RTI.push_back(childTskId);
            }
        }
    }
//    double oldFit = chrom.FitnessValue;
//    if (fabs(DcdEvl_S(chrom, true)-oldFit)> PrecisionValue){
//        cout << endl << "PEFT is wrong!";
//    }
    return chrom;
}

chromosome GnrChr_HEFT_b_S(vector<double> Rank_b) {
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnDescend(TemChrom.TskSchLst,Rank_b); //直接采样降序函数-xy
//    vector<int> ind(comConst.NumOfTsk);
//    IndexSortByValueOnAscend(ind, Rank_b);
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
//    }
    GnrMS_Evl_S(TemChrom);
    return TemChrom;
}


void InitProModel(vector<vector<double>>& PMS, chromosome BstChr) {
    vector<vector<int>> flag(comConst.NumOfTsk, vector<int> (comConst.NumOfTsk, -1));
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int k = i + 1; k < comConst.NumOfTsk; ++k) {
            flag[BstChr.TskSchLst[i]][BstChr.TskSchLst[k]] = 1;   //修改命名-xy
        }
    }
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        PMS[i][i] = 1;                            //用了相乘来计算选择概率，故定义为1
        for (int k = i+1; k < comConst.NumOfTsk; ++k) {
            if (Descendants[i].find(k) != Descendants[i].end()) {
                PMS[i][k] = 1; PMS[k][i] = 0;    //k is a descendant of i, namely i is an ancestor of k
                continue;
            }
            if (Ancestors[i].find(k) != Ancestors[i].end()) {
                PMS[i][k] = 0; PMS[k][i] = 1;    //k is an ancestor of i,
                continue;
            }
            PMS[i][k] = 0.5 + Parameter_APRE_EDA.alpha * flag[i][k]; //为与论文保存一致，修改了命名-xy
            PMS[k][i] = 1 - PMS[i][k];
//            if(PMS[i][k] == -1){   //不需要if判断了！-xy
//                PMS[i][k] = 0.5 + Parameter_APRE_EDA.theta2 * flag[i][k];
//                PMS[k][i] = 1 - PMS[i][k];
//            }
        }
    }
}

void GeneratePermutation(chromosome& Chrom,vector<vector<double>>& PMS){
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0){
            RTI.push_back(i);
        }
    }
    int TaskId = 0;
    while (!RTI.empty()){
        vector<double> TempProb,NormProb;   //格式：缩进2格，命名-xy
        double Sum = 0.0;
        for (int i = 0; i < RTI.size(); ++i) {
            double sp = 1.0;
            for (int j = 0; j < RTI.size(); ++j) {
                if(i != j){
                    //sp = sp * PMS[i][j];
                    sp *= PMS[RTI[i]][RTI[j]];
                }
            }
            TempProb.push_back(sp);
            Sum += sp;
        }
        for (int i = 0; i < TempProb.size(); ++i) {
            NormProb.push_back(TempProb[i] / Sum);
        }
        double SumProb = 0;
        double rnd = double (rand()%100) / 100;
        for (int k = 0; k < RTI.size();k++) {
            SumProb += NormProb[k];
            int CurTask = RTI[k];
            if(rnd + PrecisionValue < SumProb){
                Chrom.TskSchLst[TaskId] = CurTask;
                TaskId++;
                RTI.erase(find(RTI.begin(),RTI.end(),CurTask));
                for (int ChildId: Tasks[CurTask].children) { //-xy
                    upr[ChildId] = upr[ChildId] - 1;
                    if (upr[ChildId] == 0 ){
                        RTI.push_back(ChildId);
                    }
                }
                break;
            }
        }
    }
}

void HybridDecodingMechanism(chromosome& Chrom, vector<vector<double>>& OCT){
    chromosome Chrom_EFT = Chrom, Chrom_HEFT = Chrom, Chrom_OEFT = Chrom, Chrom_OHEFT = Chrom; //格式-xy
    vector<double >Result;
    double MS_EFT = GnrMS_EFT_S(Chrom_EFT);
    double MS_HEFT = GnrMS_HEFT_S(Chrom_HEFT);
    double MS_OEFT = GnrMS_OEFT_S(Chrom_OEFT,OCT);
    double MS_OHEFT = GnrMS_OHEFT_S(Chrom_OHEFT, OCT);
    Result.push_back(MS_EFT); Result.push_back(MS_HEFT); Result.push_back(MS_OEFT); Result.push_back(MS_OHEFT);
    std::sort(Result.begin(), Result.end());
    if (fabs(Result[0] - MS_EFT) < PrecisionValue){ //加精度控制-xy
        Chrom = Chrom_EFT;
        Chrom.Label = 0;
    } else if(fabs(Result[0] - MS_HEFT) < PrecisionValue){
        //HEFT2EFT(Chrom,Chrom_HEFT);//xinshanchu
        Chrom = Chrom_HEFT;
        Chrom.Label = 1;
    } else if(fabs(Result[0] - MS_OEFT) < PrecisionValue){
        Chrom = Chrom_OEFT;
        Chrom.Label = 2;
    } else {
        Chrom = Chrom_OHEFT;
        Chrom.Label = 3;
    }
}




void LocalIntensification(chromosome& BestChrom, vector<chromosome>& ElitePop, vector<vector<double>>& OCT){
    for (int i = 0; i < Parameter_APRE_EDA.NumOfEliteOfPop; ++i) {
        chromosome InsertChrom = PathRelinking(0,BestChrom,ElitePop[i],OCT);
        if(InsertChrom.FitnessValue > PrecisionValue){
            ElitePop[i] = InsertChrom;
            if(InsertChrom.FitnessValue + PrecisionValue < BestChrom.FitnessValue){
                BestChrom = InsertChrom;
            }
        }
    }
    sort(ElitePop.begin(),ElitePop.end(), SortPopOnFitValueByAscend); // 不确定是否要排序，应该需要-xy
    chromosome SwapChrom = PathRelinking(1,BestChrom,ElitePop[0],OCT);
    if(SwapChrom.FitnessValue > PrecisionValue){
        ElitePop[0] = SwapChrom;
        if (SwapChrom.FitnessValue + PrecisionValue < BestChrom.FitnessValue) {
            BestChrom = SwapChrom;
        }
    }
}

void UpdateProModel(vector<vector<double>>& PMS, vector<chromosome>& ElitePop){
    vector<vector<int>> count(comConst.NumOfTsk,vector<int>(comConst.NumOfTsk,0));
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j = i + 1; j < comConst.NumOfTsk; ++j) {
            for(int n = 0; n < Parameter_APRE_EDA.NumOfEliteOfPop; ++n) {
                ++count[ElitePop[n].TskSchLst[i]][ElitePop[n].TskSchLst[j]];

            }
        }
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int k = i + 1; k < comConst.NumOfTsk; ++k) {
            PMS[i][k] = (1 - Parameter_APRE_EDA.alpha) * PMS[i][k] + Parameter_APRE_EDA.alpha * (1.0 * count[i][k] / Parameter_APRE_EDA.NumOfEliteOfPop);
            PMS[k][i] = 1 - PMS[i][k];
        }
    }
}

// void UpdateProModel(vector<vector<double>>& PMS, vector<chromosome>& ElitePop){
//     // 初始化一个计数矩阵
//     vector<vector<int>> count(comConst.NumOfTsk, vector<int>(comConst.NumOfTsk, 0));
//
//     // 遍历每个精英个体
//     for (int n = 0; n < Parameter_APRE_EDA.NumOfEliteOfPop; ++n) {
//         // 获取该个体的任务调度顺序 π
//         const vector<int>& sch = ElitePop[n].TskSchLst;
//
//         // 对该个体的调度顺序，统计所有 task_i 在 task_j 前面的组合
//         for (int i = 0; i < comConst.NumOfTsk; ++i) {
//             for (int j = i + 1; j < comConst.NumOfTsk; ++j) {
//                 int t1 = sch[i]; // 排在前面的任务编号
//                 int t2 = sch[j]; // 排在后面的任务编号
//                 ++count[t1][t2]; // 统计顺序 t1 ≺ t2 出现了一次
//             }
//         }
//     }
//
//     // 更新概率模型
//     for (int i = 0; i < comConst.NumOfTsk; ++i) {
//         for (int j = 0; j < comConst.NumOfTsk; ++j) {
//             if (i == j) continue;
//             double newProb = (1.0 * count[i][j]) / Parameter_APRE_EDA.NumOfEliteOfPop;
//             PMS[i][j] = (1 - Parameter_APRE_EDA.alpha) * PMS[i][j] + Parameter_APRE_EDA.alpha * newProb;
//         }
//     }
// }



void SeletRsc_OHEFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscId, double& FinalStartTime, double& FinalEndTime, vector<vector<double>>& OCT) {
    double OHEFT = InfiniteValue * 1.0;
    for (int RscIdOfCrnTsk: Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        double TransferData = Tasks[TaskIndex].ExternalInputFileSizeSum;
        for (int ParentIndex: Tasks[TaskIndex].parents) { //calculate the ready time and transfer data of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                TransferData = TransferData + ParChildTranFileSizeSum[ParentIndex][TaskIndex];
            }
            if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]){
                ReadyTime = ch.EndTime[ParentIndex];
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc + (TransferData + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIdOfCrnTsk].bw;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk],ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + OCT[TaskIndex][RscIdOfCrnTsk] + PrecisionValue < OHEFT) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            OHEFT = EndTime + OCT[TaskIndex][RscIdOfCrnTsk];
            RscId = RscIdOfCrnTsk;
        }
    }
}

double GnrMS_Evl_S(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1, TaskId = ch.TskSchLst[i];
        double FinalEndTime = 9999999999, FinalStartTime = 0;
        SeletRsc_EFT_S(ch,ITL,TaskId,RscId,FinalStartTime,FinalEndTime);  //Find the resource that can finish the task earliest
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime);              //update ITL
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

double GnrMS_EFT_S(chromosome& ch) { //无需定义可得时间段列表，直接用资源最近的可得时间（最后一个任务的完成时间）就可以了-xy
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<double> avlTime(comConst.NumOfRsc,0);      //资源最近的可得时间（最后一个任务的完成时间）
    double makespan = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskId = ch.TskSchLst[i], RscId = -1;
        double FinalEndTime = InfiniteValue, FinalStartTime = 0;
        SeletRsc_EFT_1_S(ch, avlTime, TaskId, RscId, FinalStartTime, FinalEndTime);
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        avlTime[RscId] = FinalEndTime;
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

double GnrMS_HEFT_S(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double>> ITL;                                   //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskId = ch.TskSchLst[i], RscId = -1;
        double FinalEndTime = InfiniteValue, FinalStartTime = 0;
        SeletRsc_EFT_S(ch, ITL, TaskId, RscId, FinalStartTime, FinalEndTime);
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId], FinalStartTime, FinalEndTime);     //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

double GnrMS_OEFT_S(chromosome& ch, vector<vector<double>>& OCT) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<double> avlTime(comConst.NumOfRsc, 0);
    double makespan = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskId = ch.TskSchLst[i], RscId = -1;
        double FinalEndTime = InfiniteValue, FinalStartTime = 0;
        SeletRsc_OEFT_S(ch, avlTime, TaskId, RscId, FinalStartTime, FinalEndTime, OCT);
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        avlTime[RscId] = FinalEndTime;
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

double GnrMS_OHEFT_S(chromosome& ch, vector<vector<double>>& OCT) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                                   //the idle time-slot lists  for all resources
    double makespan = -1;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(InfiniteValue * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskId = ch.TskSchLst[i], RscId = -1;
        double FinalEndTime = -1, FinalStartTime = -1;
        SeletRsc_OHEFT_S(ch, ITL, TaskId, RscId, FinalStartTime, FinalEndTime, OCT);
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId], FinalStartTime, FinalEndTime);     //{update ITL}
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

void HEFT2EFT(chromosome& Chrom, chromosome& Chrom_HEFT){
    vector<int>InsertSet;
    vector<int>FllowTaskSet;
    GnrInsertSetAndFollowTask(InsertSet,FllowTaskSet,Chrom_HEFT);
//    vector<int>InsertSet2;
//    vector<int>FllowTaskSet2;
//    GnrInsertSetAndFollowTask2(InsertSet2,FllowTaskSet2,Chrom_HEFT);
//    for (int i = 0; i < InsertSet.size(); i++) {
//        if (InsertSet[i] != InsertSet2[i] || FllowTaskSet[i] != FllowTaskSet2[i]) {
//            cout << endl << "Insertset is wrong1";
//            break;
//        }
//    }
    if (InsertSet.empty()){
        Chrom = Chrom_HEFT;
        Chrom.Label = 0;
    } else if (transfer(InsertSet,FllowTaskSet,Chrom_HEFT)){
        Chrom = Chrom_HEFT;
        Chrom.Label = 0;
    } else {
        if (produce(Chrom_HEFT,Chrom)){  //-xy
            Chrom.Label = 0;
        } else {
            Chrom = Chrom_HEFT;
            Chrom.Label = 1;
        }
    }
}

bool produce(chromosome& chrom_HEFT, chromosome& chrom){
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    list<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    chromosome temChrom;
    IntChr(temChrom);
    vector<double> avlTime(comConst.NumOfRsc, 0);
    vector<int> ind(comConst.NumOfTsk);
    IndexSortByValueOnAscend(ind, chrom_HEFT.StartTime); //按开始时间升序排序任务

    vector<list<int>>TskLstInRsc(comConst.NumOfRsc); // 二维编码
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TskLstInRsc[chrom_HEFT.RscAlcLst[ind[i]]].push_back(ind[i]);
    }

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        temChrom.TskSchLst[i] = -1;
    }
    temChrom.FitnessValue = 0;

    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) { //从队首（即集合A中）取一个就绪任务-xy
            if (TskLstInRsc[j].size() == 0 || find(RTI.begin(),RTI.end(),TskLstInRsc[j].front()) == RTI.end()) {
                continue;
            }
            int TskId = TskLstInRsc[j].front(), RscId = -1;
            double FnlEndTime = InfiniteValue, FnlStrTime = 0;
            SeletRsc_EFT_1_S(temChrom,avlTime,TskId,RscId, FnlStrTime, FnlEndTime);
            if ( j==RscId && fabs(FnlStrTime-chrom_HEFT.StartTime[TskId])<PrecisionValue && fabs(FnlEndTime-chrom_HEFT.EndTime[TskId])<PrecisionValue ) {
                temChrom.RscAlcLst[TskId] = j;
                temChrom.TskSchLst[i] = TskId;
                temChrom.StartTime[TskId] = FnlStrTime;
                temChrom.EndTime[TskId] = FnlEndTime;
                temChrom.FitnessValue = XY_MAX(temChrom.FitnessValue,FnlEndTime);
                avlTime[j] = FnlEndTime;
                TskLstInRsc[j].erase(TskLstInRsc[j].begin());
                for (int ChildId: Tasks[TskId].children) {
                    upr[ChildId] = upr[ChildId] - 1;
                    if (upr[ChildId] == 0)  RTI.push_back(ChildId);
                }
                break;
            }
        }
        if (temChrom.TskSchLst[i] == -1) {
//            cout << endl <<"fail!";
            return false;
        }
    }
    chrom = temChrom;
//    cout << endl <<"success!";
    return true;
}

bool transfer(vector<int>& InsertSet,vector<int >& FollowTaskSet, chromosome& Chrom_HEFT){
    //try
    if (InsertSet.empty()) {
        return true;
    }
    // test
    if (InsertSet.size() != FollowTaskSet.size()){
        cout << "It is wrong!!" << endl;
    }
    //论文的原始算法实现-xy
    vector<int> tskLst = Chrom_HEFT.TskSchLst;
    for (int i = 0; i < InsertSet.size(); ++i) {
//        cout<<endl<<"i="<<i;
        int insTsk = InsertSet[i], flwTsk = FollowTaskSet[i];
        for (int k = comConst.NumOfTsk-1; k >=0; --k) { //insert insTsk in the front of flwTsk
            if (tskLst[k] == insTsk) {
                for (int p = k-1; p >=0; --p) {
                    tskLst[p+1] = tskLst[p];
                    if (tskLst[p] == flwTsk) {
                        tskLst[p] = insTsk;
                        break;
                    }
                }
                break;
            }
        }
        if (isValid(tskLst)) {
            chromosome Chrom_Tem;
            IntChr(Chrom_Tem);
            Chrom_Tem.TskSchLst = tskLst;
            GnrMS_EFT_S(Chrom_Tem);
            if (Chrom_Tem.FitnessValue > Chrom_HEFT.FitnessValue + PrecisionValue) {
                return false;
            } else if (fabs(Chrom_Tem.FitnessValue - Chrom_HEFT.FitnessValue) < PrecisionValue) {
//            cout << endl<< "1: find equal!" ;
                Chrom_HEFT = Chrom_Tem;
                return true;
            } else {
//            cout << endl<< "2: find smaller!";
                vector<int> InsertSet_new;
                vector<int> FollowTaskSet_new;
                GnrMS_HEFT_S(Chrom_Tem);
                GnrInsertSetAndFollowTask(InsertSet_new, FollowTaskSet_new, Chrom_Tem);
                transfer(InsertSet_new, FollowTaskSet_new, Chrom_Tem);
//                if (InsertSet_new.empty()){
//                    Chrom_HEFT = Chrom_Tem;
//                    return true;
//                } else {
//                    transfer(InsertSet_new, FollowTaskSet_new, Chrom_Tem);
//                }

                if (transfer(InsertSet_new, FollowTaskSet_new, Chrom_Tem)) {
                    Chrom_HEFT = Chrom_Tem;
                    return true;
                } else {
                    return false;
                }


            }
        } else {
            return false;
        }
    }
}

bool isValid(vector<int >& TskList){ //修改命名-xy
//    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
//    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        upr[i] = Tasks[i].parents.size();
//        if (upr[i] == 0){
//            RTI.push_back(i);
//        }
//    }
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        if(std::find(RTI.begin(), RTI.end(),TskList[i]) == RTI.end()){
//            return false;
//        }
//        RTI.erase(find(RTI.begin(),RTI.end(),TskList[i]));
//        for (int j = 0; j < Tasks[TskList[i]].children.size(); ++j) {
//            int ChildId = Tasks[TskList[i]].children[j];
//            upr[ChildId] = upr[ChildId] - 1;
//            if (upr[ChildId] == 0 ){
//                RTI.push_back(ChildId);
//            }
//        }
//    }
//    return true;
    //直接根据定义判断,重新实现-xy
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int k = i+1; k < comConst.NumOfTsk; ++k) {
            if (find(Tasks[TskList[k]].children.begin(),Tasks[TskList[k]].children.end(),TskList[i]) != Tasks[TskList[k]].children.end()) {
                return  false;
            }
        }
    }
    return true;
}

void GnrInsertSetAndFollowTask(vector<int>& InsertSet, vector<int>& FollowTaskSet ,chromosome& ch){
    //重新实现-xy
    vector<int> Flag(comConst.NumOfTsk,0); //Flag[i]标记任务i是否是插入任务，0：否；1：是 -xy
    for (int i = 1; i < comConst.NumOfTsk; ++i) {  //标记任务i是否为插入任务；
        for (int j = 0; j < i; ++j) {
            if ( ch.StartTime[ch.TskSchLst[j]]>ch.StartTime[ch.TskSchLst[i]]+PrecisionValue && ch.RscAlcLst[ch.TskSchLst[j]]==ch.RscAlcLst[ch.TskSchLst[i]] ) {
                Flag[ch.TskSchLst[i]] = 1;
                break;
            }
        }
    }
    vector<int> ind(comConst.NumOfTsk);
    IndexSortByValueOnAscend(ind, ch.StartTime);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (Flag[ind[i]] == 1) {
            InsertSet.push_back(ind[i]);
            for (int k = i+1; k < comConst.NumOfTsk; ++k) {
//                if (ch.RscAlcLst[ind[i]] == ch.RscAlcLst[ind[k]]) {  //一个插入任务的被插任务（following task）为同一个资源上紧接着该插入任务执行的任务-xy
//                    FollowTaskSet.push_back(ind[k]);
//                    break;
//                }
                if (ch.RscAlcLst[ind[i]] == ch.RscAlcLst[ind[k]] && Flag[ind[k]] == 0) { //一个插入任务的被插任务（following task）为同一个资源上该插入任务之后执行的第一个非插入任务-xy
                    FollowTaskSet.push_back(ind[k]);
                    break;
                }
            }
        }
    }

//    struct TOS {
//        int tskId, ordId, flg;
//        double startTime;
//    };
//    vector<vector<TOS>> TskLstOfRsc( comConst.NumOfRsc );
//    for (int i = 0; i < comConst.NumOfTsk; ++i) { //把任务及其开始时间添加到各个资源的任务列表中去
//        int tskId = ch.TskSchLst[i];
//        TskLstOfRsc[ch.RscAlcLst[tskId]].push_back({tskId, i, 0, ch.StartTime[tskId]});
//    }
//    for (int j = 0; j < comConst.NumOfRsc; ++j) {
//        for (int i = 0; i < TskLstOfRsc[j].size(); ++i) {
//            for (int k = 1; k < TskLstOfRsc[j].size() - i; ++k) { //按开始时间升序排序
//                if (TskLstOfRsc[j][k-1].startTime - PrecisionValue > TskLstOfRsc[j][k].startTime) {
//                    TOS tem = TskLstOfRsc[j][k-1];
//                    TskLstOfRsc[j][k-1] = TskLstOfRsc[j][k]; TskLstOfRsc[j][k] = tem;
//                }
//            }
//        }
//    }
//
//    for (int j = 0; j < comConst.NumOfRsc; ++j) {
//        for (int p = TskLstOfRsc[j].size()-2; p >= 0 ; --p) {
//            for (int q = TskLstOfRsc[j].size()-1; q > p; --q ) {
//                if (TskLstOfRsc[j][p].ordId > TskLstOfRsc[j][q].ordId) {
//                    InsertSet.push_back(TskLstOfRsc[j][p].tskId);
//                    TskLstOfRsc[j][p].flg = 1;
//                    for (int k = p+1; k<TskLstOfRsc.size(); ++k) {
//                        if (TskLstOfRsc[j][k].flg == 0){
//                            FollowTaskSet.push_back(TskLstOfRsc[j][k].tskId);
//                            break;
//                        }
//                    }
//                    break;
//                }
//            }
//        }
//    }
}

chromosome PathRelinking(int Y, chromosome& BestChrom, chromosome& Chrom_Elite, vector<vector<double>>& OCT){
    vector<chromosome > NewPop;
    chromosome I = Chrom_Elite;
    if(Y == 0){
        NewPop = PRInsert(Chrom_Elite,BestChrom);
    } else {
        NewPop = PRSwap(Chrom_Elite,BestChrom);
    }
/*    int count = 0;
    for(chromosome CurChrom : NewPop){
        if(!isValid(CurChrom.TskSchLst)){
            count++;
        }
    }*/
    for(chromosome CurChrom: NewPop){
        if(Chrom_Elite.Label == 0){
            GnrMS_EFT_S(CurChrom);
        } else if(Chrom_Elite.Label == 1){
            GnrMS_HEFT_S(CurChrom);
        } else if(Chrom_Elite.Label == 2){
            GnrMS_OEFT_S(CurChrom, OCT);
        } else {
            GnrMS_OHEFT_S(CurChrom, OCT);
        }
        if(CurChrom.FitnessValue + PrecisionValue < I.FitnessValue){
            I = CurChrom;
        }
    }
    if(I.FitnessValue + PrecisionValue < Chrom_Elite.FitnessValue){
        return I;
    } else {
        chromosome Chrom_NULL; // 没数据的个体
        Chrom_NULL.FitnessValue = 0;
        return Chrom_NULL;
    }
}

vector<chromosome> PRInsert(chromosome& Chrom_R, chromosome& Chrom_S){
    chromosome TemCH; // 新生成的排列
    vector<chromosome > TemPop; // 新排列的集合
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Chrom_R.TskSchLst[i] != Chrom_S.TskSchLst[i]){
            TemCH = Chrom_R;
            for (int j = i + 1; j < comConst.NumOfTsk; ++j) {
                if (Chrom_R.TskSchLst[j] == Chrom_S.TskSchLst[i]){
                    TemCH.TskSchLst[i] = Chrom_S.TskSchLst[i];
                    for (int k = i + 1; k <= j; ++k) {
                        TemCH.TskSchLst[k] = Chrom_R.TskSchLst[k - 1];
                    }
                    TemPop.push_back(TemCH);
                    Chrom_R = TemCH; //此句原算法中是没有的，但由论文中插入和交换操作后调度顺序仍是合法的证明可知其必须有，否则可能会产生非法调度顺序-xy
                    break;
                }
            }
        }
    }
    return TemPop;
}

vector<chromosome> PRSwap(chromosome& Chrom_R, chromosome& Chrom_S){
    chromosome TemCH;
    vector<chromosome > TemPop;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Chrom_R.TskSchLst[i] != Chrom_S.TskSchLst[i]){
            TemCH = Chrom_R;
            for (int j = i + 1; j < comConst.NumOfTsk; ++j) {
                if (Chrom_R.TskSchLst[j] == Chrom_S.TskSchLst[i]){
                    for (int k = j; k > i; --k) {
                        TemCH.TskSchLst[k] = Chrom_R.TskSchLst[k - 1];
                        TemCH.TskSchLst[k - 1] = Chrom_R.TskSchLst[k];
                        TemPop.push_back(TemCH);
                        Chrom_R = TemCH; //此句原算法中是没有的，但由论文中插入和交换操作后调度顺序仍是合法的证明可知其必须有，否则可能会产生非法调度顺序-xy
                    }
                    break;
                }
            }
        }
    }
    return TemPop;
}

void SeletRsc_EFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime) {
    for (int RscIdOfCrnTsk: Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        double TransferData = Tasks[TaskIndex].ExternalInputFileSizeSum;
        for (int ParentIndex: Tasks[TaskIndex].parents) { //calculate the ready time and transfer data of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                TransferData = TransferData + ParChildTranFileSizeSum[ParentIndex][TaskIndex];
            }
            if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]){
                ReadyTime = ch.EndTime[ParentIndex];
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc + (TransferData + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIdOfCrnTsk].bw;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk],ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = RscIdOfCrnTsk;
        }
    }
}

void SeletRsc_EFT_1_S(chromosome& ch, vector<double>& avlTime, int& TaskIndex, int& RscId, double& FinalStartTime, double& FinalEndTime) {
    for (int RscIdOfCrnTsk : Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        double TransferData = Tasks[TaskIndex].ExternalInputFileSizeSum;
        for (int ParentIndex: Tasks[TaskIndex].parents) { //calculate the ready time and transfer data of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                TransferData = TransferData + ParChildTranFileSizeSum[ParentIndex][TaskIndex];
            }
            if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]){
                ReadyTime = ch.EndTime[ParentIndex];
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc + (TransferData + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIdOfCrnTsk].bw;
        double StartTime = XY_MAX(avlTime[RscIdOfCrnTsk],ReadyTime); //-xy
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscId = RscIdOfCrnTsk;
        }
    }
}

void SeletRsc_OEFT_S(chromosome& ch, vector<double>& avlTime, int& TaskIndex, int& RscId, double& FinalStartTime, double& FinalEndTime, vector<vector<double>>& OCT) {
    double OEFT = InfiniteValue * 1.0;
    for (int RscIdOfCrnTsk: Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        double TransferData = Tasks[TaskIndex].ExternalInputFileSizeSum;
        for (int ParentIndex: Tasks[TaskIndex].parents) { //calculate the ready time and transfer data of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                TransferData = TransferData + ParChildTranFileSizeSum[ParentIndex][TaskIndex];
            }
            if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]){
                ReadyTime = ch.EndTime[ParentIndex];
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc + (TransferData + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIdOfCrnTsk].bw;
        double StartTime = XY_MAX(avlTime[RscIdOfCrnTsk],ReadyTime);
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + OCT[TaskIndex][RscIdOfCrnTsk] + PrecisionValue < OEFT) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            OEFT = EndTime + OCT[TaskIndex][RscIdOfCrnTsk];
            RscId = RscIdOfCrnTsk;
        }
    }
}

//void HEFT2EFT(chromosome& Chrom, chromosome& Chrom_HEFT){
//    vector<int>InsertSet;
//    vector<int>FllowTaskSet;
//    GnrInsertSetAndFollowTask(InsertSet,FllowTaskSet,Chrom_HEFT);
////    vector<int>InsertSet2;
////    vector<int>FllowTaskSet2;
////    GnrInsertSetAndFollowTask2(InsertSet2,FllowTaskSet2,Chrom_HEFT);
////    for (int i = 0; i < InsertSet.size(); i++) {
////        if (InsertSet[i] != InsertSet2[i] || FllowTaskSet[i] != FllowTaskSet2[i]) {
////            cout << endl << "Insertset is wrong1";
////            break;
////        }
////    }
//    if (InsertSet.empty()){
//        Chrom = Chrom_HEFT;
//        Chrom.Label = 0;
//    } else if (transfer(InsertSet,FllowTaskSet,Chrom_HEFT)){
//        Chrom = Chrom_HEFT;
//        Chrom.Label = 0;
//    } else {
//        if (produce(Chrom_HEFT,Chrom)){  //-xy
//            Chrom.Label = 0;
//        } else {
//            Chrom = Chrom_HEFT;
//            Chrom.Label = 1;
//        }
//    }
//}
///NEW
void IntPop(vector<vector<chromosome>>& Pops,int& stg){
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);                                //calculate the average execution time of tasks
    C_Cal_Average(cc);                                //calculate the average transfer time among tasks
    Calculate_Rank_b(Rank_b,cc, ww);          //calculate the rank of tasks -w
    chromosome chrom_HEFT = GnrChr_HEFT(Rank_b);         //generate a chromosome according to HEFT
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> TemSubPop;                       //use data structure "set" for removing the same chromosome and sorting(temSubPop)
        TemSubPop.insert(chrom_HEFT);
        int nmb = 0;                                     //record the number of attempts
        while(TemSubPop.size() < Parameter_TMGA.NumOfChrInSubPop) {
            ++nmb;
            if(nmb < 2 * Parameter_TMGA.NumOfChrInSubPop) {
                chromosome TemChrom = GnrChr_Lvl_EFT();  //generate a chromosome whose SS is generated randomly based on level and MS is generated based on EFT
                TemSubPop.insert(TemChrom);
            } else if (nmb < 4 * Parameter_TMGA.NumOfChrInSubPop){//-w
                chromosome TemChrom = GnrChr_TS_EFT();   //generate a chromosome whose SS is generated randomly based on TS and MS is generated based on EFT
                TemSubPop.insert(TemChrom);
            } else {
                chromosome TemChrom = GnrChr_TS_Rnd();   //generate a chromosome whose SS is generated random-ly based on TS and MS is generated randomly
                TemSubPop.insert(TemChrom);
            }
            if( nmb >= 4 * Parameter_TMGA.NumOfChrInSubPop){
                stg = 2;
            }
        }
        vector<chromosome> NewSubPop;
        NewSubPop.assign(TemSubPop.begin(),TemSubPop.end());
        Pops.push_back(NewSubPop);
    }
}

chromosome GnrChr_HEFT(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_Lvl_EFT(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_Lvl();
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_TS_EFT(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_TS();
    GnrMS_Evl(TemChrom);
    return TemChrom;
}

chromosome GnrChr_TS_Rnd(){
    chromosome TemChrom;
    IntChr(TemChrom);
    TemChrom.TskSchLst = GnrSS_TS();
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RandInt = rand() % Tasks[i].ElgRsc.size();
        TemChrom.RscAlcLst[i] = Tasks[i].ElgRsc[RandInt];
    }
    DcdEvl(TemChrom,true);
    return TemChrom;
}

//void IntChr(chromosome& chrom) {
//    chrom.TskSchLst.resize(comConst.NumOfTsk);
//    chrom.RscAlcLst.resize(comConst.NumOfTsk);
//    chrom.EndTime.resize(comConst.NumOfTsk);
//}
//{generate a topological sort randomly}
//vector<int> GnrSS_TS() {
//    vector<int> SS;
//    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
//    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        upr[i] = Tasks[i].parents.size();
//    }
//
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        if (upr[i] == 0) RTI.push_back(i);
//    }
//    while (!RTI.empty()) {
//        int RandVec = rand() % RTI.size();
//        int v = RTI[RandVec];
//        vector<int>::iterator iter = RTI.begin() + RandVec;
//        RTI.erase(iter);
//        for (int i = 0; i < Tasks[v].children.size(); ++i) {
//            --upr[Tasks[v].children[i]];
//            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
//        }
//        SS.push_back(v);
//    }
//    return SS;
//}

void ChrExc(vector<vector<chromosome>>& Pops){
    //{form the elite populaiton EltPop }
    set<chromosome> EltPop;
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        EltPop.insert(Pops[m].begin(),Pops[m].end());
    }
    set<chromosome>::iterator iter = EltPop.begin();
    advance(iter,Parameter_TMGA.NumOfEliteOfPop);
    EltPop.erase(iter,EltPop.end());
    //{select the best N different chromosomes form the subpopulaiton and EltPop to form new subpopulation  }
    for (int m = 0; m < Parameter_TMGA.NumOfSubPop; ++m) {
        set<chromosome> NewSubPop ;
        NewSubPop.insert(Pops[m].begin(),Pops[m].end());
        NewSubPop.insert(EltPop.begin(),EltPop.end());
        set<chromosome>::iterator it = NewSubPop.begin();
        advance(it,Parameter_TMGA.NumOfChrInSubPop);
        Pops[m].assign(NewSubPop.begin(),it);
    }
}

chromosome_R2GA GnrChr_HEFT_To_R2GA(vector<double> Rank_b) {
    int n = comConst.NumOfTsk;
    int m = comConst.NumOfRsc;
    chromosome_R2GA TemChrom;
    TemChrom.TskSchPart.resize(n);
    TemChrom.RscAlcPart.resize(n);
    TemChrom.genes.resize(2 * n);

    vector<int> ind(n);
    IndexSort(ind, Rank_b);
    vector<int> uprank_tasks(n);
//    for (int i = 0; i < n; ++i) {
//        uprank_tasks[i] = ind[n - i - 1];
//    }

    chromosome HEFT_Discrete = GnrChr_HEFT_b(Rank_b);/////二，这里是３５
////////////////////////////////////////////////////////////////////////////////
// 1. 获取任务总数
    int sz = static_cast<int>(HEFT_Discrete.StartTime.size());

// 2. 创建一个索引数组 [0, 1, 2, ..., sz-1]，代表原始的任务编号位置
    std::vector<int> indices(sz);
    std::iota(indices.begin(), indices.end(), 0);

// 3. 对索引进行排序：如果索引 i 对应的开始时间比 j 小，则 i 排在前面
    std::sort(indices.begin(), indices.end(), [&HEFT_Discrete](int i, int j) {
        return HEFT_Discrete.StartTime[i] < HEFT_Discrete.StartTime[j];
    });

// 4. 按照排好序的索引，重新生成 TskSchLst
// 此时 indices[0] 就是 StartTime 最小的那个任务在原数组中的下标
    for (int k = 0; k < sz; ++k) {
        HEFT_Discrete.TskSchLst[k] = indices[k];
    }

    DcdEvl__R2GA(HEFT_Discrete, true);//三，这里是３５／／／／／／
///////////////////////////////////////////////ｃ/
    for (int i = 0; i < n; ++i) {
        uprank_tasks[i] = HEFT_Discrete.TskSchLst[i];
    }
/////////////////////ｃ//////////////////////ｃ//////////////////////ｃ//////////////////////ｃ/
    vector<int> in_degree(n);
    for (int i = 0; i < n; ++i) {
        in_degree[i] = Tasks[i].parents.size();
    }

    vector<int> candidate_tasks;
    for (int i = 0; i < n; ++i) {
        if (in_degree[i] == 0) {
            candidate_tasks.push_back(i);
        }
    }

    for (int j = 0; j < n; ++j) {
        sort(candidate_tasks.begin(), candidate_tasks.end());
        int task_id = uprank_tasks[j];
        auto it = find(candidate_tasks.begin(), candidate_tasks.end(), task_id);
        if (it == candidate_tasks.end()) {
            task_id = candidate_tasks.front();
            it = candidate_tasks.begin();
        }

        int task_pos = static_cast<int>(distance(candidate_tasks.begin(), it));
        int candidate_count = static_cast<int>(candidate_tasks.size());
        TemChrom.TskSchPart[j] = static_cast<double>(task_pos) / candidate_count +
                                 static_cast<double>(rand()) / RAND_MAX / candidate_count;

        candidate_tasks.erase(it);
        for (int child : Tasks[task_id].children) {
            --in_degree[child];
            if (in_degree[child] == 0) {
                candidate_tasks.push_back(child);
            }
        }

        int rsc_id = HEFT_Discrete.RscAlcLst[task_id];
        TemChrom.RscAlcPart[j] = static_cast<double>(rsc_id) / m +
                                 static_cast<double>(rand()) / RAND_MAX / m;
    }

    TemChrom.genes.clear();
    TemChrom.genes.insert(TemChrom.genes.end(), TemChrom.TskSchPart.begin(), TemChrom.TskSchPart.end());
    TemChrom.genes.insert(TemChrom.genes.end(), TemChrom.RscAlcPart.begin(), TemChrom.RscAlcPart.end());
    Decode_R2GA(TemChrom);
    return TemChrom;
}

// 将实数基因解码为合法的任务调度顺序和资源分配
void Decode_R2GA(chromosome_R2GA& chrom) {
    int n = comConst.NumOfTsk;
    chrom.TskSchLst.clear();
    chrom.RscAlcLst.assign(n, -1);
    chrom.StartTime.assign(n, 0.0);
    chrom.EndTime.assign(n, 0.0);

    chromosome TemChrom;
    TemChrom.TskSchLst.resize(n);
    TemChrom.RscAlcLst.resize(n);
    TemChrom.StartTime.resize(n);
    TemChrom.EndTime.resize(n);

    vector<int> in_degree(n);
    for (int i = 0; i < n; ++i) {
        in_degree[i] = Tasks[i].parents.size();
    }

    vector<int> candidate_tasks;
    for (int i = 0; i < n; ++i) {
        if (in_degree[i] == 0) {
            candidate_tasks.push_back(i);
        }
    }

    for (int j = 0; j < n; ++j) {
        sort(candidate_tasks.begin(), candidate_tasks.end());
        int size = static_cast<int>(candidate_tasks.size());
        int task_index = static_cast<int>(chrom.genes[j] * size);
        if (task_index >= size) task_index = size - 1;

        int selected_task = candidate_tasks[task_index];
        chrom.TskSchLst.push_back(selected_task);
        TemChrom.TskSchLst[j] = selected_task;

        candidate_tasks.erase(candidate_tasks.begin() + task_index);
        for (int child : Tasks[selected_task].children) {
            --in_degree[child];
            if (in_degree[child] == 0) {
                candidate_tasks.push_back(child);
            }
        }

        int rsc_index = static_cast<int>(chrom.genes[n + j] * comConst.NumOfRsc);
        if (rsc_index >= comConst.NumOfRsc) rsc_index = comConst.NumOfRsc - 1;
        chrom.RscAlcLst[selected_task] = rsc_index;
        TemChrom.RscAlcLst[selected_task] = rsc_index;
    }

    // DcdEvl(TemChrom, true);
    DcdEvl__R2GA(TemChrom, true);
    chrom.FitnessValue = TemChrom.FitnessValue;
    chrom.StartTime = TemChrom.StartTime;
    chrom.EndTime = TemChrom.EndTime;
}

double DcdEvl__R2GA(chromosome& ch, bool IsFrw) {
    double makespan = 0;

    // 非插入策略不需要追踪空闲时间槽，因此可以移除 ITL 的定义和初始化
    // vector<set<double> > ITL;
    // for (int j = 0; j < comConst.NumOfRsc; ++j) {
    //     set<double> a;
    //     a.insert(0.0); a.insert(InfiniteValue * 1.0);
    //     ITL.push_back(a);
    // }

    // 新增一个向量来记录每个资源上最后一个任务的结束时间
    vector<double> RscLastEndTime(comConst.NumOfRsc, 0.0);

    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the resource (Rsc) allocated to the task
        double ReadyTime = 0;

        if(IsFrw) {                              //forward-loading
            for (int j = 0; j < Tasks[TaskIndex].parents.size(); ++j) {
                int ParentTask = Tasks[TaskIndex].parents[j];
                int ParentRsc = ch.RscAlcLst[ParentTask];
                double fft = ch.EndTime[ParentTask];
                if(RscIndex != ParentRsc) {
                    fft += ParChildTranFileSizeSum[ParentTask][TaskIndex] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ParentRsc].bw));
                }
                if (ReadyTime < fft) {
                    ReadyTime = fft;
                }
            }
        } else {                                 //backward-loading
            for (int j = 0; j < Tasks[TaskIndex].children.size(); ++j) {
                int ChildTask = Tasks[TaskIndex].children[j];
                int ChildRsc = ch.RscAlcLst[ChildTask];
                double fft = ch.EndTime[ChildTask];
                if(RscIndex != ChildRsc) {
                    fft += ParChildTranFileSizeSum[TaskIndex][ChildTask] / VALUE * 8 / (XY_MIN(Rscs[RscIndex].bw, Rscs[ChildRsc].bw));
                }
                if (ReadyTime < fft) {
                    ReadyTime = fft;
                }
            }
        }

        double ExecutionTime = Tasks[TaskIndex].length / Rscs[RscIndex].pc;

        //  替换 FindIdleTimeSlot 函数
        // 非插入策略：任务的开始时间是其“准备就绪时间”和“资源可用时间”中的较大值
        // RscLastEndTime[RscIndex] 记录了该资源上最后一个任务的结束时间，即资源的可用时间
        ch.StartTime[TaskIndex] = max(ReadyTime, RscLastEndTime[RscIndex]);

        ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExecutionTime;

        if (makespan < ch.EndTime[TaskIndex]) {
            makespan = ch.EndTime[TaskIndex];
        }

        //替换 UpdateITL 函数
        // 更新该资源的最后结束时间，为下一个分配到此资源的任务做准备
        RscLastEndTime[RscIndex] = ch.EndTime[TaskIndex];
    }
    ch.FitnessValue = makespan;
    return ch.FitnessValue;
}