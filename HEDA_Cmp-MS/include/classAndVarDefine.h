#ifndef CSTCHANGE_CLASSDEFINE_H

#include "math.h"
#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>

#define CSTCHANGE_CLASSDEFINE_H
using namespace std;
//{file}
class vfile {
public:
    string FileName;     //file name
    int source;          //the source of file, -1:from the shared server; i: from task i
    double size;         //the size of file
};

//{task} delete new
class Task {
public:
    double length;
    vector<int> ElgRsc;
    vector<int> parents;
    vector<int> children;
    vector<vfile> IFile;
    vector<vfile> OFile;
    double IFileSizeSum = 0.0;
    double OrignalInputFileSizeSum = 0.0;
    double ExternalInputFileSizeSum = 0.0;
    double OFileSizeSum = 0.0;
};
//end new
////{resources}
class Resource {
public:
    vector<int> ElgTsk;     //the set of tasks which it is eligible to perform
    double pc, bw;          //processing capacity, bandwidth

    Resource(int id, double pc, double bw) {
        this->pc = pc;
        this->bw = bw;
    }
};

class chromosome {
public:
    vector<int> RscAlcLst;       //resources allocation, task-to-resource mapping, (Match String)
    vector<int> TskSchLst;       //task scheduling order (Scheduling String)  -xy4
    vector<double> Code_RK;
    vector<double> EndTime;      //the finish time of task
    vector<double> StartTime;
    vector<double> RscAlcPart;
    vector<double> TskSchPart;
    vector<double> VTskSchPart;
    vector<double> VRscAlcPart;
    double pc;
    double pm;
    int Label;
    double FitnessValue;         //Fitness (makespan) -xy4
    //{sort chromosome and remove the same chromosome according to fitness value}
    bool operator<(const chromosome &otherChromosome)const {
        return this->FitnessValue + 1e-6 < otherChromosome.FitnessValue;
    }
//    //{sort chromosome according to fitness value and remove the same chromosome according to code }
//    bool operator<(const chromosome &otherChromosome)const {
//        int flag = -1;
//        for (int i = 0; i < this->TaskOrderList.size(); i++) {
//            if (this->VMAllocationList[i] != otherChromosome.VMAllocationList[i]) {
//                flag = i;
//                break;
//            }
//            if (this->TaskOrderList[i] != otherChromosome.TaskOrderList[i]) {
//                flag = i;
//                break;
//            }
//        }
//        if (flag == -1) {
//            return false;
//        } else {
//            if(fabs(this->FitnessValue-otherChromosome.FitnessValue)<1e-6) {
//                if(this->VMAllocationList[flag] != otherChromosome.VMAllocationList[flag]){
//                    return this->VMAllocationList[flag]<otherChromosome.VMAllocationList[flag];
//                }else {
//                    return this->TaskOrderList[flag]<otherChromosome.TaskOrderList[flag];
//                }
//            }else {
//                return this->FitnessValue < otherChromosome.FitnessValue;
//            }
//        }
//    }
};

class chromosome_R2GA{
public:
    vector<double> genes;
    vector<double> TskSchPart;
    vector<double> RscAlcPart;

    vector<int> TskSchLst;
    vector<int> RscAlcLst;

    double FitnessValue;
    vector<double> StartTime;
    vector<double> EndTime;

    double pc;
    double pm;
    int Label;

    bool operator<(const chromosome_R2GA &other) const {
        return this->FitnessValue + 1e-6 < other.FitnessValue;
    }

    bool operator==(const chromosome_R2GA &other) const {
        return fabs(this->FitnessValue - other.FitnessValue) < 1e-6;
    }
};
class Paramet_CGA {
public:
    int NumOfChromPerPop;
    double CrossoverRate;      //crossover rate(CrossoverRate)
    double MutationRate;       //mutation rate (MutationRate)
};

class Paramet_HGA {
public:
    int NumOfChromPerPop;
    double EliteRate;         //crossover rate(CrossoverRate)
    double MutationRate;      //mutation rate (MutationRate)
};

class Paramet_LWSGA {
public:
    int NumOfChromPerPop;
    double CrossoverRate;     //crossover rate(CrossoverRate)
};

class Paramet_MOELS {
public:
    int NumOfChromPerPop;
    float CrossoverRate;
    float MutationRate;
};

class Paramet_HPSO {
public:
    int NumOfChromPerPop;
    double InertiaWeight;
    double c1;
    double c2;
};

class Paramet_TSEDA {
public:
    int NumOfChromPerPop;
    double theta;
//    double theta2;
    double eta;
    int NumOfEliteOfPop;      //the number of elite from population for update probability
    int NumOfImproveOfPop;
    double RunTimeRatioOfStg1;
};

class Paramet_HEDA {
public:
    int NumOfChromPerPop;
    double theta1;
    double theta2;
    double eta;
    double prp;
    double RunTimeRatioOfStg1;
};


class Paramet_ADBRKGA {
public:
    int NumOfChromPerPop;     //the number of chromosomes in each population
    double alpha;
    double beta;
    double BiasesRate;
    double ImmigrationRate;
    double ImprovementRate;
};

class ComConst {
public:
    int NumOfTsk;             //the number of Tasks
    int NumOfRsc;             //the number of resources
};
///new
class Paramet_APRE_EDA {
public:
    int NumOfChormPerPop;
    int NumOfEliteOfPop;      //the number of elite from population for update probability
    double alpha;
};

class Paramet_TMGA {
public:
    int NumOfSubPop;          //the number of subpopulations
    int NumOfChrInSubPop;     //the number of chromosomes in each subpopulation
    int NumOfChrImp;          //the number of chromosomes need to be improved in each subpopulation
    int interval;             //the interval for exchange
    int NumOfEliteOfPop;      //the number of elite from populaiotns for exchange
    double MutationRate_TMGA;
    int TrmThresholdOfStg1;   //(TrmThresholdOfStg1)int
};

class Paramet_R2GA {
public:
    int NumOfChromPerPop;      // 种群规模 (Population size)
    double EliteRate;          // 精英保留比例 (Elite rate)
    int NumOfElite;            // 精英个体数量 (Number of elites)
    double MutationRate;       // 变异概率 (Mutation rate)
    double CrossoverRate;

};




//end new


extern vector<Task> Tasks;
//extern vector<Task> OriginalTasks;
//extern vector<int> Oid;
//extern vector<int> Nid;
extern vector<vector<int> > TskLstInLvl; //task list (set) in each level
extern vector<int> LevelIdOfTask;        //the level of task
extern vector<vector<double> > ParChildTranFileSizeSum;
extern vector<Resource> Rscs;
extern double MinBW ;
extern vector<double> MaxLd;
//extern vector<chromosome> population;
//extern vector<vector<chromosome> > populations;
extern vector<set<int> > Descendants;
extern vector<set<int> > Ancestors;
extern ComConst comConst;
extern double ModelScale;
extern vector<vector<chromosome>> populations;
extern Paramet_CGA Parameter_CGA;
extern Paramet_HGA Parameter_HGA;
extern Paramet_LWSGA Parameter_LWSGA;
extern Paramet_MOELS Parameter_MOELS;
extern Paramet_HPSO Parameter_HPSO;
extern Paramet_TSEDA Parameter_TSEDA;
extern Paramet_HEDA Parameter_HEDA;
extern Paramet_ADBRKGA Parameter_ADBRKGA;
//new
extern Paramet_APRE_EDA Parameter_APRE_EDA;
extern Paramet_TMGA Parameter_TMGA;
extern Paramet_R2GA Parameter_R2GA;
//extern vector<vector<double>> InitPheromone_S;
#endif //CSTCHANGE_CLASSDEFINE_H
//ndefined reference to Parameter APRE EDA

