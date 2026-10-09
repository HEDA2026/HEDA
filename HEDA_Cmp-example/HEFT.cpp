
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "HEFT.h"
double runHEFT(string XmlFile, string RscAlcFile, double& SchTime) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);

    CalculateLevelList();
    vector<double> ww(comConst.NumOfTsk, 0.0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0.0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);

    vector<double> Rank_b(comConst.NumOfTsk, 0.0);
    Calculate_Rank_b(Rank_b,cc,ww);


    // chromosome Chrom_HEFT_b = GnrChr_HEFT_b(Rank_b);

    chromosome Chrom_HEFT_b = GnrChr_HEFT(Rank_b);

//    vector<double> Rank_t(comConst.NumOfTsk, 0.0);
//    Calculate_Rank_t(Rank_t,cc,ww);
//    chromosome Chrom_HEFT_t = GnrChr_HEFT_t(Rank_t);

    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;

    ClearALL();
    //
    // //////////////
    cout<<endl<<"HEFT"<<endl;
    // cout<<endl;
    for (int i=0;i<comConst.NumOfTsk;++i) {
        cout<<Chrom_HEFT_b.EndTime[i]<<" ";
    }
    cout<<endl;
    for (int i=0;i<comConst.NumOfTsk;++i) {
        cout<<Chrom_HEFT_b.RscAlcLst[i]+1<<" ";
    }
    cout<<endl;
    cout<<Chrom_HEFT_b.FitnessValue<<endl;
    // ////////////////
    return Chrom_HEFT_b.FitnessValue;
}