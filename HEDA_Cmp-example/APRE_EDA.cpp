//
// Created by Vince on 2022/7/4.
//
#include <common.h>
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
//#include "GenerateAChrom.cpp"
#include "APRE_EDA.h"

double runAPRE_EDA(string XmlFile, string RscAlcFile,double& SchTime, int& iteration) {
    double EndFitness=63;
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_APRE_EDA();///config.cpp//////////////////////
    CalculateLevelList();
    CalculateDescendants();
    CalculateAncestors();

    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average_S(ww);//任务平均运行时间//
    Calculate_Rank_b_S(Rank_b, ww);
    vector<vector<double>> OCT(comConst.NumOfTsk,vector<double>(comConst.NumOfRsc,0));  //定义为局部变量了-xy
    Calculate_OCT_S(OCT);  // 计算共享数据模式下的OCT值-xy//////////
    chromosome Chrom_PEFT = GnrChr_PEFT_S(OCT);
    chromosome Chrom_HEFT_b = GnrChr_HEFT_b_S(Rank_b);
    chromosome BestChrom;
    if (Chrom_PEFT.FitnessValue + PrecisionValue < Chrom_HEFT_b.FitnessValue) {
        BestChrom = Chrom_PEFT;
    } else {
        BestChrom = Chrom_HEFT_b;
    }

    //////新
    vector<double> Rank_b2(comConst.NumOfTsk,  0.0);
    W_Cal_Average_S(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b( Rank_b2,cc,ww);
    chromosome Chrom_HEFT_b2 = GnrChr_HEFT_b(Rank_b2);

    if(BestChrom.FitnessValue + PrecisionValue > Chrom_HEFT_b2.FitnessValue) {
        BestChrom = Chrom_HEFT_b2;
    }

    // 概率模型初始化
    vector<vector<double>> PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, -1));
    InitProModel(PMS,BestChrom); //-xy生成一个染色体generateAChrom//////
    vector<chromosome> population(Parameter_APRE_EDA.NumOfChormPerPop);
    vector<chromosome> ElitePop(Parameter_APRE_EDA.NumOfEliteOfPop);//classAndVarious.h////////
    double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;



    // while (RunTime < SchTime + PrecisionValue){
    while (1){
        ++iteration;
        for ( int n = 0; n < Parameter_APRE_EDA.NumOfChormPerPop; ++n) {
            IntChr(population[n]);
            GeneratePermutation(population[n],PMS);
            HybridDecodingMechanism(population[n], OCT);
        }

        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        if(population[0].FitnessValue + PrecisionValue < BestChrom.FitnessValue){
            BestChrom = population[0];

        }
///////////////////////////
        if (population[0].FitnessValue == EndFitness + 1e-2) {
            SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
            ClearALL();
            /////////////////
            cout<<endl<<"APRE_EDA"<<endl;
            for (int i=0;i<comConst.NumOfTsk;++i) {
                cout<<population[0].EndTime[i]<<" ";
            }
            cout<<endl;
            for (int i=0;i<comConst.NumOfTsk;++i) {
                cout<<population[0].RscAlcLst[i]+1<<" ";
            }
            cout<<endl;
            cout<<SchTime<<" "<<population[0].FitnessValue<<endl;
            // ////////////////
            return population[0].FitnessValue;
        }

//////////////////////////////////////
        for (int i = 0; i < Parameter_APRE_EDA.NumOfEliteOfPop; ++i) {
            ElitePop[i] = population[i]; //精英种群应该为前Q个-xy
        }
        LocalIntensification(BestChrom,ElitePop,OCT);
        UpdateProModel(PMS,ElitePop);
//        ++iteration;
        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    }

    // //////////////
    cout<<endl<<"APRE_EDA"<<endl;
    // cout<<endl;
    for (int i=0;i<comConst.NumOfTsk;++i) {
        cout<<BestChrom.EndTime[i]<<" ";
    }
    cout<<endl;
    for (int i=0;i<comConst.NumOfTsk;++i) {
        cout<<BestChrom.RscAlcLst[i]+1<<" ";
    }
    cout<<endl;
    cout<<BestChrom.FitnessValue<<endl;
    // ////////////////

    return BestChrom.FitnessValue;
}
//parameter
//ConfigParameter APRE_EDA