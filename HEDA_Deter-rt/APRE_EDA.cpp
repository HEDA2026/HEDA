//
// Created by Vince on 2022/7/4.
//
#include <common.h>
#include "tools.hpp"
#include "config.h"
#include <cfloat>
#include "GenerateAChrom.h"
//#include "GenerateAChrom.cpp"
#include "APRE_EDA.h"

double runAPRE_EDA(string XmlFile, string RscAlcFile,double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_APRE_EDA();///config.cpp//////////////////////
    CalculateLevelList();
    CalculateDescendants();
    CalculateAncestors();

    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average_S(ww);
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

    // 插入终止条件相关变量 ===
    double bestFitness = DBL_MAX;
    int terminationNum = ceil(800 * ModelScale/sqrt(comConst.NumOfTsk)/Parameter_APRE_EDA.NumOfChormPerPop);
    int NumOfNoImpGen = 0;




    double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    while (1){
        vector<chromosome> NewPopulation(Parameter_APRE_EDA.NumOfChormPerPop);
        for ( int n = 0; n < Parameter_APRE_EDA.NumOfChormPerPop; ++n) {
            IntChr(population[n]);
            GeneratePermutation(population[n],PMS);
            HybridDecodingMechanism(population[n], OCT);

            // GnrMS_Evl(population[n]);//增加算法2-2
            // DcdEvl(population[n], true);//增加算法2-1

        }
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);


        for (int i = 0; i < Parameter_APRE_EDA.NumOfEliteOfPop; ++i) {
            ElitePop[i] = population[i]; //精英种群应该为前Q个-xy
        }
        LocalIntensification(BestChrom,ElitePop,OCT);
        UpdateProModel(PMS,ElitePop);
        ++iteration;
//        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        ++NumOfNoImpGen;
        if( population[0].FitnessValue + PrecisionValue < bestFitness ){
            bestFitness = population[0].FitnessValue;
            NumOfNoImpGen = 0;
        }else {
            if( NumOfNoImpGen == terminationNum ){
                SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
                break;
            }
        }
    }
    return BestChrom.FitnessValue;
}
//parameter
//ConfigParameter APRE_EDA

