////
//// Created by xieyi on 2022/5/29.
////
//
//#include "HEDA.h"
//#include "tools.hpp"
//#include "config.h"
//#include "GenerateAChrom.h"
//
//double runHEDA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
//    double RunTime = 0;                    //the variable for recording the scheduling(running) time
//
//    clock_t start = clock();
//
//    ReadFile(XmlFile, RscAlcFile);
//    ConfigParameter_HEDA();
//    CalculateLevelList();
//    CalculateDescendants();
//    CalculateAncestors();
//
//    vector<double> ww(comConst.NumOfTsk, 0);
//    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
//    vector<double> Rank_b(comConst.NumOfTsk, 0);
//    W_Cal_Average(ww);
//    C_Cal_Average(cc);
//    Calculate_Rank_b(Rank_b, cc, ww);
//
//    double MaxRank_b = Rank_b[0];
//    for (int i = 1; i < comConst.NumOfTsk; ++i) {
//        if (MaxRank_b + PrecisionValue < Rank_b[i]) {
//            MaxRank_b = Rank_b[i];
//        }
//    }
//
//    //初始化任务调度顺序概率模型
//    vector<int> NumOfAncestors(comConst.NumOfTsk);
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        NumOfAncestors[i] = Ancestors[i].size();
//    }
//    //递归计算，子孙任务的数量
//    vector<int> NumOfDescendants(comConst.NumOfTsk);
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        NumOfDescendants[i] = Descendants[i].size();
//    }
//    //总任务数-子孙任务数量
//    vector<int> NumOfNonDescendants(comConst.NumOfTsk);
//    for (int i = 0; i < comConst.NumOfTsk; ++i) {
//        NumOfNonDescendants[i] = comConst.NumOfTsk - NumOfDescendants[i];
//    }
//    int flag=0;
//    int flag2=0;
//    while(1) {
//        chromosome TemChrom1=GnrChr_TS_Rnd();
//        // chromosome TemChrom2=TemChrom1;
//        // IFBDI(TemChrom2);
//        chromosome TemChrom3=TemChrom1;
//        // TemChrom3.TskSchLst={1,4,3,6,2,8,5,7};//{0,3,2,5,1,7,4,6};
//        // TemChrom3.RscAlcLst={3,2,3,2,3,2,1,2};//{2,1,2,1,2,1,0,1};
//        TemChrom3.TskSchLst={0,3,2,5,1,7,4,6};
//        TemChrom3.RscAlcLst={2,1,2,1,2,1,0,1};
//
//        // TemChrom3.TskSchLst={6,3,21,4,0,1,26,5,11,8,24,29,25,31,9,22,28,16,7,2,33,10,30,14,17,15,27,34,13,19,20,12,32,35,23,18,36,37,43,40,42,38,39,45,44,41,46,47,48,49};
//        // TemChrom3.RscAlcLst={2,3,5,4,4,5,4,4,3,4,2,1,4,4,4,3,1,0,4,2,5,1,3,5,4,3,4,5,2,4,5,2,3,3,1,2,0,1,4,2,0,5,3,4,4,3,3,0,3,4};
//
//        // TemChrom3.TskSchLst={0,1,2,4,5,11,24,9,3,20,8,19,18,42,40,85,6,10,7,14,22,17,46,36,73,93,15,31,63,82,64,86,37,128,76,75,139,29,60,59,16,30,61,123,34,12,69,41,83,26,84,148,54,94,33,68,67,192,131,62,126,39,32,140,125,190,79,38,80,143,77,78,25,141,207,35,71,52,23,105,81,146,21,43,106,88,124,45,91,147,44,149,90,213,50,13,27,127,89,191,66,65,28,129,150,193,154,72,257,212,214,278,49,100,157,110,145,211,256,120,320,74,384,137,70,144,195,87,135,133,197,199,261,263,254,318,151,188,56,136,252,316,209,134,205,273,132,198,130,196,262,48,189,97,51,152,99,187,98,251,161,163,225,216,102,253,215,153,138,109,217,174,317,55,53,112,107,58,381,47,315,108,379,142,172,96,443,119,236,95,159,169,255,111,204,201,101,319,321,184,165,383,210,227,385,208,449,104,206,103,203,167,274,272,338,270,231,202,268,114,113,334,178,176,380,240,445,271,269,335,507,259,323,387,194,121,122,185,92,223,117,160,155,221,219,158,118,285,156,182,222,177,246,280,265,344,183,267,170,57,331,115,116,180,238,242,162,168,302,164,232,226,171,327,275,337,339,300,364,403,401,258,467,322,325,234,389,391,296,186,455,233,333,250,314,397,295,248,312,376,304,166,298,336,362,200,228,264,290,277,266,341,281,287,332,179,330,279,181,224,396,343,243,451,218,282,407,346,328,394,345,460,409,473,329,393,276,398,340,342,402,406,244,404,306,470,260,288,368,366,245,430,326,247,309,310,308,372,428,324,388,370,434,374,249,311,410,436,468,498,283,373,378,347,411,475,440,426,442,173,235,229,175,237,241,299,382,239,444,504,301,363,230,289,351,292,220,386,284,307,297,509,294,371,361,446,508,448,303,405,510,469,360,424,408,305,365,450,359,453,369,435,348,313,433,472,488,286,432,377,400,350,375,390,441,356,454,466,439,505,447,395,492,506,354,494,293,496,358,438,500,357,490,349,423,422,502,486,399,413,465,418,291,355,412,352,416,459,414,503,462,480,471,420,392,484,367,437,464,421,431,478,499,429,457,353,474,415,482,427,497,491,461,417,485,493,456,452,479,458,476,425,487,489,501,419,481,463,477,483,495,511};
//        // TemChrom3.RscAlcLst={0,0,4,4,0,4,0,3,5,5,2,2,0,2,1,5,0,3,1,2,0,5,0,0,1,0,0,1,1,4,5,3,0,5,5,2,3,2,5,0,1,3,2,3,3,2,4,3,2,4,1,1,4,5,1,2,1,2,0,1,1,5,4,1,0,4,4,4,2,1,2,1,0,1,3,0,4,2,2,4,4,0,1,1,4,1,1,4,0,0,0,0,3,1,1,4,5,3,5,5,0,4,0,2,5,1,1,5,5,5,0,0,1,4,3,5,4,3,3,2,0,4,5,1,1,3,1,4,1,4,4,1,4,5,3,0,0,0,0,1,0,3,2,0,1,4,0,5,1,4,3,0,4,4,0,4,5,0,2,0,0,3,2,4,4,5,3,3,5,2,3,2,3,2,4,3,4,0,5,5,3,5,5,0,0,1,3,3,3,2,1,1,1,4,0,2,4,4,5,4,4,1,3,4,2,1,4,1,5,3,2,1,4,1,5,1,1,1,3,4,5,5,1,1,0,0,1,1,5,1,2,1,1,2,4,4,1,5,4,4,1,2,0,1,3,5,0,5,5,3,5,1,2,0,4,2,2,1,2,5,3,0,0,1,3,1,4,4,5,4,4,5,4,0,2,4,4,5,0,3,4,0,3,3,3,1,3,1,1,4,1,5,2,4,2,0,1,3,3,3,3,4,1,5,0,4,4,5,3,0,3,2,3,5,3,4,1,3,1,3,5,2,0,5,5,4,1,0,5,4,4,1,3,4,0,0,0,2,0,4,5,3,2,4,0,5,0,3,4,5,0,0,3,4,4,0,4,0,3,3,5,1,1,1,1,3,2,5,2,2,3,0,3,1,4,3,1,4,4,4,0,4,5,5,1,5,2,1,1,5,5,5,4,0,4,3,2,1,1,2,4,4,0,3,3,4,5,2,3,3,1,4,2,3,4,4,5,4,4,2,5,1,2,1,2,3,0,5,0,5,0,5,2,1,2,1,5,2,5,2,0,2,3,1,2,4,4,5,4,0,3,1,2,3,0,0,5,2,2,3,1,5,3,2,3,5,4,0,1,3,1,5,4,0,3,1,4,2,4,0,4,5,2,4,5,1,4,4,0,3,4,3,1,5,5,5,2,0,0,0,3,5,1,2,0,3,5,1,3,3,4,1};
//
//        // TemChrom3.TskSchLst={0,22,26,25,28,27,17,6,18,9,5,1,12,16,4,19,21,24,10,29,54,46,37,8,7,23,11,56,3,36,14,35,13,50,39,52,41,49,40,38,34,42,53,32,45,33,31,20,47,51,44,55,2,15,48,43,30,57,61,82,81,63,65,76,75,78,58,71,79,67,62,60,84,69,66,104,74,92,101,97,88,72,105,68,86,80,94,106,107,83,95,109,59,64,91,89,108,77,93,102,103,87,90,98,73,99,85,70,96,110,130,133,115,116,134,132,114,123,128,120,119,113,127,129,124,100,111,112,122,117,118,125,135,153,139,149,157,137,126,121,136,144,147,151,141,145,156,152,148,154,158,143,150,159,176,180,174,146,170,173,168,160,167,182,198,140,164,131,195,177,138,161,192,175,171,196,183,179,181,190,202,162,193,199,197,204,155,216,189,172,184,210,201,203,222,220,169,224,191,223,213,194,215,219,217,212,211,163,214,186,178,142,166,188,209,205,218,225,244,207,238,235,232,230,185,233,234,229,165,237,236,231,187,243,208,227,239,242,200,221,241,228,206,240,226,245,247,262,246,257,264,248,249,267,253,255,259,250,277,268,263,275,258,256,260,252,280,270,266,276,274,254,281,251,269,271,272,273,265,278,282,291,293,284,297,261,286,279,294,298,285,290,289,288,283,295,299,302,314,305,301,304,306,311,300,307,313,315,321,317,287,310,320,329,322,326,319,316,330,309,335,340,331,325,336,333,344,303,324,334,353,348,349,339,343,328,296,338,342,356,347,312,352,292,318,351,308,346,323,337,350,332,345,357,365,363,368,361,355,367,359,360,358,362,364,327,341,369,378,376,374,370,373,375,379,371,372,380,389,354,388,385,384,381,366,377,383,386,387,390,396,398,395,397,382,393,391,394,392,399,404,406,400,402,405,403,407,401,410,409,408,413,411,412,414,417,419,416,415,420,421,425,422,426,418,424,423,428,427,429,430,432,431,433};
//        // TemChrom3.RscAlcLst={1,5,3,5,4,5,2,5,5,3,1,3,0,3,2,3,2,1,2,1,5,5,4,5,4,4,4,2,3,5,3,2,3,1,2,2,2,2,0,4,4,0,1,5,5,4,4,0,4,5,2,1,2,1,2,1,5,4,3,5,4,4,0,3,1,1,0,5,0,3,0,1,1,2,5,5,5,3,2,5,1,5,0,0,3,0,3,4,5,0,4,5,0,0,0,0,3,3,4,3,0,5,2,5,1,5,3,2,4,5,2,3,2,5,4,1,1,5,5,3,5,0,3,3,5,0,1,2,1,3,4,4,2,5,0,0,4,2,3,4,4,4,0,1,3,5,4,0,0,0,5,2,0,4,4,0,1,2,0,2,3,3,0,1,1,1,3,4,5,3,0,2,2,3,1,0,3,4,4,3,4,4,1,2,2,4,2,1,0,1,5,5,2,4,1,1,4,0,3,5,4,4,4,1,5,4,5,1,2,2,3,1,3,4,5,2,3,0,0,1,4,4,0,2,4,0,5,5,4,2,1,5,2,3,5,2,0,4,5,4,1,0,1,4,3,4,0,5,4,2,2,4,3,1,4,2,3,2,2,0,4,3,3,3,4,1,2,4,3,1,3,2,1,1,2,3,0,0,2,2,4,0,0,1,1,4,1,4,2,1,5,4,2,5,1,5,2,0,2,4,4,3,5,0,2,3,5,1,2,4,5,0,4,1,2,2,4,5,4,2,0,1,0,0,0,4,4,4,1,5,1,1,4,4,0,1,3,0,2,1,5,2,0,4,4,5,3,5,0,0,3,1,3,0,3,0,1,0,0,1,3,3,3,4,0,1,2,2,5,1,2,2,5,1,5,1,4,1,2,5,1,0,2,5,4,5,5,3,2,4,1,3,5,0,5,4,3,2,5,4,4,3,2,1,1,2,4,2,2,0,3,3,5,4,0,2,2,3,3,5,4,4,3,3,0,4,1,3,3,2,4,4,5,4};
//
//        DcdEvl(TemChrom3,true);
//        chromosome TemChrom4=TemChrom3;
//        LBCAI(TemChrom3);
//        // if (TemChrom1.FitnessValue>TemChrom2.FitnessValue+ 1e-2) {
//        //     flag=1;
//        //     cout<<"origin"<<endl;
//        //     for (int i=0;i<comConst.NumOfTsk;++i) {
//        //         cout<<TemChrom1.TskSchLst[i]+1<<" ";
//        //     }
//        //     cout<<endl;
//        //     for (int i=0;i<comConst.NumOfTsk;++i) {
//        //         cout<<TemChrom1.EndTime[i]<<" ";
//        //     }
//        //     cout<<endl;
//        //     for (int i=0;i<comConst.NumOfTsk;++i) {
//        //         cout<<TemChrom1.RscAlcLst[i]+1<<" ";
//        //     }
//        //     cout<<endl;
//        //     cout<<SchTime<<" "<<TemChrom1.FitnessValue<<endl;
//        //
//        // }
//        if (TemChrom4.FitnessValue>TemChrom3.FitnessValue+ 1e-2) {
//            flag2++;
//            // if (TemChrom3.FitnessValue < TemChrom3.FitnessValue + 1e-2) {
//            cout<<"origin"<<endl;
//            for (int i=0;i<comConst.NumOfTsk;++i) {
//                cout<<TemChrom4.TskSchLst[i]+1<<" ";
//            }
//            cout<<endl;
//            for (int i=0;i<comConst.NumOfTsk;++i) {
//                cout<<TemChrom4.EndTime[i]<<" ";
//            }
//            cout<<endl;
//            for (int i=0;i<comConst.NumOfTsk;++i) {
//                cout<<TemChrom4.RscAlcLst[i]+1<<" ";
//            }
//            cout<<endl;
//            cout<<SchTime<<" "<<TemChrom4.FitnessValue<<endl;
//            //
//                cout<<"LBCAI"<<endl;
//                for (int i=0;i<comConst.NumOfTsk;++i) {
//                    cout<<TemChrom3.TskSchLst[i]+1<<" ";
//                }
//                cout<<endl;
//                for (int i=0;i<comConst.NumOfTsk;++i) {
//                    cout<<TemChrom3.EndTime[i]<<" ";
//                }
//                cout<<endl;
//                for (int i=0;i<comConst.NumOfTsk;++i) {
//                    cout<<TemChrom3.RscAlcLst[i]+1<<" ";
//                }
//                cout<<endl;
//                // cout<<SchTime<<" "<<TemChrom3.FitnessValue<<endl;
//                // return Chrom_lb.FitnessValue;
//            // }
//        }
//        if (flag2==1) {
//            // cout<<SchTime<<" qq"<<TemChrom3.FitnessValue<<endl;
//            RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
//            SchTime = RunTime;
//            cout<<SchTime<<" "<<TemChrom3.FitnessValue<<endl;
//            break;
//        }
//    }
//    return 0;
//}



 //
 // Created by xieyi on 2022/5/29.
 //

 #include "HEDA.h"
 #include "tools.hpp"
 #include "config.h"
 #include "GenerateAChrom.h"

 double runHEDA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
     double EndFitness=54;
     double RunTime = 0;
     clock_t start = clock();
     ReadFile(XmlFile, RscAlcFile);
     ConfigParameter_HEDA();
     CalculateLevelList();
     CalculateDescendants();
     CalculateAncestors();

     vector<double> ww(comConst.NumOfTsk, 0);
     vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
     vector<double> Rank_b(comConst.NumOfTsk, 0);
     W_Cal_Average(ww);
     C_Cal_Average(cc);
     Calculate_Rank_b(Rank_b, cc, ww);

     double MaxRank_b = Rank_b[0];
     for (int i = 1; i < comConst.NumOfTsk; ++i) {
         if (MaxRank_b + PrecisionValue < Rank_b[i]) {
             MaxRank_b = Rank_b[i];
         }
     }

     //初始化任务调度顺序概率模型
     vector<int> NumOfAncestors(comConst.NumOfTsk);
     for (int i = 0; i < comConst.NumOfTsk; ++i) {
         NumOfAncestors[i] = Ancestors[i].size();
     }
     //递归计算，子孙任务的数量
     vector<int> NumOfDescendants(comConst.NumOfTsk);
     for (int i = 0; i < comConst.NumOfTsk; ++i) {
         NumOfDescendants[i] = Descendants[i].size();
     }
     //总任务数-子孙任务数量
     vector<int> NumOfNonDescendants(comConst.NumOfTsk);
     for (int i = 0; i < comConst.NumOfTsk; ++i) {
         NumOfNonDescendants[i] = comConst.NumOfTsk - NumOfDescendants[i];
     }

     vector<vector<double> > PMR(comConst.NumOfTsk, vector<double>(comConst.NumOfRsc, 0));
     InitProModelOfResAlc(PMR);
     vector<vector<double> > PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
     InitProModelOfTskSch(PMS, NumOfAncestors, NumOfNonDescendants, Rank_b); //PMS[i][k] represents the probability that the k-th scheduled task is task i

     vector<chromosome> Population(Parameter_HEDA.NumOfChromPerPop);
     chromosome Chrom_gb = GnrChr_HEFT_b(Rank_b);

 //     chromosome TemChrom1 = GnrChr_Lvl_EFT();
 //     chromosome TemChrom2 = GnrChr_TS_EFT();
 //     chromosome TemChrom3 = GnrChr_TS_Rnd();
 //     if (Chrom_gb.FitnessValue>=TemChrom1.FitnessValue) {
 //          Chrom_gb=TemChrom1 ;
 //     }
 //     if (Chrom_gb.FitnessValue>=TemChrom2.FitnessValue) {
 //         Chrom_gb=TemChrom2 ;
 //     }
 //     if (Chrom_gb.FitnessValue>=TemChrom3.FitnessValue) {
 //         Chrom_gb=TemChrom3 ;
 //     }
  //////////////////////////////


     // ///////////////////////
     // if (Chrom_gb.FitnessValue <= EndFitness + 1e-2) {
     //     SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
     //     ClearALL();
     //     /////////////////
     //     cout<<endl<<"HEDA"<<endl;
     //     for (int i=0;i<comConst.NumOfTsk;++i) {
     //         cout<<Chrom_gb.TskSchLst[i]+1<<" ";
     //     }
     //     cout<<endl;
     //     for (int i=0;i<comConst.NumOfTsk;++i) {
     //         cout<<Chrom_gb.EndTime[i]<<" ";
     //     }
     //     cout<<endl;
     //     for (int i=0;i<comConst.NumOfTsk;++i) {
     //         cout<<Chrom_gb.RscAlcLst[i]+1<<" ";
     //     }
     //     cout<<endl;
     //     cout<<SchTime<<" "<<Chrom_gb.FitnessValue<<endl;
     //     // ////////////////
     //     return Chrom_gb.FitnessValue;
     // }
     // ////////////////////////////////
     vector<double> eta_TSO(comConst.NumOfTsk);
     RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;

     // IntPop1(populations);

     while (RunTime + PrecisionValue < Parameter_HEDA.RunTimeRatioOfStg1 * SchTime){

       //////////////////////
         #pragma omp parallel for
         for (int i = 0; i < comConst.NumOfTsk; ++i) {
             eta_TSO[i] = pow(Rank_b[i]/MaxRank_b, (1-RunTime/SchTime) * Parameter_HEDA.eta);
         }
         #pragma omp parallel for
         for ( int n = 0; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
             Population[n] = GnrTskLstOfChr_prp(PMS, eta_TSO);
             GnrMS_Evl(Population[n]);
         }

         for (int n = 0; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
             if (Population[n].FitnessValue + PrecisionValue < Chrom_gb.FitnessValue)
                 Chrom_gb = Population[n];
         }
     // ///////////////////////
     //     if (Chrom_gb.FitnessValue <= EndFitness + 1e-2) {
     //         SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
     //         ClearALL();
     //         /////////////////
     //         cout<<endl<<"HEDA"<<endl;
     //         for (int i=0;i<comConst.NumOfTsk;++i) {
     //             cout<<Chrom_gb.TskSchLst[i]+1<<" ";
     //         }
     //         cout<<endl;
     //         for (int i=0;i<comConst.NumOfTsk;++i) {
     //             cout<<Chrom_gb.EndTime[i]<<" ";
     //         }
     //         cout<<endl;
     //         for (int i=0;i<comConst.NumOfTsk;++i) {
     //             cout<<Chrom_gb.RscAlcLst[i]+1<<" ";
     //         }
     //         cout<<endl;
     //         cout<<SchTime<<" "<<Chrom_gb.FitnessValue<<endl;
     //         // ////////////////
     //         return Chrom_gb.FitnessValue;
     //     }
     // ////////////////////////////////
         UpdatePMR(PMR, Chrom_gb); UpdatePMS(PMS, Chrom_gb);

         ++iteration;
         RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
     }
     // while (RunTime + PrecisionValue < SchTime){
     while (1){
         RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
         if ( RunTime >= SchTime ) {
             SchTime = RunTime;
             break;
         }
         #pragma omp parallel for
         for (int i = 0; i < comConst.NumOfTsk; ++i) {
             eta_TSO[i] = pow(Rank_b[i]/MaxRank_b, (1-RunTime/SchTime) * Parameter_HEDA.eta);
         }
         #pragma omp parallel for
         for (int n = 0; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
             Population[n] = GnrTskLstOfChr_prp(PMS, eta_TSO);
             GnrRscLstOfChr(Population[n],PMR);
             DcdEvl(Population[n],true);
         }

         chromosome Chrom_lb = Population[0];
         for (int n = 1; n < Parameter_HEDA.NumOfChromPerPop; ++n) {
             if (Population[n].FitnessValue + PrecisionValue < Chrom_lb.FitnessValue)
                 Chrom_lb = Population[n];
         }
 //

         if (Chrom_gb.FitnessValue<=Chrom_lb.FitnessValue) {
             Chrom_lb=Chrom_gb;
         }
         /////////////////
         cout<<"origin"<<endl;
         for (int i=0;i<comConst.NumOfTsk;++i) {
             cout<<Chrom_gb.TskSchLst[i]+1<<" ";
         }
         cout<<endl;
         for (int i=0;i<comConst.NumOfTsk;++i) {
             cout<<Chrom_gb.EndTime[i]<<" ";
         }
         cout<<endl;
         for (int i=0;i<comConst.NumOfTsk;++i) {
             cout<<Chrom_gb.RscAlcLst[i]+1<<" ";
         }
         cout<<endl;
         cout<<SchTime<<" "<<Chrom_gb.FitnessValue<<endl;
         chromosome TemChrom_lb = Chrom_gb;
         // //////////////////////////
         IFBDI(Chrom_lb);
         if (Chrom_lb.FitnessValue < Chrom_gb.FitnessValue + 1e-2) {
             //SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
             ////////
             cout<<"IFBDI"<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<Chrom_lb.TskSchLst[i]+1<<" ";
             }

             cout<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<Chrom_lb.EndTime[i]<<" ";
             }
             cout<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<Chrom_lb.RscAlcLst[i]+1<<" ";
             }
             cout<<endl;
             cout<<SchTime<<" "<<Chrom_lb.FitnessValue<<endl;

             //////////////////////
             //ClearALL();
             return Chrom_lb.FitnessValue;
         }

         LBCAI(TemChrom_lb);
         if (TemChrom_lb.FitnessValue < Chrom_gb.FitnessValue + 1e-2) {
             //SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
             ////////
             cout<<"LBCAI"<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<TemChrom_lb.TskSchLst[i]+1<<" ";
             }
             cout<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<TemChrom_lb.EndTime[i]<<" ";
             }
             cout<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<TemChrom_lb.RscAlcLst[i]+1<<" ";
             }
             cout<<endl;
             cout<<SchTime<<" "<<TemChrom_lb.FitnessValue<<endl;

             //////////////////////
             //ClearALL();
             return Chrom_lb.FitnessValue;
         }



         if (Chrom_lb.FitnessValue + PrecisionValue < Chrom_gb.FitnessValue) {
             Chrom_gb = Chrom_lb;

         }
         // ///////////////////////
         // if (Chrom_gb.FitnessValue <= EndFitness + 1e-2) {
         //     SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
         //     ClearALL();
         //     /////////////////
         //     cout<<endl<<"HEDA"<<endl;
         //     for (int i=0;i<comConst.NumOfTsk;++i) {
         //         cout<<Chrom_gb.TskSchLst[i]+1<<" ";
         //     }
         //     cout<<endl;
         //     for (int i=0;i<comConst.NumOfTsk;++i) {
         //         cout<<Chrom_gb.EndTime[i]<<" ";
         //     }
         //     cout<<endl;
         //     for (int i=0;i<comConst.NumOfTsk;++i) {
         //         cout<<Chrom_gb.RscAlcLst[i]+1<<" ";
         //     }
         //     cout<<endl;
         //     cout<<SchTime<<" "<<Chrom_gb.FitnessValue<<endl;
         //     // ////////////////
         //     return Chrom_gb.FitnessValue;
         // }
         // ////////////////////////////////
         UpdatePMR(PMR, Chrom_gb); UpdatePMS(PMS, Chrom_gb);

         ++iteration;
         RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
     }
     SchTime = RunTime;

     // //////////////
     cout<<endl<<"HEDA"<<endl;
     for (int i=0;i<comConst.NumOfTsk;++i) {
         cout<<Chrom_gb.TskSchLst[i]+1<<" ";
     }
     cout<<endl;
     for (int i=0;i<comConst.NumOfTsk;++i) {
         cout<<Chrom_gb.EndTime[i]<<" ";
     }
     cout<<endl;
     for (int i=0;i<comConst.NumOfTsk;++i) {
         cout<<Chrom_gb.RscAlcLst[i]+1<<" ";
     }
     cout<<endl;
     cout<<Chrom_gb.FitnessValue<<endl;
     // ////////////////

     return Chrom_gb.FitnessValue;
 }