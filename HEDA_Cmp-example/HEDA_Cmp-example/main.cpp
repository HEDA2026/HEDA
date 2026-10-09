#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "HEFT.h"
#include "HGA.h"
#include "LWSGA.h"
#include "HPSO.h"
#include "HEDA.h"
#include "APRE_EDA.h"
#include "TMGA.h"
#include "R2GA.hpp"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm} -xy4
    map<string, double> SchTime;
    SchTime["Montage25_1.0"] = 3.162;
    SchTime["Montage50_1.0"] = 11.154;
    SchTime["Montage100_1.0"] = 57.030;
//    SchTime["CyberShake30_1.0"] = 4.639;
    SchTime["CyberShake30_1.0"] = 4.639;
    SchTime["CyberShake50_1.0"] = 12.749;
    SchTime["CyberShake100_1.0"] = 77.858;

    SchTime["CyberShake8_0.7"] = 5.00;
    SchTime["CyberShake8_1.0"] = 2.00;

    SchTime["Epigenomics24_1.0"] = 2.474;
    SchTime["Epigenomics47_1.0"] = 7.615;
    SchTime["Epigenomics100_1.0"] = 30.426;

    SchTime["Ligo30_1.0"] = 3.390;
    SchTime["Ligo50_1.0"] = 7.677;
    SchTime["Ligo100_1.0"] = 32.455;

    SchTime["gaussian9_1.0"] = 0.745;
    SchTime["gaussian209_1.0"] = 102.464;
    SchTime["gaussian324_1.0"] = 236.955;
    //SchTime["gaussian434_1.0"] = 120.0;

    SchTime["fft40_1.0"] = 6.074;
    SchTime["fft96_1.0"] = 30.894;
    SchTime["fft224_1.0"] = 140.820;
    //SchTime["fft512_1.0"] = 160.0;
    //SchTime["FlashAttention1716_1.0"] = 6.820;

    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(0);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;  //NumOfTask, RscAvlRatio
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

           double HEFT_SchTime  = 0 ;
           double HEFT_Result = runHEFT(XmlFile, RscAlcFile, HEFT_SchTime);
           ClearALL();
//
             double HGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
           // double HGA_SchTime =0;
             int HGA_Iteration = 0;
             double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
             ClearALL();
////
           double LWSGA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
           // double LWSGA_SchTime = 0;
             int LWSGA_Iteration = 0;
             double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
             ClearALL();
// // // //
             double HPSO_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
           // double HPSO_SchTime=0.01;
             int HPSO_Iteration = 0;
             double HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
             ClearALL();
// // // ///新加的 APRE——DEA
            double APRE_EDA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
          // double APRE_EDA_SchTime  = 0;
            int APRE_EDA_Iteration = 0;
            double APRE_EDA_Result = runAPRE_EDA(XmlFile, RscAlcFile, APRE_EDA_SchTime, APRE_EDA_Iteration);
            ClearALL();
//
//              double TMGA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//              int TMGA_Iteration = 0;
//              double TMGA_Result = runTMGA(XmlFile, RscAlcFile, TMGA_SchTime, TMGA_Iteration);
//              ClearALL();
// /////////
         double R2GA_SchTime  =  SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        double R2GA_SchTime  = 0;
         int R2GA_Iteration = 0;
         double R2GA_Result = runR2GA(XmlFile, RscAlcFile, R2GA_SchTime, R2GA_Iteration);
         ClearALL();

        double HEDA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
         // double HEDA_SchTime =0.001;
         int HEDA_Iteration = 0;
         double HEDA_Result = runHEDA(XmlFile, RscAlcFile, HEDA_SchTime, HEDA_Iteration);
         ClearALL();



        //results are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(5);
        // outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "<< endl
//                 << HEFT_Result << " " << HEFT_SchTime << " "<< endl
//                 << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "<< endl
        //         << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "<< endl
        //         << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "<< endl
        //         << APRE_EDA_Result  << " " << APRE_EDA_SchTime << " " << APRE_EDA_Iteration  << " "<<endl
        //         << TMGA_Result  << " " << TMGA_SchTime  << " " << TMGA_Iteration  << " "<<endl
        //         << HEDA_Result  << " " << HEDA_SchTime << " " << HEDA_Iteration  << " "
        //         << endl<<endl;
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                 << HEFT_Result << " " << HEFT_SchTime << " "
                 << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                 << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
                 << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "
                 << APRE_EDA_Result  << " " << APRE_EDA_SchTime << " " << APRE_EDA_Iteration  << " "
//                 << TMGA_Result  << " " << TMGA_SchTime  << " " << TMGA_Iteration  << " "
                << R2GA_Result  << " " << R2GA_SchTime  << " " << R2GA_Iteration  << " "
                << HEDA_Result  << " " << HEDA_SchTime << " " << HEDA_Iteration  << " "
                << endl;
        //outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " " << HEFT_Result << " " << HGA_Result << " " << LWSGA_Result << " " << HPSO_Result<<" " << APRE_EDA_Result<<" " << HEDA_Result
               //<< endl;
        // outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " " << HEFT_Result
        //         << endl;
        outfile.close();
        DeleteFirstLineInFile("../fileList.txt"); //delete the first line in the file
    } while (1);
    //return 0;
}
