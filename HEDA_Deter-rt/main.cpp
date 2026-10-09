#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "HEFT.h"
#include "HGA.h"
#include "LWSGA.h"
#include "CGA.h"
#include "HPSO.h"
#include "APRE_EDA.h"
#include "TMGA.h"
#include "R2GA.h"


using namespace std;

int main() {
    srand((int) time(0));
    //{clear "result"}
    ofstream outfile("../result.txt", ios::out);
    outfile.close();
    string Model, NumOfTask, RscAvlRatio;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
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
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        cout <<endl<< Model << " " << NumOfTask << " " << RscAvlRatio << " ";

         double HGA_SchTime  = 0;
         int HGA_Iteration = 0;
         double HGA_Result = runHGA(XmlFile, RscAlcFile, HGA_SchTime, HGA_Iteration);
         ClearALL();

         double LWSGA_SchTime  = 0;
         int LWSGA_Iteration = 0;
         double LWSGA_Result = runLWSGA(XmlFile, RscAlcFile, LWSGA_SchTime, LWSGA_Iteration);
         ClearALL();
//
//        double CGA_SchTime  = 0;
//        int CGA_Iteration = 0;
//        double CGA_Result = runCGA(XmlFile, RscAlcFile, CGA_SchTime, CGA_Iteration);
 //        ClearALL();

         double HPSO_SchTime  = 0;
         int HPSO_Iteration = 0;
         double HPSO_Result = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
         ClearALL();

         double APRE_EDA_SchTime  = 0;
         int APRE_EDA_Iteration = 0;
         double APRE_EDA_Result = runAPRE_EDA(XmlFile, RscAlcFile, APRE_EDA_SchTime, APRE_EDA_Iteration);
         ClearALL();

        // double TMGA_SchTime  =  0;
        // int TMGA_Iteration = 0;
        // double TMGA_Result = runTMGA(XmlFile, RscAlcFile, TMGA_SchTime, TMGA_Iteration);
        // ClearALL();

        double R2GA_SchTime = 0;
        int R2GA_Iteration = 0;
        double R2GA_Result = runR2GA(XmlFile, RscAlcFile, R2GA_SchTime, R2GA_Iteration);
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
        //         << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "<< endl
        //         << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "<< endl
        //      //   << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "<< endl
        //         << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "<< endl
        //         << APRE_EDA_Result << " " << APRE_EDA_SchTime << " " << APRE_EDA_Iteration << " "<< endl
        //         << endl;
        outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                 << HGA_Result << " " << HGA_SchTime << " " << HGA_Iteration << " "
                 << LWSGA_Result << " " << LWSGA_SchTime << " " << LWSGA_Iteration << " "
              //  << CGA_Result << " " << CGA_SchTime << " " << CGA_Iteration << " "<< endl
                 << HPSO_Result << " " << HPSO_SchTime << " " << HPSO_Iteration << " "
                 << APRE_EDA_Result << " " << APRE_EDA_SchTime << " " << APRE_EDA_Iteration << " "
                // << TMGA_Result  << " " << TMGA_SchTime  << " " << TMGA_Iteration  << " "
                <<  R2GA_Result  << " " << R2GA_SchTime  << " " << R2GA_Iteration  << " "
                << endl;
        // outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " " << HGA_SchTime << " " << LWSGA_SchTime << " " << HPSO_SchTime << " " << APRE_EDA_SchTime
        //         << endl;
        outfile.close();
        //delete the first line in the file
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}
