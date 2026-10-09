#include <fstream>
#include <sstream>
#include "common.h"
#include "config.h"
#include "HEDA.h"

using namespace std;

int main() {
    srand((int) time(0));

    //{set the runtime (termination) time of the algorithm} -xy4
    map<string, double> SchTime;

    // SchTime["Montage25_1.0"]=3.267;
    // SchTime["Montage50_1.0"]=11.922;
    // SchTime["Montage100_1.0"]=62.800;
    //
    // SchTime["CyberShake30_1.0"]=4.868;
    // SchTime["CyberShake50_1.0"]=13.698;
    // SchTime["CyberShake100_1.0"]=84.084;
    //
    // SchTime["Epigenomics24_1.0"]=3.048;
    // SchTime["Epigenomics47_1.0"]=8.433;
    // SchTime["Epigenomics100_1.0"]=34.750;
    //
    // SchTime["Ligo30_1.0"]=3.977;
    // SchTime["Ligo50_1.0"]=9.264;
    // SchTime["Ligo100_1.0"]=32.074;
    //
    // SchTime["gaussian9_1.0"]=0.891;
    // SchTime["gaussian209_1.0"]=116.487;
    // SchTime["gaussian324_1.0"]=276.288;
    // //SchTime["gaussian434_1.0"]=;
    //
    // SchTime["fft40_1.0"]=7.342;
    // SchTime["fft96_1.0"]=31.663;
    // SchTime["fft224_1.0"]=158.779;
    //SchTime["fft512_1.0"] = 160.0;

    SchTime["Montage25_1.0"] = 5.274;
    SchTime["Montage50_1.0"] = 19.545;
    SchTime["Montage100_1.0"] =101.995;

    SchTime["CyberShake30_1.0"] = 7.997;
    SchTime["CyberShake50_1.0"] = 22.401;
    SchTime["CyberShake100_1.0"] = 138.458;

    SchTime["Epigenomics24_1.0"] = 4.311;
    SchTime["Epigenomics47_1.0"] = 13.189;
    SchTime["Epigenomics100_1.0"] = 54.792;

    SchTime["Ligo30_1.0"] = 6.465;
    SchTime["Ligo50_1.0"] = 14.682;
    SchTime["Ligo100_1.0"] = 32.455;

    SchTime["gaussian9_1.0"] = 1.363;
    SchTime["gaussian209_1.0"] = 180.709;
    SchTime["gaussian324_1.0"] = 419.334;
    //SchTime["gaussian434_1.0"] = 120.0;

    SchTime["fft40_1.0"] = 11.608;
    SchTime["fft96_1.0"] = 51.266;
    SchTime["fft224_1.0"] = 255.423;



    //

    string strLine;
    ifstream iFile("../ExpParaSet.txt");
    if (!iFile) {
        cout << "filelist open failed!\n";
        exit(1);
    }
    while(getline(iFile,strLine)){
        istringstream is(strLine);
        Orthogonal TemOrthogonal;
        is >> TemOrthogonal.PopSizeFactor >> TemOrthogonal.theta1 >> TemOrthogonal.theta2 >> TemOrthogonal.eta
           >> TemOrthogonal.prp >> TemOrthogonal.RunTimeRatioOfStg1;
        orthogonal.push_back(TemOrthogonal);
    }
    iFile.close();

    string Model, NumOfTask, RscAvlRatio;
    do {
        string strLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << "filelist open failed!\n";
            exit(1);
        }
        getline(iFile, strLine);
        if (strLine.size() < 1) {
            cout << "Empty input file(fileList)" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile;
        string RscAlcFile;
        istringstream is(strLine);
        is >> Model >> NumOfTask >> RscAvlRatio;
        XmlFile = Model + "_" + NumOfTask + "_0.xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        int index =0;
        for(Orthogonal TemOrthogonal: orthogonal){
            ++index;
            ofstream outfile("../OrthResultOutput/result"+to_string(index)+".txt", ios::app);
            if (!outfile) {
                cout << "Open the result file failure...\n";
                exit(0);
            }
            outfile.setf(ios::fixed, ios::floatfield);
            outfile.precision(3);
            cout <<endl<< "Parameter" + to_string(index) << " " << Model << " " << NumOfTask << " " << RscAvlRatio << " ";
            for (int times = 0; times < 10; ++times) {
                double HEDA_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
                int HEDA_Iteration = 0;
                double HEDA_Result = runHEDA(XmlFile, RscAlcFile, TemOrthogonal, HEDA_SchTime, HEDA_Iteration);
            //
                // cout<<time<<endl;
                //
                ClearALL();
                outfile << Model << " " << NumOfTask << " " << RscAvlRatio << " "
                        << HEDA_Result << " " << HEDA_SchTime << " " << HEDA_Iteration
                        << endl;
            }
            outfile.close();
        }
        DeleteFirstLineInFile("../fileList.txt");
    } while (1);
}

