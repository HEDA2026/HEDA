
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runLWSGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    /////add
    double EndFitness = 62;
    //////
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);     //read model information-w
    ConfigParameter_LWSGA();           //set the parameter values
    CalculateLevelList();              //calculate the levels of tasks
    vector<chromosome> population(Parameter_LWSGA.NumOfChromPerPop);
//    #pragma omp parallel for
    for ( int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n ) {
        chromosome chrom;
        IntChr(chrom);
        chrom.TskSchLst = GnrSS_Lvl();
        for ( int i = 0; i < comConst.NumOfTsk; ++i ) {
            chrom.RscAlcLst[i] = Tasks[i].ElgRsc[rand() % Tasks[i].ElgRsc.size()];
        }
        population[n] = chrom;
    }
//    #pragma omp parallel for
    for ( int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n ) {
        DcdEvl(population[n], true);                                             //decoding
        ///////
         if (population[n].FitnessValue <= EndFitness + 1e-2) {
             SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
             ClearALL();
             /////////////////
             cout<<endl<<"LWSGA"<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<population[n].EndTime[i]<<" ";
             }
             cout<<endl;
             for (int i=0;i<comConst.NumOfTsk;++i) {
                 cout<<population[n].RscAlcLst[i]+1<<" ";
             }
             cout<<endl;
             cout<<SchTime<<" "<<population[n].FitnessValue<<endl;
             // ////////////////
             return population[n].FitnessValue;
         }
        ///////////
    }
    sort(population.begin(), population.end(), SortPopOnFitValueByAscend);       //sorting

    while (1) {
        ++iteration;
        //{terminate the algorithm according to running time}
        // double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
        // if ( RunTime >=2 * SchTime ) {
        //     SchTime = RunTime;
        //     break;
        // }

        // if (population[0].FitnessValue <= EndFitness + 1e-2) {
        //
        //     return population[0].FitnessValue;
        // }


        vector<chromosome> NewPopulation(Parameter_LWSGA.NumOfChromPerPop) ;
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; n += 2 ) {
            int parent1 = -1;
            int parent2 = -1;
            SelectionTournament(parent1, parent2,Parameter_LWSGA.NumOfChromPerPop); //select two chromosomes using the tournament method
            double rand = RandomDouble(0, 1);
            chromosome TemChromosome1 = population[parent1];
            chromosome TemChromosome2 = population[parent2];
            if ( rand < Parameter_LWSGA.CrossoverRate ) {
                Crossover_LWSGA(TemChromosome1, TemChromosome2);                    //crossover
            } else {
                Mutation_LWSGA(TemChromosome1);                                     //mutation
                Mutation_LWSGA(TemChromosome2);
            }
            NewPopulation[n] = TemChromosome1;
            NewPopulation[n+1] = TemChromosome2;
        }
        #pragma omp parallel for
        for ( int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n ) {
            DcdEvl(NewPopulation[n], true);
        }

        for (int n = 0; n < Parameter_LWSGA.NumOfChromPerPop; ++n) {
            if (NewPopulation[n].FitnessValue <= EndFitness + 1e-2) {
                SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
                ClearALL();
                /////////////////
                cout<<endl<<"LWSGA"<<endl;
                for (int i=0;i<comConst.NumOfTsk;++i) {
                    cout<<NewPopulation[n].EndTime[i]<<" ";
                }
                cout<<endl;
                for (int i=0;i<comConst.NumOfTsk;++i) {
                    cout<<NewPopulation[n].RscAlcLst[i]+1<<" ";
                }
                cout<<endl;
                cout<<SchTime<<" "<<NewPopulation[n].FitnessValue<<endl;
                // ////////////////
                return NewPopulation[n].FitnessValue;
            }
        }

        //{generate the next population}
        population.insert(population.end(),NewPopulation.begin(),NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        population.resize(Parameter_LWSGA.NumOfChromPerPop);

    }

    ////////////////
    cout<<endl<<"LWSGA"<<endl;
    for (int i=0;i<comConst.NumOfTsk;++i) {
        cout<<population[0].EndTime[i]<<" ";
    }
    cout<<endl;
    for (int i=0;i<comConst.NumOfTsk;++i) {
        cout<<population[0].RscAlcLst[i]+1<<" ";
    }
    cout<<endl;
    cout<<population[0].FitnessValue<<endl;
    // ////////////////

    return population[0].FitnessValue;
}


