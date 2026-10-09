//
// R2GA translated from Source-Code-main/R2GA/SCHEDULER/task/geneticscheduler.py
//
#include "config.h"
#include <algorithm>
#include <ctime>
#include <vector>
#include "GenerateAChrom.h"
#include "tools.hpp"

namespace {
const int R2GA_POP_SIZE = 500;

inline double Rand01() {
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}
}

double runR2GA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();

    const int n = comConst.NumOfTsk;
    iteration = 0;

    vector<double> ww(n, 0.0);
    vector<vector<double>> cc(n, vector<double>(n, 0.0));
    vector<double> Rank_b(n, 0.0);
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b, cc, ww);

    // Python init_population(): the first chromosome is the HEFT/uprank-based one.
    vector<chromosome_R2GA> population(R2GA_POP_SIZE);
    population[0] = GnrChr_HEFT_To_R2GA(Rank_b);

    // Python init_population(): remaining chromosomes are random real-coded genes.
    for (int i = 1; i < R2GA_POP_SIZE; ++i) {
        population[i].genes.resize(2 * n);
        for (int j = 0; j < 2 * n; ++j) {
            population[i].genes[j] = Rand01();
        }
        Decode_R2GA(population[i]);
    }
    sort(population.begin(), population.end(), CompareR2GA);

    const int elite_count = R2GA_POP_SIZE / 2;

    while (true) {
        double current_run_time = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
        if (current_run_time >= SchTime) {
            SchTime = current_run_time;
            break;
        }

        ++iteration;

        // Python select(): keep the best half after sorting by makespan.
        sort(population.begin(), population.end(), CompareR2GA);
        vector<chromosome_R2GA> next_population;
        next_population.reserve(R2GA_POP_SIZE);
        for (int i = 0; i < elite_count; ++i) {
            next_population.push_back(population[i]);
        }

        // Python crossover(): pair population[0], population[1], ... from the elite half.
        vector<chromosome_R2GA> offspring;
        offspring.reserve(elite_count);
        for (int i = 0; i < elite_count - 1; i += 2) {
            chromosome_R2GA child1 = population[i];
            chromosome_R2GA child2 = population[i + 1];

            int crossover_point1 = 1 + rand() % (n - 2);
            int crossover_point2 = n + 1 + rand() % (n - 2);
            for (int j = crossover_point1; j < crossover_point2; ++j) {
                swap(child1.genes[j], child2.genes[j]);
            }

            offspring.push_back(child1);
            offspring.push_back(child2);
        }

        // Python mutate(): each offspring mutates one task gene and one resource gene.
        for (size_t i = 0; i < offspring.size(); ++i) {
            int pos1 = rand() % n;
            int pos2 = n + rand() % n;
            offspring[i].genes[pos1] = Rand01();
            offspring[i].genes[pos2] = Rand01();
            Decode_R2GA(offspring[i]);
            next_population.push_back(offspring[i]);
        }

        population.swap(next_population);
        sort(population.begin(), population.end(), CompareR2GA);
    }

    return population[0].FitnessValue;
}
