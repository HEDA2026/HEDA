
#ifndef CSTCHANGE_TOOLS_HPP
#include "common.h"
#define CSTCHANGE_TOOLS_HPP

void CalculateLevelList();
void CalculateDescendants();
void CalculateAncestors();
bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b);
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
void IndexSortByValueOnAscend(vector<int>& ind, vector<int>& value);
void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& fitness);
void IndexSortByValueOnDescend(vector<int>& ind, vector<int>& value);
void IndexSortByValueOnDescend(vector<int>& ind, vector<double>& fitness);
double RandomDouble(int start, int end);
double RandomDouble2(int start, int end);
int CalculateParNum(int taskIndex, vector<int>& markList);
int CalculateSonNum(int taskIndex, vector<int>& markList);

void IndexSort(vector<int>& ind, vector<int>& value);
void IndexSort(vector<int>& ind, vector<double>& fitness);
void FndLvl(chromosome& a ,vector<vector<int>>& islvl);
// bool CompareR2GA(const chromosome& a, const chromosome& b) ;
bool CompareR2GA(const chromosome_R2GA& a, const chromosome_R2GA& b);
#endif //CSTCHANGE_TOOLS_HPP
