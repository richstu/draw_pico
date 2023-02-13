#include <vector>

std::vector<float> eff_elhltpt23(float pt, float abseta);
std::vector<float> eff_dielleg12(float pt, float abseta);
std::vector<float> eff_muhltpt17(float pt, float abseta);
std::vector<float> eff_dimuleg8(float pt, float abseta);
std::vector<float> eff_isoel3235(float pt, float abseta);
std::vector<float> eff_isomu2427(float pt, float abseta);

std::vector<float> effsig_elhltpt23(float pt, float abseta);
std::vector<float> effsig_dielleg12(float pt, float abseta);
std::vector<float> effsig_muhltpt17(float pt, float abseta);
std::vector<float> effsig_dimuleg8(float pt, float abseta);
std::vector<float> effsig_isoel3235(float pt, float abseta);
std::vector<float> effsig_isomu2427(float pt, float abseta);

float eff_hig19014_ele23(float pt, float abseta);
float eff_hig19014_ele12(float pt, float abseta);
float eff_hig19014_mu17(float pt, float abseta);
float eff_hig19014_mu8(float pt, float abseta);
