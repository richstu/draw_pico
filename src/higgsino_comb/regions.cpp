#include <algorithm>
#include <stdlib.h>
#include <regex>
#include <string>
#include "higgsino/hig_utilities.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_opt.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino_comb/vardef.hpp"
#include "higgsino_comb/regions.hpp"

using namespace std;

// function to get a NamedFunc of cuts, given a list of cuts and a cut to remove (for N-1 plots)
const string get_cuts(vector<pair<string, NamedFunc>> cut_list, string removed_cut){
  stringstream cuts;
  cuts << "njet>=2"; // basic cut for both 4b and bbgg
  for (auto const selection : cut_list){
    if (selection.first == removed_cut) { continue; }
    else { cuts << "&&" << selection.second; }
  }
  return cuts.str();
}

const string get_cuts_Nmk(vector<pair<string, NamedFunc>> cut_list, vector<string> removed_cuts){
  stringstream cuts;
  cuts << "njet>=2"; // basic cut for both 4b and bbgg
  for (auto const selection : cut_list){
    bool skip = false;
    for (auto const cut : removed_cuts){
      if (selection.first == cut) { skip = true; }
    }
    if (skip) { continue; }
    else { cuts << "&&" << selection.second; }
  }
  return cuts.str();
}


// specific regions for 4b analysis
namespace regions_4b {

//  const NamedFunc nb_res("nb_res", [](const Baby &b) -> NamedFunc::ScalarType{ return vardef::num_b.GetScalar(b) >= 2; });
//  const NamedFunc nb_inv("nb_inv", [](const Baby &b) -> NamedFunc::ScalarType{ return vardef::num_b.GetScalar(b) < 2; });
/*  const NamedFunc dphi_res("dphi_res", [](const Baby &b) -> NamedFunc::ScalarType{
    double cut;
    for (int i = 0; i < min(b.njet(),4); i++){
      cut = (i>2) ? 0.5 : 0.3;
      if (vardef::dphi_vec.GetVector(b)[i] < cut) {return false;}
    } 
    return true;
  });
*/

  vector<pair<string, NamedFunc>> cuts_4b_res = {
    {"nvl", "nvlep==0"},
    {"ntk", "ntk==0"},
    {"met", "met>150"},
    {"njet", "(njet==4 || njet==5)"},
    {"hig_am", "hig_df_cand_am[0]<200"},
    {"hig_dm", "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
  //  {"nb", nb_res},
  //  {"dphi", dphi_res},
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR1 = {
    {"nlep", "nlep==1"},
    {"lep_pt", "lep_pt[0]>50"},
    {"mt", "mt<100"},
    {"njet", "(njet==4 || njet==5)"},
    {"hig_am", "hig_df_cand_am[0]<200"},
    {"hig_dm", "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
  //  {"nb", nb_res},
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR2 = {
    {"nlep", "nlep==2 && nll>0"},
    {"met", "met<50"},
    {"lep_pt", "lep_pt[0]>50"},
    {"mll", "ll_m[0]>80 && ll_m[0]<100"},
    {"njet", "(njet==4 || njet==5)"},
    {"hig_am", "hig_df_cand_am[0]<200"},
    {"hig_dm", "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
    {"nb", "(nbdfm==1 || nbdfm==0)"},
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR3 = {
    {"nvl", "nvlep==0"},
    {"ntk", "ntk==0"},
    {"met", "met>250"},
    {"njet", "(njet==4 || njet==5)"},
    {"hig_am", "hig_df_cand_am[0]<200"},
    {"hig_dm", "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
    {"nb", "(nbdfm==1 || nbdfm==0)"},
  //  {"dphi", dphi_res},
   
  };

  const string res_baseline = get_cuts(cuts_4b_res, "");
  const string res_nm_higam = get_cuts(cuts_4b_res, "hig_am");
  const string res_nm_met_higam = get_cuts_Nmk(cuts_4b_res, {"hig_am", "hig_dm"});

  const string CR1_el = get_cuts(cuts_4b_CR1, "") + "&& nel>0";
  const string CR1_mu = get_cuts(cuts_4b_CR1, "") + "&& nmu>0";

  const string CR2_el = get_cuts(cuts_4b_CR2, "") + "&& nel>0";
  const string CR2_mu = get_cuts(cuts_4b_CR2, "") + "&& nmu>0";
  
  const string CR3 = get_cuts(cuts_4b_CR3, "");

}

