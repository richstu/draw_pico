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

// function to get cuts for N-k plots
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

  const NamedFunc sig_decay_4b("sig_decay_4b", [](const Baby &b) -> NamedFunc::ScalarType{
    std::string file = *b.FileNames().begin();
    if (file.find("TChiHH") != file.npos){
      int n_trueb = 0;
      for (size_t idx = 0; idx < b.mc_id()->size(); ++idx){
        if ((b.mc_id()->at(idx) == 5 || b.mc_id()->at(idx) == -5) && b.mc_mom()->at(idx) == 25){
          n_trueb += 1;
        }
      }
      return (n_trueb == 4);
    }
    else {return true;}
  });


  const NamedFunc dphi_res("dphi_res", [](const Baby &b) -> NamedFunc::ScalarType{
    double cut;
    for (int i = 0; i < min(b.njet(),4); i++){
      cut = (i>1) ? 0.5 : 0.3;
      if (vardef::dphi_vec.GetVector(b)[i] < cut) {return false;}
    } 
    return true;
  });


  vector<pair<string, NamedFunc>> cuts_4b_res = {
    {"nvl", 	  "nvlep==0"},
    {"ntk", 	  "ntk==0"},
    {"met", 	  "met>150"},
    {"njet", 	  "(njet==4 || njet==5)"},
    {"nb", 	  "( (nbdft==2 && nbdfm==2) || (nbdft>=2 && nbdfm==3 && nbdfl==3) || (nbdft>=2 && nbdfm>=3 && nbdfl>=4) )"},
    {"fakemet",   "(met/mht)<2 && (met/met_calo)<2"},
    {"hig_am", 	  "hig_df_cand_am[0]<200"},
    {"hig_dm", 	  "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
  };

  vector<pair<string, NamedFunc>> cuts_4b_boo = {
    {"met",    "met>300"},
    {"ht",     "ht>600"},
    {"nak8",   "nfjet>1"},
    {"ak8_pt", "fjet_pt[0]>300 && fjet_pt[1]>300"},
    {"mj_sd",  "fjet_msoftdrop[0]>60 && fjet_msoftdrop[0]<260 && fjet_msoftdrop[1]>60 && fjet_msoftdrop[1]<260"},
//    {"mj_pnet", "fjet_pnet_m[0]>60 && fjet_pnet_m[0]<260 && fjet_pnet_m[1]>60 && fjet_pnet_m[1]<260"}
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR1 = {
    {"nlep", 	  "nlep==1"},
    {"lep_pt", 	  "lep_pt[0]>50"},
    {"mt", 	  "mt<100"},
    {"njet", 	  "(njet==4 || njet==5)"},
    {"nb", 	  "( (nbdft==2 && nbdfm==2) || (nbdft>=2 && nbdfm==3 && nbdfl==3) || (nbdft>=2 && nbdfm>=3 && nbdfl>=4) )"},
    {"hig_am", 	  "hig_df_cand_am[0]<200"},
    {"hig_dm", 	  "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR2 = {
    {"nlep", 	  "nlep==2 && nll>0"},
    {"met",	  "met<50"},
    {"lep_pt",	  "lep_pt[0]>50"},
    {"mll",	  "ll_m[0]>80 && ll_m[0]<100"},
    {"njet",	  "(njet==4 || njet==5)"},
    {"hig_am",	  "hig_df_cand_am[0]<200"},
    {"hig_dm",	  "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
    {"nb", 	  "(nbdfm==1 || nbdfm==0)"},
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR3 = {
    {"nvl", 	  "nvlep==0"},
    {"ntk", 	  "ntk==0"},
    {"met", 	  "met>250"},
    {"njet", 	  "(njet==4 || njet==5)"},
    {"hig_am", 	  "hig_df_cand_am[0]<200"},
    {"hig_dm", 	  "hig_df_cand_dm[0]<40"},
    {"hig_drmax", "hig_df_cand_drmax[0]<2.2"},
    {"nb", 	  "(nbdfm==1 || nbdfm==0)"},
  };

  vector<pair<string, NamedFunc>> cuts_4b_CR1b = {
    {"met",    "met>300"},
    {"ht",     "ht>600"},
    {"nak8",   "nfjet>1"},
    {"ak8_pt", "fjet_pt[0]>300 && fjet_pt[1]>300"},
    {"mj_sd",  "fjet_msoftdrop[0]>60 && fjet_msoftdrop[0]<260 && fjet_msoftdrop[1]>60 && fjet_msoftdrop[1]<260"},
    {"nlep",   "nlep==1"},
    {"lep_pt", "lep_pt[0]>50"},
    {"mt",     "mt<100"},
//    {"mj_pnet", "fjet_pnet_m[0]>60 && fjet_pnet_m[0]<260 && fjet_pnet_m[1]>60 && fjet_pnet_m[1]<260"}
  };


  const NamedFunc res_baseline = sig_decay_4b && get_cuts(cuts_4b_res, "") && dphi_res;
  const NamedFunc res_nm_met = sig_decay_4b && get_cuts(cuts_4b_res, "met") && dphi_res;
  const NamedFunc res_nm_njet = sig_decay_4b && "njet>=4" && get_cuts(cuts_4b_res, "njet") && dphi_res;
  const NamedFunc res_nm_nb = sig_decay_4b && get_cuts(cuts_4b_res, "nb") && dphi_res;
  const NamedFUnc res_nm_fakemet = sig_decay_4b && get_cuts(cuts_4b_res, "fakemet");
  const NamedFunc res_nm_higam = sig_decay_4b && get_cuts(cuts_4b_res, "hig_am") && dphi_res;
  const NamedFunc res_nm_higdm = sig_decay_4b && get_cuts(cuts_4b_res, "hig_dm") && dphi_res;
  const NamedFunc res_nm_higdrmax = sig_decay_4b && get_cuts(cuts_4b_res, "hig_drmax") && dphi_res;
  const NamedFunc res_nm_met_higam = sig_decay_4b && get_cuts_Nmk(cuts_4b_res, {"hig_am", "met"}) && dphi_res;
  const NamedFunc res_test = get_cuts(cuts_4b_res, "") && dphi_res;

  const NamedFunc boo_skimcuts = sig_decay_4b && !res_baseline && "nvlep==0&&ntk==0" && "met/mht<2" && "met/met_calo<2" && "met>150" && "njet>=2" && dphi_res;
  const NamedFunc boo_baseline = boo_skimcuts && get_cuts(cuts_4b_boo,"");
  const NamedFunc boo_nm_met = boo_skimcuts && get_cuts(cuts_4b_boo, "met");
  const NamedFunc boo_nm_ht = boo_skimcuts && get_cuts(cuts_4b_boo, "ht");
  const NamedFunc boo_nm_ak8pt = boo_skimcuts && get_cuts(cuts_4b_boo, "ak8_pt");
  const NamedFunc boo_nm_msd = boo_skimcuts && get_cuts(cuts_4b_boo, "mj_sd");

  const NamedFunc CR1_el = get_cuts(cuts_4b_CR1, "") && "nel>0";
  const NamedFunc CR1_mu = get_cuts(cuts_4b_CR1, "") && "nmu>0";
  const NamedFunc CR2_el = get_cuts(cuts_4b_CR2, "") && "nel>0";
  const NamedFunc CR2_mu = get_cuts(cuts_4b_CR2, "") && "nmu>0";
  const NamedFunc CR3 = get_cuts(cuts_4b_CR3, "") && !dphi_res;
  const NamedFunc CR1b_el = get_cuts(cuts_4b_CR1b, "") && "nel>0";
  const NamedFunc CR1b_mu = get_cuts(cuts_4b_CR1b, "") && "nmu>0";

}

