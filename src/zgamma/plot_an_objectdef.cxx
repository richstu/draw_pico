#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/categorization_utilities.hpp"
#include "zgamma/controlregion_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;
using namespace CRUtilities;


//--This is a skeleton file meant to be used to easily make plotting files
NamedFunc photon_pt_tr("photon_pt_tr",[](const Baby &b) -> NamedFunc::ScalarType{ TLorentzVector ph_tr = AssignGamma(b,true); return ph_tr.Pt(); });
NamedFunc photon_eta_tr("photon_eta_tr",[](const Baby &b) -> NamedFunc::ScalarType{ TLorentzVector ph_tr = AssignGamma(b,true); return ph_tr.Eta(); });
//NamedFunc Nphtrue_g1("photon_eta_tr",[](const Baby &b) -> NamedFunc::ScalarType{ return AssignGamma(b,true).Pt() > 0.5;  });

NamedFunc l1_pt_tr("l1_pt_tr",[](const Baby &b) -> NamedFunc::ScalarType{   
  double pt1 = AssignL1(b,true).Pt();
  double pt2 = AssignL2(b,true).Pt();
  return pt1 > pt2 ? pt1 : pt2; 
});

NamedFunc l2_pt_tr("l2_pt_tr",[](const Baby &b) -> NamedFunc::ScalarType{   
  double pt1 = AssignL1(b,true).Pt();
  double pt2 = AssignL2(b,true).Pt();
  return pt1 < pt2 ? pt1 : pt2; 
});

NamedFunc lep_pt_tr("lep_pt_tr",[](const Baby &b) -> NamedFunc::VectorType{   
  double pt1 = AssignL1(b,true).Pt();
  double pt2 = AssignL2(b,true).Pt();
  return {pt1,pt2}; 
});

bool lep_tr_pass(const Baby &b, double pt_req, double eta_req){
  TLorentzVector l1 = AssignL1(b,true);
  TLorentzVector l2 = AssignL2(b,true);
  return l1.Pt() > pt_req && l2.Pt() > pt_req && fabs(l1.Eta()) < eta_req && fabs(l2.Eta()) < eta_req; 
}

bool ph_tr_pass(const Baby &b, double pt_req = 15, double eta_req = 2.5){
  TLorentzVector ph = AssignGamma(b,true);
  return ph.Pt() > pt_req && fabs(ph.Eta()) < eta_req && (fabs(ph.Eta()) < 1.4442 || fabs(ph.Eta()) > 1.566) ; 
}

NamedFunc mu_tr_pass("mu_tr_pass",[](const Baby &b) -> NamedFunc::ScalarType{ return lep_tr_pass(b,5,2.4);});
NamedFunc el_tr_pass("el_tr_pass",[](const Baby &b) -> NamedFunc::ScalarType{ return lep_tr_pass(b,7,2.5);});
NamedFunc photon_tr_pass("photon_tr_pass",[](const Baby &b) -> NamedFunc::ScalarType{ return ph_tr_pass(b);});


NamedFunc l1_eta_tr("l1_eta_tr",[](const Baby &b) -> NamedFunc::ScalarType{   
  TLorentzVector l1 = AssignL1(b,true);
  TLorentzVector l2 = AssignL2(b,true);
  return l1.Pt() > l2.Pt() ? l1.Eta() : l2.Eta(); 
});

NamedFunc l2_eta_tr("l2_eta_tr",[](const Baby &b) -> NamedFunc::ScalarType{   
  TLorentzVector l1 = AssignL1(b,true);
  TLorentzVector l2 = AssignL2(b,true);
  return l1.Pt() < l2.Pt() ? l1.Eta() : l2.Eta(); 
});

NamedFunc lep_eta_tr("lep_eta_tr",[](const Baby &b) -> NamedFunc::VectorType{   
  double l1eta = AssignL1(b,true).Eta();
  double l2eta = AssignL2(b,true).Eta();
  return {l1eta,l2eta}; 
});
/*
int Zchild(const Baby &b){
  for(size_t i = 0; i < b.mc_id()->size(); i++){
    int mcid = abs(b.mc_id()->at(i));
    int mcmom= b.mc_mom()->at(i);
    if( !(mcid == 11 || mcid == 13) ){ continue;}
    if(mcmom == 23) { return mcid; }
  }
  return -1;
}
*/
//Get the truth photon index
int ph_trmatch(const Baby &b){
  //Get truth photon index
  int idx_ph = 0;
  for(size_t i = 0; i < b.mc_id()->size(); i++){
      int mcid = abs(b.mc_id()->at(i));
      int mcmom= b.mc_mom()->at(i);
      if( !(mcid == 22) ){ continue;}
      if(mcmom == 25) { idx_ph = i; break; }
  }

  double min_pt = 9999; double min_dr = 9999;
  int idx_trmatch = -1;

  //Loop through reco photons and check reco photons against the truth photon
  for(size_t idx = 0; idx < b.photon_pt()->size(); idx++){
    double diff_rel_pt = fabs(b.mc_pt()->at(idx_ph) - b.photon_pt()->at(idx))/(b.mc_pt()->at(idx_ph));
    double dr = dR( b.photon_eta()->at(idx), b.mc_eta()->at(idx_ph), b.photon_phi()->at(idx), b.mc_phi()->at(idx_ph) );
    
    //Check that this is a closer photon
    if( diff_rel_pt > min_pt || dr > min_dr){continue;}
    min_pt = diff_rel_pt;
    min_dr = dr;
    idx_trmatch = idx;
  }

  //Require at least within 5 GeV and 0.2 in dR
  if(min_pt > 0.2 || min_dr > 0.2){ return -1;}
  return idx_trmatch;
}

//NamedFunc Zee("Zee",[](const Baby &b) -> NamedFunc::ScalarType{ return Zchild(b)==11; });
//NamedFunc Zmumu("Zmumu",[](const Baby &b) -> NamedFunc::ScalarType{ return Zchild(b)==13; });

//Functions to get idmva for a truth matched photon
NamedFunc has_trmatch("has_trmatch",[](const Baby &b) -> NamedFunc::ScalarType{ return ph_trmatch(b) != -1; });
NamedFunc ph_trmatch_idmva("ph_trmatch_idmva",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_idmva() -> at(ph_trmatch(b)); });
NamedFunc ph_trmatch_wp80("ph_trmatch_idmva",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_id80() -> at(ph_trmatch(b)); });


//NamedFunc that returns a vector of the filter idx 
NamedFunc filter_vector("filter_vector",[](const Baby &b) -> NamedFunc::VectorType{
  //List of filters to use
  vector<bool> filters = {b.pass_goodv(), b.pass_cschalo_tight(), b.pass_hbhe(), b.pass_hbheiso(), b.pass_ecaldeadcell(), 
                          b.pass_badpfmu(), b.pass_badpfmudz(), b.pass_hfnoisyhits(), b.pass_eebadsc(), b.pass_badcalib(), b.pass()};
  

  //Add to the vector fail_filter_idx if event does not pass the filter
  vector<double> fail_filter_idx = {};
  for(unsigned int idx = 0; idx < filters.size(); idx++){ 
    if(!filters[idx]){ fail_filter_idx.push_back(static_cast<double>(idx) + 0.5 ); } //Just want to make sure it ends up in the right bin
  }
  fail_filter_idx.push_back(filters.size()); 
  
  //Return vector to add to plot
  return fail_filter_idx; 
});




int main() {
  gErrorIgnoreLevel = 8000;
  Palette colors("txt/colors.txt","default");

  vector<PlotOpt> ops = CRUtilities::an_plotting_options("shapes"); ops[0].Bottom(PlotOptTypes::BottomType::off).Title(PlotOptTypes::TitleType::simulation);
  vector<PlotOpt> ops_shape = CRUtilities::an_plotting_options("shapes"); ops_shape[0].Title(PlotOptTypes::TitleType::simulation);
  vector<PlotOpt> ops_2D = an_plotting_options("2D"); ops_2D[0].YAxis(YAxisType::log).LogMinimum(0.01).LogMaximum(10);
  vector<PlotOpt> ops_1D = CRUtilities::an_plotting_options("1D");
  vector<PlotOpt> ops_root = CRUtilities::an_plotting_options("1D root");
  
  string lumi_r23 = CRUtilities::R23Lumi;
  string lumi_r3  = CRUtilities::R3Lumi;
  string lumi_r2  = CRUtilities::R2Lumi;

  //Run2
  string samples_file = CRUtilities::SamplesFile;
  //vector<shared_ptr<Process>> procs_r2vsr3 = ZgUtilities::procs_with_sig_scale(samples_file, "SignalRun2vsRun3", 1); 
  //vector<shared_ptr<Process>> procs_svb    = ZgUtilities::procs_with_sig_scale(samples_file, "SignalBackgroundAllYears", 1); 
  vector<shared_ptr<Process>> procs      = ZgUtilities::procs_with_sig_scale(samples_file, "SignalAllYearsSplitUnskimmed", 1);  
  vector<shared_ptr<Process>> procs_data = ZgUtilities::procs_with_sig_scale(samples_file, "DataAllYears", 1); 
  vector<shared_ptr<Process>> procs_r2 = ZgUtilities::procs_with_sig_scale(samples_file, "SignalRun2Split", 1); 
  vector<shared_ptr<Process>> procs_r3 = ZgUtilities::procs_with_sig_scale(samples_file, "SignalRun3Split", 1); 
  //vector<shared_ptr<Process>> procs_usbkg  = ZgUtilities::procs_with_sig_scale(samples_file, "SignalAllYearsBkg", 1); 
 
  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc baseline = (Zee || Zmumu) && "weight/w_lumi < 10.0";
  
  //This block of code handles all of the plot making
  PlotMaker pm;

  NamedFunc selection = baseline;
 
  pm.Push<Hist1D>(Axis(80,    0,   80, photon_pt_tr,     "p_{T}(#gamma) [GeV]",       {15}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_photon_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, photon_eta_tr,    "#eta(#gamma)",        {-2.4,2.4}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_photon_eta").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(50, -1.0,  1.0, ph_trmatch_idmva, "#gamma IDMVA",         {0.14,0.42}), selection && has_trmatch, procs_r2, ops).Weight(wgt).Tag("ShortName:an_objectdef_photon_reco_idmva_r2").LuminosityTag(lumi_r2);
  pm.Push<Hist1D>(Axis(50, -1.0,  1.0, ph_trmatch_idmva, "#gamma IDMVA",         {0.203,0.42}), selection && has_trmatch, procs_r3, ops).Weight(wgt).Tag("ShortName:an_objectdef_photon_reco_idmva_r3").LuminosityTag(lumi_r3);
  
  //electrons
  selection = baseline && Zee;
  pm.Push<Hist1D>(Axis(80,    0,  120, l1_pt_tr,  "p_{T}(e_{1}) [GeV]", {7}),        selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_e1_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, l1_eta_tr, "#eta(e_{1})",        {-2.5,2.5}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_e1_eta").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80,    0,  120, l2_pt_tr,  "p_{T}(e_{2}) [GeV]", {7}),        selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_e2_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, l2_eta_tr, "#eta(e_{2})",        {-2.5,2.5}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_e2_eta").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80,    0,  120, lep_pt_tr,  "p_{T}(e) [GeV]",    {7}),        selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_el_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, lep_eta_tr, "#eta(e)",           {-2.5,2.5}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_el_eta").LuminosityTag(lumi_r23);

  //muons
  selection = baseline && Zmumu; 
  pm.Push<Hist1D>(Axis(80,    0,  120, l1_pt_tr,   "p_{T}(#mu_{1}) [GeV]", {5}),        selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_mu1_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, l1_eta_tr,  "#eta(#mu_{1})",        {-2.4,2.4}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_mu1_eta").LuminosityTag(lumi_r23); 
  pm.Push<Hist1D>(Axis(80,    0,  120, l2_pt_tr,   "p_{T}(#mu_{2}) [GeV]", {5}),        selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_mu2_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, l2_eta_tr,  "#eta(#mu_{2})",        {-2.4,2.4}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_mu2_eta").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80,    0,  120, lep_pt_tr,  "p_{T}(#mu) [GeV]",     {5}),        selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_mu_pt").LuminosityTag(lumi_r23);
  pm.Push<Hist1D>(Axis(80, -4.0,  4.0, lep_eta_tr, "#eta(#mu)",            {-2.4,2.4}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_objectdef_mu_eta").LuminosityTag(lumi_r23);

  //weight

  selection = (Zee || Zmumu);

  pm.Push<Hist1D>(Axis(50, 0, 2.5, "weight/w_lumi", "w_{total}/w_{lumi}", {}), selection, procs_r2, ops).Weight(wgt).Tag("ShortName:an_objectdef_weight_run2").LuminosityTag(lumi_r2);
  pm.Push<Hist1D>(Axis(50, 0, 2.5, "weight/w_lumi", "w_{total}/w_{lumi}", {}), selection, procs_r3, ops).Weight(wgt).Tag("ShortName:an_objectdef_weight_run3").LuminosityTag(lumi_r3);

  pm.Push<Hist1D>(Axis(50, 0, 2.5, "weight/w_lumi", "w_{total}/w_{lumi}", {}), ZgFunctions::tightened_baseline, procs_r2, ops).Weight(wgt).Tag("ShortName:an_objectdef_weight_bs_run2").LuminosityTag(lumi_r2);
  pm.Push<Hist1D>(Axis(50, 0, 2.5, "weight/w_lumi", "w_{total}/w_{lumi}", {}), ZgFunctions::tightened_baseline, procs_r3, ops).Weight(wgt).Tag("ShortName:an_objectdef_weight_bs_run3").LuminosityTag(lumi_r3);

 
  vector<NamedFunc> weight_list  = {"w_lep", "w_photon", "w_phshape", "w_nnlo", "w_bhig", "w_bhig_df", "w_btag", "w_btag_df", "w_pu", "w_prefire", "w_trig"};

  vector<string> weight_labels  = {"w_lep", "w_photon", "w_phshape", "w_nnlo", "w_bhig", "w_bhig_df", "w_btag", "w_btag_df", "w_pu", "w_prefire", "w_trig"};
  
  vector<string> weight_names = {"w_{e,ID} #cdot w_{e,Imini} #cdot w_{#mu,ID} #cdot w_{#mu,Imini}", "w_{#gamma,ID}", "w_{#gamma,shape}", "w_{NNLO}", 
                                 "w_{b tag}", "w_{b tag deepflav}", "w_{b higg}", "w_{b higg deepflav}", "w_{PU}", "w_{prefire}", "w_{trig}"};

  //For all the weights plot them individually, some may be 0 or unused
  for(unsigned int idx = 0; idx < weight_list.size(); idx++){
    NamedFunc single_weight = weight_list[idx];
    string    weight_name   = weight_names[idx]; 
    string    weight_label  = weight_labels[idx]; 

    //Plot all the weights. Pruned the ones that don't matter (at least to my knowledge)!
    pm.Push<Hist1D>(Axis(50, 0, 2.5, single_weight, weight_name, {}), selection, procs_r2, ops).Weight(wgt).Tag("ShortName:an_objectdef_" + weight_label + "_run2").LuminosityTag(lumi_r2);
    pm.Push<Hist1D>(Axis(50, 0, 2.5, single_weight, weight_name, {}), selection, procs_r3, ops).Weight(wgt).Tag("ShortName:an_objectdef_" + weight_label + "_run3").LuminosityTag(lumi_r3);
  }


  //Truth level object requirements & selections
  vector<TableRow> truth_object_reqs = {};
  truth_object_reqs.push_back(TableRow("Event has $\\PZ \\rightarrow \\Pep\\Pem$ or $\\PZ \\rightarrow \\PGmp\\PGmm$", selection, 0, 0, wgt));
  truth_object_reqs.push_back(TableRow("$\\pt(\\ell)$, $|\\eta(\\ell)|$ requirements",   selection && ((mu_tr_pass && Zmumu) || (Zee && el_tr_pass)), 0, 0, wgt));
  truth_object_reqs.push_back(TableRow("$\\pt(\\PGm) > 5$ GeV, $|\\eta(\\PGm)| < 2.4$",  selection && (Zmumu && mu_tr_pass), 1, 0, wgt));
  truth_object_reqs.push_back(TableRow("$\\pt(\\Pe) > 7$ GeV, $|\\eta(\\Pe)| < 2.5$ ",   selection && (Zee && el_tr_pass), 0, 1, wgt));
  truth_object_reqs.push_back(TableRow("$\\pt(\\PGg) > 15$ GeV, $|\\eta(\\PGg)| < 2.5$", selection && ((mu_tr_pass && Zmumu) || (Zee && el_tr_pass)) && photon_tr_pass, 0, 0, wgt));
  truth_object_reqs.push_back(TableRow("Has truth matched $\\PGg$",                      selection && ((mu_tr_pass && Zmumu) || (Zee && el_tr_pass)) && photon_tr_pass && has_trmatch, 0, 0, wgt));
  truth_object_reqs.push_back(TableRow("Truth matched $\\PGg$ passes ID WP80",           selection && ((mu_tr_pass && Zmumu) || (Zee && el_tr_pass)) && photon_tr_pass && has_trmatch && ph_trmatch_wp80, 0, 0, wgt));
  pm.Push<Table>("an_object_definitions_tr_selections",  truth_object_reqs,  procs, 0, 1, 0, 0, 0, 1).LuminosityTag("200").Precision(1);

  //Plot impact of filters
  int nfilters = 12;

  //Baseline selection, but with the pass variable removed
  selection = "nllphoton>0" && pass_trigs_and_pt && "photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] > 80 && ll_m[0] < 100 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && ll_m[0] + llphoton_m[0] > 185";

  pm.Push<Hist1D>(Axis(nfilters, 0.00, nfilters, filter_vector, "Events failing filters", {}), selection, procs_data, ops_root).Weight(wgt).Tag("ShortName:an_objectdef_ev_failing_").LuminosityTag(lumi_r23);

  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.MakePlots(1); 
}


