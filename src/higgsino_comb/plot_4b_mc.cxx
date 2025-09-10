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
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Reader.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/sample_loader.hpp"
#include "higgsino_comb/vardef.hpp"
#include "higgsino_comb/regions.hpp"
#include "higgsino_comb/weights.hpp"

using namespace std;
using namespace PlotOptTypes;
/*

 NamedFunc w_run2("w_run2",[](const Baby &b) -> NamedFunc::ScalarType{
   double w_year = 1;
   if(b.SampleTypeString().Contains("-")) w_year = 1;
   else{
     //bool isdiph = false;
     //bool isQCD2 = false;
     bool isttWToLNu = false;
     bool isSMS = false;
     std::string file = *b.FileNames().begin();
     //if (file.find("DiPhotonJetsBox") != file.npos) isdiph = true;
     //if (file.find("QCD_Pt-40") != file.npos) isQCD2 = true;
     if (file.find("TTWJetsToLNu") != file.npos) isttWToLNu = true;
     if (file.find("TChiHH") != file.npos) isSMS = true;
     if (b.SampleTypeString()=="2016") {
         w_year = 19.5;
     } else if (b.SampleTypeString()=="2016APV"){
         w_year = 16.8;
     } else if (b.SampleTypeString()=="2017"){
         w_year = 41.48;
     } else if (b.SampleTypeString()=="2018"){
         w_year = 59.83;
     }
     if (isSMS) w_year = 137.61*vardef::xsec_ratio.GetScalar(b);
     w_year = w_year*b.w_lumi();// *b.w_trig();
     //if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); / *w_year = w_year*b.mc_weight();* / //w_lumi alr    eady knows the sign
     //if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
     if (isttWToLNu) w_year *= 0.2163/66680; // n2p adelie_v4 had bug where ttWJetsToLNu picked up WJetsToLNu xsec
   }
   return w_year; 

   });
*/

 NamedFunc w_run2("w_run2", [](const Baby &b) -> NamedFunc::ScalarType{ return weights::w_lumi_hb.GetScalar(b); });
// NamedFunc w_lumi("w_lumi", [](const Baby &b) -> NamedFunc::ScalarType{ return weights::w_lumi_hb.GetScalar(b); });

 NamedFunc met = "met";

 NamedFunc hig_am = "hig_df_cand_am[0]";
 NamedFunc hig_dm = "hig_df_cand_dm[0]";
 NamedFunc hig_drmax = "hig_df_cand_drmax[0]";

 NamedFunc jets_hh_pt("jets_hh_pt", [](const Baby &b) -> NamedFunc::VectorType{ return vardef::jets_hh_pt.GetVector(b); });
 NamedFunc jets_hh_eta("jets_hh_eta", [](const Baby &b) -> NamedFunc::VectorType{ return vardef::jets_hh_eta.GetVector(b); });
 NamedFunc jets_hh_phi("jets_hh_phi", [](const Baby &b) -> NamedFunc::VectorType{ return vardef::jets_hh_phi.GetVector(b); });
 NamedFunc deepflav_hh_bscore("deepflav_hh_bscore", [](const Baby &b) -> NamedFunc::VectorType{ return vardef::deepflav_hh_bscore.GetVector(b);});

 NamedFunc num_b("num_b", [](const Baby &b) -> NamedFunc::ScalarType{ return vardef::num_b.GetScalar(b);});
 NamedFunc njet = "njet";

 NamedFunc ht = "ht";
 NamedFunc num_ak8 = "nfjet";
 

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  // triggers
//  NamedFunc pass_filters("pass_filters", [](const Baby &b) -> NamedFunc::ScalarType{
//    if (!b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_goodv() || !b.pass_cschalo_tight() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet() || !b.pass_low_neutral_jet() || !b.pass_htratio_dphi_tight() || !b.pass_ecalnoisejet() ){ return false; }
//    return true;
//  });

  NamedFunc pass_filters("pass_filters", [](const Baby &b) -> NamedFunc::ScalarType{ return filters::pass_filters.GetScalar(b); });

  //Declares the samples

  string mcfolder = "/net/cms18/cms18r0/pico/NanoAODv9/higgsino_4b_adelie_v3/2016/mc/unskimmed/";
  string new_mcfolder = "/net/cms37/data1/mhussain/HH-MET/nano2pico_apr2025/nano2pico/test_sigprod_higgsino/signal/unskimmed/";

  auto procs_500 = Process::MakeShared<Baby_pico>("m_nlsp=500,m_lsp=1", Process::Type::signal, kGreen,{mcfolder+"pico_SMS-TChiHH_mChi-500_mLSP-0*.root"}, "1");
  auto procs_500_new = Process::MakeShared<Baby_pico>("m_nlsp=500,m_lsp=1, FastSim 2017", Process::Type::signal, kMagenta,{new_mcfolder+"pico_SMS-TChiHH_mChi-500_mLSP-1*.root"}, "1");
  auto procs_800_new = Process::MakeShared<Baby_pico>("m_nlsp=800,m_lsp=1, FastSim 2017", Process::Type::signal, kBlue,{new_mcfolder+"pico_SMS-TChiHH_mChi-800_mLSP-1*.root"}, "1");


  vector<shared_ptr<Process>> procs_sig_res = {procs_500};
  vector<shared_ptr<Process>> procs_new = {procs_500_new, procs_800_new};

  PlotOpt log_lumi("txt/plot_styles.txt","Std1D");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::both)
          .YTitleOffset(1.)
          .LogMinimum(0.001)
          //.LogMaximum(10000)
          //.LinMaximum(20000)
          .AutoYAxis(true) //false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .CanvasWidth(900)
          .FileExtensions({"pdf"});
          //.Bottom(BottomType::ratio); //_cut_upper); //Added this to find best cuts. remove this later 

  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .LogMinimum(0.1)
                                     //.LogMaximum(2000)
                                     .Overflow(OverflowType::both)
                                     .CanvasWidth(600)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::linear) // originally YAxisType::log
                                     .Title(TitleType::info)
				     .FileExtensions({"pdf"})};


  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);//.Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_lumi.Stack(StackType::data_norm)};
  vector<PlotOpt> ops_log = {log_lumi().Stack(StackType::data_norm)};
  vector<PlotOpt> ops_shapes = {lin_lumi.Stack(StackType::shapes)};
  vector<PlotOpt> ops_2D = {bkg_hist};
  

//  NamedFunc test_regions = regions_4b::res_4b;
  vector<pair<string, NamedFunc>> cuts_4b_res = regions_4b::cuts_4b_res;
  
  const NamedFunc sig_decay_4b("sig_decay_4b", [](const Baby &b) -> NamedFunc::ScalarType{ return regions_4b::sig_decay_4b.GetScalar(b); });
  const NamedFunc dphi_res("dphi_res", [](const Baby &b) -> NamedFunc::ScalarType{ return regions_4b::dphi_res.GetScalar(b); });
  NamedFunc nb2 = "nbdft==2 && nbdfm==2";
  NamedFunc nb3 = "nbdft>=2 && nbdfm==3 && nbdfl==3";
  NamedFunc nb4 = "nbdft>=2 && nbdfm>=3 && nbdfl>=4";
  NamedFunc nb_resolved = nb2 || nb3 || nb4;

  NamedFunc res_baseline_m_higcuts = regions_4b::res_nm_higcuts;
  const NamedFunc res_baseline = regions_4b::res_baseline;

  PlotMaker pm;
  
  pm.Push<Hist1D>(Axis(40,0,1200,met, "p_{T}^{miss} [GeV]", {}), sig_decay_4b, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_met_truthsel_wrun2");
  pm.Push<Hist1D>(Axis(40,0,500, jets_hh_pt, "Jet p_{T} [GeV]", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_updateddphi_jets_hh_pt_resbaseline_m_higcuts_nofilters_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_eta, "Jet #eta", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_jets_hh_eta_resbaseline_m_higcuts_nofilters_wrun2");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14, jets_hh_phi, "Jet #phi", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_jets_hh_phi_resbaseline_m_higcuts_nofilters_wrun2");
  pm.Push<Hist1D>(Axis(20,0,1, deepflav_hh_bscore, "Jet DeepFlavor score", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_jets_hh_deepflavscore_resbaseline_m_higcuts_nofilters_wrun2");
  pm.Push<Hist1D>(Axis(40,0,400, hig_am,  "<m_{bb}> (DeepFlav) [GeV]", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_hig_am_resbaseline_m_higcuts_nofilters_wrun2");
  pm.Push<Hist1D>(Axis(25,0,200, hig_dm,  "#Deltam_{HH} (DeepFlav) [GeV]", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_hig_dm_resbaseline_m_higcuts_nofilters_wrun2");
  pm.Push<Hist1D>(Axis(20,0,4, hig_drmax,  "#DeltaR_{max} (DeepFlav)", {}), res_baseline_m_higcuts, procs_new, ops).Weight(w_run2).Tag("ShortName:procsnew_dpnew_hig_drmax_resbaseline_m_higcuts_nofilters_wrun2");



/*  pm.Push<Hist1D>(Axis(50,0,500,jets_hh_pt, "Jet p_{T} [GeV]", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:jets_hh_pt_baseline");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14,jets_hh_eta, "Jet #eta", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:jets_hh_eta_baseline");
  pm.Push<Hist1D>(Axis(20,-3.14,3.14,jets_hh_phi, "Jet #phi", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:jets_hh_phi_baseline");
  pm.Push<Hist1D>(Axis(20,0,1,deepflav_hh_bscore, "Jet DeepFlavor score", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:deepflav_hh_bscore_baseline");
  pm.Push<Hist1D>(Axis(20,0,400,hig_am, "<m_{bb}> [GeV]", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:higam_baseline");
  pm.Push<Hist1D>(Axis(20,0,200,hig_dm, "#Deltam_{HH} [GeV]", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:higdm_baseline");
  pm.Push<Hist1D>(Axis(20,0,4,hig_drmax, "#DeltaR_{max}", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:higdrmax_baseline");
*/

  pm.min_print_ = true;
  pm.MakePlots(1);
  //pm.MakePlots(137.6);
  
}


