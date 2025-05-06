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

using namespace std;
using namespace PlotOptTypes;


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
     if (isSMS) w_year = 137.61;//*vardef::xsec_ratio.GetScalar(b); //since we have only 2016 v7 MC
     w_year = w_year*b.w_lumi();//*b.w_trig();
     //if(isdiph) w_year *= b.mc_weight()>0 ? b.mc_weight():-b.mc_weight(); /*w_year = w_year*b.mc_weight();*/ //w_lumi alr    eady knows the sign
     //if(isQCD2) w_year *= -0.1131001131; //due to it being unable to find cross section in nano2pico
     if (isttWToLNu) w_year *= 0.2163/66680; // n2p adelie_v4 had bug where ttWJetsToLNu picked up WJetsToLNu xsec
   }
   return w_year; 

   });


 NamedFunc met = "met";

 NamedFunc hig_am = "hig_df_cand_am[0]";
 NamedFunc hig_dm = "hig_df_cand_dm[0]";
 NamedFunc hig_drmax = "hig_df_cand_drmax[0]";

 const NamedFunc dphi_vec("dphi_vec", [](const Baby &b) -> NamedFunc::VectorType{
   vector<double> dphi;
   for (int i=0; i < min(b.njet(), 4); i++){
     dphi.push_back( fabs(b.jet_phi()->at(b.jet_ordered_pt_indices()->at(i)) - b.met_phi()) );
   }
   return dphi;
 });

 NamedFunc jets_hh_pt("jets_hh_pt", [](const Baby &b) -> NamedFunc::VectorType{ return vardef::jets_hh_pt.GetVector(b); });

 NamedFunc num_b("num_b", [](const Baby &b) -> NamedFunc::ScalarType{ return vardef::num_b.GetScalar(b);});

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

  NamedFunc type("rtype",  [](const Baby &b) -> NamedFunc::ScalarType{ return b.type(); });

  // triggers
  NamedFunc pass_filters("pass_filters", [](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_goodv() || !b.pass_cschalo_tight() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet() || !b.pass_low_neutral_jet() || !b.pass_htratio_dphi_tight() || !b.pass_ecalnoisejet() ){ return false; }
    return true;
  });

  //Declares the samples

  string mcfolder = "/net/cms18/cms18r0/pico/NanoAODv9/higgsino_4b_adelie_v3/2016/mc/unskimmed/";

  auto procs_150 = Process::MakeShared<Baby_pico>("m_nlsp=150,m_lsp=1", Process::Type::signal, kRed,{mcfolder+"pico_SMS-TChiHH_mChi-150_mLSP-0*.root"}, "1");
  auto procs_200 = Process::MakeShared<Baby_pico>("m_nlsp=200,m_lsp=1", Process::Type::signal, kViolet,{mcfolder+"pico_SMS-TChiHH_mChi-200_mLSP-0*.root"}, "1");
  auto procs_300 = Process::MakeShared<Baby_pico>("m_nlsp=300,m_lsp=1", Process::Type::signal, kCyan,{mcfolder+"pico_SMS-TChiHH_mChi-300_mLSP-0*.root"}, "1");
  auto procs_500 = Process::MakeShared<Baby_pico>("m_nlsp=500,m_lsp=1", Process::Type::signal, kGreen,{mcfolder+"pico_SMS-TChiHH_mChi-500_mLSP-0*.root"}, "1");
  auto procs_600 = Process::MakeShared<Baby_pico>("m_nlsp=600,m_lsp=1", Process::Type::signal, kGreen,{mcfolder+"pico_SMS-TChiHH_mChi-600_mLSP-0*.root"}, "1");
  auto procs_800 = Process::MakeShared<Baby_pico>("m_nlsp=800,m_lsp=1", Process::Type::signal, kOrange,{mcfolder+"pico_SMS-TChiHH_mChi-800_mLSP-0*.root"}, "1");
  auto procs_900 = Process::MakeShared<Baby_pico>("m_nlsp=900,m_lsp=1", Process::Type::signal, kBlue,{mcfolder+"pico_SMS-TChiHH_mChi-900_mLSP-0*.root"}, "1");
  auto procs_1000 = Process::MakeShared<Baby_pico>("m_nlsp=1000,m_lsp=1", Process::Type::signal, kRed+3,{mcfolder+"pico_SMS-TChiHH_mChi-1000_mLSP-0*.root"}, "1");
  auto procs_500_300 = Process::MakeShared<Baby_pico>("m_nlsp=500,m_lsp=300", Process::Type::signal, kViolet,{mcfolder+"pico_SMS-TChiHH_mChi-500_mLSP-300_*.root"}, "1");


  vector<shared_ptr<Process>> procs_sig_res = {procs_200, procs_300, procs_500, procs_800, procs_1000};
  vector<shared_ptr<Process>> procs_sig_boo = {procs_200, procs_300, procs_500, procs_800, procs_1000};

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
				     .FileExtensions({"pdf", "root"})};


  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);//.Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_lumi.Stack(StackType::data_norm)};
  vector<PlotOpt> ops_log = {log_lumi().Stack(StackType::data_norm)};
  vector<PlotOpt> ops_shapes = {lin_lumi.Stack(StackType::shapes)};
  vector<PlotOpt> ops_2D = {bkg_hist};
  

//  NamedFunc test_regions = regions_4b::res_4b;

  const NamedFunc res_baseline = regions_4b::res_baseline;

  PlotMaker pm;
//  pm.Push<Hist1D>(Axis(4,0,4,num_b, "N_{b}", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:testplot");
  pm.Push<Hist1D>(Axis(30,0,300,met, "p_{T}^{miss}", {}), res_baseline, procs_sig_res, ops).Weight(w_run2).Tag("ShortName:testplot");

  pm.min_print_ = true;
  pm.MakePlots(1);
  //pm.MakePlots(137.6);
  
}


