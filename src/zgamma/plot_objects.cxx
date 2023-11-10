//This script generates some tables and plots for ground-up H->Zgamma studies using unskimmed samples
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to plot
//
//possible string options: 
//

#include "core/test.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TError.h"

#include "core/baby.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"
#include "core/utilities.hpp"
#include "zgamma/apply_zg_trigeffs.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using NamedFuncUtilities::FilterNamedFunc;
using NamedFuncUtilities::MapNamedFunc;
using NamedFuncUtilities::MultiReduceNamedFunc;
using NamedFuncUtilities::ReduceNamedFunc;
using NamedFuncUtilities::reduce_max;
using NamedFuncUtilities::reduce_maxfirst;
using NamedFuncUtilities::reduce_sublead;
using NamedFuncUtilities::reduce_subleadfirst;
using NamedFuncUtilities::reduce_sum;
using PlotOptTypes::BottomType;
using PlotOptTypes::OverflowType;
using PlotOptTypes::StackType;
using PlotOptTypes::TitleType;
using PlotOptTypes::YAxisType;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;
using ZgFunctions::w_years;
using ZgFunctions::mc_mommom;
using ZgFunctions::mc_isFirstCopy;
using ZgFunctions::lead_el_pt;
using ZgFunctions::sublead_el_pt;
using ZgFunctions::lead_mu_pt;
using ZgFunctions::sublead_mu_pt;

double map_abs(double x) {
  if (x>0) return x;
  return -1.0*x;
}

int main(){
  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------
  
  NamedFunc w_sigx50("w_sigx50",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() <= 205000)
      return 50.0;
    return 1.0;
  });

  const NamedFunc w_ewkzgamma("w_ewkzgamma",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.FirstFileName().find("ZGamma2JToGamma2L2J_EWK") != std::string::npos) {
      float xs = 0.1145*1000.0; //fb
      float year_nevts = 230000.0; //2016
      if (abs(b.SampleType())==2016) year_nevts = 500000.0;
      if (abs(b.SampleType())==2017) year_nevts = 498000.0;
      return xs/year_nevts;
    }
    return b.w_lumi();
  });

  const NamedFunc year("year", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString()=="2016APV"||
        b.SampleTypeString()=="-2016APV")
      return 2016;
    else if (b.SampleTypeString()=="2016"||
        b.SampleTypeString()=="-2016")
      return 2016.5;
    else if (b.SampleTypeString()=="2017"||
        b.SampleTypeString()=="-2017")
      return 2017;
    else if (b.SampleTypeString()=="2018"||
        b.SampleTypeString()=="-2018")
      return 2018;
    else if (b.SampleTypeString()=="2022"||
        b.SampleTypeString()=="-2022")
      return 2022;
    else if (b.SampleTypeString()=="2022EE"||
        b.SampleTypeString()=="-2022EE")
      return 2022.5;
    else if (b.SampleTypeString()=="2023"||
        b.SampleTypeString()=="-2023")
      return 2023;
    //else if (b.SampleTypeString()=="2023BPix"||
    //    b.SampleTypeString()=="-2023BPix")
      return 2023.5;
  });

  const NamedFunc mc_mu_pt = FilterNamedFunc("mc_pt","(mc_id==13||mc_id==-13)"&&mc_isFirstCopy);
  const NamedFunc mc_el_pt = FilterNamedFunc("mc_pt","(mc_id==11||mc_id==-11)"&&mc_isFirstCopy);
  const NamedFunc mc_photon_pt = FilterNamedFunc("mc_pt","mc_id==22"&&mc_isFirstCopy);
  const NamedFunc mc_mu_eta = FilterNamedFunc("mc_eta","(mc_id==13||mc_id==-13)"&&mc_isFirstCopy);
  const NamedFunc mc_el_eta = FilterNamedFunc("mc_eta","(mc_id==11||mc_id==-11)"&&mc_isFirstCopy);
  const NamedFunc mc_photon_eta = FilterNamedFunc("mc_eta","mc_id==22"&&mc_isFirstCopy);

  const NamedFunc nmc_mu = ReduceNamedFunc("(mc_id==13||mc_id==-13)"&&mc_isFirstCopy, reduce_sum).Name("nmc_mu");
  const NamedFunc nmc_el = ReduceNamedFunc("(mc_id==11||mc_id==-11)"&&mc_isFirstCopy, reduce_sum).Name("nmc_el");
  const NamedFunc nmc_lep = nmc_mu+nmc_el;
  const NamedFunc nmc_photon = ReduceNamedFunc("mc_id==22"&&mc_isFirstCopy, reduce_sum).Name("nmc_photon");

  const NamedFunc mc_lead_mu_pt = ReduceNamedFunc(mc_mu_pt, reduce_max)
      .Name("mc_lead_mu_pt");
  const NamedFunc mc_lead_el_pt = ReduceNamedFunc(mc_el_pt, reduce_max)
      .Name("mc_lead_el_pt");
  const NamedFunc mc_sublead_mu_pt = ReduceNamedFunc(mc_mu_pt, reduce_sublead)
      .Name("mc_sublead_mu_pt");
  const NamedFunc mc_sublead_el_pt = ReduceNamedFunc(mc_el_pt, reduce_sublead)
      .Name("mc_sublead_el_pt");
  const NamedFunc mc_lead_photon_pt = ReduceNamedFunc(mc_photon_pt, reduce_max)
      .Name("mc_lead_photon_pt");

  const NamedFunc mc_lead_mu_eta = MultiReduceNamedFunc(
      {mc_mu_pt, mc_mu_eta}, reduce_maxfirst).Name("mc_lead_mu_eta");
  const NamedFunc mc_lead_el_eta = MultiReduceNamedFunc(
      {mc_el_pt, mc_el_eta}, reduce_maxfirst).Name("mc_lead_el_eta");
  const NamedFunc mc_sublead_mu_eta = MultiReduceNamedFunc(
      {mc_mu_pt, mc_mu_eta}, reduce_subleadfirst).Name("mc_sublead_mu_eta");
  const NamedFunc mc_sublead_el_eta = MultiReduceNamedFunc(
      {mc_el_pt, mc_el_eta}, reduce_subleadfirst).Name("mc_sublead_el_eta");
  const NamedFunc mc_lead_photon_eta = MultiReduceNamedFunc(
      {mc_photon_pt, mc_photon_eta}, reduce_maxfirst).Name("mc_lead_photon_eta");

  const NamedFunc has_z_to_mumu = ReduceNamedFunc("(mc_id==13||mc_id==-13)&&mc_mom==23"&&mc_mommom==25,
      reduce_max).Name("has_z_to_mumu");
  const NamedFunc has_z_to_ee = ReduceNamedFunc("(mc_id==11||mc_id==-11)&&mc_mom==23"&&mc_mommom==25,
      reduce_max).Name("has_z_to_ee");

  const NamedFunc nmc_mu_acc = ReduceNamedFunc(
      "(mc_id==13||mc_id==-13&&mc_pt>5&&mc_eta>-2.4&&mc_eta<2.4)"&&mc_isFirstCopy, 
      reduce_sum).Name("nmc_mu_acc");
  const NamedFunc nmc_el_acc = ReduceNamedFunc(
      "(mc_id==11||mc_id==-11&&mc_pt>7&&mc_eta>-2.5&&mc_eta<2.5)"&&mc_isFirstCopy, 
      reduce_sum).Name("nmc_el_acc");
  const NamedFunc nmc_lep_acc = nmc_mu_acc+nmc_el_acc;
  const NamedFunc nmc_photon_acc = ReduceNamedFunc(
      "(mc_id==22&&mc_pt>15&&mc_eta>-2.5&&mc_eta<2.5)"&&mc_isFirstCopy, 
      reduce_sum).Name("nmc_photon_acc");

  //measure of quality of 2 muons
  //X = no id, L = loose, M = medium, T = tight
  //0=XX, 1=XL, 2=XM, 3=XT, 4=LL, 5=LM, 6=LT, 7=MM, 8=MT, 9=TT
  const NamedFunc mu_quality("mu_quality",[](const Baby &b) -> NamedFunc::ScalarType{
    int nmut = 0, nmum = 0, nmul = 0;
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_tightid()->at(imu)) nmut++;
      else if (b.mu_mediumid()->at(imu)) nmum++;
      else if (b.mu_id()->at(imu)) nmul++;
    }
    if (nmut>=2) return 9;
    else if (nmum>=1 && nmut==1) return 8;
    else if (nmum>=2) return 7;
    else if (nmul>=1 && nmut==1) return 6;
    else if (nmul>=1 && nmum==1) return 5;
    else if (nmul>=2) return 4;
    else if (nmut==1) return 3;
    else if (nmum==1) return 2;
    else if (nmul==1) return 1;
    return 0;
  });

  //measure of quality of 2 electrons
  //X = no id, L = loose, M = WP90, T = WP80
  //0=XX, 1=XL, 2=XM, 3=XT, 4=LL, 5=LM, 6=LT, 7=MM, 8=MT, 9=TT
  const NamedFunc el_quality("el_quality",[](const Baby &b) -> NamedFunc::ScalarType{
    int nelt = 0, nelm = 0, nell = 0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_id80()->at(iel)) nelt++;
      else if (b.el_id90()->at(iel)) nelm++;
      else if (b.el_idLoose()->at(iel)) nell++;
    }
    if (nelt>=2) return 9;
    else if (nelm>=1 && nelt==1) return 8;
    else if (nelm>=2) return 7;
    else if (nell>=1 && nelt==1) return 6;
    else if (nell>=1 && nelm==1) return 5;
    else if (nell>=2) return 4;
    else if (nelt==1) return 3;
    else if (nelm==1) return 2;
    else if (nell==1) return 1;
    return 0;
  });

  //run2 H->Zgamma photon id
  const NamedFunc photon_idvl("photon_idvl",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> photon_id;
    for (unsigned iph = 0; iph < b.photon_eta()->size(); iph++) {
      float abs_eta = fabs(b.photon_eta()->at(iph));
      float mva = b.photon_idmva()->at(iph);
      if ((abs_eta < 1.4442 && mva > -0.4) || (abs_eta > 1.566 && abs_eta < 2.5 && mva > -0.58))
        photon_id.push_back(1.0);
      else
        photon_id.push_back(0.0);
    }
    return photon_id;
  });

  NamedFunc zg_baseline_nolep = NamedFunc("nlep>=2 && nphoton>=1 && (ll_m[0]>50) && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)").Name("baseline_no_leptonptcut");

  const NamedFunc ptloose_photon_sig = NamedFunc("photon_elveto&&photon_drmin>0.4"&&photon_idvl)
      .Name("ptloose_photon_sig");
  const NamedFunc nptloose_photon = ReduceNamedFunc(ptloose_photon_sig,reduce_sum)
      .Name("nptloose_photon");
  const NamedFunc lead_photon_pt = ReduceNamedFunc(FilterNamedFunc("photon_pt",ptloose_photon_sig),
      reduce_max).Name("lead_photon_pt");

  const NamedFunc isoloose_mu_sig = NamedFunc("mu_sip3d<4&&(mu_id||(mu_pt>200&&mu_highptid))")
      .Name("isoloose_mu_sig");
  const NamedFunc nisoloose_mu = ReduceNamedFunc(isoloose_mu_sig,reduce_sum)
      .Name("nisoloose_mu");
  const NamedFunc isoloose_mu_pt = FilterNamedFunc("mu_pt",isoloose_mu_sig)
      .Name("isoloose_mu_pt");
  const NamedFunc isoloose_mu_reliso = FilterNamedFunc("mu_reliso",isoloose_mu_sig)
      .Name("isoloose_mu_reliso");

  const NamedFunc drloose_photon_sig = NamedFunc("photon_elveto"&&photon_idvl)
      .Name("drloose_photon_sig");
  const NamedFunc ndrloose_photon = ReduceNamedFunc(drloose_photon_sig,reduce_sum)
      .Name("ndrloose_photon");
  const NamedFunc drloose_photon_pt = FilterNamedFunc("photon_pt",drloose_photon_sig)
      .Name("drloose_photon_pt");
  const NamedFunc drloose_photon_drmin = FilterNamedFunc("photon_drmin",drloose_photon_sig)
      .Name("drloose_photon_drmin");

  const NamedFunc lead_mu_reliso = MultiReduceNamedFunc({isoloose_mu_pt,isoloose_mu_reliso},reduce_maxfirst)
      .Name("lead_mu_reliso");
  const NamedFunc sublead_mu_reliso = MultiReduceNamedFunc({isoloose_mu_pt,isoloose_mu_reliso},reduce_subleadfirst)
      .Name("sublead_mu_reliso");
  const NamedFunc lead_photon_drmin = MultiReduceNamedFunc({drloose_photon_pt,drloose_photon_drmin},reduce_maxfirst)
      .Name("lead_photon_drmin");

  const NamedFunc sig_mu_pt = FilterNamedFunc("mu_pt","mu_sig").Name("sig_mu_pt");
  const NamedFunc sig_mu_absdz = MapNamedFunc(FilterNamedFunc("mu_dz","mu_sig"),
      map_abs).Name("sig_mu_absdz");
  const NamedFunc sig_mu_absdxy = MapNamedFunc(FilterNamedFunc("mu_dxy","mu_sig"),
      map_abs).Name("sig_mu_absdxy");
  const NamedFunc lead_mu_absdz = MultiReduceNamedFunc({sig_mu_pt,sig_mu_absdz},
      reduce_maxfirst).Name("lead_mu_absdz");
  const NamedFunc sublead_mu_absdz = MultiReduceNamedFunc({sig_mu_pt,sig_mu_absdz},
      reduce_subleadfirst).Name("sublead_mu_absdz");
  const NamedFunc lead_mu_absdxy = MultiReduceNamedFunc({sig_mu_pt,sig_mu_absdxy},
      reduce_maxfirst).Name("lead_mu_absdxy");
  const NamedFunc sublead_mu_absdxy = MultiReduceNamedFunc({sig_mu_pt,sig_mu_absdxy},
      reduce_subleadfirst).Name("sublead_mu_absdxy");
  const NamedFunc sig_el_pt = FilterNamedFunc("el_pt","el_sig").Name("sig_el_pt");
  const NamedFunc sig_el_absdz = MapNamedFunc(FilterNamedFunc("el_dz","el_sig"),
      map_abs).Name("sig_el_absdz");
  const NamedFunc sig_el_absdxy = MapNamedFunc(FilterNamedFunc("el_dxy","el_sig"),
      map_abs).Name("sig_el_absdxy");
  const NamedFunc lead_el_absdz = MultiReduceNamedFunc({sig_el_pt,sig_el_absdz},
      reduce_maxfirst).Name("lead_el_absdz");
  const NamedFunc sublead_el_absdz = MultiReduceNamedFunc({sig_el_pt,sig_el_absdz},
      reduce_subleadfirst).Name("sublead_el_absdz");
  const NamedFunc lead_el_absdxy = MultiReduceNamedFunc({sig_el_pt,sig_el_absdxy},
      reduce_maxfirst).Name("lead_el_absdxy");
  const NamedFunc sublead_el_absdxy = MultiReduceNamedFunc({sig_el_pt,sig_el_absdxy},
      reduce_subleadfirst).Name("sublead_el_absdxy");

  const NamedFunc nmu_med = ReduceNamedFunc("mu_sig&&mu_mediumid",reduce_sum).Name("nmu_med");
  const NamedFunc nel_med = ReduceNamedFunc("el_sig&&el_id90",reduce_sum).Name("nel_med");
  const NamedFunc nlep_med = nmu_med+nel_med;

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  //Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;

  string prod_folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/");
  //std::set<int> years = {2016,2017,2018};
  std::set<std::string> years = {"2016APV","2016","2017","2018"};
  std::string lumi_string = "138";
  //std::set<int> years = {2016};
  //std::string lumi_string = "36";

  bool multiply_sig = true;
  std::string sig_name = "gg#rightarrow H#rightarrow Z#gamma";
  NamedFunc weight = NamedFunc(w_ewkzgamma*w_years).Name("weight");
  NamedFunc weight_noscale = NamedFunc(w_ewkzgamma*w_years).Name("weight");
  NamedFunc w_earlyrun3 = NamedFunc(199.9/137.6);
  if (multiply_sig) {
    sig_name = "gg#rightarrow H#rightarrow Z#gamma (x50)";
    weight = w_ewkzgamma*w_sigx50*w_years;
    weight.Name("weight_sigx50");
  }
  Palette mc_colors("txt/colors_zgamma.txt","default");

  std::cout << "Loading samples:\n";
  //auto proc_smzg      = Process::MakeShared<Baby_pico>("Z+#gamma", back, mc_colors("zgtollg"),
  //                         attach_folder(prod_folder,years,"mc/unskimmed",{"*ZGToLLG*"}), "1");
  //auto proc_dy        = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, mc_colors("dyjets"),
  //                         attach_folder(prod_folder,years,"mc/unskimmed",{"*DYJets*"}), "stitch_dy");
  //auto proc_tt        = Process::MakeShared<Baby_pico>("tt", back, mc_colors("tt"),
  //                         attach_folder(prod_folder,years,"mc/unskimmed",{"*TTTo2L2Nu*"}), "1");
  //auto proc_tt1l      = Process::MakeShared<Baby_pico>("Fake lepton", back, mc_colors("fakelep"),
  //                         attach_folder(prod_folder,years,"mc/unskimmed",{"*TTG*"}), 
  //                         "(ntruel+ntrumu)==1");
  ////auto proc_hzg       = Process::MakeShared<Baby_pico>(sig_name, sig, kRed,
  ////                         attach_folder(prod_folder,years,"mc/unskimmed",
  ////                         {"*HToZG*M-125*","*HToZG_M125*"}), "1");
  //auto proc_hzg       = Process::MakeShared<Baby_pico>(sig_name, sig, kRed,
  //                         attach_folder(prod_folder,years,"mc/unskimmed",
  //                         {"*GluGluHToZG*M-125*"}), "1");
  //vector<shared_ptr<Process>> procs = {proc_dy, proc_tt, proc_smzg, proc_tt1l, proc_hzg};
  vector<shared_ptr<Process>> procs = {};

  //vector<shared_ptr<Process>> procs_all = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMoreTrigs");
  vector<shared_ptr<Process>> procs_all = {};

  auto proc_gghzg = Process::MakeShared<Baby_pico>("gg#rightarrow H#rightarrow Z#gamma", sig, kRed,
                        attach_folder(prod_folder,years,"mc/unskimmed",
                        {"*GluGluHToZG*M-125*"}), "1");
  auto proc_qqhzg = Process::MakeShared<Baby_pico>("VBF H#rightarrow Z#gamma", sig, kGreen,
                        attach_folder(prod_folder,years,"mc/unskimmed",
                        {"*VBFHToZG*M-125*"}), "1");
  auto proc_whzg  = Process::MakeShared<Baby_pico>("WH#rightarrow Z#gamma", sig, kBlue,
                        attach_folder(prod_folder,years,"mc/unskimmed",
                        {"*WplusH_HToZG*M-125*","*WminusH_HToZG*M-125*"}), "1");
  auto proc_zhzg  = Process::MakeShared<Baby_pico>("ZH#righarrow Z#gamma", sig, kBlue,
                        attach_folder(prod_folder,years,"mc/unskimmed",
                        {"*ZH_HToZG*M-125*"}), "1");
  auto proc_tthzg = Process::MakeShared<Baby_pico>("ttH#rightarrow Z#gamma", sig, kBlue,
                        attach_folder(prod_folder,years,"mc/unskimmed",
                        {"*ttHToZG_M125*"}), "1");
  vector<shared_ptr<Process>> procs_signal = {proc_gghzg, proc_qqhzg, proc_whzg, proc_zhzg, proc_tthzg};

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.Title(TitleType::info)
          .Stack(StackType::signal_overlay)
          .Overflow(OverflowType::none)
          .Bottom(BottomType::off)
          .UseCMYK(false)
          .LegendColumns(1)
          .ShowBackgroundError(false)
          .FileExtensions({"pdf"});
  PlotOpt shapes_norm("txt/plot_styles.txt", "CMSPaper");
  shapes_norm.Title(TitleType::info)   
             .Bottom(BottomType::off)
             .YAxis(YAxisType::linear)
             .Stack(StackType::shapes)
             .LegendColumns(3);
  vector<PlotOpt> ops = {lin_lumi, shapes_norm};
  vector<PlotOpt> ops_linlumi = {PlotOpt("txt/plot_styles.txt","LinLumi"),PlotOpt("txt/plot_styles.txt","Shapes")}; 

  //std::string lumi_string = "41.5";

  //std::string base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/signal/";
  //std::string base_folder_ul = "/net/cms17/cms17r0/pico/NanoAODv2/nano/2017/";
  //std::string base_folder_ulv9 = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2016/";
  //std::string base_folder_v7data = "/net/cms17/cms17r0/pico/NanoAODv7/nano/2017/data/";

  //std::vector<std::shared_ptr<Process>> procs;
  //procs.push_back(Process::MakeShared<Baby_nano>("HToZG", 
  //    Process::Type::signal, kBlack,
  //    {base_folder+"*GluGluHToZG_ZToLL*"},"1"));

  ////2017 triggers
  //NamedFunc hlt_el_trigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  //NamedFunc hlt_mu_trigger = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8";
  //NamedFunc hlt_single_el_trigger = "HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  //NamedFunc hlt_single_mu_trigger = "HLT_IsoMu27||HLT_IsoMu24";
  ////
  //NamedFunc hlt_el_trigger_plus = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele35_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG";
  //NamedFunc hlt_mu_trigger_plus = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24";
  //NamedFunc hlt_mu_trigger_plusplus = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8||HLT_IsoMu27||HLT_IsoMu24||HLT_Mu17_Photon30_IsoCaloId";
  //NamedFunc l1_el_trigger = "L1_SingleEG24||L1_SingleEG34er2p1||L1_SingleIsoEG24er2p1||L1_SingleIsoEG24||L1_DoubleEG_18_17||L1_DoubleEG_22_10||L1_DoubleEG_LooseIso23_10||L1_DoubleEG_LooseIso24_10";
  //NamedFunc l1_mu_trigger = "L1_DoubleMu_12_5||L1_DoubleMu_15_7_SQ||L1_DoubleMu_15_7_SQ_Mass_Min4";

  //std::vector<std::shared_ptr<Process>> full_procs;
  ////full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  ////    Process::Type::background, TColor::GetColor("#ffb400"),
  ////    {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-"
  ////    "pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000*"},"1"));
  ////full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
  ////    Process::Type::background, TColor::GetColor("#16bac5"),
  ////    {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_*"},"1"));
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  //    Process::Type::background, TColor::GetColor("#ffb400"),
  //    {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__035CF0D8-DDC5-6D46-B172-3BDC98CC624D*",
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__223A6F1C-C8C0-9F43-94C3-233AD51ABD32*",
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__29EF0133-B16A-E747-9E02-89BC6682B669*", //temp?
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__25118058-0F30-8A41-B82B-C7EAFEE55379*", //temp?
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__24950148-3722-B84D-9050-C996F5205F0A*"},"1"));
  //full_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+#gamma", 
  //    Process::Type::background, TColor::GetColor("#16bac5"),
  //    {base_folder_ul+"mc/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__230000__03E636EE-BAB0-CA42-BFB1-28E772F2D53E*"},"1"));
  //full_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
  //    Process::Type::signal, TColor::GetColor("#ff0000"),
  //    {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  //std::vector<std::shared_ptr<Process>> dy_procs;
  //dy_procs.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  //    Process::Type::background, TColor::GetColor("#ffb400"),
  //    {base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__035CF0D8-DDC5-6D46-B172-3BDC98CC624D*",
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__223A6F1C-C8C0-9F43-94C3-233AD51ABD32*",
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__29EF0133-B16A-E747-9E02-89BC6682B669*", //temp?
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__25118058-0F30-8A41-B82B-C7EAFEE55379*", //temp?
  //    base_folder_ul+"mc/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv2__106X_mc2017_realistic_v8-v1__260000__24950148-3722-B84D-9050-C996F5205F0A*"},"1"));

  //std::vector<std::shared_ptr<Process>> dy_procs_2016;
  //dy_procs_2016.push_back(Process::MakeShared<Baby_nano>("Z/#gamma*+jets", 
  //    Process::Type::background, TColor::GetColor("#ffb400"),
  //    {base_folder_ulv9+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0082C29D-E74C-024A-BE9B-97B29EE7A4A2*",
  //    base_folder_ulv9+"mc/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL16NanoAODv9__106X_mcRun2_asymptotic_v17-v1__30000__0718C107-8960-6B44-B96A-C60D53D52A95*"},"1"));

  //std::vector<std::shared_ptr<Process>> signal_trigger_procs;
  //signal_trigger_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (Dilepton Triggers)", 
  //    Process::Type::signal, TColor::GetColor("#ff0000"),
  //    {base_folder_ul+"signal/GluGluHToZG*"},hlt_el_trigger||hlt_mu_trigger));
  //signal_trigger_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (Only Single Lepton Triggers)", 
  //    Process::Type::signal, TColor::GetColor("#0000ff"),
  //    {base_folder_ul+"signal/GluGluHToZG*"},(hlt_single_el_trigger||hlt_single_mu_trigger)&&!(hlt_el_trigger||hlt_mu_trigger)));
  //
  //std::vector<std::shared_ptr<Process>> signal_procs;
  //signal_procs.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG (x100)", 
  //    Process::Type::background, TColor::GetColor("#ff0000"),
  //    {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  //std::vector<std::shared_ptr<Process>> signal_procs_noscale;
  //signal_procs_noscale.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
  //    Process::Type::signal, TColor::GetColor("#ff0000"),
  //    {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  //std::vector<std::shared_ptr<Process>> signal_procs_noscale_2d;
  //signal_procs_noscale_2d.push_back(Process::MakeShared<Baby_nano>("gg#rightarrow H#rightarrow ZG", 
  //    Process::Type::background, TColor::GetColor("#ff0000"),
  //    {base_folder_ul+"signal/GluGluHToZG*"},"1"));

  //std::vector<std::shared_ptr<Process>> procs_onelep_data;
  //procs_onelep_data.push_back(Process::MakeShared<Baby_nano>("Data", 
  //    Process::Type::data, kBlack,
  //    //{base_folder_v7data+"SingleElectron__Run2017B*",
  //    //base_folder_v7data+"SingleMuon__Run2017B*"},"1"));
  //    {base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__014FD4FF-A181-2843-B548-BB3198F16E88*",
  //    base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__08B5DC81-D780-A54B-80FE-C94EBD267ACA*",
  //    base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__0ACA0341-E365-2544-99F0-FBA4C92C0301*",
  //    base_folder_v7data+"SingleElectron__Run2017B__02Apr2020-v1__230000__0C083B57-67B2-FF43-AF2E-56C850A094D6*",
  //    base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__0C07D411-3FE0-974E-A5A4-25BB763E4038*",
  //    base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__0D1DA690-909B-F74F-8DB7-48D6F786F8D8*",
  //    base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__12EF443F-9CD8-7842-9563-2BBD627A28FB*",
  //    base_folder_v7data+"SingleMuon__Run2017B__02Apr2020-v1__230000__1A14183C-9C11-E947-9D71-6B89994127D9*"},"1"));
  //    
  //std::vector<std::shared_ptr<Process>> procs_met_data;
  //procs_met_data.push_back(Process::MakeShared<Baby_nano>("Data", 
  //    Process::Type::data, kBlack,
  //    {base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__00DCCA4E-F5F1-F84D-A6EC-2956ACAB6E02*",
  //    base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__0B4DB0F0-168B-544A-B5F9-A5B3ACFE7F1E*",
  //    base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__0F13F71C-9A13-5641-A5C1-4955D60DE44A*",
  //    base_folder_v7data+"MET__Run2017B__02Apr2020-v1__230000__1781BB70-1AD7-2F49-8BF6-A6DFE7C58167*",
  //    base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__034FA7B4-A09F-0D47-B752-2897BE4FE43F*",
  //    base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__0D3E9423-8FB8-9548-B63E-B94913999718*",
  //    base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__0EE78B5A-1800-0341-A33E-29EFDCAB3FAA*",
  //    base_folder_v7data+"MET__Run2017F__02Apr2020-v1__30000__1450C885-647B-074B-BAEB-84763C291B92*"},"1"));

  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  //ops should include shapes

  //new stuff
  bool plot_truth = false;
  bool plot_objects = false;
  bool plot_signal_acceptance = true;

  if (plot_signal_acceptance) {

    pm.Push<Table>("zgacc_objects", vector<TableRow>{
      TableRow("No Selection, $Z\\to\\mu\\mu$", 
          has_z_to_mumu,0,0,weight_noscale*w_earlyrun3),
      TableRow("Objects in acceptance", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1,0,0,weight_noscale*w_earlyrun3),
      TableRow("Leptons reco'd", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1&&"nmu>=2",0,0,weight_noscale*w_earlyrun3),
      TableRow("Objects reco'd", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1&&"nmu>=2&&nphoton>=1",0,0,weight_noscale*w_earlyrun3),
      TableRow("Trigger", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1&&"nmu>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton),
          0,0,weight_noscale*w_earlyrun3),
      TableRow("1 lepton medium ID", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1&&"nmu>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton)&&nmu_med>=1,
          0,0,weight_noscale*w_earlyrun3),
      TableRow("Lead lepton pt>20", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1&&"nmu>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton)&&lead_mu_pt>=20,
          0,0,weight_noscale*w_earlyrun3),
      TableRow("Sublead lepton pt>10", 
          has_z_to_mumu&&nmc_mu_acc>=2&&nmc_photon_acc>=1&&"nmu>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton)&&sublead_mu_pt>=10,
          0,0,weight_noscale*w_earlyrun3),
      TableRow("\\hline No Selection, $Z\\to ee$", 
          has_z_to_ee,0,0,weight_noscale*w_earlyrun3),
      TableRow("Objects in acceptance", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1,0,0,weight_noscale*w_earlyrun3),
      TableRow("Leptons reco'd", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1&&"nel>=2",0,0,weight_noscale*w_earlyrun3),
      TableRow("Objects reco'd", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1&&"nel>=2&&nphoton>=1",0,0,weight_noscale*w_earlyrun3),
      TableRow("Trigger", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1&&"nel>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton),
          0,0,weight_noscale*w_earlyrun3),
      TableRow("1 lepton medium ID", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1&&"nel>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton)&&nel_med>=1,
          0,0,weight_noscale*w_earlyrun3),
      TableRow("Lead lepton pt>20", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1&&"nel>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton)&&lead_el_pt>=20,
          0,0,weight_noscale*w_earlyrun3),
      TableRow("Sublead lepton pt>10", 
          has_z_to_ee&&nmc_el_acc>=2&&nmc_photon_acc>=1&&"nel>=2&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton)&&sublead_el_pt>=10,
          0,0,weight_noscale*w_earlyrun3),
      TableRow("\\hline $Z\\to\\mu\\mu$ and three leptons", 
          has_z_to_mumu&&nmc_lep>=3,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects in acceptance", 
          has_z_to_mumu&&nmc_lep>=3&&nmc_lep_acc>=3&&nmc_photon_acc>=1,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects reco'd", 
          has_z_to_mumu&&nmc_lep>=3&&nmc_lep_acc>=3&&nmc_photon_acc>=1&&"nlep>=3&&nphoton>=1",0,0,weight_noscale*w_earlyrun3),
      TableRow("Trigger", 
          has_z_to_mumu&&nmc_lep>=3&&nmc_lep_acc>=3&&nmc_photon_acc>=1&&"nlep>=3&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton),0,0,weight_noscale*w_earlyrun3),
      TableRow("\\hline $Z\\to ee$ and three leptons", 
          has_z_to_ee&&nmc_lep>=3,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects in acceptance", 
          has_z_to_ee&&nmc_lep>=3&&nmc_lep_acc>=3&&nmc_photon_acc>=1,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects reco'd", 
          has_z_to_ee&&nmc_lep>=3&&nmc_lep_acc>=3&&nmc_photon_acc>=1&&"nlep>=3&&nphoton>=1",0,0,weight_noscale*w_earlyrun3),
      TableRow("Trigger", 
          has_z_to_ee&&nmc_lep>=3&&nmc_lep_acc>=3&&nmc_photon_acc>=1&&"nlep>=3&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton),0,0,weight_noscale*w_earlyrun3),
      TableRow("\\hline $Z\\to\\mu\\mu$ and four leptons", 
          has_z_to_mumu&&nmc_lep>=4,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects in acceptance", 
          has_z_to_mumu&&nmc_lep>=4&&nmc_lep_acc>=4&&nmc_photon_acc>=1,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects reco'd", 
          has_z_to_mumu&&nmc_lep>=4&&nmc_lep_acc>=4&&nmc_photon_acc>=1&&"nlep>=4&&nphoton>=1",0,0,weight_noscale*w_earlyrun3),
      TableRow("Trigger", 
          has_z_to_mumu&&nmc_lep>=4&&nmc_lep_acc>=4&&nmc_photon_acc>=1&&"nlep>=4&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton),0,0,weight_noscale*w_earlyrun3),
      TableRow("\\hline $Z\\to ee$ and four leptons", 
          has_z_to_ee&&nmc_lep>=4,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects in acceptance", 
          has_z_to_ee&&nmc_lep>=4&&nmc_lep_acc>=4&&nmc_photon_acc>=1,0,0,weight_noscale*w_earlyrun3),
      TableRow("Obects reco'd", 
          has_z_to_ee&&nmc_lep>=4&&nmc_lep_acc>=4&&nmc_photon_acc>=1&&"nlep>=4&&nphoton>=1",0,0,weight_noscale*w_earlyrun3),
      TableRow("Trigger", 
          has_z_to_ee&&nmc_lep>=4&&nmc_lep_acc>=4&&nmc_photon_acc>=1&&"nlep>=4&&nphoton>=1"
          &&(ZgFunctions::HLT_pass_singlelepton||ZgFunctions::HLT_pass_dilepton),0,0,weight_noscale*w_earlyrun3),
    },procs_signal,false,true,false,false,false,true).LuminosityTag(lumi_string).Precision(2);
  }

  if (plot_objects) {
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_mu_pt, "Lead muon p_{T} [GeV]", {}),
        zg_baseline_nolep&&"ll_lepid[0]==13",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, sublead_mu_pt, "Sublead muon p_{T} [GeV]", {}),
        zg_baseline_nolep&&"ll_lepid[0]==13",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_el_pt, "Lead electron p_{T} [GeV]", {}),
        zg_baseline_nolep&&"ll_lepid[0]==11",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, sublead_el_pt, "Sublead electron p_{T} [GeV]", {}),
        zg_baseline_nolep&&"ll_lepid[0]==11",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, lead_mu_absdz, "Lead muon |d_{z}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==13",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, sublead_mu_absdz, "Sublead muon |d_{z}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==13",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, lead_mu_absdxy, "Lead muon |d_{xy}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==13",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, sublead_mu_absdxy, "Sublead muon |d_{xy}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==13",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, lead_el_absdz, "Lead electron |d_{z}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==11",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, sublead_el_absdz, "Sublead electron |d_{z}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==11",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, lead_el_absdxy, "Lead electron |d_{xy}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==11",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(60, 0.0, 0.5, sublead_el_absdxy, "Sublead electron |d_{xy}|", {}),
        zg_baseline_nolep&&"ll_lepid[0]==11",
        procs_all, ops_linlumi).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Table>("zgacc_objects", vector<TableRow>{
      TableRow("Baseline and Z->ee", 
          zg_baseline_nolep&&"ll_lepid[0]==11",0,0,weight_noscale),
      TableRow("Electron pt > (25,15) or (35/30,any)", 
          zg_baseline_nolep&&"ll_lepid[0]==11"&&((lead_el_pt>25&&sublead_el_pt>15)||(lead_el_pt>35||(lead_el_pt>30&&year<2017))),0,0,weight_noscale),
      TableRow("Baseline and Z->\\mu\\mu", 
          zg_baseline_nolep&&"ll_lepid[0]==13",0,0,weight_noscale),
      TableRow("Muon pt > (20,10) or (25,any)", 
          zg_baseline_nolep&&"ll_lepid[0]==13"&&((lead_mu_pt>20&&sublead_mu_pt>10)||(lead_mu_pt>25)),0,0,weight_noscale),
    },procs_all,false,true,false,false,false,true).LuminosityTag(lumi_string).Precision(2);
  }

  if (plot_truth) {
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, mc_lead_mu_pt, "Generator lead muon p_{T} [GeV]", {}),
        nmc_mu>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, mc_sublead_mu_pt, "Generator sublead muon p_{T} [GeV]", {}),
        nmc_mu>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, mc_lead_el_pt, "Generator lead electron p_{T} [GeV]", {}),
        nmc_el>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, mc_sublead_el_pt, "Generator sublead electron p_{T} [GeV]", {}),
        nmc_el>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, mc_lead_photon_pt, "Generator photon p_{T} [GeV]", {}),
        nmc_photon>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);

    pm.Push<Hist1D>(Axis(30, -4.0, 4.0, mc_mu_eta, "Generator muon #eta", {}),
        nmc_mu>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 5.0, drloose_photon_drmin, "Photon minimum #Delta R_{#gamma,l}", {}),
        ndrloose_photon>0.0&&"nlep>0",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);

    pm.Push<Hist1D>(Axis(30, -4.0, 4.0, mc_lead_mu_eta, "Generator lead muon #eta", {}),
        nmc_mu>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, -4.0, 4.0, mc_sublead_mu_eta, "Generator sublead muon #eta", {}),
        nmc_mu>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, -4.0, 4.0, mc_lead_el_eta, "Generator lead electron #eta", {}),
        nmc_el>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, -4.0, 4.0, mc_sublead_el_eta, "Generator sublead electron #eta", {}),
        nmc_el>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, -4.0, 4.0, mc_lead_photon_eta, "Generator photon #eta", {}),
        nmc_photon>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    //bins - xx xl xm xt ll lm lt mm mt tt
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, mu_quality, "Muon quality", {}),
        nmc_mu>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(10, -0.5, 9.5, el_quality, "Electron quality", {}),
        nmc_el>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    //reco objects
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_mu_pt, "Lead muon p_{T} [GeV]", {}),
        "nmu>0",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, sublead_mu_pt, "Sublead muon p_{T} [GeV]", {}),
        "nmu>1",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_el_pt, "Lead electron p_{T} [GeV]", {}),
        "nel>0",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, sublead_el_pt, "Sublead electron p_{T} [GeV]", {}),
        "nel>1",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 0.5, lead_mu_reliso, "Lead muon I_{rel}", {}),
        nisoloose_mu>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 0.5, sublead_mu_reliso, "Sublead muon I_{rel}", {}),
        nisoloose_mu>1,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_photon_pt, "Lead photon p_{T} [GeV]", {}),
        nptloose_photon>0.0,
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 5.0, lead_photon_drmin, "Lead photon minimum #Delta R_{#gamma,l}", {}),
        ndrloose_photon>0.0&&"nlep>0",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    //eliminate fsr
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_mu_pt, "Lead muon p_{T} [GeV]", {}),
        "nmu>1&&ll_m[0]>81&&ll_m[0]<101",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, sublead_mu_pt, "Sublead muon p_{T} [GeV]", {}),
        "nmu>1&&ll_m[0]>81&&ll_m[0]<101",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, lead_el_pt, "Lead electron p_{T} [GeV]", {}),
        "nel>1&&ll_m[0]>81&&ll_m[0]<101",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 100.0, sublead_el_pt, "Sublead electron p_{T} [GeV]", {}),
        "nel>1&&ll_m[0]>81&&ll_m[0]<101",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
    pm.Push<Hist1D>(Axis(30, 0.0, 5.0, lead_photon_drmin, "Lead photon minimum #Delta R_{#gamma,l}", {}),
        ndrloose_photon>0.0&&"nll>0&&ll_m[0]>81&&ll_m[0]<91",
        procs, ops).Weight(weight).Tag("zgacc")
        .LuminosityTag(lumi_string);
  }
  // /new stuff

  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Muon p_{T}", {}),
  //    (Truth_zdecay_pdgid==13)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==13||GenPart_pdgId==-13",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genmuonpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(40, -4.0, 4.0, "GenPart_eta", "Gen Muon #eta", {}),
  //    (Truth_zdecay_pdgid==13)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==13||GenPart_pdgId==-13",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genmuoneta")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Lead Muon p_{T}", {}),
  //    (Truth_zdecay_pdgid==13)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==13||GenPart_pdgId==-13"&&GenPart_isLeadLepton,
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genleadmuonpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Sublead Muon p_{T}", {}),
  //    (Truth_zdecay_pdgid==13)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==13||GenPart_pdgId==-13"&&!GenPart_isLeadLepton,
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_gensublmuonpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Electron p_{T}", {}),
  //    (Truth_zdecay_pdgid==11)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==11||GenPart_pdgId==-11",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genelectronpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(40, -4.0, 4.0, "GenPart_eta", "Gen Electron #eta", {}),
  //    (Truth_zdecay_pdgid==11)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==11||GenPart_pdgId==-11",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genelectroneta")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Lead Electron p_{T}", {}),
  //    (Truth_zdecay_pdgid==11)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==11||GenPart_pdgId==-11"&&GenPart_isLeadLepton,
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genleadelectronpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Sublead Electron p_{T}", {}),
  //    (Truth_zdecay_pdgid==11)&&(GenPart_pdgIdMother==23)&&"GenPart_pdgId==11||GenPart_pdgId==-11"&&!GenPart_isLeadLepton,
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_gensublelectronpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Photon p_{T}", {}),
  //    (GenPart_pdgIdMother==25)&&"GenPart_pdgId==22",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genphotonpt")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(40, -4.0, 4.0, "GenPart_eta", "Gen Photon #eta", {}),
  //    (GenPart_pdgIdMother==25)&&"GenPart_pdgId==22",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genphotoneta")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Muon p_{T}", {}),
  //    (GenPart_pdgIdMother==23)&&"GenPart_pdgId==13||GenPart_pdgId==-13&&(GenPart_eta>-2.4&&GenPart_eta<2.4)",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genmuonpt_etacut")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Electron p_{T}", {}),
  //    (GenPart_pdgIdMother==23)&&"GenPart_pdgId==11||GenPart_pdgId==-11&&(GenPart_eta>-2.5&&GenPart_eta<2.5)",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genelectronpt_etacut")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Hist1D>(Axis(30, 0.0, 100.0, "GenPart_pt", "Gen Photon p_{T}", {}),
  //    (GenPart_pdgIdMother==25)&&"GenPart_pdgId==22&&(GenPart_eta>-2.5&&GenPart_eta<2.5)",
  //    signal_procs_noscale, plt_lin).Weight(weight_noscale)
  //    .Tag("FixName:acc_genphotonpt_etacut")
  //    .LuminosityTag(lumi_string);
  //pm.Push<Table>("zgaccept_"+options.year_string, vector<TableRow>{
  //  TableRow("\\hline $Z\\rightarrow e^{+}e^{-}$ decays", 
  //      (Truth_zdecay_pdgid==11),0,0,weight_noscale),
  //  TableRow("Electrons in acceptance", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAccept,0,0,weight_noscale),
  //  TableRow("Photons in acceptance", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAccept&&Truth_PhotonAccept,0,0,weight_noscale),
  //  TableRow("Electrons in scale factor acceptance", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAccept,0,0,weight_noscale),
  //  TableRow("Photons in scale factor acceptance", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt,0,0,weight_noscale),
  //  TableRow("Electrons reconstructed", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco,0,0,weight_noscale),
  //  TableRow("Photons reconstructed", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco,0,0,weight_noscale),
  //  TableRow("L1 trigger", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco&&l1_el_trigger,0,0,weight_noscale),
  //  TableRow("HLT Trigger", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale),
  //  TableRow("HLT Trigger (Data)", 
  //      (Truth_zdecay_pdgid==11)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco&&l1_el_trigger&&hlt_el_trigger,0,0,weight_noscale*ScaleFactor_Triggers),

  //  TableRow("\\hline $Z\\rightarrow \\mu^{+}\\mu^{-}$ decays", 
  //      (Truth_zdecay_pdgid==13),0,0,weight_noscale),
  //  TableRow("Muons in acceptance", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAccept,0,0,weight_noscale),
  //  TableRow("Photons in acceptance", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAccept&&Truth_PhotonAccept,0,0,weight_noscale),
  //  TableRow("Muons in scale factor acceptance", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAccept,0,0,weight_noscale),
  //  TableRow("Photons in scale factor acceptance", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt,0,0,weight_noscale),
  //  TableRow("Muons reconstructed", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco,0,0,weight_noscale),
  //  TableRow("Photons reconstructed", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco,0,0,weight_noscale),
  //  TableRow("L1 trigger", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco&&l1_mu_trigger,0,0,weight_noscale),
  //  TableRow("HLT Trigger", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale),
  //  TableRow("HLT Trigger (Data)", 
  //      (Truth_zdecay_pdgid==13)&&Truth_LeptonAcceptHighPt&&Truth_PhotonAcceptHighPt&&Truth_LeptonReco&&Truth_PhotonReco&&l1_mu_trigger&&hlt_mu_trigger,0,0,weight_noscale*ScaleFactor_Triggers),
  //},signal_procs_noscale,false,true,false,false,false,true).LuminosityTag(lumi_string).Precision(2);
  
  pm.min_print_ = true;
  pm.MakePlots(1.);

}
