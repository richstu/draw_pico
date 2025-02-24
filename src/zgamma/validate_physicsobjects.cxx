/**
 * Script to validate physics objects between pico productions
 */

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TLorentzVector.h"

#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/plot_maker.hpp"
#include "core/process.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

#include "core/baby.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"

using PlotOptTypes::OverflowType;
using std::invalid_argument;
using std::shared_ptr;
using std::string;
using std::set;
using std::vector;
using ZgFunctions::lead_jet_pt;
using ZgFunctions::lead_jet_eta;
using ZgFunctions::lead_jet_phi;
using ZgFunctions::lead_jet_m;
using ZgFunctions::sublead_jet_pt;
using ZgFunctions::sublead_jet_eta;
using ZgFunctions::sublead_jet_phi;
using ZgFunctions::sublead_jet_m;
using ZgFunctions::w_years;
using ZgFunctions::photon_relpterr;
using ZgUtilities::SetProcessesBackground;
using ZgUtilities::ZgSampleLoader;

int main() {

  //---------------------------------------------------------------------------
  //                                       settings
  //---------------------------------------------------------------------------

  string production = "pinnaclesv0";
  string years = "2016";

  //---------------------------------------------------------------------------
  //                                    initialization
  //---------------------------------------------------------------------------

  gErrorIgnoreLevel = 6000;

  //set lumi string
  string lumi_string = "138";
  set<string> year_set;
  set<string> data_dirs;
  if (years == "2016APV") {
    lumi_string = "20";
    year_set = {"2016APV"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9UCSB2/htozgamma_pinnacles_v0/2016APV/"};
  }
  else if (years == "2016") {
    lumi_string = "17";
    year_set = {"2016"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/"};
  }
  else if (years == "2017") {
    lumi_string = "41";
    year_set = {"2017"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/"};
  }
  else if (years == "2018") {
    lumi_string = "60";
    year_set = {"2018"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/"};
  }
  else if (years == "2022") {
    lumi_string = "8";
    year_set = {"2022"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/"};
  }
  else if (years == "2022EE") {
    lumi_string = "27";
    year_set = {"2022EE"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/"};
  }
  else if (years == "2023") {
    lumi_string = "18";
    year_set = {"2023"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/"};
  }
  else if (years == "2023BPix") {
    lumi_string = "10";
    year_set = {"2023BPix"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/"};
  }
  else if (years == "Run2") {
    lumi_string = "138";
    year_set = {"2016APV","2016","2017","2018"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv9UCSB2/htozgamma_pinnacles_v0/2016APV/",
                 "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2016/",
                 "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2017/",
                 "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_pinnacles_v0/2018/"};
  }
  else if (years == "Run3") {
    lumi_string = "62";
    year_set = {"2022","2022EE","2023","2023BPix"};
    data_dirs = {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022/",
                 "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2022EE/",
                 "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023/",
                 "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_pinnacles_v0/2023BPix/"};
  }
  else {
    throw invalid_argument("Unknown year");
  }

  //load processes
  vector<shared_ptr<Process>> procs;
  vector<shared_ptr<Process>> procs_mc;
  vector<shared_ptr<Process>> procs_data;
  if (production == "pinnaclesv0") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("DATA_PRODUCTION_YEARS",data_dirs)
        .LoadSamples("txt/samples_zgamma.txt","llcr");
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("DATA_PRODUCTION_YEARS",data_dirs)
        .LoadSamples("txt/samples_zgamma.txt","llcr",
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("DATA_PRODUCTION_YEARS",data_dirs)
        .LoadSamples("txt/samples_zgamma.txt","llcr",
        {Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else if (production == "lassenv0") {
    set<string> production_names = {
        "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/",
        "/net/cms11/cms11r0/pico/NanoAODv11p9/htozgamma_lassen_v0/"};
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("PRODUCITON_NAME",production_names)
        .LoadSamples("txt/samples_zgamma.txt","llcr");
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("PRODUCITON_NAME",production_names)
        .LoadSamples("txt/samples_zgamma.txt","llcr",
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .SetMacro("PRODUCITON_NAME",production_names)
        .LoadSamples("txt/samples_zgamma.txt","llcr",{Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else if (production == "kingscanyonv1") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv1.txt","llcr");
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv1.txt","llcr",
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv1.txt","llcr",
        {Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else if (production == "kingscanyonv0") {
    procs = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt","llcr");
    procs_mc = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt","llcr",
        {Process::Type::background});
    procs_data = ZgSampleLoader().SetMacro("YEARS",year_set)
        .LoadSamples("txt/samples_zgamma_kingscanyonv0.txt","llcr",
        {Process::Type::data});
    SetProcessesBackground(procs_data);
  }
  else {
    throw invalid_argument("unknown production");
  }

  //load plot opts
  vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumiData")
                         .RatioMinimum(0.6).RatioMaximum(1.4)};
  vector<PlotOpt> ops_root = {PlotOpt("txt/plot_styles.txt","LinLumiDataRoot")};
  vector<PlotOpt> ops_mc = {PlotOpt("txt/plot_styles.txt","LinLumi")}; //for debugging
  vector<PlotOpt> ops_log = {PlotOpt("txt/plot_styles.txt","LogLumiDataRoot")};
  vector<PlotOpt> ops2d = {PlotOpt("txt/plot_styles.txt","LinLumi2DRoot")};

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------

  NamedFunc zee_ctrlregion = NamedFunc("nel>=2&&nphoton==0&&pass&&ll_m[0]>80&&ll_m[0]<100").Name("Ztoee_CR");
  NamedFunc zmm_ctrlregion = NamedFunc("nmu>=2&&nphoton==0&&pass&&ll_m[0]>80&&ll_m[0]<100").Name("Ztomumu_CR");
  NamedFunc zeey_ctrlregion = NamedFunc("nel>=2&&nphoton==1&&pass&&llphoton_m[0]>80&&llphoton_m[0]<100").Name("Ztoeegamma_CR");
  NamedFunc zmmy_ctrlregion = NamedFunc("nmu>=2&&nphoton==1&&pass&&llphoton_m[0]>80&&llphoton_m[0]<100").Name("Ztomumugamma_CR");
  NamedFunc zeeyjj_ctrlregion = NamedFunc(zeey_ctrlregion&&"njet>=2")
                                .Name("Ztoeegamma_jj_CR");
  NamedFunc zmmyjj_ctrlregion = NamedFunc(zmmy_ctrlregion&&"njet>=2")
                                .Name("Ztomumugamma_jj_CR");

  //in case MET looks weird
  NamedFunc pass_ra2b_filters = NamedFunc(
      "pass_muon_jet&&pass_low_neutral_jet&&pass_htratio_dphi_tight&&pass_ecalnoisejet")
      .Name("RA2b_filters");

  //alias to avoid naming confusion
  NamedFunc llphoton_coscaptheta0 = NamedFunc("llphoton_cosTheta[0]").Name("llphoton_coscaptheta0");

  //returns ll mass calculated with corrected muon pts
  const NamedFunc ll_corrected_m("ll_corrected_m",
      [](const Baby &b) -> NamedFunc::ScalarType{
    //only recalculate for Muons
    if (b.ll_lepid()->at(0)==11) return b.ll_m()->at(0);
    TLorentzVector m1, m2;
    m1.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i1()->at(0)),
                    b.mu_eta()->at(b.ll_i1()->at(0)),
                    b.mu_phi()->at(b.ll_i1()->at(0)),
                    0.106);
    m2.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i2()->at(0)),
                    b.mu_eta()->at(b.ll_i2()->at(0)),
                    b.mu_phi()->at(b.ll_i2()->at(0)),
                    0.106);
    return (m1+m2).M();
  });

  //returns llgamma mass calculated with corrected muon pts
  const NamedFunc llphoton_corrected_m("llphoton_corrected_m",
      [](const Baby &b) -> NamedFunc::ScalarType{
    //only recalculate for Muons
    if (b.ll_lepid()->at(0)==11) return b.llphoton_m()->at(0);
    TLorentzVector m1, m2, ph;
    m1.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i1()->at(0)),
                    b.mu_eta()->at(b.ll_i1()->at(0)),
                    b.mu_phi()->at(b.ll_i1()->at(0)),
                    0.106);
    m2.SetPtEtaPhiM(b.mu_corrected_pt()->at(b.ll_i2()->at(0)),
                    b.mu_eta()->at(b.ll_i2()->at(0)),
                    b.mu_phi()->at(b.ll_i2()->at(0)),
                    0.106);
    ph.SetPtEtaPhiM(b.photon_pt()->at(0),
                    b.photon_eta()->at(0),
                    b.photon_phi()->at(0),
                    0.0);
    return (m1+m2+ph).M();
  });

  //Z pt reweighting for DY amcatnlo sample
  const NamedFunc w_z_pt("w_z_pt",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    float z_pt = 0.0;
    if (b.type() >= 6000 && b.type() < 7000 && b.nll()>0)
      z_pt = b.ll_pt()->at(0);
    else if (b.type() >= 17000 && b.type() < 18000 && b.nllphoton()>0)
      z_pt = b.llphoton_pt()->at(0);
    else
      return 1;
    if (z_pt<2.0)                      {return 0.751271;}
    else if (z_pt>2.0 && z_pt<4.0)     {return 0.793257;}
    else if (z_pt>4.0 && z_pt<6.0)     {return 0.891502;}
    else if (z_pt>6.0 && z_pt<8.0)     {return 0.997759;}
    else if (z_pt>8.0 && z_pt<10.0)    {return 1.06628;}
    else if (z_pt>10.0 && z_pt<12.0)   {return 1.09597;}
    else if (z_pt>12.0 && z_pt<14.0)   {return 1.10345;}
    else if (z_pt>14.0 && z_pt<16.0)   {return 1.09802;}
    else if (z_pt>16.0 && z_pt<18.0)   {return 1.08284;}
    else if (z_pt>18.0 && z_pt<20.0)   {return 1.0655;}
    else if (z_pt>20.0 && z_pt<25.0)   {return 1.03228;}
    else if (z_pt>25.0 && z_pt<30.0)   {return 0.987943;}
    else if (z_pt>30.0 && z_pt<35.0)   {return 0.956658;}
    else if (z_pt>35.0 && z_pt<40.0)   {return 0.944815;}
    else if (z_pt>40.0 && z_pt<50.0)   {return 0.954598;}
    else if (z_pt>50.0 && z_pt<60.0)   {return 0.971563;}
    else if (z_pt>60.0 && z_pt<80.0)   {return 0.975438;}
    else if (z_pt>80.0 && z_pt<100.0)  {return 0.974402;}
    else if (z_pt>100.0 && z_pt<200.0) {return 0.966944;}
    return 0.947066;
  });

  const NamedFunc w_lumi_years("w_lumi*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_years("w_lumi*w_lep*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_trig_years("w_lumi*w_lep*w_trig*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_trig_pu_years("w_lumi*w_lep*w_trig*w_pu*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*b.w_pu()*w_years.GetScalar(b);
  });

  const NamedFunc w_lumi_lep_trig_pu_prefire_years("w_lumi*w_lep*w_trig*w_pu*w_prefire*w_years",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.w_lumi()*b.w_lep()*b.w_trig()*b.w_pu()*b.w_prefire()*w_years.GetScalar(b);
  });

  const NamedFunc weight_syslep_up("weight_syslep_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_lep()*b.sys_lep()->at(0)*w_years.GetScalar(b);
  });

  const NamedFunc weight_syslep_dn("weight_syslep_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_lep()*b.sys_lep()->at(1)*w_years.GetScalar(b);
  });

  const NamedFunc weight_systrig_up("weight_systrig_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_trig()*b.sys_trig()->at(0)*w_years.GetScalar(b);
  });

  const NamedFunc weight_systrig_dn("weight_systrig_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_trig()*b.sys_trig()->at(1)*w_years.GetScalar(b);
  });

  const NamedFunc weight_sysprefire_up("weight_sysprefire_up",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_prefire()*b.sys_prefire()->at(0)*w_years.GetScalar(b);
  });

  const NamedFunc weight_sysprefire_dn("weight_sysprefire_dn",
      [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()<1000) return 1.0;
    return b.weight()/b.w_prefire()*b.sys_prefire()->at(1)*w_years.GetScalar(b);
  });
  
  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------
  
  string dijet_balance_str = "|#sum_{ll,#gamma,j1,j2} #vec{p}_{T}|/(#sum_{ll,#gamma,j1,j2} |#vec{p}_{T}|)";
  dijet_balance_str = "Dijet balance (vector p_{T} sum over scalar p_{T} sum)";

  PlotMaker pm;
  pm.multithreaded_ = true;
  pm.min_print_ = true;

  if (false) {
    vector<double> ll_reweight_bins = {0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,80.0,100.0,200.0,300.0};
    pm.Push<Hist1D>(Axis(ll_reweight_bins, "ll_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops_root).Weight("weight"*w_years).Tag("zgvalidate_"+production+"_"+years);
  }

  //for (unsigned iwgt = 0; iwgt < 7; iwgt++) {
  if (true) {
    unsigned iwgt = 0;

    //set weight
    NamedFunc weight = "weight"*w_years*w_z_pt;
    if (iwgt==1) weight = w_lumi_years;
    if (iwgt==2) weight = w_lumi_lep_years;
    if (iwgt==3) weight = w_lumi_lep_trig_years;
    if (iwgt==4) weight = "weight"*w_years*w_z_pt;
    //if (iwgt==4) weight = w_lumi_lep_trig_pu_prefire_years;
    //if (iwgt==6) weight = w_lumi_lep_trig_pu_years;
    //if (iwgt==7) weight = weight_sysprefire_up*w_years*w_z_pt;
    //if (iwgt==8) weight = weight_sysprefire_dn*w_years*w_z_pt;
    //if (iwgt==7) weight = weight_syslep_up*w_years*w_z_pt;
    //if (iwgt==8) weight = weight_syslep_dn*w_years*w_z_pt;
    //if (iwgt==9) weight = weight_systrig_up*w_years*w_z_pt;
    //if (iwgt==10) weight = weight_systrig_dn*w_years*w_z_pt;

    //Z->ee plots (10 total)
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "ll_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-6.0,6.0, "ll_eta[0]", "Z candidate #eta", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "ll_phi[0]", "Z candidate #phi", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,80,100, "ll_m[0]", "Z candidate m [GeV]", {}), 
        zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[0]", "Lead electron I_{rel}", {}), 
    //    zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[1]", "Sublead electron I_{rel}", {}), 
    //    zee_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[0]", "Lead electron d_{z} [cm]", {}), 
    //    zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[1]", "Sublead electron d_{z} [cm]", {}), 
    //    zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[0]", "Lead electron d_{xy} [cm]", {}), 
    //    zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[1]", "Sublead electron d_{xy} [cm]", {}), 
    //    zee_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
    //    zee_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
    //    zee_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
    //    zee_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
    //    zee_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");

    //Z->mumu plots (10 total)
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "ll_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-6.0,6.0, "ll_eta[0]", "Z candidate #eta", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "ll_phi[0]", "Z candidate #phi", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,80,100, "ll_m[0]", "Z candidate m [GeV]", {}), 
        zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[0]", "Lead muon corrected p_{T} [GeV]", {}), 
    //    zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[1]", "Sublead muon corrected p_{T} [GeV]", {}), 
    //    zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[0]", "Lead muon I_{rel}", {}), 
    //    zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[1]", "Sublead muon I_{rel}", {}), 
    //    zmm_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[0]", "Lead muon d_{z} [cm]", {}), 
    //    zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[1]", "Sublead muon d_{z} [cm]", {}), 
    //    zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[0]", "Lead muon d_{xy} [cm]", {}), 
    //    zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[1]", "Sublead muon d_{xy} [cm]", {}), 
    //    zmm_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
    //    zmm_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
    //    zmm_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
    //    zmm_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
    //    zmm_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist1D>(Axis(25,71,111, ll_corrected_m, "Z candidate corrected m [GeV]", {}), 
    //    "nmu==2&&nphoton==0", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z->eey plots (13+12=25 total)
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[0]", "Lead electron p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "el_pt[1]", "Sublead electron p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,80.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "llphoton_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-6.0,6.0, "llphoton_eta[0]", "Z candidate #eta", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "llphoton_phi[0]", "Z candidate #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,80,100, "llphoton_m[0]", "Z candidate m [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //advanced
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.2, photon_relpterr, "Photon #sigma_{E}/E", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,4.0, "photon_drmin[0]", "Photon #Delta R_{min}(l)", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,5.0, "photon_drmax[0]", "Photon #Delta R_{max}(l)", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        "met>50"&&zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfm", "N_{b medium}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "njet", "N_{jet}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(3,1.5,4.5, "nlep", "N_{lep}", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "llphoton_costheta[0]", "ll#gamma cos(#theta)", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, llphoton_coscaptheta0, "ll#gamma cos(#Theta)", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "llphoton_psi[0]", "ll#gamma #phi", {}), 
        zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[0]", "Lead electron I_{rel}", {}), 
    //    zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "el_reliso[1]", "Sublead electron I_{rel}", {}), 
    //    zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[0]", "Lead electron d_{z} [cm]", {}), 
    //    zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "el_dz[1]", "Sublead electron d_{z} [cm]", {}), 
    //    zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[0]", "Lead electron d_{xy} [cm]", {}), 
    //    zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "el_dxy[1]", "Sublead electron d_{xy} [cm]", {}), 
    //    zeey_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
    //    zeey_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
    //    zeey_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[0]", "Lead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[0]", "Lead electron #phi", {}), 
    //    zeey_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "el_eta[1]", "Sublead electron #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "el_phi[1]", "Sublead electron #phi", {}), 
    //    zeey_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
    //    zeey_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
    //    zeey_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "jet_pt", "Jet p_{T} [GeV]", {}), 
    //    zeey_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
    //    zeey_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
    //    zeey_ctrlregion&&"njet>0", procs, ops).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
    //    zeey_ctrlregion&&"jet_isgood", procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
    //    zeey_ctrlregion&&"jet_isgood", procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfl", "N_{b loose}", {}), 
    //    zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdft", "N_{b tight}", {}), 
    //    zeey_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
    //    zeey_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z->mumuy plots (13+12=25 total)
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[0]", "Lead muon p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_pt[1]", "Sublead muon p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,80.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "llphoton_pt[0]", "Z candidate p_{T} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-6.0,6.0, "llphoton_eta[0]", "Z candidate #eta", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "llphoton_phi[0]", "Z candidate #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,80,100, "llphoton_m[0]", "Z candidate m [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //advanced
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,0.2, photon_relpterr, "Photon #sigma_{E}/E", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,4.0, "photon_drmin[0]", "Photon #Delta R_{min}(l)", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,5.0, "photon_drmax[0]", "Photon #Delta R_{max}(l)", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
        "met>50"&&zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "njet", "N_{jet}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfm", "N_{b medium}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(3,1.5,4.5, "nlep", "N_{lep}", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, "llphoton_costheta[0]", "ll#gamma cos(#theta)", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-1.0,1.0, llphoton_coscaptheta0, "ll#gamma cos(#Theta)", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "llphoton_psi[0]", "ll#gamma #phi", {}), 
        zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[0]", "Lead muon corrected p_{T} [GeV]", {}), 
    //    zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "mu_corrected_pt[1]", "Sublead muon corrected p_{T} [GeV]", {}), 
    //    zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[0]", "Lead muon I_{rel}", {}), 
    //    zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,0.5, "mu_reliso[1]", "Sublead muon I_{rel}", {}), 
    //    zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[0]", "Lead muon d_{z} [cm]", {}), 
    //    zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-1.0,1.0, "mu_dz[1]", "Sublead muon d_{z} [cm]", {}), 
    //    zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[0]", "Lead muon d_{xy} [cm]", {}), 
    //    zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-0.5,0.5, "mu_dxy[1]", "Sublead muon d_{xy} [cm]", {}), 
    //    zmmy_ctrlregion, procs, ops_log).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
    //    zmmy_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
    //    zmmy_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[0]", "Lead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[0]", "Lead muon #phi", {}), 
    //    zmmy_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "mu_eta[1]", "Sublead muon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "mu_phi[1]", "Sublead muon #phi", {}), 
    //    zmmy_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
    //    zmmy_ctrlregion, procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-2.5,2.5, "photon_eta[0]", "Photon #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "photon_phi[0]", "Photon #phi", {}), 
    //    zmmy_ctrlregion, procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "jet_pt", "Jet p_{T} [GeV]", {}), 
    //    zmmy_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
    //    zmmy_ctrlregion&&"jet_isgood", procs, ops).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
    //    zmmy_ctrlregion&&"njet>0", procs, ops).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
    //    zmmy_ctrlregion&&"jet_isgood", procs_mc, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_mc");
    //pm.Push<Hist2D>(Axis(25,-5.0,5.0, "jet_eta", "Jet #eta", {}), 
    //    Axis(25,-3.1416,3.1416, "jet_phi", "Jet #phi", {}), 
    //    zmmy_ctrlregion&&"jet_isgood", procs_data, ops2d).Weight(weight)
    //    .Tag("zgvalidate_"+production+"_"+years+"_data");
    //pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdfl", "N_{b loose}", {}), 
    //    zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nbdft", "N_{b tight}", {}), 
    //    zmmy_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,0.0,100.0, "met", "p_{T}^{miss} [GeV]", {}), 
    //    zmmy_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, "met_phi", "p_{T}^{miss} #phi", {}), 
    //    zmmy_ctrlregion&&"pass"&&pass_ra2b_filters, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    //pm.Push<Hist1D>(Axis(25,71,111, llphoton_corrected_m, "Z candidate corrected m [GeV]", {}), 
    //    "nmu==2&&nphoton==1", procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z to eegamma + dijet CR (15 plots)
    pm.Push<Hist1D>(Axis(25,30.0,200.0, lead_jet_pt, "Lead jet p_{T} [GeV]", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,30.0,100.0, sublead_jet_pt, "Sublead jet p_{T} [GeV]", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, lead_jet_eta, "Lead jet #eta", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, sublead_jet_eta, "Sublead jet #eta", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, lead_jet_phi, "Lead jet #phi", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, sublead_jet_phi, "Sublead jet #phi", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,50.0, lead_jet_m, "Lead jet mass [GeV]", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,30.0, sublead_jet_m, "Sublead jet mass [GeV]", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,3.1416, "dijet_dphi", "Dijet #Delta #phi", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,8.0, "dijet_deta", "Dijet #Delta #eta", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,6.0, "photon_jet1_dr[0]", "#Delta R(#gamma, jet 1)", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,6.0, "photon_jet2_dr[0]", "#Delta R(#gamma, jet 2)", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,3.1416, "llphoton_dijet_dphi[0]", "#Delta #phi(ll#gamma, jj)", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,5.0, "photon_zeppenfeld[0]", "|#eta_{#gamma}-(#eta_{j1}+#eta_{j2})/2|", {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,1.0, "llphoton_dijet_balance[0]", dijet_balance_str, {}), 
        zeeyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);

    //Z to mumugamma + dijet CR (15 plots)
    pm.Push<Hist1D>(Axis(25,30.0,200.0, lead_jet_pt, "Lead jet p_{T} [GeV]", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,30.0,100.0, sublead_jet_pt, "Sublead jet p_{T} [GeV]", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, lead_jet_eta, "Lead jet #eta", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-5.0,5.0, sublead_jet_eta, "Sublead jet #eta", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, lead_jet_phi, "Lead jet #phi", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,-3.1416,3.1416, sublead_jet_phi, "Sublead jet #phi", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,50.0, lead_jet_m, "Lead jet mass [GeV]", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,30.0, sublead_jet_m, "Sublead jet mass [GeV]", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,3.1416, "dijet_dphi", "Dijet #Delta #phi", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,8.0, "dijet_deta", "Dijet #Delta #eta", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,6.0, "photon_jet1_dr[0]", "#Delta R(#gamma, jet 1)", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,6.0, "photon_jet2_dr[0]", "#Delta R(#gamma, jet 2)", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,3.1416, "llphoton_dijet_dphi[0]", "#Delta #phi(ll#gamma, jj)", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,5.0, "photon_zeppenfeld[0]", "|#eta_{#gamma}-(#eta_{j1}+#eta_{j2})/2|", {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
    pm.Push<Hist1D>(Axis(25,0.0,1.0, "llphoton_dijet_balance[0]", dijet_balance_str, {}), 
        zmmyjj_ctrlregion, procs, ops).Weight(weight).Tag("zgvalidate_"+production+"_"+years);
  }

  pm.SetLuminosityTag(lumi_string).MakePlots(1.0,"validation");

  return 0;
}
