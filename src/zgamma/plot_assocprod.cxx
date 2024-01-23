#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TChain.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "TLorentzVector.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "zgamma/nano_functions.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace NamedFuncUtilities;
using namespace ZgFunctions;
using namespace ZgUtilities;

int main() {

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------

  //setup
  gErrorIgnoreLevel = 6000;

  std::vector<std::shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMC");
  std::vector<std::shared_ptr<Process>> procs_tth = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMC");
  std::vector<std::shared_ptr<Process>> procs_vh = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMC");
  std::vector<std::shared_ptr<Process>> procs_vbf = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMC");
  std::vector<std::shared_ptr<Process>> proc_zg_tth = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","ttHtoZG");
  std::vector<std::shared_ptr<Process>> proc_zg_vh = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","VHtoZG");
  std::vector<std::shared_ptr<Process>> proc_zg_vbf = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","VBFHtoZG");
  std::vector<std::shared_ptr<Process>> procs_split = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllSplit");
  procs_tth.insert(procs_tth.end(),proc_zg_tth.begin(),proc_zg_tth.end());
  procs_vh.insert(procs_vh.end(),proc_zg_vh.begin(),proc_zg_vh.end());
  procs_vbf.insert(procs_vbf.end(),proc_zg_vbf.begin(),proc_zg_vbf.end());

  std::vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumi").Overflow(OverflowType::none)}; 
  std::vector<PlotOpt> ops_sorb = {PlotOpt("txt/plot_styles.txt","LinSorb"),PlotOpt("txt/plot_styles.txt","LinSorbUpper"),PlotOpt("txt/plot_styles.txt","LinLumi"),PlotOpt("txt/plot_styles.txt","Shapes")}; 
  std::vector<PlotOpt> ops_log = {PlotOpt("txt/plot_styles.txt","LogLumi")}; 
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","Shapes"), PlotOpt("txt/plot_styles.txt","LinLumi")}; 
  std::vector<PlotOpt> ops_2d = {PlotOpt("txt/plot_styles.txt","LogLumi2D")}; 

  //------------------------------------------------------------------------------------
  //                 Temporary processes until fatjets added to pico
  //------------------------------------------------------------------------------------
  
  std::string nano_folder = "/net/cms17/cms17r0/pico/NanoAODv9/nano/2017/mc/"; 
  Palette colors("txt/colors_zgamma.txt","default");
  Process::Type back =  Process::Type::background;
  //Process::Type data =  Process::Type::data;
  Process::Type sig =  Process::Type::signal;

  //DY, ZG, tt, ttG, VBS ZG, ttX, multiboson, signal
  //note: can't use ZZG since only leptonic decays
  std::vector<std::string> file_names = {"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__2430000*",
                                         "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__2*",
                                         "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v1__2510000__3*",
                                         "*TTGJets*",
                                         "EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__3*",
                                         "ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__1*",
                                         "ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8__RunIISummer20UL17NanoAODv9__106X_mc2017_realistic_v9-v2__2520000*",
                                         "WW_Tune*","WZ_Tune*","ZZ_Tune*"};
  std::vector<std::string> fname_nowc = {"DYJetsToLL","ZGToLLG_01J_5f_Tune",
                                         "TTTo2L2Nu","TTGJets","EWKZ2Jets",
                                         "ttZJets","ttWJets_Tune","WW_Tune",
                                         "WZ_Tune","ZZ_Tune"};
  std::vector<std::string> samp_names = {"Z+Fake Photon","Z+#gamma","t#bar{t}",
                                         "t#bar{t}#gamma","VBS Z+#gamma","t#bar{t}Z",
                                         "t#bar{t}W","WW","WZ","ZZ"};
  std::vector<std::string> samp_color = {"dyjets","zgtollg","tt","ttg","vbszg",
                                         "ttx","ttx","vv","vv","vv"};
  std::vector<float>       sample_xss = {6077.22,55.48,87.339,0.342,6.215,0.5407,
                                         0.4611,118.7,51.11,16.91};
  std::vector<int>         samp_sizes = {};

  std::cout << "Adding NanoAOD samples and counting entries\n";
  std::vector<std::shared_ptr<Process>> procs_nano;
  TChain * entries_chain;
  for (unsigned isample = 0; isample < file_names.size(); isample++) {
    procs_nano.push_back(Process::MakeShared<Baby_nano>(samp_names[isample], back, 
                         colors(samp_color[isample]),{nano_folder+file_names[isample]},
                         "1"));
                         //NanoFunctions::HLT_pass_dilepton&&NanoFunctions::stitch));
    entries_chain = new TChain("Events");
    entries_chain->Add((nano_folder+file_names[isample]).c_str());
    samp_sizes.push_back(entries_chain->GetEntries());
    std::cout << fname_nowc[isample] << " has " << samp_sizes.back() << " entries\n";
    delete entries_chain;
  }

  procs_nano.push_back(Process::MakeShared<Baby_nano>("associated H#rightarrow Z#gamma", sig,
                       kRed,{nano_folder+"WplusH_HToZG*M-125*",nano_folder+"WminusH_HToZG*M-125*",
                       nano_folder+"ZH_HToZG*M-125*",nano_folder+"*HToZG_M125*"},
                       NanoFunctions::HLT_pass_dilepton&&NanoFunctions::stitch));

  const NamedFunc weight_nano("weight_nano",[fname_nowc,samp_sizes,sample_xss](const Baby &b) -> NamedFunc::ScalarType{
    float nevs = 1;
    float xsec = 0;
    for (unsigned isample = 0; isample < fname_nowc.size(); isample++) {
      if (b.FirstFileName().find(fname_nowc[isample]) != std::string::npos) {
        nevs = samp_sizes[isample];
        xsec = sample_xss[isample];
      }
    }
    if (b.FirstFileName().find("HToZG") != std::string::npos) {
      //determine signal production mode from truth
      bool found_top = false;
      bool found_wplus = false;
      bool found_wminus = false;
      bool found_z = false;
      int vbf_jet = 0;
      for (unsigned imc = 0; imc < b.GenPart_pdgId()->size(); imc++) {
        if (b.GenPart_pdgId()->at(imc)==6) found_top=true;
        if (b.GenPart_pdgId()->at(imc)==-6) found_top=true;
        if (b.GenPart_pdgId()->at(imc)==24) found_wplus=true;
        if (b.GenPart_pdgId()->at(imc)==-24) found_wminus=true;
        if (b.GenPart_pdgId()->at(imc)==23) found_z=true;
        //note this VBF selection criteria accidentally selects O(1%) of ggF events
        if (b.GenPart_pdgId()->at(imc)>=5 &&
            b.GenPart_pdgId()->at(imc)<=5 &&
            (b.GenPart_statusFlags()->at(imc) & 0x1)==0x1 &&
            ((b.GenPart_statusFlags()->at(imc) & 0x1000)!=0 || imc<2)) vbf_jet++;
      }
      if (found_top)         {xsec = 0.00077738; nevs = 200000;}
      else if (found_wplus)  {xsec = 0.0012739;  nevs = 299983;}
      else if (found_wminus) {xsec = 0.00080789; nevs = 299988;}
      else if (found_z)      {xsec = 0.001355;   nevs = 279885;}
      else if (vbf_jet>=2)   {xsec = 0.00058543; nevs = 200000;}
      else                   {xsec = 0.0074997;  nevs = 400000;}
    }
    //todo add signal, do the same TODO
    float w_year = 0;
    if (b.SampleType()==2016) w_year = 36.32264; 
    else if (b.SampleType()==2017) w_year = 41.52756;
    else if (b.SampleType()==2018) w_year = 59.67377;
    return w_year*xsec*1000.0/nevs*340.0/41.52756; //currently only using 2017
  });
  
  const NamedFunc tighter_baseline_nano = NamedFunc(NanoFunctions::zg_baseline&&
      NanoFunctions::Lead_SignalPhoton_mvaID>0.5&&NanoFunctions::ZCandidate_mass>80&&
      NanoFunctions::ZCandidate_mass<100).Name("new_baseline");

  //lead fatjet pt
  const NamedFunc Lead_FatJet_pt = ReduceNamedFunc("FatJet_pt", reduce_max)
      .Name("Lead_FatJet_pt");

  //lead w tag score
  const NamedFunc Lead_FatJet_deepTag_WvsQCD = MultiReduceNamedFunc(
      {"FatJet_pt","FatJet_deepTag_WvsQCD"}, reduce_maxfirst)
      .Name("Lead_FatJet_deepTag_WvsQCD");

  //lead jet softdrop mass
  const NamedFunc Lead_FatJet_msoftdrop = MultiReduceNamedFunc({"FatJet_pt","FatJet_msoftdrop"}, 
      reduce_maxfirst).Name("Lead_FatJet_msoftdrop");

  //max jet softdrop mass
  const NamedFunc Max_FatJet_msoftdrop = ReduceNamedFunc("FatJet_msoftdrop",reduce_max)
      .Name("Max_FatJet_msoftdrop");

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------

  //increase signal weight by 20x
  NamedFunc w_sigx20("w_sigx20",[](const Baby &b) -> NamedFunc::ScalarType{
    double w_lumi = 1.0;
    if(b.type() >= 200000 && b.type() <= 205000)
      w_lumi *= 20.0;
    return w_lumi;
  });

  //associated lepton pdg id
  const NamedFunc assoc_lep_pdgid("assoc_lep_pdgid",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_pt = 999;
    bool z_is_el = false;
    float pdgid = 0;
    unsigned int z_lep_i1 = static_cast<unsigned int>(b.ll_i1()->at(0));
    unsigned int z_lep_i2 = static_cast<unsigned int>(b.ll_i2()->at(0));
    if (b.ll_lepid()->at(0)==11)
      z_is_el = true;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
      if (b.el_pt()->at(iel)<min_pt) {
        min_pt = b.el_pt()->at(iel);
        pdgid = 11;
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
      if (b.mu_pt()->at(imu)<min_pt) {
        min_pt = b.mu_pt()->at(imu);
        pdgid = 13;
      }
    }
    return pdgid;
  });

  //associated lepton ID quality: 0 = no id, 1 = loose id, 2 = medium id, 3 = tight id
  const NamedFunc assoc_lep_id("assoc_lep_id",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_pt = 999;
    bool z_is_el = false;
    float quality = 0;
    unsigned int z_lep_i1 = static_cast<unsigned int>(b.ll_i1()->at(0));
    unsigned int z_lep_i2 = static_cast<unsigned int>(b.ll_i2()->at(0));
    if (b.ll_lepid()->at(0)==11)
      z_is_el = true;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
      if (b.el_pt()->at(iel)<min_pt) {
        min_pt = b.el_pt()->at(iel);
        if (b.el_id80()->at(iel)) quality = 3;
        else if (b.el_id90()->at(iel)) quality = 2;
        else if (b.el_idLoose()->at(iel)) quality = 1;
        else quality = 0;
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
      if (b.mu_pt()->at(imu)<min_pt) {
        min_pt = b.mu_pt()->at(imu);
        if (b.mu_tightid()->at(imu)) quality = 3;
        else if (b.mu_mediumid()->at(imu) 
            || (b.mu_highptid()->at(imu) && b.mu_pt()->at(imu)>200)) quality = 2;
        else if (b.mu_id()->at(imu)) quality = 1;
        else quality = 0;
      }
    }
    return quality;
  });

  //associated lepton mini Iso
  const NamedFunc assoc_lep_miniso("assoc_lep_miniso",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_pt = 999;
    bool z_is_el = false;
    float miniso = 0;
    unsigned int z_lep_i1 = static_cast<unsigned int>(b.ll_i1()->at(0));
    unsigned int z_lep_i2 = static_cast<unsigned int>(b.ll_i2()->at(0));
    if (b.ll_lepid()->at(0)==11)
      z_is_el = true;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
      if (b.el_pt()->at(iel)<min_pt) {
        min_pt = b.el_pt()->at(iel);
        miniso = b.el_miniso()->at(iel);
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
      if (b.mu_pt()->at(imu)<min_pt) {
        min_pt = b.mu_pt()->at(imu);
        miniso = b.mu_miniso()->at(imu);
      }
    }
    return miniso;
  });

  //associated lepton pt
  const NamedFunc assoc_lep_pt("assoc_lep_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_pt = 999;
    bool z_is_el = false;
    unsigned int z_lep_i1 = static_cast<unsigned int>(b.ll_i1()->at(0));
    unsigned int z_lep_i2 = static_cast<unsigned int>(b.ll_i2()->at(0));
    if (b.ll_lepid()->at(0)==11)
      z_is_el = true;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
      if (b.el_pt()->at(iel)<min_pt) {
        min_pt = b.el_pt()->at(iel);
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
      if (b.mu_pt()->at(imu)<min_pt) {
        min_pt = b.mu_pt()->at(imu);
      }
    }
    return min_pt;
  });

  //associated lepton mt
  const NamedFunc assoc_lep_mt("assoc_lep_mt",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nlep()<3) return 0;
    float max_pt = 0;
    float max_phi = 0;
    bool z_is_el = false;
    unsigned int z_lep_i1 = static_cast<unsigned int>(b.ll_i1()->at(0));
    unsigned int z_lep_i2 = static_cast<unsigned int>(b.ll_i2()->at(0));
    if (b.ll_lepid()->at(0)==11)
      z_is_el = true;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (z_is_el && (iel==z_lep_i1 || iel==z_lep_i2)) continue;
      if (b.el_sig()->at(iel) && b.el_pt()->at(iel)>max_pt) {
        max_pt = b.el_pt()->at(iel);
        max_phi = b.el_phi()->at(iel);
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (!z_is_el && (imu==z_lep_i1 || imu==z_lep_i2)) continue;
      if (b.mu_sig()->at(imu) && b.mu_pt()->at(imu)>max_pt) {
        max_pt = b.mu_pt()->at(imu);
        max_phi = b.mu_phi()->at(imu);
      }
    }
    return sqrt(2.0*b.met()*max_pt*(1.0-cos(b.met_phi()-max_phi)));
  });

  //min lepton pt
  const NamedFunc min_lep_pt("min_lep_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_pt = 999;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel) && b.el_pt()->at(iel)<min_pt) {
        min_pt = b.el_pt()->at(iel);
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu) && b.mu_pt()->at(imu)<min_pt) {
        min_pt = b.mu_pt()->at(imu);
      }
    }
    return min_pt;
  });

  //min lepton ID quality: 0 = no id, 1 = loose id, 2 = medium id, 3 = tight id
  const NamedFunc min_lep_id("min_lep_id",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_quality = 3;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel)) {
        float quality = 0;
        if (b.el_id80()->at(iel)) quality = 3;
        else if (b.el_id90()->at(iel)) quality = 2;
        else if (b.el_idLoose()->at(iel)) quality = 1;
        if (quality < min_quality)
          min_quality = quality;
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu)) {
        float quality = 0;
        if (b.mu_tightid()->at(imu)) quality = 3;
        else if (b.mu_mediumid()->at(imu) 
            || (b.mu_highptid()->at(imu) && b.mu_pt()->at(imu)>200)) quality = 2;
        else if (b.mu_id()->at(imu)) quality = 1;
        if (quality < min_quality)
          min_quality = quality;
      }
    }
    return min_quality;
  });

  //max lepton miniso
  const NamedFunc max_lep_miniso("max_lep_miniso",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_miniso = 0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (b.el_sig()->at(iel) && b.el_miniso()->at(iel)>max_miniso) {
        max_miniso = b.el_miniso()->at(iel);
      }
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (b.mu_sig()->at(imu) && b.mu_miniso()->at(imu)>max_miniso) {
        max_miniso = b.mu_miniso()->at(imu);
      }
    }
    return max_miniso;
  });

  //top mass estimate
  const NamedFunc closest_top_mass("closest_top_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    //first find 2 highest b-tag score jets
    float lead_btag_score = -1;
    float subl_btag_score = -1;
    unsigned bidx1 = -1;
    unsigned bidx2 = -1;
    TLorentzVector j1, j2, w;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (!b.jet_isgood()->at(ijet)) continue;
      float this_btag_score = b.jet_deepflav()->at(ijet);
      if (this_btag_score > lead_btag_score) {
        subl_btag_score = lead_btag_score;
        lead_btag_score = this_btag_score;
        bidx2 = bidx1;
        bidx1 = ijet;
      }
      else if (this_btag_score > subl_btag_score) {
        subl_btag_score = this_btag_score;
        bidx2 = ijet;
      }
    }
    //from other jets, find pair closest to W
    float closest_w = -999;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (!b.jet_isgood()->at(ijet1)) continue;
      if (ijet1==bidx1 || ijet1==bidx2) continue;
      j1.SetPtEtaPhiM(b.jet_pt()->at(ijet1), b.jet_eta()->at(ijet1), 
                      b.jet_phi()->at(ijet1), b.jet_m()->at(ijet1));
      for (unsigned ijet2 = ijet1+1; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (!b.jet_isgood()->at(ijet2)) continue;
        if (ijet2==bidx1 || ijet2==bidx2) continue;
        j2.SetPtEtaPhiM(b.jet_pt()->at(ijet2), b.jet_eta()->at(ijet2), 
                        b.jet_phi()->at(ijet2), b.jet_m()->at(ijet2));
        float this_mass = (j1+j2).M();
        if (fabs(this_mass-80.4) < fabs(closest_w-80.4)) {
          closest_w = this_mass;
          w = j1+j2;
        }
      }
    }
    //refit W??
    //combine W with one of the b-tagged jets and take value closest to top
    j1.SetPtEtaPhiM(b.jet_pt()->at(bidx1), b.jet_eta()->at(bidx1), 
                    b.jet_phi()->at(bidx1), b.jet_m()->at(bidx1));
    j2.SetPtEtaPhiM(b.jet_pt()->at(bidx2), b.jet_eta()->at(bidx2), 
                    b.jet_phi()->at(bidx2), b.jet_m()->at(bidx2));
    float top1_m = (j1+w).M();
    float top2_m = (j2+w).M();
    if (fabs(top1_m-172.7)<fabs(top2_m-172.7)) return top1_m;
    return top2_m;
  });

  //dr top mass estimate
  const NamedFunc dr_top_mass("dr_top_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    //first find 2 highest b-tag score jets
    float lead_btag_score = -1;
    float subl_btag_score = -1;
    unsigned bidx1 = -1;
    unsigned bidx2 = -1;
    TLorentzVector j1, j2, w;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (!b.jet_isgood()->at(ijet)) continue;
      float this_btag_score = b.jet_deepflav()->at(ijet);
      if (this_btag_score > lead_btag_score) {
        subl_btag_score = lead_btag_score;
        lead_btag_score = this_btag_score;
        bidx2 = bidx1;
        bidx1 = ijet;
      }
      else if (this_btag_score > subl_btag_score) {
        subl_btag_score = this_btag_score;
        bidx2 = ijet;
      }
    }
    //from other jets, find pair closest to W
    float closest_w = -999;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (!b.jet_isgood()->at(ijet1)) continue;
      if (ijet1==bidx1 || ijet1==bidx2) continue;
      j1.SetPtEtaPhiM(b.jet_pt()->at(ijet1), b.jet_eta()->at(ijet1), 
                      b.jet_phi()->at(ijet1), b.jet_m()->at(ijet1));
      for (unsigned ijet2 = ijet1+1; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (!b.jet_isgood()->at(ijet2)) continue;
        if (ijet2==bidx1 || ijet2==bidx2) continue;
        j2.SetPtEtaPhiM(b.jet_pt()->at(ijet2), b.jet_eta()->at(ijet2), 
                        b.jet_phi()->at(ijet2), b.jet_m()->at(ijet2));
        float this_mass = (j1+j2).M();
        if (fabs(this_mass-80.4) < fabs(closest_w-80.4)) {
          closest_w = this_mass;
          w = j1+j2;
        }
      }
    }
    //combine W with one of the b-tagged jets and take value closest to top
    j1.SetPtEtaPhiM(b.jet_pt()->at(bidx1), b.jet_eta()->at(bidx1), 
                    b.jet_phi()->at(bidx1), b.jet_m()->at(bidx1));
    j2.SetPtEtaPhiM(b.jet_pt()->at(bidx2), b.jet_eta()->at(bidx2), 
                    b.jet_phi()->at(bidx2), b.jet_m()->at(bidx2));
    if (w.DeltaR(j1)<w.DeltaR(j2))
      return (j1+w).M();
    return (j2+w).M();
  });

  //top mass estimate
  const NamedFunc closest_top_mass_constrained("closest_top_mass_constrained",[](const Baby &b) -> NamedFunc::ScalarType{
    //first find 2 highest b-tag score jets
    float lead_btag_score = -1;
    float subl_btag_score = -1;
    unsigned bidx1 = -1;
    unsigned bidx2 = -1;
    TLorentzVector j1, j2, w;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (!b.jet_isgood()->at(ijet)) continue;
      float this_btag_score = b.jet_deepflav()->at(ijet);
      if (this_btag_score > lead_btag_score) {
        subl_btag_score = lead_btag_score;
        lead_btag_score = this_btag_score;
        bidx2 = bidx1;
        bidx1 = ijet;
      }
      else if (this_btag_score > subl_btag_score) {
        subl_btag_score = this_btag_score;
        bidx2 = ijet;
      }
    }
    //from other jets, find pair closest to W
    float closest_w = -999;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (!b.jet_isgood()->at(ijet1)) continue;
      if (ijet1==bidx1 || ijet1==bidx2) continue;
      j1.SetPtEtaPhiM(b.jet_pt()->at(ijet1), b.jet_eta()->at(ijet1), 
                      b.jet_phi()->at(ijet1), b.jet_m()->at(ijet1));
      for (unsigned ijet2 = ijet1+1; ijet2 < b.jet_pt()->size(); ijet2++) {
        if (!b.jet_isgood()->at(ijet2)) continue;
        if (ijet2==bidx1 || ijet2==bidx2) continue;
        j2.SetPtEtaPhiM(b.jet_pt()->at(ijet2), b.jet_eta()->at(ijet2), 
                        b.jet_phi()->at(ijet2), b.jet_m()->at(ijet2));
        float this_mass = (j1+j2).M();
        if (fabs(this_mass-80.4) < fabs(closest_w-80.4)) {
          closest_w = this_mass;
          w = j1+j2;
        }
      }
    }
    //refit W
    w.SetPtEtaPhiM(w.Pt(),w.Eta(),w.Phi(),80.4);
    //combine W with one of the b-tagged jets and take value closest to top
    j1.SetPtEtaPhiM(b.jet_pt()->at(bidx1), b.jet_eta()->at(bidx1), 
                    b.jet_phi()->at(bidx1), b.jet_m()->at(bidx1));
    j2.SetPtEtaPhiM(b.jet_pt()->at(bidx2), b.jet_eta()->at(bidx2), 
                    b.jet_phi()->at(bidx2), b.jet_m()->at(bidx2));
    float top1_m = (j1+w).M();
    float top2_m = (j2+w).M();
    if (fabs(top1_m-172.7)<fabs(top2_m-172.7)) return top1_m;
    return top2_m;
  });

  //full dr top mass estimate
  const NamedFunc full_dr_top_mass("full_dr_top_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    //first find 2 highest b-tag score jets
    float lead_btag_score = -1;
    unsigned bidx1 = -1;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (!b.jet_isgood()->at(ijet)) continue;
      float this_btag_score = b.jet_deepflav()->at(ijet);
      if (this_btag_score > lead_btag_score) {
        lead_btag_score = this_btag_score;
        bidx1 = ijet;
      }
    }
    //get closest 2 other jets
    float minm_dr = 999;
    float next_dr = 999;
    unsigned widx1 = -1;
    unsigned widx2 = -1;
    for (unsigned ijet1 = 0; ijet1 < b.jet_pt()->size(); ijet1++) {
      if (!b.jet_isgood()->at(ijet1)) continue;
      if (ijet1==bidx1) continue;
      float this_dr = deltaR(b.jet_eta()->at(bidx1), b.jet_phi()->at(bidx1),
                             b.jet_eta()->at(ijet1), b.jet_phi()->at(ijet1));
      if (this_dr < minm_dr) {
        next_dr = minm_dr;
        minm_dr = this_dr;
        widx2 = widx1;
        widx1 = ijet1;
      }
      else if (this_dr < next_dr) {
        next_dr = this_dr;
        widx2 = ijet1;
      }
    }
    TLorentzVector j1, j2, j3;
    j1.SetPtEtaPhiM(b.jet_pt()->at(widx1), b.jet_eta()->at(widx1), 
                    b.jet_phi()->at(widx1), b.jet_m()->at(widx1));
    j2.SetPtEtaPhiM(b.jet_pt()->at(widx2), b.jet_eta()->at(widx2), 
                    b.jet_phi()->at(widx2), b.jet_m()->at(widx2));
    j3.SetPtEtaPhiM(b.jet_pt()->at(bidx1), b.jet_eta()->at(bidx1), 
                    b.jet_phi()->at(bidx1), b.jet_m()->at(bidx1));
    return (j1+j2+j3).M();
  });

  //minimum dphi between a jet and ptmiss
  const NamedFunc min_dphi_jet_met("min_dphi_jet_met",[](const Baby &b) -> NamedFunc::ScalarType{
    float min_dphi = 3.1415;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (!b.jet_isgood()->at(ijet)) continue;
      float this_dphi = deltaPhi(b.jet_phi()->at(ijet),b.met_phi());
      if (this_dphi < min_dphi) min_dphi = this_dphi;
    }
    return min_dphi;
  });

  //maximum dphi between a jet and ptmiss
  const NamedFunc max_dphi_jet_met("max_dphi_jet_met",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_dphi = 0;
    for (unsigned ijet = 0; ijet < b.jet_pt()->size(); ijet++) {
      if (!b.jet_isgood()->at(ijet)) continue;
      float this_dphi = deltaPhi(b.jet_phi()->at(ijet),b.met_phi());
      if (this_dphi > max_dphi) max_dphi = this_dphi;
    }
    return max_dphi;
  });

  //Dphi between H and MET
  const NamedFunc dphi_h_met("dphi_h_met",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.llphoton_phi()->size()==0) return 0;
    return deltaPhi(b.llphoton_phi()->at(0),b.met_phi());
  });

  const NamedFunc tighter_baseline = NamedFunc(zg_baseline&&"photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100").Name("new_baseline");

  //fix weight for Zgamma sample
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

  const NamedFunc january_baseline = NamedFunc("(nel>=2||nmu>=2) && nphoton>=1 && photon_id80[0] && (ll_m[0]>81 && ll_m[0]<101) && (ll_charge[0]==0) && ((photon_pt[0]/llphoton_m[0])>15.0/110.0) && ((llphoton_m[0]+ll_m[0])>185)" && trig_plateau_cuts).Name("baseline");

  //fat jet signal branch including DR 0.8 separation from photons/leptons
  const NamedFunc fjet_sig("fjet_sig",[](const Baby &b) -> NamedFunc::VectorType{
    vector<double> isgood;
    for (unsigned ifjet = 0; ifjet < b.fjet_pt()->size(); ifjet++) {
      if (b.fjet_pt()->at(ifjet) < 200) {
        isgood.push_back(0);
        continue;
      }
      for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
        if (b.el_sig()->at(iel)) {
          if (deltaR(b.fjet_eta()->at(ifjet), b.fjet_phi()->at(ifjet),
                     b.el_eta()->at(iel), b.el_phi()->at(iel))<0.8) {
            isgood.push_back(0);
            continue;
          }
        }
      }
      for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
        if (b.mu_sig()->at(imu)) {
          if (deltaR(b.fjet_eta()->at(ifjet), b.fjet_phi()->at(ifjet),
                     b.mu_eta()->at(imu), b.mu_phi()->at(imu))<0.8) {
            isgood.push_back(0);
            continue;
          }
        }
      }
      for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
        if (b.photon_sig()->at(iph)) {
          if (deltaR(b.fjet_eta()->at(ifjet), b.fjet_phi()->at(ifjet),
                     b.photon_eta()->at(iph), b.photon_phi()->at(iph))<0.8) {
            isgood.push_back(0);
            continue;
          }
        }
      }
      isgood.push_back(1);
    }
    return isgood;
  });

  const NamedFunc sig_fjet_pt = FilterNamedFunc("fjet_pt",fjet_sig).Name("sig_fjet_pt");
  const NamedFunc sig_fjet_pnet_wtag = FilterNamedFunc("fjet_pnet_wtag",fjet_sig).Name("sig_fjet_pnet_wtag");
  const NamedFunc sig_fjet_pnet_ttag = FilterNamedFunc("fjet_pnet_ttag",fjet_sig).Name("sig_fjet_pnet_ttag");

  const NamedFunc lead_sig_fjet_pt = ReduceNamedFunc(sig_fjet_pt,reduce_max).Name("lead_sig_fjet_pt");
  const NamedFunc lead_sig_fjet_pnet_wtag = MultiReduceNamedFunc(
      {sig_fjet_pt,sig_fjet_pnet_wtag},reduce_maxfirst).Name("lead_sig_fjet_pnet_wtag");
  const NamedFunc lead_sig_fjet_pnet_ttag = MultiReduceNamedFunc(
      {sig_fjet_pt,sig_fjet_pnet_ttag},reduce_maxfirst).Name("lead_sig_fjet_pnet_ttag");

  const NamedFunc nsig_fjet = ReduceNamedFunc(fjet_sig,reduce_sum).Name("nsig_fjet");

  //NamedFuncs
  
  //NamedFunc mc_higgs_pt("mc_higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==25) {
  //      return b.mc_pt()->at(imc);
  //    }
  //  }
  //  return -999;
  //});

  //NamedFunc mc_photon_pt("mc_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
  //      return b.mc_pt()->at(imc);
  //    }
  //  }
  //  return -999;
  //});

  //NamedFunc mc_deltar_zg("mc_deltar_zg",[](const Baby &b) -> NamedFunc::ScalarType{
  //  TVector3 ph;
  //  TVector3 Z;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
  //      ph.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
  //    }
  //    if (b.mc_id()->at(imc)==23 && b.mc_mom()->at(imc)==25) {
  //      Z.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
  //    }
  //  }
  //  return ph.DeltaR(Z);
  //});

  //NamedFunc mc_llg_in_acceptance("mc_llg_in_acceptance",[](const Baby &b) -> NamedFunc::ScalarType{
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
  //      if (!(b.mc_pt()->at(imc)>15 && abs(b.mc_eta()->at(imc))<2.5))
  //        return 0;
  //    }
  //    if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
  //      int mom_idx = b.mc_momidx()->at(imc);
  //      if (mom_idx != -1) {
  //        if (b.mc_mom()->at(mom_idx)==25) {
  //          float eta_cut = 2.5;
  //          float pt_cut = 15;
  //          if (abs(b.mc_id()->at(imc))==13) {
  //            eta_cut = 2.4;
  //            pt_cut = 10;
  //          }
  //          if (!(b.mc_pt()->at(imc)>pt_cut && abs(b.mc_eta()->at(imc))<eta_cut))
  //            return 0;
  //        }
  //      }
  //    }
  //  }
  //  return 1;
  //});

  //NamedFunc mc_llg_reco("mc_llg_reco",[](const Baby &b) -> NamedFunc::ScalarType{
  //  TVector3 l1, l2, ph;
  //  TVector3 reco_obj;
  //  int z_decay_pdgid = 11;
  //  int lep_idx = 0;
  //  bool found_ph(false), found_l1(false), found_l2(false);
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
  //      ph.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
  //    }
  //    if ((abs(b.mc_id()->at(imc))==11 || abs(b.mc_id()->at(imc))==13)&&b.mc_mom()->at(imc)==23) {
  //      int mom_idx = b.mc_momidx()->at(imc);
  //      if (mom_idx != -1) {
  //        if (b.mc_mom()->at(mom_idx)==25) {
  //          z_decay_pdgid = abs(b.mc_id()->at(imc));
  //          if (lep_idx==0) {
  //            l1.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
  //            lep_idx = 1;
  //          }
  //          else {
  //            l2.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
  //          }
  //        }
  //      }
  //    }
  //  }
  //  for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
  //    reco_obj.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
  //    if (reco_obj.DeltaR(ph)<0.2) {
  //      found_ph = true;
  //    }
  //  }
  //  if (z_decay_pdgid==13) {
  //    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
  //      reco_obj.SetPtEtaPhi(b.mu_pt()->at(imu),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
  //      if (reco_obj.DeltaR(l1)<0.2) {
  //        found_l1 = true;
  //      }
  //      else if (reco_obj.DeltaR(l2)<0.2) {
  //        found_l2 = true;
  //      }
  //    }
  //  }
  //  else {
  //    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
  //      reco_obj.SetPtEtaPhi(b.el_pt()->at(iel),b.el_eta()->at(iel),b.el_phi()->at(iel));
  //      if (reco_obj.DeltaR(l1)<0.2) {
  //        found_l1 = true;
  //      }
  //      else if (reco_obj.DeltaR(l2)<0.2) {
  //        found_l2 = true;
  //      }
  //    }
  //  }
  //  if (found_ph&&found_l1&&found_l2) return 1;
  //  return 0;
  //});

  //NamedFunc sigphoton_mindr("sigphoton_mindr",[](const Baby &b) -> NamedFunc::ScalarType{
  //  for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
  //    if (b.photon_pt()->at(iph)<15) continue;
  //    if (abs(b.photon_eta()->at(iph))>2.5) continue;
  //    if (abs(b.photon_eta()->at(iph))<1.4442 && b.photon_idmva()->at(iph) < -0.4) continue;
  //    if (abs(b.photon_eta()->at(iph))>1.4442 && abs(b.photon_eta()->at(iph))<1.566) continue;
  //    if (abs(b.photon_eta()->at(iph))>1.566 && b.photon_idmva()->at(iph) < -0.58) continue;
  //    if (!b.photon_elveto()->at(iph)) continue;
  //    //use first signal photon
  //    TVector3 ph, lep;
  //    ph.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
  //    float mindr = 999;
  //    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
  //      if (b.el_sig()->at(iel)) {
  //        lep.SetPtEtaPhi(b.el_pt()->at(iel),b.el_eta()->at(iel),b.el_phi()->at(iel));
  //        if (ph.DeltaR(lep)<mindr) mindr = ph.DeltaR(lep);
  //      }
  //    }
  //    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
  //      if (b.mu_sig()->at(imu)) {
  //        lep.SetPtEtaPhi(b.mu_pt()->at(imu),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
  //        if (ph.DeltaR(lep)<mindr) mindr = ph.DeltaR(lep);
  //      }
  //    }
  //    return mindr;
  //  }
  //  return 999;
  //});

  //NamedFunc coslowertheta("costheta",[](const Baby &b) -> NamedFunc::ScalarType{
  //  return cos_theta(b);
  //});

  //NamedFunc coscaptheta("coscaptheta",[](const Baby &b) -> NamedFunc::ScalarType{
  //  return cos_Theta(b);
  //});

  //NamedFunc llgphi("phi",[](const Baby &b) -> NamedFunc::ScalarType{
  //  return Getphi(b);
  //});

  //NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
  //  double w_lumi = b.w_lumi();
  //  //fix mis-weighted assoc. prod
  //  if(b.type() >= 200200 && b.type() <= 200500)
  //    w_lumi /= 0.100974;
  //  //fix mis-weighted ttW and WWZ
  //  if (b.FirstFileName().find("ttW") != std::string::npos)
  //    w_lumi = 610.5/27662138;
  //  if (b.FirstFileName().find("WWZ") != std::string::npos)
  //    w_lumi = 354.0/10116400;
  //  return w_lumi;
  //});

  ////multiplicative factor for signal
  //NamedFunc sigx50("sigx50",[](const Baby &b) -> NamedFunc::ScalarType{
  //  double w_lumi = 1.0;
  //  if(b.type() >= 200000 && b.type() <= 200500)
  //    w_lumi *= 50.0;
  //  return w_lumi;
  //});

  //NamedFunc sigx500("sigx500",[](const Baby &b) -> NamedFunc::ScalarType{
  //  double w_lumi = 1.0;
  //  if(b.type() >= 200000 && b.type() <= 200500)
  //    w_lumi *= 500.0;
  //  return w_lumi;
  //});

  //NamedFunc higgs_pt("higgs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if(b.type() >= 200000 && b.type() <= 205000)
  //    return b.w_lumi();
  //  return b.w_lumi();
  //});


  //NamedFunc truth_mllg("truth_mllg",[](const Baby &b) -> NamedFunc::ScalarType{
  //  //calculate truth mllg for dy (pre-stitch) events with ISR
  //  int lep_num = 0;
  //  float max_ph_pt = 0;
  //  TLorentzVector lep1, lep2, ph;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (abs(b.mc_id()->at(imc))==11) {
  //      if (lep_num==0)
  //        lep1.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.000511);
  //      else
  //        lep2.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.000511);
  //      lep_num++;
  //    }
  //    if (abs(b.mc_id()->at(imc))==13) {
  //      if (lep_num==0)
  //        lep1.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.106);
  //      else
  //        lep2.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.106);
  //      lep_num++;
  //    }
  //    if (b.mc_id()->at(imc)==22) {
  //      if (b.mc_pt()->at(imc) > max_ph_pt) {
  //        max_ph_pt = b.mc_pt()->at(imc);
  //        ph.SetPtEtaPhiM(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc),0.0);
  //      }
  //    }
  //  }
  //  if (lep_num<2 || max_ph_pt <= 0) return -999;
  //  return (lep1+lep2+ph).M();
  //});

  //NamedFunc truth_photon_pt("truth_photon_pt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  float max_ph_pt = 0;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22) {
  //      if (b.mc_pt()->at(imc) > max_ph_pt) {
  //        max_ph_pt = b.mc_pt()->at(imc);
  //      }
  //    }
  //  }
  //  return max_ph_pt;
  //});

  //NamedFunc truth_photon_abseta("truth_photon_abseta",[](const Baby &b) -> NamedFunc::ScalarType{
  //  float max_ph_pt = 0;
  //  float ph_abseta = 0;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22) {
  //      if (b.mc_pt()->at(imc) > max_ph_pt) {
  //        max_ph_pt = b.mc_pt()->at(imc);
  //        ph_abseta = abs(b.mc_eta()->at(imc)); 
  //      }
  //    }
  //  }
  //  return ph_abseta;
  //});

  //NamedFunc stitch_zg("stitch_zg",[](const Baby &b) -> NamedFunc::ScalarType{
  //  if (b.type()==6200) return b.stitch_dy();
  //  return 1;
  //});

  //NamedFunc truth_lepton_minpt("truth_lepton_minpt",[](const Baby &b) -> NamedFunc::ScalarType{
  //  float min_lep_pt = 999;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
  //      if (b.mc_pt()->at(imc) < min_lep_pt) {
  //        min_lep_pt = b.mc_pt()->at(imc);
  //      }
  //    }
  //  }
  //  return min_lep_pt;
  //});

  //NamedFunc truth_lepton_eta_acc("truth_lepton_eta_acc",[](const Baby &b) -> NamedFunc::ScalarType{
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
  //      if (abs(b.mc_eta()->at(imc))>2.4)
  //        return 0.0;
  //    }
  //  }
  //  return 1.0;
  //});

  //NamedFunc mc_photon_drmin("mc_photon_drmin",[](const Baby &b) -> NamedFunc::ScalarType{
  //  float ph_eta = -999;
  //  float ph_phi = 0;
  //  float min_dr = 999;
  //  float dr = 999;
  //  bool found_ph = false;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (b.mc_id()->at(imc)==22 && b.mc_mom()->at(imc)==25) {
  //      ph_eta = b.mc_eta()->at(imc);
  //      ph_phi = b.mc_phi()->at(imc);
  //      found_ph = true;
  //    }
  //  }
  //  if (!found_ph) return 999;
  //  for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
  //    if (abs(b.mc_id()->at(imc))==11||abs(b.mc_id()->at(imc))==13) {
  //      dr = deltaR(ph_eta, ph_phi, b.mc_eta()->at(imc), b.mc_phi()->at(imc));
  //      if (dr < min_dr)
  //        min_dr = dr;
  //    }
  //  }
  //  return min_dr;
  //});

  //NamedFunc photon_sig_nodr("photon_sig_nodr",[](const Baby &b) -> NamedFunc::VectorType{
  //  std::vector<double> photon_sig_nodr_;
  //  for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
  //    float abseta = abs(b.photon_eta()->at(iph));
  //    float idmva = b.photon_idmva()->at(iph);
  //    float pt = b.photon_pt()->at(iph);
  //    if (pt>15 && ((abseta<1.4442&&idmva>-0.4)||(abseta>1.566&&abseta<2.5&&idmva>-0.58)))
  //      photon_sig_nodr_.push_back(1);
  //    else
  //      photon_sig_nodr_.push_back(0);
  //  }
  //  return photon_sig_nodr_;
  //});

  //NamedFunc nphoton_nodr("nphoton_nodr",[photon_sig_nodr](const Baby &b) -> NamedFunc::ScalarType{
  //  std::vector<double> photon_sig_nodr_ = photon_sig_nodr.GetVector(b);
  //  int nph = 0;
  //  for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
  //    if (photon_sig_nodr_[iph]>0)
  //      nph++;
  //  }
  //  return nph;
  //});

  //NamedFunc true_photon_recod("true_photon_recod",[photon_sig_nodr](const Baby &b) -> NamedFunc::ScalarType{
  //  std::vector<double> photon_sig_nodr_ = photon_sig_nodr.GetVector(b);
  //  for (unsigned iph = 0; iph < b.photon_pt()->size(); iph++) {
  //    if (photon_sig_nodr_[iph]>0) {
  //      if (b.photon_pflavor()->at(iph)==1)
  //        return 1;
  //    }
  //  }
  //  return 0;
  //});

  //NamedFunc leps_recod("leps_recod",[](const Baby &b) -> NamedFunc::ScalarType{
  //  int reco_num = 0;
  //  float true_lep_eta = 999;
  //  float true_lep_phi = 0;
  //  float true_antilep_eta = 999;
  //  float true_antilep_phi = 0;
  //  for (unsigned imc = 0; imc < b.mc_pt()->size(); imc++) {
  //    int pdgid = b.mc_id()->at(imc);
  //    if ((pdgid == 11 || pdgid == 13) && b.mc_mom()->at(imc)==23) {
  //      true_lep_eta = b.mc_eta()->at(imc);
  //      true_lep_phi= b.mc_phi()->at(imc);
  //    }
  //    if ((pdgid == -11 || pdgid == -13) && b.mc_mom()->at(imc)==23) {
  //      true_antilep_eta = b.mc_eta()->at(imc);
  //      true_antilep_phi= b.mc_phi()->at(imc);
  //    }
  //  }
  //  for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
  //    if (b.el_sig()->at(iel)) {
  //      float dr_lep = deltaR(b.el_eta()->at(iel), b.el_phi()->at(iel), true_lep_eta, true_lep_phi);
  //      float dr_antilep = deltaR(b.el_eta()->at(iel), b.el_phi()->at(iel), true_antilep_eta, true_antilep_phi);
  //      if (dr_lep < 0.3 || dr_antilep < 0.3)
  //        reco_num++;
  //    }
  //  }
  //  for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
  //    if (b.mu_sig()->at(imu)) {
  //      float dr_lep = deltaR(b.mu_eta()->at(imu), b.mu_phi()->at(imu), true_lep_eta, true_lep_phi);
  //      float dr_antilep = deltaR(b.mu_eta()->at(imu), b.mu_phi()->at(imu), true_antilep_eta, true_antilep_phi);
  //      if (dr_lep < 0.3 || dr_antilep < 0.3)
  //        reco_num++;
  //    }
  //  }
  //  return reco_num;
  //});

  NamedFunc higgs_window = "llphoton_m[0]>122&&llphoton_m[0]<128";

  //------------------------------------------------------------------------------------
  //                                        BDTs
  //------------------------------------------------------------------------------------

  //bool plot_sig_accept = false;
  //bool plot_other_objs = true;
  bool plot_tth_leptonic = false;
  bool plot_tth_hadronic = false;
  bool plot_vh_leptonic = false;
  bool plot_vh_met = false;
  bool plot_vh_hadronic = false;
  bool plot_vbf = false;
  bool plot_yield_tables = false;
  bool plot_vh_hadronic_pico = true;
  bool plot_mllgs = false;

  MVAWrapper tthlep_bdt_reader("tthlep_bdt");
  NamedFunc tthlep_bdt_score = 1.0;
  MVAWrapper tthhad_bdt_reader("tthhad_bdt");
  NamedFunc tthhad_bdt_score = 1.0;
  MVAWrapper vhlep_bdt_reader("vhlep_bdt");
  NamedFunc vhlep_bdt_score = 1.0;
  MVAWrapper vhmet_bdt_reader("vhmet_bdt");
  NamedFunc vhmet_bdt_score = 1.0;
  MVAWrapper kin_bdt_reader("kinematic_bdt");
  NamedFunc kin_bdt_score = 1.0;
  if (plot_tth_leptonic) {
    tthlep_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
    tthlep_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
    tthlep_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
    tthlep_bdt_reader.SetVariable("min_lepton_quality",min_lep_id);
    tthlep_bdt_reader.SetVariable("assoc_lepton_quality",assoc_lep_id);
    tthlep_bdt_reader.SetVariable("max_lepton_miniso",max_lep_miniso);
    tthlep_bdt_reader.SetVariable("assoc_lepton_miniso",assoc_lep_miniso);
    tthlep_bdt_reader.SetVariable("min_lepton_pt",min_lep_pt);
    tthlep_bdt_reader.SetVariable("assoc_lepton_pt",assoc_lep_pt);
    tthlep_bdt_reader.SetVariable("assoc_lepton_pdgid",assoc_lep_pdgid);
    tthlep_bdt_reader.SetVariable("ht","ht");
    tthlep_bdt_reader.SetVariable("met","met");
    tthlep_bdt_reader.SetVariable("njet_f","njet");
    tthlep_bdt_reader.SetVariable("nbdfl_f","nbdfl");
    tthlep_bdt_reader.SetVariable("nbdfm_f","nbdfm");
    tthlep_bdt_reader.SetVariable("nbdft_f","nbdft");
    tthlep_bdt_reader.SetVariable("nlep_f","nlep");
    tthlep_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_tthlep_bdt_BDT.weights.xml");
    tthlep_bdt_score = tthlep_bdt_reader.GetDiscriminant();
  }
  if (plot_tth_hadronic) {
    tthhad_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
    tthhad_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
    tthhad_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
    tthhad_bdt_reader.SetVariable("ht","ht");
    tthhad_bdt_reader.SetVariable("met","met");
    tthhad_bdt_reader.SetVariable("njet_f","njet");
    tthhad_bdt_reader.SetVariable("nbdfl_f","nbdfl");
    tthhad_bdt_reader.SetVariable("nbdfm_f","nbdfm");
    tthhad_bdt_reader.SetVariable("nbdft_f","nbdft");
    tthhad_bdt_reader.SetVariable("top_mass",closest_top_mass);
    tthhad_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_tthhad_bdt_BDT.weights.xml");
    tthhad_bdt_score = tthhad_bdt_reader.GetDiscriminant();
  }
  if (plot_vh_leptonic) {
    vhlep_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
    vhlep_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
    vhlep_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
    vhlep_bdt_reader.SetVariable("min_lepton_quality",min_lep_id);
    vhlep_bdt_reader.SetVariable("max_lepton_miniso",max_lep_miniso);
    vhlep_bdt_reader.SetVariable("min_lepton_pt",min_lep_pt);
    vhlep_bdt_reader.SetVariable("assoc_lepton_quality",assoc_lep_id);
    vhlep_bdt_reader.SetVariable("assoc_lepton_miniso",assoc_lep_miniso);
    vhlep_bdt_reader.SetVariable("assoc_lepton_pt",assoc_lep_pt);
    vhlep_bdt_reader.SetVariable("assoc_lepton_pdgid",assoc_lep_pdgid);
    vhlep_bdt_reader.SetVariable("assoc_lepton_mt",assoc_lep_mt);
    vhlep_bdt_reader.SetVariable("ht","ht");
    vhlep_bdt_reader.SetVariable("met","met");
    vhlep_bdt_reader.SetVariable("njet_f","njet");
    vhlep_bdt_reader.SetVariable("nlep_f","nlep");
    vhlep_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_vhlep_bdt_BDT.weights.xml");
    vhlep_bdt_score = vhlep_bdt_reader.GetDiscriminant();
  }
  if (plot_vh_met) {
    vhmet_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
    vhmet_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
    vhmet_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
    vhmet_bdt_reader.SetVariable("dphi_h_met",dphi_h_met);
    vhmet_bdt_reader.SetVariable("met","met");
    vhmet_bdt_reader.SetVariable("njet_f","njet");
    vhmet_bdt_reader.SetVariable("nlep_f","nlep");
    vhmet_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_vhmet_bdt_BDT.weights.xml");
    vhmet_bdt_score = vhmet_bdt_reader.GetDiscriminant();
  }
  if (plot_mllgs) {
    kin_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
    kin_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
    kin_bdt_reader.SetVariable("max_dR",photon_drmax);
    kin_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
    kin_bdt_reader.SetVariable("cosTheta","llphoton_cosTheta[0]");
    kin_bdt_reader.SetVariable("costheta","llphoton_costheta[0]");
    kin_bdt_reader.SetVariable("phi","llphoton_psi[0]");
    kin_bdt_reader.SetVariable("photon_res","photon_pterr[0]/photon_pt[0]");
    kin_bdt_reader.SetVariable("photon_rapidity","photon_eta[0]");
    kin_bdt_reader.SetVariable("l1_rapidity",lead_lepton_eta);
    kin_bdt_reader.SetVariable("l2_rapidity",sublead_lepton_eta);
    kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_idmvacut_run2_BDT.weights.xml");
    kin_bdt_score = kin_bdt_reader.GetDiscriminant();
  }
  
  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  pm.multithreaded_ = true;

  //move acceptance studies to another file
  //if (plot_sig_accept) {
  //  pm.Push<Hist1D>(Axis(40,0,4.0, mc_photon_drmin, "Truth Min #Delta R_{#gamma, l}", {}), 
  //      "1", procs_sig_full, ops).Weight(wgt*sigx50).Tag("zgtth");
  //  pm.Push<Table>("zgtth__reco", vector<TableRow>{
  //    TableRow("Inclusive", 
  //        "1",0,0,wgt),
  //    TableRow("1 Leptons Reconstructed", 
  //        leps_recod>=1,0,0,wgt),
  //    TableRow("2 Leptons Reconstructed", 
  //        leps_recod>=2,0,0,wgt),
  //    TableRow("Photon reconstructed (no dr)", 
  //        leps_recod>=2&&true_photon_recod,0,0,wgt),
  //    TableRow("Photon reconstructed (dr)", 
  //        leps_recod>=2&&"nphoton>0"&&true_photon_recod,0,0,wgt),
  //    TableRow("Photon reconstructed in slot 0 (no dr)", 
  //        leps_recod>=2&&nphoton_nodr>0.0&&"photon_pflavor[0]==1",0,0,wgt),
  //    TableRow("Photon reconstructed in slot 0 (dr)", 
  //        leps_recod>=2&&"nphoton>0&&photon_pflavor[0]==1",0,0,wgt),
  //    TableRow("Triggers", 
  //        leps_recod>=2&&"nphoton>0"&&true_photon_recod&&trigs,0,0,wgt),
  //    TableRow("b-Jet", 
  //        leps_recod>=2&&"nphoton>0"&&true_photon_recod&&trigs&&"nbdfl>=1",0,0,wgt),
  //    TableRow("Lepton", 
  //        leps_recod>=2&&"nlep>=3&&nphoton>0"&&true_photon_recod&&trigs,0,0,wgt),
  //    TableRow("Lepton+b-Jet", 
  //        leps_recod>=2&&"nlep>=3&&nphoton>0"&&true_photon_recod&&trigs&&"nbdfl>=1",0,0,wgt),
  //  },procs_sig_full,false,true,false,false,false,true).Precision(3);
  //}

  //if (plot_other_objs) {

  //  //to determine b-tagging
  //  pm.Push<Table>("zgtth__btag", vector<TableRow>{
  //    TableRow("Inclusive", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130",0,0,wgt),
  //    TableRow("1 Loose", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=1",0,0,wgt),
  //    TableRow("1 Medium", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=1",0,0,wgt),
  //    TableRow("1 Tight", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Loose", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2",0,0,wgt),
  //    TableRow("Loose+Medium", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdfm>=1",0,0,wgt),
  //    TableRow("Loose+Tight", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Medium", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2",0,0,wgt),
  //    TableRow("Medium+Tight", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Tight", 
  //        "nlep>=2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=2",0,0,wgt),

  //    TableRow("2 L", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130",0,0,wgt),
  //    TableRow("1 Loose", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=1",0,0,wgt),
  //    TableRow("1 Medium", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=1",0,0,wgt),
  //    TableRow("1 Tight", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Loose", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2",0,0,wgt),
  //    TableRow("Loose+Medium", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdfm>=1",0,0,wgt),
  //    TableRow("Loose+Tight", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Medium", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2",0,0,wgt),
  //    TableRow("Medium+Tight", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Tight", 
  //        "nlep==2&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=2",0,0,wgt),

  //    TableRow("3 L", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130",0,0,wgt),
  //    TableRow("1 Loose", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=1",0,0,wgt),
  //    TableRow("1 Medium", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=1",0,0,wgt),
  //    TableRow("1 Tight", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Loose", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2",0,0,wgt),
  //    TableRow("Loose+Medium", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdfm>=1",0,0,wgt),
  //    TableRow("Loose+Tight", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfl>=2&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Medium", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2",0,0,wgt),
  //    TableRow("Medium+Tight", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdfm>=2&&nbdft>=1",0,0,wgt),
  //    TableRow("2 Tight", 
  //        "nlep>=3&&nphoton>=1&&llphoton_m[0]>120&&llphoton_m[0]<130&&nbdft>=2",0,0,wgt),
  //  },procs,false,true,false,false,false,true).Precision(3);
  //}
  
  //ttH leptonic studies
  if (plot_tth_leptonic) {
    pm.multithreaded_ = true;
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1"&&max_lep_miniso<0.1, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1&&njet>=3"&&max_lep_miniso<0.1, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1&&(njet>=3||nlep>=4)"&&max_lep_miniso<0.1, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1&&njet>=2"&&((("ht>=200||njet>=4"||llphoton_rel_pt>1.0)&&min_lep_id>=2.0&&max_lep_miniso<0.15&&min_lep_pt>12.0)||("nlep>=4&&ht>120"&&max_lep_miniso<0.2)), procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1"&&max_lep_miniso<0.15, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1"&&min_lep_id>=2.0, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl>=1"&&tthlep_bdt_score>0.22, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,-1.0,1.0, tthlep_bdt_score, "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,-1.0,1.0, tthlep_bdt_score, "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl>=1&&llphoton_m[0]>120&&llphoton_m[0]<130", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.0, llphoton_rel_pt, "Higgs candidate p_{T}/m", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.0, llphoton_rel_pt, "Higgs candidate p_{T}/m", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1"&&min_lep_id>=2.0&&max_lep_miniso<0.15, procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, min_lep_id, "Minimum lepton ID quality", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, min_lep_id, "Minimum lepton ID quality", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, max_lep_miniso, "Maximum lepton I_{mini}", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, max_lep_miniso, "Maximum lepton I_{mini}", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, max_lep_miniso, "Maximum lepton I_{mini}", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1"&&min_lep_id>=2.0, procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,50.0, min_lep_pt, "Minimum lepton p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,50.0, min_lep_pt, "Minimum lepton p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,0.0,300.0, "ht", "H_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,0.0,300.0, "ht", "H_{T} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,0.0,300.0, "ht", "H_{T} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1"&&min_lep_id>=2.0&&max_lep_miniso<0.15, procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(7,-0.5,6.5, "njet", "N_{J}", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(7,-0.5,6.5, "njet", "N_{J}", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(7,-0.5,6.5, "njet", "N_{J}", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1"&&min_lep_id>=2.0&&max_lep_miniso<0.15, procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, "nbdfl", "N_{bl}", {}), 
        tighter_baseline&&"nlep==3&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, "nbdfl", "N_{bl}", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
  }

  //ttH hadronic studies
  if (plot_tth_hadronic) {
    pm.multithreaded_ = true;
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=2", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=2", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=2&&nbdfm>=1", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=2&&nbdfm>=1&&njet>=5", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=2&&njet>=5", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=5", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5"&&closest_top_mass>130&&closest_top_mass<220, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=3"&&tthhad_bdt_score>0.14, procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_drmin[0]<1.2&&met<120", 
        procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.0, "photon_drmin[0]", "Minimum #Delta R(l,#gamma)", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_idmva[0]>0.6", 
        procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.0, llphoton_rel_pt, "p_{Tll#gamma}/m_{ll#gamma}", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_idmva[0]>0.6", 
        procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,300, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_idmva[0]>0.6", 
        procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=2&&nbdfm>=1&&njet>=6", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=2&&njet>=6", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=6", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=6", procs_tth, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, "nbdfl", "N_{bl}", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, "nbdfm", "N_{bm}", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(10,-0.5,9.5, "njet", "N_{j}", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(50,0.0,500.0, "ht", "H_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(50,0.0,500.0, "ht", "H_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(10,0.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        tighter_baseline&&"nlep==2&&nbdfl>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, closest_top_mass, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=4", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, closest_top_mass, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, dr_top_mass, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=4", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, dr_top_mass, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, full_dr_top_mass, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=4", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, full_dr_top_mass, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, closest_top_mass_constrained, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=4", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,400.0, closest_top_mass_constrained, "Top candidate m [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5", procs_tth, ops_sorb).Weight("w_lumi"*w_years*w_run3*w_sigx20).Tag("zgassoc");
  }

  if (plot_vh_leptonic) {
    pm.multithreaded_ = true;
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pdgid==11, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pdgid==13, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=4&&nbdfl==0", procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0"&&assoc_lep_pt>20.0&&assoc_lep_miniso<0.2&&assoc_lep_id>1.0, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0"&&assoc_lep_pt>20.0&&assoc_lep_miniso<0.2&&assoc_lep_id>1.0, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0"&&assoc_lep_pt>20.0&&assoc_lep_miniso<0.2&&assoc_lep_id>1.0, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0"&&((assoc_lep_pt>20.0&&assoc_lep_miniso<0.2&&assoc_lep_id>1.0)||("nlep>=4")), procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl==0"&&vhlep_bdt_score>0.16, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0&&photon_drmin[0]<1.5"&&min_lep_pt>12.0&&min_lep_id>=2.0&&max_lep_miniso<0.2&&("met>30"||assoc_lep_mt>50), procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,0.0,1.0, max_lep_miniso, "Maximum lepton I_{mini}", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(20,0.0,50.0, min_lep_pt, "Minimum lepton p_{T} [GeV]", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    //pm.Push<Hist1D>(Axis(4,-0.5,3.5, min_lep_id, "Minimum lepton ID quality", {}), 
    //    tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, assoc_lep_miniso, "Associated lepton I_{mini}", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,50.0, assoc_lep_pt, "Associated lepton p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, assoc_lep_id, "Associated lepton ID quality", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl==0", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, assoc_lep_miniso, "Associated lepton I_{mini}", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0&&assoc_lep_pdgid==11, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,50.0, assoc_lep_pt, "Associated lepton p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pdgid==11, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, assoc_lep_id, "Associated lepton ID quality", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_pdgid==11, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, assoc_lep_miniso, "Associated lepton I_{mini}", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0&&assoc_lep_pdgid==13, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,50.0, assoc_lep_pt, "Associated lepton p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pdgid==13, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, assoc_lep_id, "Associated lepton ID quality", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_pdgid==13, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,0.0,300.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0&&assoc_lep_pdgid==11, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,0.0,300.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0&&assoc_lep_pdgid==13, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,100.0, assoc_lep_mt, "m_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0&&assoc_lep_pdgid==11, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,100.0, assoc_lep_mt, "m_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0&&assoc_lep_pdgid==13, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,80.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(10,0.0,1.0, "photon_idmva[0]", "Photon IDMVA [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0"&&assoc_lep_pt>20&&assoc_lep_id>2.0, procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
  }

  if (plot_vh_met) {
    pm.multithreaded_ = false;
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>80", procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>110", procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>110"&&dphi_h_met>2.3, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>110"&&dphi_h_met>2.3&&"llphoton_pt[0]>75", procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>105"&&dphi_h_met>2.1&&"llphoton_pt[0]>70", procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>90"&&vhmet_bdt_score>0.12, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>95&&photon_idmva[0]>0.55&&photon_drmin[0]<1.5"&&llphoton_rel_pt>0.7&&dphi_h_met>1.5, procs_vh, ops).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,80.0,320.0, "met", "p_{T}^{miss} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>80", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.1415, min_dphi_jet_met, "Min #Delta #phi(jet,p_{T}^{miss}) [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>80", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.1415, max_dphi_jet_met, "Max #Delta #phi(jet,p_{T}^{miss}) [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>80", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,3.1415, dphi_h_met, "#Delta #phi(H,p_{T}^{miss}) [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>110", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,80.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>110", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, "llphoton_pt[0]", "Higgs candidate p_{T} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>110", procs_vh, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
  }

  if (plot_vh_hadronic) {
    pm.Push<Hist1D>(Axis(20,100.0,160.0, NanoFunctions::HiggsCandidate_mass, "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&
        Lead_FatJet_pt>200&&Lead_FatJet_deepTag_WvsQCD>0.9&&Lead_FatJet_msoftdrop>60&&Lead_FatJet_msoftdrop<110,
        procs_nano, ops).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(15,200.0,500.0, Lead_FatJet_pt, "Lead AK8 Jet p_{T} [GeV]", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&Lead_FatJet_pt>200,
        procs_nano, ops_sorb).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, Lead_FatJet_msoftdrop, "Lead AK8 jet m_{softdrop} [GeV]", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&Lead_FatJet_pt>200,
        procs_nano, ops_log).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,200.0, Max_FatJet_msoftdrop, "Maximum AK8 jet m_{softdrop} [GeV]", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&Lead_FatJet_pt>200,
        procs_nano, ops_log).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, Lead_FatJet_deepTag_WvsQCD, "Lead AK8 jet W-tag score", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&Lead_FatJet_pt>200,
        procs_nano, ops_sorb).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, Lead_FatJet_deepTag_WvsQCD, "Lead AK8 jet W-tag score", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&Lead_FatJet_pt>200,
        procs_nano, ops_log).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, Lead_FatJet_deepTag_WvsQCD, "Lead AK8 jet W-tag score", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&
        Lead_FatJet_pt>200&&Lead_FatJet_msoftdrop>60&&Lead_FatJet_msoftdrop<110,
        procs_nano, ops_sorb).Weight(weight_nano).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, Lead_FatJet_deepTag_WvsQCD, "Lead AK8 jet W-tag score", {}), 
        tighter_baseline_nano&&NanoFunctions::nSignalLepton==2&&NanoFunctions::nJet_bdfm==0.0&&"MET_pt<110"&&
        Lead_FatJet_pt>200&&Lead_FatJet_msoftdrop>60&&Lead_FatJet_msoftdrop<110,
        procs_nano, ops_log).Weight(weight_nano).Tag("zgassoc");
  }

  if (plot_vbf) {
    pm.multithreaded_ = false;
    pm.Push<Hist1D>(Axis(50,0.0,1000.0, "dijet_m", "m_{jj} [GeV]", {}), 
        tighter_baseline&&"njet>=2", procs_vbf, ops_sorb).Weight("w_lumi"*w_years*w_run3).Tag("zgassoc");
  }

  if (plot_yield_tables) { 
    pm.Push<Table>("mc_table", vector<TableRow>{
      TableRow("ttH leptonic selection weighted", 
          tighter_baseline&&"nlep>=3&&nbdfl>=1&&njet>=2"&&((("ht>=200||njet>=4"||llphoton_rel_pt>1.0)&&min_lep_id>=2.0&&max_lep_miniso<0.15&&min_lep_pt>12.0)||("nlep>=4&&ht>120"&&max_lep_miniso<0.2))&&"llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,w_ewkzgamma*w_years),
      TableRow("ttH leptonic selection unweighted", 
          tighter_baseline&&"nlep>=3&&nbdfl>=1&&njet>=2"&&((("ht>=200||njet>=4"||llphoton_rel_pt>1.0)&&min_lep_id>=2.0&&max_lep_miniso<0.15&&min_lep_pt>12.0)||("nlep>=4&&ht>120"&&max_lep_miniso<0.2))&&"llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,"1"),
      TableRow("ttH hadronic selection weighted", 
          tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_drmin[0]<1.2&&met<120"&&"llphoton_m[0]>100&&llphoton_m[0]<150", 
          0,0,w_ewkzgamma*w_years),
      TableRow("ttH hadronic selection unweighted", 
          tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_drmin[0]<1.2&&met<120"&&"llphoton_m[0]>100&&llphoton_m[0]<150", 
          0,0,"1"),
      TableRow("VH leptonic selection weighted", 
          tighter_baseline&&"nlep>=3&&nbdfl==0&&photon_drmin[0]<1.5"&&min_lep_pt>12.0&&min_lep_id>=2.0&&max_lep_miniso<0.2&&("met>30"||assoc_lep_mt>50)&&"llphoton_m[0]>100&&llphoton_m[0]<150", 
          0,0,w_ewkzgamma*w_years),
      TableRow("VH leptonic selection unweighted", 
          tighter_baseline&&"nlep>=3&&nbdfl==0&&photon_drmin[0]<1.5"&&min_lep_pt>12.0&&min_lep_id>=2.0&&max_lep_miniso<0.2&&("met>30"||assoc_lep_mt>50)&&"llphoton_m[0]>100&&llphoton_m[0]<150", 
          0,0,"1"),
      TableRow("VH MET selection weighted", 
          tighter_baseline&&"nlep==2&&nbdfm==0&&met>95&&photon_idmva[0]>0.55&&photon_drmin[0]<1.5"&&llphoton_rel_pt>0.7&&dphi_h_met>1.5&&"llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,w_ewkzgamma*w_years),
      TableRow("VH MET selection unweighted", 
          tighter_baseline&&"nlep==2&&nbdfm==0&&met>95&&photon_idmva[0]>0.55&&photon_drmin[0]<1.5"&&llphoton_rel_pt>0.7&&dphi_h_met>1.5&&"llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,"1"),
      TableRow("Simple VBF selection weighted", 
          tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&dijet_m>750&&llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,w_ewkzgamma*w_years),
      TableRow("Simple VBF selection unweighted", 
          tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&dijet_m>750&&llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,"1"),
      TableRow("Simple ggF selection weighted", 
          tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&dijet_m<750&&llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,w_ewkzgamma*w_years),
      TableRow("Simple ggF selection unweighted", 
          tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&dijet_m<750&&llphoton_m[0]>100&&llphoton_m[0]<150",
          0,0,"1"),
      TableRow("\\hline Baseline selection", 
          tighter_baseline,
          0,0,w_ewkzgamma*w_years),
      TableRow("Baseline selection unweighted", 
          tighter_baseline,
          0,0,"1"),
    },procs_split,false,true,false,false,false,true).Precision(3);
  }

  if (plot_mllgs) {
    pm.multithreaded_ = false;
    //tth leptonic
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep>=3&&nbdfl>=1&&njet>=2"&&((("ht>=200||njet>=4"||llphoton_rel_pt>1.0)&&min_lep_id>=2.0&&max_lep_miniso<0.15&&min_lep_pt>12.0)||("nlep>=4&&ht>120"&&max_lep_miniso<0.2)), 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //tth hadronic
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm>=1&&njet>=5&&photon_drmin[0]<1.2&&met<120", 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //vh 4l
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==4&&nbdfl==0"&&min_lep_pt>12.0&&min_lep_id>=2.0&&max_lep_miniso<0.2, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //vh 3l
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==3&&nbdfl==0&&photon_drmin[0]<1.5"&&min_lep_pt>12.0&&min_lep_id>=2.0&&max_lep_miniso<0.2&&("met>30"||assoc_lep_mt>50), 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //vh met
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>95&&photon_idmva[0]>0.55&&photon_drmin[0]<1.5"&&llphoton_rel_pt>0.7&&dphi_h_met>1.5, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //vh/bbh 2b
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&njet<=4&&nbdfl>=2&&nbdfm>=1&&met<95", 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //vbf categories
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m>=750"
        &&kin_bdt_score<0.16, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m>=750"
        &&kin_bdt_score>0.16&&kin_bdt_score<0.42, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m>=750"
        &&kin_bdt_score>0.42, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    //ggf categories
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m<750"
        &&kin_bdt_score<-0.12, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m<750"
        &&kin_bdt_score>-0.12&&kin_bdt_score<0.1, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m<750"
        &&kin_bdt_score>0.1&&kin_bdt_score<0.36, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(15,100.0,160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        tighter_baseline&&"nlep==2&&nbdfm==0&&nbdfl<2&&met<95&&dijet_m<750"
        &&kin_bdt_score>0.36, 
        procs, ops).Weight(w_ewkzgamma*w_years*w_run3).Tag("zgassoc");
  }

  if (plot_vh_hadronic_pico) {
    pm.multithreaded_ = false;
    pm.Push<Hist1D>(Axis(4,-0.5,3.5, nsig_fjet, "N_{AK8 Jet}", {}), 
        january_baseline&&"nbm==0"&&higgs_window,
        procs, ops_shapes).Weight("weight"*w_years).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,200.0,500.0, lead_sig_fjet_pt, "Fat jet p_{T} [GeV]", {}), 
        january_baseline&&"nbm==0"&&higgs_window&&(nsig_fjet>=1),
        procs, ops_shapes).Weight("weight"*w_years).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, lead_sig_fjet_pnet_wtag, "Fat jet W-tag score", {}), 
        january_baseline&&"nbm==0"&&higgs_window&&(nsig_fjet>=1),
        procs, ops_shapes).Weight("weight"*w_years).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,1.0, lead_sig_fjet_pnet_ttag, "Fat jet top-tag score", {}), 
        january_baseline&&"nbm==0"&&higgs_window&&(nsig_fjet>=1),
        procs, ops_shapes).Weight("weight"*w_years).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,4.0, "llphoton_pt[0]/llphoton_m[0]", "p_{Tll#gamma}/m_{ll#gamma}", {}), 
        january_baseline&&"nbm==0"&&higgs_window&&(nsig_fjet>=1),
        procs, ops_shapes).Weight("weight"*w_years).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(20,0.0,4.0, "llphoton_pt[0]/llphoton_m[0]", "p_{Tll#gamma}/m_{ll#gamma}", {}), 
        january_baseline&&"nbm==0"&&higgs_window&&(nsig_fjet>=1)&&lead_sig_fjet_pnet_wtag>0.7,
        procs, ops).Weight("weight"*w_years).Tag("zgassoc");
    pm.Push<Hist1D>(Axis(30,100,160, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        january_baseline&&"nbm==0"&&(nsig_fjet>=1)&&lead_sig_fjet_pnet_wtag>0.7&&"(llphoton_pt[0]/llphoton_m[0])>1.5",
        procs, ops).Weight("weight"*w_years).Tag("zgassoc");
  }

  pm.min_print_ = true;
  pm.SetLuminosityTag("200").MakePlots(1.0);

  return 0;
}
