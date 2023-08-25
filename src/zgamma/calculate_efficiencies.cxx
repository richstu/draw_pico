/**
 * File used to generate various efficiency measurements and scale factors for 
 * the H->Zgamma analysis
 */

#include <cmath>
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
#include "TFile.h"
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
#include "core/process.hpp"
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

  //std::vector<std::shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllMoreTrigs");
  //reshuffle DY and ZG into real and fake photons
  string prod_folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v3/");
  std::set<int> years = {2016,2017,2018};
  Process::Type back =  Process::Type::background;
  Palette mc_colors("txt/colors_zgamma.txt","default");
  const NamedFunc trig_and_stitch = (HLT_pass_dilepton||HLT_pass_singlelepton)&&stitch;

  auto proc_smzg      = Process::MakeShared<Baby_pico>("Z+Real Photon", back, mc_colors("zgtollg"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*DYJets*","*ZGToLLG*"}), 
                           trig_and_stitch&&"photon_pflavor[0]==1");
  auto proc_dy        = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, mc_colors("dyjets"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*DYJets*"}), 
                           trig_and_stitch&&"photon_pflavor[0]!=1");
  auto proc_tt        = Process::MakeShared<Baby_pico>("tt", back, mc_colors("tt"),
                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*TTTo2L2Nu*"}), 
                           trig_and_stitch);
  auto proc_dy_unskimmed = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, mc_colors("dyjets"),
                           attach_folder(prod_folder,years,"mc/unskimmed",{"*DYJets*"}),"1");

  vector<shared_ptr<Process>> procs = {proc_smzg, proc_tt, proc_dy};
  vector<shared_ptr<Process>> procs_dy = {proc_dy_unskimmed};

  std::vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumi").Overflow(OverflowType::none)}; 
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","LinLumi").Overflow(OverflowType::none),
                                     PlotOpt("txt/plot_styles.txt","Shapes").Overflow(OverflowType::none)}; 
  std::vector<PlotOpt> ops_twodim = {PlotOpt("txt/plot_styles.txt","Eff2D").Overflow(OverflowType::none)}; 

  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------


  //max lepton dxy
  const NamedFunc max_lep_dxy("max_lep_dxy",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_lep_dxy_ = 0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (fabs(b.el_dxy()->at(iel)) > max_lep_dxy_)
        max_lep_dxy_ = fabs(b.el_dxy()->at(iel));
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (fabs(b.mu_dxy()->at(imu)) > max_lep_dxy_)
        max_lep_dxy_ = fabs(b.mu_dxy()->at(imu));
    }
    return max_lep_dxy_;
  });

  //max lepton dz
  const NamedFunc max_lep_dz("max_lep_dz",[](const Baby &b) -> NamedFunc::ScalarType{
    float max_lep_dz_ = 0;
    for (unsigned iel = 0; iel < b.el_pt()->size(); iel++) {
      if (!b.el_sig()->at(iel)) continue;
      if (fabs(b.el_dz()->at(iel)) > max_lep_dz_)
        max_lep_dz_ = fabs(b.el_dz()->at(iel));
    }
    for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
      if (!b.mu_sig()->at(imu)) continue;
      if (fabs(b.mu_dz()->at(imu)) > max_lep_dz_)
        max_lep_dz_ = fabs(b.mu_dz()->at(imu));
    }
    return max_lep_dz_;
  });

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

  //llphoton z momentum
  const NamedFunc llphoton_pz("llphoton_pz",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.llphoton_m()->size()==0) return 0;
    return b.llphoton_pt()->at(0)*sinh(b.llphoton_eta()->at(0));
  });

  //llphoton abs cosTheta
  const NamedFunc llphoton_abscosTheta("llphoton_abscosTheta",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.llphoton_m()->size()==0) return 0;
    return fabs(b.llphoton_cosTheta()->at(0));
  });

  //reconstructed pt for truth muons; if no reco muon, returns truth pt
  //N.B. reco efficiency is ~1 and SFs are binned in reco quantities. 
  //Furthermore, for muons reco and truth quantities are nearly identical anyway
  const NamedFunc mc_muon_reco_pt("mc_muon_reco_pt",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> reco_pt;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==13) { //is muon
        float min_dr = 999;
        unsigned mu_idx = 999;
        for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
          float dr = deltaR(b.mc_eta()->at(imc),b.mc_phi()->at(imc),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
          if (dr < 0.2 && dr < min_dr) {
            mu_idx = imu;
            min_dr = dr;
          }
        }
        if (mu_idx != 999)
          reco_pt.push_back(b.mu_pt()->at(mu_idx));
        else
          reco_pt.push_back(b.mc_pt()->at(imc));
      }
    }
    return reco_pt;
  });

  //reconstructed eta for truth muons; if no reco muon, returns truth eta
  const NamedFunc mc_muon_reco_abseta("mc_muon_reco_abseta",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> reco_abseta;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==13) { //is muon
        float min_dr = 999;
        unsigned mu_idx = 999;
        for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
          float dr = deltaR(b.mc_eta()->at(imc),b.mc_phi()->at(imc),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
          if (dr < 0.2 && dr < min_dr) {
            mu_idx = imu;
            min_dr = dr;
          }
        }
        if (mu_idx != 999)
          reco_abseta.push_back(fabs(b.mu_eta()->at(mu_idx)));
        else
          reco_abseta.push_back(fabs(b.mc_eta()->at(imc)));
      }
    }
    return reco_abseta;
  });

  //if MC muon is reconstrcuted and passes analysis cuts
  const NamedFunc mc_muon_reco_sig("mc_muon_reco_sig",[](const Baby &b) -> NamedFunc::VectorType{
    std::vector<double> reco_sig;
    for (unsigned imc = 0; imc < b.mc_id()->size(); imc++) {
      if (abs(b.mc_id()->at(imc))==13) { //is muon
        float min_dr = 999;
        float sig = 0;
        for (unsigned imu = 0; imu < b.mu_pt()->size(); imu++) {
          float dr = deltaR(b.mc_eta()->at(imc),b.mc_phi()->at(imc),b.mu_eta()->at(imu),b.mu_phi()->at(imu));
          if (dr < 0.2 && dr < min_dr) {
            min_dr = dr;
            if (b.mu_sig()->at(imu)) sig = 1;
            else sig = 0;
          }
        }
        reco_sig.push_back(sig);
      }
    }
    return reco_sig;
  });

  const NamedFunc tighter_baseline = NamedFunc(zg_baseline&&"photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100").Name("new_baseline");

  const NamedFunc zg_el_cuts = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)";
  const NamedFunc zg_mu_cuts = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)";
  NamedFunc photon_sf_baseline = NamedFunc("nlep==2 && nphoton==1 && (photon_drmin[0]>0.4)" && (zg_el_cuts || zg_mu_cuts)).Name("photon_sf_region");
  
  const NamedFunc phmva = "photon_idmva[0]";
  const NamedFunc phpt = "photon_pt[0]";
  const NamedFunc phabseta("phabseta",[](const Baby &b) -> NamedFunc::ScalarType{
    return fabs(b.photon_eta()->at(0));
  });

  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  pm.multithreaded_ = true;

  bool calc_photon_sf = false;
  bool calc_muon_eff = true;

  const std::vector<double> muon_abseta_bins = {0.0, 0.9, 1.2, 2.1, 2.4};
  const std::vector<double> muon_pt_bins = {5.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,500.0};

  if (calc_muon_eff) {
    //caculate muon MC Reco+ID+Iso efficiency as a function of pt and eta
    pm.Push<Hist2D>(Axis(muon_abseta_bins, mc_muon_reco_abseta, "Reconstructed |#eta|", {}), 
        Axis(muon_pt_bins, mc_muon_reco_pt, "Reconstructed p_{T} [GeV]", {}), 
        "1", procs_dy, ops_twodim).Weight("1").Tag("FixName:muon_eff_den");
    pm.Push<Hist2D>(Axis(muon_abseta_bins, mc_muon_reco_abseta, "Reconstructed |#eta|", {}), 
        Axis(muon_pt_bins, mc_muon_reco_pt, "Reconstructed p_{T} [GeV]", {}), 
        mc_muon_reco_sig, procs_dy, ops_twodim).Weight("1").Tag("FixName:muon_eff_num");
  }

  if (calc_photon_sf) {
    //mu mu gamma proof of concept
    //inclusive
    pm.Push<Hist1D>(Axis(30,70.0,100.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
        photon_sf_baseline, procs, ops).Weight(w_ewkzgamma*w_years).Tag("zgsf");
    pm.Push<Hist1D>(Axis(50,-1.0,1.0, "photon_idmva[0]", "Photon IDMVA", {}), 
        photon_sf_baseline, procs, ops).Weight(w_ewkzgamma*w_years).Tag("zgsf");
    pm.Push<Hist1D>(Axis(50,15.0,100.0, "photon_pt[0]", "Photon p_{T} [GeV]", {}), 
        photon_sf_baseline, procs, ops).Weight(w_ewkzgamma*w_years).Tag("zgsf");
    pm.Push<Hist1D>(Axis(50,-2.5,2.5, "photon_eta[0]", "Photon |#eta|", {}), 
        photon_sf_baseline, procs, ops).Weight(w_ewkzgamma*w_years).Tag("zgsf");

    //exclusive
    const std::vector<float> idmva_bins = {-1.0,0.2,0.4,0.6,0.7,0.8,0.9,1.0};
    const std::vector<float> pt_bins = {15.0,20.0,35.0,50.0,200.0}; 
    //not sure if we have stats for last two bins or if it even matters
    const std::vector<float> abseta_bins = {0.0,0.8,1.4442,2.0,2.5}; //we veto the gap
    for (unsigned imva = 0; imva < (idmva_bins.size()-1); imva++) {
      pm.Push<Hist1D>(Axis(30,70.0,100.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
          photon_sf_baseline&&(phmva>idmva_bins[imva])&&(phmva<idmva_bins[imva+1]),
          procs, ops).Weight(w_ewkzgamma*w_years).Tag("zgsf");
      for (unsigned ipt = 0; ipt < (pt_bins.size()-1); ipt++) {
        for (unsigned ieta = 0; ieta < (abseta_bins.size()-1); ieta++) {
          pm.Push<Hist1D>(Axis(30,70.0,100.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), 
              photon_sf_baseline&&(phmva>idmva_bins[imva])&&(phmva<idmva_bins[imva+1])&&
              (phpt>pt_bins[ipt])&&(phpt<pt_bins[ipt+1])&&
              (phabseta>abseta_bins[ieta])&&(phabseta<abseta_bins[ieta+1]),
              procs, ops).Weight(w_ewkzgamma*w_years).Tag("FixName:zgsf__llphoton_m0__photon_sf_region_idmvabin"+std::to_string(imva)+"_ptbin"+std::to_string(ipt)+"_etabin"+std::to_string(ieta));
        }
      }
    }
  }

  pm.min_print_ = true;
  pm.SetLuminosityTag("138").MakePlots(1.0);

  //------------------------------------------------------------------------------------
  //                                   postprocessing
  //------------------------------------------------------------------------------------
  if (calc_muon_eff) {
    Hist2D* muon_eff_den_hist = static_cast<Hist2D*>(pm.GetFigure("muon_eff_den").get());
    Hist2D* muon_eff_num_hist = static_cast<Hist2D*>(pm.GetFigure("muon_eff_num").get());
    TH2D muon_eff_den = muon_eff_den_hist->GetBkgHist(true);
    TH2D muon_eff_num = muon_eff_num_hist->GetBkgHist(true);
    muon_eff_num.Divide(&muon_eff_den);
    TCanvas c;
    muon_eff_num.Draw();
    c.Print("plots/muon_efficiency.pdf");
    TFile* out_file = TFile::Open("muon_efficiency.root","RECREATE");
    out_file->WriteObject(&muon_eff_num, "muon_efficiency");
    out_file->Write();
    out_file->Close();
  }

  return 0;
}
