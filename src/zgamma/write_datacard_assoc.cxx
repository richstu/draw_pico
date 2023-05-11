/**
 * This script generates a histograms needed to build a H->Zgamma datacard
 */

#include "core/test.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"

#include "core/baby.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/mva_wrapper.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/table.hpp"
#include "core/utilities.hpp"
#include "zgamma/apply_zg_trigeffs.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

using namespace ZgFunctions;
using namespace std;
using namespace PlotOptTypes;

//int main(int argc, char *argv[]){
int main(){
  //------------------------------------------------------------------------------------
  //                                   NamedFuncs
  //------------------------------------------------------------------------------------
  
  const NamedFunc w_sigx20("w_sigx20",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 200000 && b.type() <= 205000)
      return 20.0;
    return 1.0;
  });

  const NamedFunc ggf_only = NamedFunc("nlep==2&&njet<2&&nbdfm==0&&met<80").Name("ggf_only");

  const NamedFunc trig = HLT_pass_dilepton; //||HLT_pass_singlelepton;

  MVAWrapper kin_bdt_reader("kinematic_bdt");
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
  //kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_idmvacut_run2_BDT.weights.xml");
  kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_idmvacut_run2_BDT.weights.xml");
  //kin_bdt_reader.BookMVA("/homes/abarzdukas/ZGamma/MVA_Training/IDMVA_Training/test/dataset/weights/TMVAClassification_TightIDMVA_mll_lDR_mreW_BDTG.weights.xml");
  NamedFunc bdt_score = kin_bdt_reader.GetDiscriminant();

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

  //Dphi between H and MET
  const NamedFunc dphi_h_met("dphi_h_met",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.llphoton_phi()->size()==0) return 0;
    return deltaPhi(b.llphoton_phi()->at(0),b.met_phi());
  });

  const NamedFunc tighter_baseline = NamedFunc(zg_baseline&&"photon_idmva[0]>0.5&&ll_m[0]>80&&ll_m[0]<100").Name("new_baseline");

  const NamedFunc tthlep_selections = NamedFunc("njet>=2"&&((("ht>=200||njet>=4"||llphoton_rel_pt>1.0)&&min_lep_id>=2.0&&max_lep_miniso<0.15&&min_lep_pt>12.0)||("nlep>=4&&ht>120"&&max_lep_miniso<0.2))).Name("ttH_leptonic_selections");

  const NamedFunc vhlep_selections = NamedFunc("photon_drmin[0]<1.5"&&min_lep_pt>12.0&&min_lep_id>=2.0&&max_lep_miniso<0.2&&("met>30"||assoc_lep_mt>50)).Name("VH_leptonic_selections");

  const NamedFunc tthhad_selections = NamedFunc("njet>=5&&photon_drmin[0]<1.2&&met<120").Name("ttH_hadronic_selections");

  const NamedFunc vhmet_selections =  NamedFunc("photon_idmva[0]>0.55&&photon_drmin[0]<1.5"&&llphoton_rel_pt>0.7&&dphi_h_met>1.5).Name("VH_MET_selections");

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------

  gErrorIgnoreLevel = 6000;
  std::string lumi_string = "340";
  NamedFunc weight = "w_lumi"*w_years*w_run3;

  //auto proc_pseudodata  = Process::MakeShared<Baby_pico>("pseudodata", data, kBlack,
  //                           attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                           {"*ZGToLLG*","*DYJets*","*TTTo2L2Nu*","*HToZG*M-125*","*HToZG_M125*"}), 
  //                           trig&&stitch);

  //need to get the histograms out, so no SampleLoader for now
  //
  std::vector<std::shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","All");

  std::vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","LinLumi").Overflow(OverflowType::none)}; 
  std::vector<PlotOpt> ops_2d = {PlotOpt("txt/plot_styles.txt","Scatter")};
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","Shapes")};

  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------
  
  PlotMaker pm;

  //BDT discriminant
  //pm.Push<Hist1D>(Axis(120, -1.1, 1.1, bdt_score, "BDT Discriminant", {}),
  //      zg_baseline&&"photon_idmva[0]>0.5",
  //      procs, ops_shapes).Weight(weight).Tag("FixName:zgdatacard_bdt_ver8");
  //pm.Push<Hist2D>(Axis(20, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
  //      Axis(120, -1.1, 1.1, bdt_score, "BDT Discriminant", {}),
  //      zg_baseline&&"photon_idmva[0]>0.5",
  //      procs_bak, ops_2d).Weight(weight).Tag("FixName:zgdatacard_bdt_mll_corr_ver8");
  //histograms for datacard
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep>=3&&nbdfl>=1"&&tthlep_selections,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin0_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep>=3&&nbdfl==0"&&vhlep_selections,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin1_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm>=1"&&tthhad_selections,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin2_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met>=95"&&vhmet_selections,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin3_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&(njet>=2&&dijet_m>=750)"
        &&bdt_score<0.16,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin4_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&(njet>=2&&dijet_m>=750)"
        &&bdt_score>0.16&&bdt_score<0.42,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin5_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&(njet>=2&&dijet_m>=750)"
        &&bdt_score>0.42,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin6_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&!(njet>=2&&dijet_m>=750)"
        &&bdt_score<-0.12,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin7_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&!(njet>=2&&dijet_m>=750)"
        &&bdt_score>-0.12&&bdt_score<0.1,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin8_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&!(njet>=2&&dijet_m>=750)"
        &&bdt_score>0.1&&bdt_score<0.36,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin9_ver9");
  pm.Push<Hist1D>(Axis(30, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        tighter_baseline&&"nlep==2&&nbdfm==0&&met<95&&!(njet>=2&&dijet_m>=750)"
        &&bdt_score>0.36,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin10_ver9");
  
  pm.multithreaded_ = false; //BDTs only handle single threads
  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_string).MakePlots(1.);

  //unique_ptr's were confusing ROOT so back to ordinary pointers it is
  TFile* out_file = TFile::Open("test_templatefit_ver9.root","RECREATE");
  std::vector<float> dat_yields;
  std::vector<float> sig_yields;
  std::vector<float> bak_yields;
  int nbins = 11;
  for (int ibin = 0; ibin < nbins; ibin++) {
    std::string bin_num = std::to_string(ibin);
    Hist1D* metnb_hist1d = static_cast<Hist1D*>(pm.GetFigure("zgdatacardhist_bdtbin"+bin_num+"_ver9").get());
    //note: the histograms from GetComponent->scaled_hist_ are stacked!
    //this means that the proc_dy TH1D has ZG and TTbar added to it
    TH1D* sig_hist = new TH1D((static_cast<Hist1D::SingleHist1D*>(metnb_hist1d->GetComponent(procs[7].get())))->scaled_hist_);
    TH1D* bak_hist = new TH1D((static_cast<Hist1D::SingleHist1D*>(metnb_hist1d->GetComponent(procs[0].get())))->scaled_hist_);
    if (ibin==0) //ttH lep stats are terrible
      bak_hist->Smooth(20);
    else
      bak_hist->Smooth(12);
    TH1D* dat_hist = static_cast<TH1D*>(bak_hist->Clone("dat_hist"));
    dat_hist->Add(sig_hist);
    sig_yields.push_back(sig_hist->Integral());
    bak_yields.push_back(bak_hist->Integral());
    dat_yields.push_back(dat_hist->Integral());
    TDirectory * bin0_dir = out_file->mkdir(("bin"+bin_num).c_str());
    bin0_dir->WriteObject(sig_hist, "sig");
    bin0_dir->WriteObject(bak_hist, "bkg");
    bin0_dir->WriteObject(dat_hist, "data_obs");
    delete sig_hist;
    delete bak_hist;
    delete dat_hist;
  }
  out_file->Write();

  ofstream datacard_file;
  datacard_file.open("test_datacard_ver9.txt",ios::out);
  datacard_file << "max  1  number of categories\n";
  datacard_file << "jmax 1  number of samples minus one\n";
  datacard_file << "kmax *  number of nuisance parameters\n\n";
  datacard_file << "shapes * * test_templatefit_ver9.root $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n\n";
  datacard_file << std::right << std::setw(26) << "bin";
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << "bin"+std::to_string(ibin);
  datacard_file << "\n";
  datacard_file << std::right << std::setw(26) << "Observation" ;
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << dat_yields[ibin];
  datacard_file << "\n\n";
  datacard_file << std::right << std::setw(26) << "bin";
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << "bin"+std::to_string(ibin)
                  << std::right << std::setw(26) << "bin"+std::to_string(ibin);
  datacard_file << "\n";
  datacard_file << std::right << std::setw(26) << "process";
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << "sig"
                  << std::right << std::setw(26) << "bkg";
  datacard_file << "\n";
  datacard_file << std::right << std::setw(26) << "process";
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << "0"
                  << std::right << std::setw(26) << "1";
  datacard_file << "\n";
  datacard_file << std::right << std::setw(26) << "rate";
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << sig_yields[ibin]
                  << std::right << std::setw(26) << bak_yields[ibin];
  datacard_file << "\n";
  datacard_file.close();

}
