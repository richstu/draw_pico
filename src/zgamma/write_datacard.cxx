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
  kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_kinbdt_masscut_run2_BDT.weights.xml");
  //kin_bdt_reader.BookMVA("/homes/abarzdukas/ZGamma/MVA_Training/IDMVA_Training/test/dataset/weights/TMVAClassification_TightIDMVA_mll_lDR_mreW_BDTG.weights.xml");
  NamedFunc bdt_score = kin_bdt_reader.GetDiscriminant();

  MVAWrapper ph_bdt_reader("photon_bdt");
  ph_bdt_reader.SetVariable("photon_ptransverse","photon_pt[0]");
  ph_bdt_reader.SetVariable("photon_rapidity","photon_eta[0]");
  ph_bdt_reader.SetVariable("ph_irel","photon_reliso[0]");
  ph_bdt_reader.SetVariable("ph_r9","photon_r9[0]");
  ph_bdt_reader.SetVariable("ph_sieie","photon_sieie[0]");
  ph_bdt_reader.SetVariable("ph_hoe","photon_hoe[0]");
  ph_bdt_reader.SetVariable("photon_mva","photon_idmva[0]");
  ph_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_phtree_ph_BDT.weights.xml");
  NamedFunc ph_bdt_score = ph_bdt_reader.GetDiscriminant();

  MVAWrapper ph_kin_bdt_reader("photon_kinematic_bdt");
  ph_kin_bdt_reader.SetVariable("photon_newmva",ph_bdt_score);
  ph_kin_bdt_reader.SetVariable("pt_mass","llphoton_pt[0]/llphoton_m[0]");
  ph_kin_bdt_reader.SetVariable("min_dR","photon_drmin[0]");
  ph_kin_bdt_reader.SetVariable("max_dR",photon_drmax);
  ph_kin_bdt_reader.SetVariable("cosTheta","llphoton_cosTheta[0]");
  ph_kin_bdt_reader.SetVariable("costheta","llphoton_costheta[0]");
  ph_kin_bdt_reader.SetVariable("phi","llphoton_psi[0]");
  ph_kin_bdt_reader.BookMVA("/homes/oshiro/analysis/small_phys_utils/dataset/weights/shuffled_phtree_bdt_comb_BDT.weights.xml");
  NamedFunc phkin_bdt_score = ph_kin_bdt_reader.GetDiscriminant();

  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------

  gErrorIgnoreLevel = 6000;
  Process::Type back =  Process::Type::background;
  Process::Type sig =  Process::Type::signal;
  Process::Type data =  Process::Type::data;

  string prod_folder("/net/cms17/cms17r0/pico/NanoAODv9/htozgamma_deathvalley_v2/");
  std::set<int> years = {2016,2017,2018};
  std::string lumi_string = "340";

  NamedFunc weight = "w_lumi"*w_years*w_run3;
  Palette mc_colors("txt/colors_zgamma.txt","default");

  auto proc_dy          = Process::MakeShared<Baby_pico>("Z+Fake Photon", back, mc_colors("dyjets"),
                             attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*DYJets*"}), 
                             trig&&"stitch_dy");
  auto proc_tt          = Process::MakeShared<Baby_pico>("Other", back, mc_colors("tt"),
                             attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*TTTo2L2Nu*"}), 
                             trig);
  auto proc_smzg        = Process::MakeShared<Baby_pico>("Z+#gamma", back, mc_colors("zgtollg"),
                             attach_folder(prod_folder,years,"mc/merged_zgmc_llg",{"*ZGToLLG*"}), 
                             trig);
  //auto proc_bkg         = Process::MakeShared<Baby_pico>("Background", back, mc_colors("zgtollg"),
  //                          attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
  //                          {"*ZGToLLG*","*DYJets*","*TTTo2L2Nu*"}), 
  //                          trig&&stitch);
  auto proc_hzg         = Process::MakeShared<Baby_pico>("H#rightarrow Z#gamma", sig, kRed,
                             attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                             {"*HToZG*M-125*","*HToZG_M125*"}), trig);
  auto proc_pseudodata  = Process::MakeShared<Baby_pico>("pseudodata", data, kBlack,
                             attach_folder(prod_folder,years,"mc/merged_zgmc_llg",
                             {"*ZGToLLG*","*DYJets*","*TTTo2L2Nu*","*HToZG*M-125*","*HToZG_M125*"}), 
                             trig&&stitch);

  //need to get the histograms out, so no SampleLoader for now
  //vector<shared_ptr<Process>> procs = ZgUtilities::ZgSampleLoader().GetSamples("txt/zg_samples.txt","FitStudies");
  vector<shared_ptr<Process>> procs = {proc_dy, proc_tt, proc_smzg, proc_hzg, proc_pseudodata};
  vector<shared_ptr<Process>> procs_nodata = {proc_dy, proc_tt, proc_smzg, proc_hzg};
  vector<shared_ptr<Process>> procs_bak = {proc_dy, proc_tt, proc_smzg};

  std::vector<PlotOpt> ops = {PlotOpt("txt/plot_styles.txt","MassFitPlot")};
  std::vector<PlotOpt> ops_nobottom = {PlotOpt("txt/plot_styles.txt","MassFitPlot").Bottom(BottomType::off)};
  std::vector<PlotOpt> ops_2d = {PlotOpt("txt/plot_styles.txt","Scatter")};
  std::vector<PlotOpt> ops_shapes = {PlotOpt("txt/plot_styles.txt","Shapes")};

  //------------------------------------------------------------------------------------
  //                                   plots and tables
  //------------------------------------------------------------------------------------
  
  std::string version_num_str = "10";
  int nbins = 4;
  
  PlotMaker pm;

  //BDT discriminant
  //pm.Push<Hist1D>(Axis(120, -1.1, 1.1, bdt_score, "BDT Discriminant", {}),
  //      zg_baseline&&"photon_idmva[0]>0.5",
  //      procs, ops_shapes).Weight(weight).Tag("FixName:zgdatacard_bdt_ver"+version_num_str);
  //pm.Push<Hist2D>(Axis(20, 100.0, 160.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
  //      Axis(120, -1.1, 1.1, bdt_score, "BDT Discriminant", {}),
  //      zg_baseline&&"photon_idmva[0]>0.5",
  //      procs_bak, ops_2d).Weight(weight).Tag("FixName:zgdatacard_bdt_mll_corr_ver"+version_num_str);
  //histograms for datacard
  //pm.Push<Hist1D>(Axis(70, 100.0, 170.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
  //      NamedFunc(zg_baseline&&"photon_idmva[0]>0.5").Name("modified_baseline"),
  //      //zg_baseline&&ggf_only&&"photon_idmva[0]>0.5"&&bdt_score<-0.12,
  //      procs_nodata, ops_nobottom).Weight(weight).Tag("FixName:zgdatacardhist_thesis_ver"+version_num_str);
  //pm.Push<Hist1D>(Axis(70, 100.0, 170.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
  //      zg_baseline&&"photon_idmva[0]>0.5",
  //      //zg_baseline&&ggf_only&&"photon_idmva[0]>0.5"&&bdt_score<-0.12,
  //      procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin0_ver"+version_num_str);
  pm.Push<Hist1D>(Axis(70, 100.0, 170.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        zg_baseline&&bdt_score<0.02,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin0_ver"+version_num_str);
  pm.Push<Hist1D>(Axis(70, 100.0, 170.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        zg_baseline&&ggf_only&&bdt_score>0.02&&bdt_score<0.28,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin1_ver"+version_num_str);
  pm.Push<Hist1D>(Axis(70, 100.0, 170.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        zg_baseline&&ggf_only&&bdt_score>0.28&&bdt_score<0.52,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin2_ver"+version_num_str);
  pm.Push<Hist1D>(Axis(70, 100.0, 170.0, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
        zg_baseline&&ggf_only&&bdt_score>0.52,
        procs, ops).Weight(weight).Tag("FixName:zgdatacardhist_bdtbin3_ver"+version_num_str);
  
  pm.multithreaded_ = false; //BDTs only handle single threads
  pm.min_print_ = true;
  pm.SetLuminosityTag(lumi_string).MakePlots(1.);

  //unique_ptr's were confusing ROOT so back to ordinary pointers it is
  TFile* out_file = TFile::Open(("test_templatefit_ver"+version_num_str+".root").c_str(),"RECREATE");
  std::vector<float> dat_yields;
  std::vector<float> sig_yields;
  std::vector<float> bak_yields;
  for (int ibin = 0; ibin < nbins; ibin++) {
    std::string bin_num = std::to_string(ibin);
    Hist1D* metnb_hist1d = static_cast<Hist1D*>(pm.GetFigure("zgdatacardhist_bdtbin"+bin_num+"_ver"+version_num_str).get());
    //note: the histograms from GetComponent->scaled_hist_ are stacked!
    //this means that the proc_dy TH1D has ZG and TTbar added to it
    TH1D* sig_hist = new TH1D((static_cast<Hist1D::SingleHist1D*>(metnb_hist1d->GetComponent(proc_hzg.get())))->scaled_hist_);
    TH1D* bak_hist = new TH1D((static_cast<Hist1D::SingleHist1D*>(metnb_hist1d->GetComponent(proc_dy.get())))->scaled_hist_);
    bak_hist->Smooth(8);
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

  std::cout << "Writing datacard out" << std::endl;
  ofstream datacard_file;
  datacard_file.open(("test_datacard_ver"+version_num_str+".txt").c_str(),ios::out);
  datacard_file << "max  1  number of categories\n";
  datacard_file << "jmax 1  number of samples minus one\n";
  datacard_file << "kmax *  number of nuisance parameters\n\n";
  datacard_file << "shapes * * test_templatefit_ver" << version_num_str << ".root $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n\n";
  datacard_file << std::right << std::setw(26) << "bin";
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << "bin"+std::to_string(ibin);
  datacard_file << "\n";
  datacard_file << std::right << std::setw(26) << "Observation" ;
  for (int ibin = 0; ibin < nbins; ibin++)
    datacard_file << std::right << std::setw(26) << static_cast<int>(dat_yields[ibin]+0.5);
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
