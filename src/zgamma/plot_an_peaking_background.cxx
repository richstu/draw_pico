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
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;
using namespace CatUtilities;

//--This is a skeleton file meant to be used to easily make plotting files

  std::unordered_map<std::string,double>  year_lumi{{"2016",19.5},{"2016APV",16.8},{"2017",41.48},{"2018",59.83}};

int ph_sel_wout_iselidx_and_dr(const Baby &b){
  int idx_keep = -1;
  double ph_pt = -1;
  for(unsigned int idx_ph = 0; idx_ph < b.photon_pt()->size(); idx_ph++){
    double ph_pt_temp = b.photon_pt()->at(idx_ph);
    if(!(b.photon_id80() -> at(idx_ph)) || !(ph_pt_temp>15) || !(b.photon_elveto() -> at(idx_ph))){ continue;}
    if(fabs(b.photon_eta()->at(idx_ph)) > 2.5 || (fabs(b.photon_eta()->at(idx_ph)) < 1.566 && fabs(b.photon_eta()->at(idx_ph)) > 1.4442)){continue;}
    if(ph_pt_temp > ph_pt){ ph_pt = ph_pt_temp; idx_keep = static_cast<int>(idx_ph); } 
  }
  return idx_keep;
}

NamedFunc Nphold_g1("Nphold_g1",[](const Baby &b) -> NamedFunc::ScalarType{ return (ph_sel_wout_iselidx_and_dr(b)>-1); });

NamedFunc photon_drmin_old("photon_drmin_old",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_drmin() -> at( ph_sel_wout_iselidx_and_dr(b) ); });
NamedFunc photon_pt_old("photon_pt_old",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt() -> at( ph_sel_wout_iselidx_and_dr(b) ); });
NamedFunc photon_idmva_old("photon_idmva_old",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_idmva() -> at( ph_sel_wout_iselidx_and_dr(b) ); });
NamedFunc photon_eta_old("photon_eta_old",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_eta() -> at( ph_sel_wout_iselidx_and_dr(b) ); });

NamedFunc llphoton_m_old("llphoton_m_old",[](const Baby &b) -> NamedFunc::ScalarType{
  unsigned int ph_idx = ph_sel_wout_iselidx_and_dr(b);
  TLorentzVector gamma; gamma.SetPtEtaPhiM(b.photon_pt()->at(ph_idx),b.photon_eta() -> at(ph_idx),b.photon_phi()->at(ph_idx),0);
  TLorentzVector ll = AssignZ(b);
  return (ll+gamma).M(); 
});

NamedFunc ptom_sel_old("ptom_sel_old",[](const Baby &b) -> NamedFunc::ScalarType{
  unsigned int ph_idx = ph_sel_wout_iselidx_and_dr(b);
  TLorentzVector gamma; gamma.SetPtEtaPhiM(b.photon_pt()->at(ph_idx),b.photon_eta() -> at(ph_idx),b.photon_phi()->at(ph_idx),0);
  TLorentzVector ll = AssignZ(b);
  return (b.photon_pt() -> at(ph_sel_wout_iselidx_and_dr(b)) )/(ll+gamma).M(); 

});

NamedFunc mll_mllphoton_old("mll_mllphoton_old",[](const Baby &b) -> NamedFunc::ScalarType{
  unsigned int ph_idx = ph_sel_wout_iselidx_and_dr(b);
  TLorentzVector gamma; gamma.SetPtEtaPhiM(b.photon_pt()->at(ph_idx),b.photon_eta() -> at(ph_idx),b.photon_phi()->at(ph_idx),0);
  TLorentzVector ll = AssignZ(b);
  return (ll.M() + (ll+gamma).M()); 
});

int main() {
//int main(int argc, char *argv[]) {
  //Uses txt/zg_samples.txt to define needed processes

  vector<shared_ptr<Process>> procs    = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt", "AllSkimLLMoreTrigs",10);
  vector<shared_ptr<Process>> procs_mc = ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt", "AllSkimLLMoreTrigsData",10);

  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");


  //Defines the plot options for 1D and 2D plots
  PlotOpt style("txt/plot_styles.txt","Eff2D");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .YTitleOffset(1.)
                                     .Overflow(OverflowType::both)
                                     .LabelSize(0.04)
                                     .UseCMYK(true)  
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
  vector<PlotOpt> ops_2D = {bkg_hist};

  PlotOpt lin_lumi("txt/plot_styles.txt","CMSPaper");
  lin_lumi.YAxis(YAxisType::linear)
          .Stack(StackType::signal_overlay)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .Title(TitleType::info)
          .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_lumi()};

  PlotOpt an_plot("txt/plot_styles.txt","CMSPaper");
  vector<PlotOpt> cms_paper = {an_plot};

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline; 
  NamedFunc old_tight_baseline = Nphold_g1 && "nll>0" && pass_trigs_and_pt && ptom_sel_old > 15.0/110.0 && "ll_m[0] > 80 && ll_m[0] < 100" && llphoton_m_old > 100 && llphoton_m_old < 180 && mll_mllphoton_old > 185 && "pass";
 
  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "";
  NamedFunc selection = "1";
  NamedFunc not_hwin = "llphoton_m[0] < 120 || llphoton_m[0] > 130";
  NamedFunc old_not_hwin = llphoton_m_old < 120 || llphoton_m_old > 130;

  //Make the plots with only MC
  selection  = old_tight_baseline;
  pm.Push<Hist1D>(Axis(80, 100, 180, llphoton_m_old, "m_{ll#gamma} [GeV]", {}), selection, procs_mc, ops).Weight(wgt).Tag("ShortName:an_peakingbackground_wpkgbkg_llphoton_m");
  pm.Push<Hist2D>(
    Axis(50, 0,   0.5,  photon_drmin_old, "#Delta R(l,#gamma) [GeV]", {}),
    Axis(40, 100, 180,  llphoton_m_old, "m_{ll#gamma} [GeV]", {}),
    selection, procs_mc, ops_2D).Tag("ShortName:an_peakingbackground_wpkgbkg_mll_mlly");

  selection = old_tight_baseline && old_not_hwin;
  pm.Push<Hist1D>(Axis(80, 100, 180, llphoton_m_old, "m_{ll#gamma} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_peakingbackground_wpkgbkg_llphoton_m_data");
  pm.Push<Hist2D>(
    Axis(50, 0,   0.5,  photon_drmin_old, "#Delta R(l,#gamma) [GeV]", {}),
    Axis(40, 100, 180,  llphoton_m_old, "m_{ll#gamma} [GeV]", {}),
    selection, procs, ops_2D).Tag("ShortName:an_peakingbackground_wpkgbkg_drmin_mlly_data");

  selection  = tight_baseline;
  pm.Push<Hist1D>(Axis(80, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), selection, procs_mc, ops).Weight(wgt).Tag("ShortName:an_peakingbackground_removed_llphoton_m");  
  pm.Push<Hist2D>(
    Axis(50, 0,   0.5, "photon_drmin[0]", "#Delta R(l,#gamma) [GeV]", {}),
    Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
    selection, procs_mc, ops_2D).Tag("ShortName:an_peakingbackground_removed_drmin_mlly");

  selection  = tight_baseline && not_hwin;
  pm.Push<Hist1D>(Axis(80, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ShortName:an_peakingbackground_removed_llphoton_m_data");  
  pm.Push<Hist2D>(
    Axis(50, 0,   0.5, "photon_drmin[0]", "#Delta R(l,#gamma) [GeV]", {}),
    Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
    selection, procs, ops_2D).Tag("ShortName:an_peakingbackground_removed_drmin_mlly_data");


  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("137.61").MakePlots(1);

}

