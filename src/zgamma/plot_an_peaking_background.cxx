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


  int old_photon_reconstruction(const Baby &b){
    int idx_keep= -1; double max_ph_pt = 0;
    for(unsigned int idx_i = 0; idx_i < b.photon_pt() -> size(); idx_i++){
      if(!(b.photon_id80() -> at(idx_i)) || !(b.photon_pt()->at(idx_i)>15) || !(b.photon_elveto() -> at(idx_i))){ continue;}
      if(b.photon_drmin() -> at(idx_i) < 0.4){continue;}
      if(fabs(b.photon_eta()->at(idx_i)) > 2.5 || (fabs(b.photon_eta()->at(idx_i)) < 1.566 && fabs(b.photon_eta()->at(idx_i)) > 1.4442)){continue;}
      if(b.photon_pt() -> at(idx_i) > max_ph_pt){
        idx_keep = static_cast<int>(idx_i);
        max_ph_pt  = b.photon_pt() -> at(idx_i);
      }
    }
    return idx_keep;
  }

  NamedFunc mlly_oldreco("mlly_oldreco",[](const Baby &b) -> NamedFunc::ScalarType{
    int idx_ph = old_photon_reconstruction(b);
    TLorentzVector ll,gamma;
    ll.SetPtEtaPhiM(b.ll_pt() ->at(0), b.ll_eta() ->at(0), b.ll_phi() ->at(0), b.ll_m() ->at(0));
    gamma.SetPtEtaPhiM(b.photon_pt() ->at(idx_ph), b.photon_eta() ->at(idx_ph), b.photon_phi() ->at(idx_ph), 0);

    return (ll+gamma).M(); 
  });

  NamedFunc baseline_old_ph_reco("baseline_old_ph_reco",[](const Baby &b) -> NamedFunc::ScalarType{ 
    int idx_ph = old_photon_reconstruction(b);
    if(idx_ph<0){return false;}
    if(b.nlep()<2){return false;}
    TLorentzVector ll,gamma;
    ll.SetPtEtaPhiM(b.ll_pt() ->at(0), b.ll_eta() ->at(0), b.ll_phi() ->at(0), b.ll_m() ->at(0));
    gamma.SetPtEtaPhiM(b.photon_pt() ->at(idx_ph), b.photon_eta() ->at(idx_ph), b.photon_phi() ->at(idx_ph), 0);
    double mlly = (ll+gamma).M();

    if(ll.M() < 80 || ll.M() > 100){ return false; }
    if(gamma.Pt()/mlly < 15/110.0){ return false; }
    if(mlly<100 || mlly>180){return false;}
    if(mlly + ll.M() < 185){return false;}
    return true;
  });



  NamedFunc hasNanoPhoton("hasNanoPhoton",[](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt() -> size() > 0; });
  NamedFunc isPhotonUnique("isPhotonUnique", [](const Baby &b) -> NamedFunc::ScalarType{
    int ph_elidx = b.photon_elidx() ->at(0);
    if( ((ph_elidx == b.ll_i1()->at(0)) || (ph_elidx == b.ll_i2()->at(0))) && (b.ll_lepid() -> at(0) ==11) ){return false;}
    return true;
  });

   NamedFunc ph_un_nsig("ph_un_nsig", [](const Baby &b) -> NamedFunc::ScalarType{
    int ph_elidx = b.photon_elidx() ->at(0);
    if(b.type() < 200000 && ((ph_elidx == b.ll_i1()->at(0)) || (ph_elidx == b.ll_i2()->at(0))) && (b.ll_lepid() -> at(0) ==11) ){return true;}
    else if(b.type() >= 200000){ return true; }
    return false;
  });

  
  NamedFunc mlly_unique("mlly_unique", [](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector ph(0,0,0,0);
    TLorentzVector ll; ll.SetPtEtaPhiM(b.ll_pt() -> at(0), b.ll_eta()->at(0), b.ll_phi() -> at(0), b.ll_m() -> at(0) );
    int ph_elidx = -1;
    for(unsigned int idx_ph=0; idx_ph < b.photon_pt() -> size(); idx_ph++){
      ph_elidx = b.photon_elidx() -> at(idx_ph);
      if( ((ph_elidx == b.ll_i1()->at(0)) || (ph_elidx == b.ll_i2()->at(0))) && (b.ll_lepid() -> at(0) ==11) ){continue;}
      if(!(b.photon_id80() -> at(idx_ph)) || (b.photon_pt() ->at(idx_ph) < 15)){ continue; }
      ph.SetPtEtaPhiM(b.photon_pt()->at(idx_ph),b.photon_eta() ->at(idx_ph), b.photon_phi() ->at(idx_ph),0);
    }
    if(ph.Pt()/(ll+ph).M() < 15/110){return 0;}
    if(ph.Pt()<15){return 0;}
    
    return (ll+ph).M();
  });


  NamedFunc drmin_elType("drmin_elType",[](const Baby &b) -> NamedFunc::ScalarType{
    double min_dr = 100;
    int min_el = 100;
    double temp_dr;

    if(b.ll_lepid()->at(0)==13){return -2;}

    TLorentzVector ph;ph.SetPtEtaPhiM(b.photon_pt() -> at(0),b.photon_eta() -> at(0),b.photon_phi() -> at(0),0);
    TLorentzVector el_t;
    for(int idx_el = 0; idx_el < b.nel(); idx_el++){
      el_t.SetPtEtaPhiM(b.el_pt() -> at(idx_el),b.el_eta() -> at(idx_el),b.el_phi() -> at(idx_el),0);
      temp_dr = el_t.DeltaR(ph);
      min_el = temp_dr < min_dr ?  idx_el : min_el;
      min_dr = temp_dr < min_dr ?  temp_dr : min_dr;
    }
    if(min_el==100){return -2;}

    return b.el_ecal() -> at(min_el);
  });

  
  NamedFunc el_matched("el_matched",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TLorentzVector ph,l1;
    ph.SetPtEtaPhiM(b.photon_pt() -> at(0), b.photon_eta() -> at(0), b.photon_phi() -> at(0), b.photon_pt() -> at(0));
    bool is_el = false;
    for(size_t idx_mc = 0; idx_mc < b.mc_pt() -> size(); idx_mc++){
      if( abs(b.mc_id() -> at(idx_mc)) == 11 ){ continue; }
      l1.SetPtEtaPhiM(b.mc_pt() -> at(idx_mc), b.mc_eta() -> at(idx_mc), b.mc_phi() -> at(idx_mc), b.mc_mass() -> at(idx_mc) );
      if(l1.DeltaR(ph) < 0.05) { is_el = true; }
    }
    return is_el;
  });

  NamedFunc jet_matched("jet_matched",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TLorentzVector ph,l1;
    ph.SetPtEtaPhiM(b.photon_pt() -> at(0), b.photon_eta() -> at(0), b.photon_phi() -> at(0), b.photon_pt() -> at(0));
    bool is_jet = false;
    for(size_t idx_mc = 0; idx_mc < b.mc_pt() -> size(); idx_mc++){
      if( abs(b.mc_id() -> at(idx_mc)) < 10 ){ continue; }
      l1.SetPtEtaPhiM(b.mc_pt() -> at(idx_mc), b.mc_eta() -> at(idx_mc), b.mc_phi() -> at(idx_mc), b.mc_mass() -> at(idx_mc) );
      if(l1.DeltaR(ph) < 0.05) { is_jet = true; }
    }
    return is_jet;
  });

  NamedFunc dr_yj("dr_yj",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TLorentzVector ph,j1; double min_dr = 100;
    ph.SetPtEtaPhiM(b.photon_pt() -> at(0), b.photon_eta() -> at(0), b.photon_phi() -> at(0), 0);
    for(size_t idx_jet = 0; idx_jet < b.jet_pt() -> size(); idx_jet++){
      if( !(b.jet_isgood() -> at(idx_jet)) ){ continue; }
      j1.SetPtEtaPhiM(b.jet_pt() -> at(idx_jet), b.jet_eta() -> at(idx_jet), b.jet_phi() -> at(idx_jet), b.jet_m() -> at(idx_jet) );
      if(j1.DeltaR(ph) < min_dr) { min_dr = j1.DeltaR(ph); }
    }
    return min_dr;
  });

  bitset<16> mc_status;
  double dr_ytpart(const Baby &b){ 
    TLorentzVector ph,tpart; double min_dr = 100; int idx_mindr = -1;
    ph.SetPtEtaPhiM(b.photon_pt() -> at(0), b.photon_eta() -> at(0), b.photon_phi() -> at(0), 0);
    for(size_t idx_mc = 0; idx_mc < b.mc_pt() -> size(); idx_mc++){
      mc_status = b.mc_statusflag() -> at(idx_mc);
      if( !(mc_status[0] || mc_status[8]) ){ continue; }
      tpart.SetPtEtaPhiM(b.mc_pt() -> at(idx_mc), b.mc_eta() -> at(idx_mc), b.mc_phi() -> at(idx_mc), b.mc_mass() -> at(idx_mc) );
      if(tpart.DeltaR(ph) < min_dr) { min_dr = tpart.DeltaR(ph); idx_mindr = idx_mc;}
    }
    return idx_mindr;
  }
  

 NamedFunc y_drtm_flav("y_drtm_flav",[](const Baby &b) -> NamedFunc::ScalarType{ 
    return abs(b.mc_id() -> at( dr_ytpart(b) ));
  });


//Currently the skeleton has the year as an input variable, but this can be removed
//This is only used in defining processes
//int main() {
int main() {

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

  //Uses txt/zg_samples.txt to define needed processes
  vector<shared_ptr<Process>> procs =  ZgUtilities::procs_with_sig_scale("txt/samples_zgamma.txt","AllMoreTrigs",10);

  Process::Type back =  Process::Type::background;

  int skim_type = 3; 
  std::string folder="unskimmed"; std::string folder_data=folder;
  if(skim_type==1){folder="merged_zgmc_llg";folder_data="merged_zgdata_llg";}
  else if(skim_type==2){folder="skim_llg";folder_data="skim_llg";}
  else if(skim_type==3){folder="skim_ll";folder_data="skim_ll";}


  std::string bfolder("/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/");
  std::string sample_name = "amcatnloFXFX";

  std::string mc_path_2016( bfolder+"2016/mc/" + folder + "/");
  std::string mc_path_2016APV( bfolder+"2016APV/mc/" + folder + "/");
  std::string mc_path_2017( bfolder+"2017/mc/" + folder + "/");
  std::string mc_path_2018( bfolder+"2018/mc/" + folder + "/");

  std::string data_path_2016( bfolder+"2016/data/" + folder_data  + "/");
  std::string data_path_2016APV( bfolder+"2016APV/data/" + folder_data + "/");
  std::string data_path_2017( bfolder+"2017/data/" + folder_data + "/");
  std::string data_path_2018( bfolder+"2018/data/" + folder_data  + "/");

                                              
  auto proc_DYJets_El  = Process::MakeShared<Baby_pico>("DY + Fake #gamma=e",   back, TColor::GetColor("#bf00ff"), 
                      {mc_path_2016+"*DYJetsToLL*" +sample_name + "*.root",mc_path_2016APV+"*DYJetsToLL*" +sample_name + "*.root",
                       mc_path_2017+"*DYJetsToLL*" +sample_name + "*.root",mc_path_2018+"*DYJetsToLL*" +sample_name + "*.root"}, 
                       "(trig_double_el||trig_double_mu||trig_single_el||trig_single_mu)&&use_event" && hasNanoPhoton && "photon_pflavor[0]==11");

  auto proc_DYJets_ph  = Process::MakeShared<Baby_pico>("DY + Fake: #gamma=#gamma",   back, TColor::GetColor("#f4ffa1"), 
                      {mc_path_2016+"*DYJetsToLL*" +sample_name + "*.root",mc_path_2016APV+"*DYJetsToLL*" +sample_name + "*.root",
                       mc_path_2017+"*DYJetsToLL*" +sample_name + "*.root",mc_path_2018+"*DYJetsToLL*" +sample_name + "*.root"}, 
                       "(trig_double_el||trig_double_mu||trig_single_el||trig_single_mu)&&use_event" && hasNanoPhoton && "photon_pflavor[0]==1");

  auto proc_DYJets_nphel  = Process::MakeShared<Baby_pico>("DY + Fake: #gamma!=e,#gamma",   back, TColor::GetColor("#ee7600"), 
                      {mc_path_2016+"*DYJetsToLL*" +sample_name + "*.root",mc_path_2016APV+"*DYJetsToLL*" +sample_name + "*.root",
                       mc_path_2017+"*DYJetsToLL*" +sample_name + "*.root",mc_path_2018+"*DYJetsToLL*" +sample_name + "*.root"}, 
                       "(trig_double_el||trig_double_mu||trig_single_el||trig_single_mu)&&use_event" && hasNanoPhoton && "photon_pflavor[0]==0");


  vector<shared_ptr<Process>> procs_mc = ZgUtilities::ZgSampleLoader().LoadSamples("txt/samples_zgamma.txt","AllSkimLLTrigsData");
  vector<shared_ptr<Process>> procs_nosig = procs_mc;
  procs_nosig.erase(procs_nosig.end()-1);

  //  procs_mc.erase(procs_mc.begin());
  procs_mc.insert(procs_mc.begin(),proc_DYJets_nphel);
  procs_mc.insert(procs_mc.begin(),proc_DYJets_El);
  procs_mc.insert(procs_mc.begin(),proc_DYJets_ph);

  //This defines the splits based on lepton flavor that will be plotted
  vector<NamedFunc> lep = { "ll_lepid[0]==11","ll_lepid[0]==13","1"};
  vector<string> lep_lab= {"_ee", "_mumu", "_ll"};

  //These are all the selections used as a part of the baseline selection from the zg_functions.cpp file
  NamedFunc tight_baseline = ZgFunctions::tightened_baseline;

 
  //This block of code handles all of the plot making
  PlotMaker pm;
  string pltNames = "";
  NamedFunc sel_cat = "1";


  sel_cat  = tight_baseline;

  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), sel_cat, procs, ops).Weight(wgt).Tag("ShortName:an_peakingbackground_removed_llphoton_m");

  pm.Push<Hist2D>(
    Axis(70,50,120,  "ll_m[0]", "m_{ll} [GeV]", {}),
    Axis(40, 100, 180,  "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
    sel_cat, procs, ops_2D).Tag("ShortName:an_peakingbackground_mll_mlly");
  
  pm.Push<Hist2D>(
    Axis(50, 0,   0.5,  "photon_drmin[0]", "m_{ll} [GeV]", {}),
    Axis(40, 100, 180,  "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
    sel_cat, procs, ops_2D).Tag("ShortName:an_peakingbackground_drmin_mlly");


  //This block controls the verbose or non-verbose plotting option and whether or not to use multiple threads
  pm.min_print_ = true;
  pm.multithreaded_ = true;

  //Luminosity option for the plots
  pm.SetLuminosityTag("137.61").MakePlots(1);

}
