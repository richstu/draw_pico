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

//temp
#include <Math/Rotation3D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/Vector4D.h>
#include <Math/Boost.h>
#include <Math/BoostX.h>
#include <Math/BoostY.h>
#include <Math/BoostZ.h>

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
#include "zgamma/zg_utilities.hpp"
#include "zgamma/zg_functions.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;
using namespace ZgFunctions;

//------------------------------------GIST--------------------------------------//
//This plot file makes a number of plots studying the effects of jet vetos in 2018
//and Run 3. No inputs are taken. 
//------------------------------------------------------------------------------//

//NamedFuncs
                     
NamedFunc cutbitcheck(int bit){//NamedFuncs don't support bitwise operators seemingly
  return NamedFunc("comp",[bit](const Baby &b)->NamedFunc::ScalarType{
    return (bit & (b.zg_cutBitMap()))>0;
  });
}

NamedFunc catbitcheck(int bit){
  return NamedFunc("comp",[bit](const Baby &b)->NamedFunc::ScalarType{
    return (bit == (b.zg_categorizationBitMap()));
  });
}

NamedFunc getJetpT(NamedFunc idx){
  return NamedFunc("pt",[idx](const Baby &b)->NamedFunc::ScalarType{
    return b.jet_pt()->at(idx.GetScalar(b));
  });
}

NamedFunc getJetEta(NamedFunc idx){
  return NamedFunc("eta",[idx](const Baby &b)->NamedFunc::ScalarType{
    return b.jet_eta()->at(idx.GetScalar(b));
  });
}

NamedFunc getJetPhi(NamedFunc idx){
  return NamedFunc("phi",[idx](const Baby &b)->NamedFunc::ScalarType{
    return b.jet_phi()->at(idx.GetScalar(b));
  });
}


NamedFunc sideband(){//using directe.root branches
  return NamedFunc("sigblind",[](const Baby &b)->NamedFunc::ScalarType{
    return (0b000000000001 & (b.zg_cutBitMap()))!=0;
  });
}


NamedFunc Nminusk(vector<NamedFunc> cuts, unsigned int skip){
  NamedFunc total("1");
  for(unsigned int j(0); j < cuts.size(); j++){
    if(j <= skip) total = total && cuts.at(j);
  }
  return total;
}

 NamedFunc good_idx("good_idx",[](const Baby &b) -> NamedFunc::ScalarType{
    int out = -1;
    int count = 0;
    for(unsigned int i(0);i<b.jet_pt()->size();i++){
      if(b.jet_isgood()->at(i) == 1) {
        count+=1;
        if(count > 2) break;
        else if(count == 2){
          out = i;
        }
      }
    }
    return out;
});

NamedFunc good_pt("good_pt",[](const Baby &b) -> NamedFunc::ScalarType{
   return b.jet_pt()->at(good_idx.GetScalar(b));
});

NamedFunc has_jet("has_jet",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_pt()->size()>0;
});


NamedFunc has_1min_jet("has_1min_jet",[](const Baby &b) -> NamedFunc::ScalarType{
  bool out = false;
  if(b.jet_pt()->size()>0){
    for(unsigned int i(0); i<b.jet_pt()->size(); i++){
      if(b.jet_isgood_min()->at(i)==true && b.jet_pt()->at(i)>0) out = true;
      //if((!b.jet_islep()->at(i) && !b.jet_isphoton()->at(i) && (fabs(b.jet_eta()->at(i)) <= 4.7) && b.jet_id()->at(i)) && b.jet_pt()->at(i)>0) out = true;
    }
  }
  return out;
});

NamedFunc has_2min_jet("has_2min_jet",[](const Baby &b) -> NamedFunc::ScalarType{
  bool out = false;
  int count = 0;
  if(b.jet_pt()->size()>=1){
    for(unsigned int i(0); i<b.jet_pt()->size(); i++){
      if(b.jet_isgood_min()->at(i)==true && count == 1) out = true;
      else if(b.jet_isgood_min()->at(i)==true && count == 0) count+=1;
    }
  }
  return out;
});

NamedFunc etaj1("etaj1",[](const Baby &b) -> NamedFunc::ScalarType{//eta distributions without the jet horn fix
  int index = -999;
  double max = 0;
  for (unsigned int i(0); i<b.jet_pt()->size(); i++){
    if(!b.jet_islep()->at(i) && !b.jet_isphoton()->at(i) && (fabs(b.jet_eta()->at(i)) <= 4.7) && b.jet_id()->at(i)){
      if (b.jet_pt()->at(i)>=max){
        max = b.jet_pt()->at(i);
        index = i;
      }
    }
  }
  if (index == -999){
    cout<<b.jet_pt()->size()<<endl;
    for (unsigned int i(0);i<b.jet_pt()->size(); i++){
      cout<<b.jet_pt()->at(i);
    }
    cout<<b.event()<<endl;
    cout<<b.run()<<endl;
  }
  return b.jet_eta()->at(index);
});

NamedFunc etaj2("etaj2",[](const Baby &b) -> NamedFunc::ScalarType{//eta distributions without the jet horn fix
  int index1 = 0;
  int index2 = 0;
  double temp = -1;
  double max = 0;
  for (unsigned int i(0); i<b.jet_pt()->size();i++){
    if(!b.jet_islep()->at(i) && !b.jet_isphoton()->at(i) && (fabs(b.jet_eta()->at(i)) <= 4.7) && b.jet_id()->at(i)){ 
      if(b.jet_pt()->at(i)>max){
        temp = max;
        max = b.jet_pt()->at(i);
        index2 = index1;
        index1 = i;
      }else if(b.jet_pt()->at(i)>=temp){
        temp = b.jet_pt()->at(i);
        index2 = i;
      }
    }
  }
  return b.jet_eta()->at(index2);
});


NamedFunc gl_jet("gl_jet",[](const Baby &b) -> NamedFunc::ScalarType{
  bool out = false;
  if(b.jet_pt()->size()>0) out = b.jet_isgood_min()->at(0); 
  return out;                                                                                                           
});
NamedFunc has_jet2("has_jet2",[](const Baby &b) -> NamedFunc::ScalarType{
  return b.jet_pt()->size()>1;
});
NamedFunc gl_notmap("gl_notmap",[](const Baby &b) -> NamedFunc::ScalarType{
  bool out = false;
  if(b.jet_pt()->size()>0) out = (b.jet_isvetomap()->at(0)==0 && b.jet_isgood_min()->at(0)==1);
  return out;
});
NamedFunc gl_nothem("gl_nothem",[](const Baby &b) -> NamedFunc::ScalarType{
  bool out = false;
  if(b.jet_pt()->size()>0) out = (b.jet_isvetohem()->at(0)==0 && b.jet_isgood_min()->at(0)==1);
  return out;
});


NamedFunc has_clean_jet("has_clean_jet",[](const Baby &b) -> NamedFunc::ScalarType{//cleaned jets without the veto map
  int out = -1;
  for(unsigned int ijet(0); ijet<b.jet_pt()->size();ijet++){
    if(b.jet_isgood()->at(ijet)) return ijet;
  }
  return out;
});

NamedFunc has_clean_jet2("has_clean_jet2",[](const Baby &b) -> NamedFunc::ScalarType{
  int out = -1;
  if(has_clean_jet.GetScalar(b) < 0) return out;
  int count = 0;
  for(unsigned int ijet(0); ijet<b.jet_pt()->size();ijet++){
    if(count == 1 && b.jet_isgood()->at(ijet)) return ijet;
    else if(count == 0 && b.jet_isgood()->at(ijet)) count+=1;
  }
  return out;
});

NamedFunc has_clean_nvjet("has_clean_nvjet",[](const Baby &b) -> NamedFunc::ScalarType{//cleaned jets without the veto map
  int out = -1;
  for(unsigned int ijet(0); ijet<b.jet_pt()->size();ijet++){
    if(b.jet_isgood()->at(ijet) && b.jet_isvetohem()->at(ijet)==0 && b.jet_isvetomap()->at(ijet)==0) return ijet;
  }
  return out;
});

NamedFunc has_clean_nvjet2("has_clean_nvjet2",[](const Baby &b) -> NamedFunc::ScalarType{
  int out = -1;
  if(has_clean_nvjet.GetScalar(b) < 0) return out;
  int count = 0;
  for(unsigned int ijet(0); ijet<b.jet_pt()->size();ijet++){
    if(count == 1 && b.jet_isgood()->at(ijet) && b.jet_isvetohem()->at(ijet)==0 && b.jet_isvetomap()->at(ijet)==0) return ijet;
    else if(count == 0 && b.jet_isgood()->at(ijet) && b.jet_isvetohem()->at(ijet)==0 && b.jet_isvetomap()->at(ijet)==0) count+=1;
  }
  return out;
});

NamedFunc zero = "0";

NamedFunc tester("tester",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.llphoton_deta()->size()<1) cout<<b.llphoton_deta()->size()<<endl;
  return true;
});

NamedFunc hastz("hastz",[](const Baby &b) -> NamedFunc::ScalarType{
  if(b.mc_pt()->size()==0) return false;
  bool foundz = false;
  for(unsigned int imc(0); imc<b.mc_pt()->size(); imc++){
    if(b.mc_id()->at(imc)==23) foundz = true;
  }
  return foundz;
});

NamedFunc tz_phi("tz_phi",[](const Baby &b) -> NamedFunc::ScalarType{
  double phiout = -10.;
  for(unsigned int imc(0); imc<b.mc_pt()->size(); imc++){
    if(b.mc_id()->at(imc)==23) phiout = b.mc_phi()->at(imc);
  }
  return phiout;
});

int main(){
 
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");

//===========================Define Processes=============================//

  auto data_2022_asmc     = Process::MakeShared<Baby_pico>("Data 2022", Process::Type::background, kBlack, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/data/skim_ll/*e.root"},"1");
  auto data_2022EE_asmc   = Process::MakeShared<Baby_pico>("Data 2022EE", Process::Type::background, kBlack, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/data/skim_ll/*e.root"},"1");
  auto data_2023_asmc     = Process::MakeShared<Baby_pico>("Data 2023", Process::Type::background, kBlack, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/data/skim_ll/*e.root"},"1");
  auto data_2023BPix_asmc = Process::MakeShared<Baby_pico>("Data 2023BPix", Process::Type::background, kBlack, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/data/skim_ll/*e.root"},"1");

  //I think we can combine the years for this.
  //auto data_2022
  //auto data_2022EE
  //auto data_2023
  //auto data_2023BPix
  auto data_Run3          = Process::MakeShared<Baby_pico>("Data Run3", Process::Type::data, kBlack, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/data/skim_ll/*e.root"},"1");

  auto dyj_2022           = Process::MakeShared<Baby_pico>("Z+FakePhoton 2022", Process::Type::background, TColor::GetColor("#FFA90E"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV*e.root"},"1");
  auto dyj_2022EE         = Process::MakeShared<Baby_pico>("Z+FakePhoton 2022EE", Process::Type::background, TColor::GetColor("#FFA90E"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/skim_ll/*DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV*e.root"},"1");
  auto dyj_2023           = Process::MakeShared<Baby_pico>("Z+FakePhoton 2023", Process::Type::background, TColor::GetColor("#FFA90E"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/skim_ll/*DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV*e.root"},"1");
  auto dyj_2023BPix       = Process::MakeShared<Baby_pico>("Z+FakePhoton 2023BPix", Process::Type::background, TColor::GetColor("#FFA90E"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/skim_ll/*DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV*e.root"},"1");
 
  auto dyj_Run3           = Process::MakeShared<Baby_pico>("Z+FakePhoton Run 3", Process::Type::background, TColor::GetColor("#FFA90E"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV*e.root"},"1");

  auto smzy_2022          = Process::MakeShared<Baby_pico>("Z+#gamma 2022", Process::Type::background, TColor::GetColor("#3F90DA"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*DYGto2LG-1Jets*PTG-10to100*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*DYGto2LG-1Jets*PTG-100to200*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*DYGto2LG-1Jets*PTG-200to400*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*DYGto2LG-1Jets*PTG-400to600*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022/mc/skim_ll/*DYGto2LG-1Jets*PTG-600*e.root"},"1");
  auto smzy_2022EE        = Process::MakeShared<Baby_pico>("Z+#gamma 2022EE", Process::Type::background, TColor::GetColor("#3F90DA"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/skim_ll/*DYGto2LG-1Jets*PTG-10to100*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/skim_ll/*DYGto2LG-1Jets*PTG-100to200*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/skim_ll/*DYGto2LG-1Jets*PTG-200to400*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/skim_ll/*DYGto2LG-1Jets*PTG-400to600*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2022EE/mc/skim_ll/*DYGto2LG-1Jets*PTG-600*e.root"},"1");
  auto smzy_2023          = Process::MakeShared<Baby_pico>("Z+#gamma 2023", Process::Type::background, TColor::GetColor("#3F90DA"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/skim_ll/*DYGto2LG-1Jets*PTG-10to100*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/skim_ll/*DYGto2LG-1Jets*PTG-100to200*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/skim_ll/*DYGto2LG-1Jets*PTG-200to400*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/skim_ll/*DYGto2LG-1Jets*PTG-400to600*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023/mc/skim_ll/*DYGto2LG-1Jets*PTG-600*e.root"},"1");
  auto smzy_2023BPix      = Process::MakeShared<Baby_pico>("Z+#gamma 2023BPix", Process::Type::background, TColor::GetColor("#3F90DA"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/skim_ll/*DYGto2LG-1Jets*PTG-10to100*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/skim_ll/*DYGto2LG-1Jets*PTG-100to200*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/skim_ll/*DYGto2LG-1Jets*PTG-200to400*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/skim_ll/*DYGto2LG-1Jets*PTG-400to600*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/2023BPix/mc/skim_ll/*DYGto2LG-1Jets*PTG-600*e.root"},"1");

  auto smzy_Run3          = Process::MakeShared<Baby_pico>("Z+#gamma Run 3", Process::Type::background, TColor::GetColor("#3F90DA"), {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*DYGto2LG-1Jets*PTG-10to100*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*DYGto2LG-1Jets*PTG-100to200*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*DYGto2LG-1Jets*PTG-200to400*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*DYGto2LG-1Jets*PTG-400to600*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*DYGto2LG-1Jets*PTG-600*e.root"},"1");

  auto signal_Run3        = Process::MakeShared<Baby_pico>("Combined Signal Run 3 x100", Process::Type::signal, kRed, {"/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*VBFHtoZG_Zto2L_M-125_TuneCP5_withDipoleRecoil_13p6TeV*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*W*H_HtoZG_WtoAll_Zto2L_M-125_TuneCP5_13p6TeV*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*ZH_ZtoAll_HtoZGto2LG_M-125_TuneCP5_13p6TeV*e.root", "/net/cms11/cms11r0/pico/NanoAODv12/htozgamma_redwood_v0/202*/mc/skim_ll/*ttHtoZG_Zto2L_M-125_TuneCP5_13p6TeV*e.root"}, "1");

  vector<string> procs_name = {"_run3datamc_","_2022data_","_2022mc_","_2022EEdata_","_2022EEmc_","_2023data_","_2023mc_","_2023BPixdata_","_2023BPixmc_", "_run3mc_"};
  vector<vector<shared_ptr<Process>>> procs_test = {{data_Run3, dyj_Run3, smzy_Run3},{data_2022_asmc},{dyj_2022, smzy_2022},{data_2022EE_asmc},{dyj_2022EE, smzy_2022EE},{data_2023_asmc},{dyj_2023, smzy_2023},{data_2023BPix_asmc},{dyj_2023BPix, smzy_2023BPix}};
  
  //vector<string> procs_name {"_just_1", "data_as_mc"};
  //vector<vector<shared_ptr<Process>>> procs_test = {{data_Run3,dyj_Run3,smzy_Run3},{data_2023BPix_asmc}};
  //This code creates the plotting styles for the 1D and 2D histograms
  PlotOpt style_2D("txt/plot_styles.txt","LogLumi2D");
  style_2D.Title(TitleType::preliminary)
          .Stack(PlotOptTypes::StackType::data_norm);
  PlotOpt style_lin_datanorm("txt/plot_styles.txt","LogLumi");
  style_lin_datanorm.Title(TitleType::preliminary)
          .YAxis(YAxisType::linear)
          .Stack(StackType::data_norm)
          .Overflow(OverflowType::none)
          .PrintVals(true)
          .Bottom(BottomType::ratio);

  PlotOpt style_lin_sigtop = style_lin_datanorm().Stack(StackType::signal_on_top);

  vector<PlotOpt> ops_2D = {style_2D};
  vector<PlotOpt> lin_datanorm = {PlotOpt("txt/plot_styles.txt","LinLumi").Stack(StackType::data_norm).ShowBackgroundError(1).Title(TitleType::info).Bottom(BottomType::ratio)};
  vector<PlotOpt> shapes      = {PlotOpt("txt/plot_styles.txt","LinLumi").Stack(StackType::shapes).ShowBackgroundError(1).Title(TitleType::simulation)};
  vector<PlotOpt> lin_sigtop = {style_lin_sigtop};
  vector<PlotOpt> lin_sigover = {PlotOpt("txt/plot_styles.txt","LinLumi").Stack(StackType::signal_overlay).ShowBackgroundError(1).Title(TitleType::simulation).Bottom(BottomType::ratio)};

  NamedFunc wgt_std("wgt_std",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }
    //return b.w_lumi()*w_years.GetScalar(b);
    //if( b.type() >= 200000 && b.type() <= 200500 ){return b.weight()*w_years.GetScalar(b)*100;}
    if(b.type() >= 200000 && b.type() <=200500 ){return b.w_lumi()*b.w_lep()*b.w_fs_lep()*b.w_bhig()*b.w_trig()*b.w_pu()*b.w_prefire()*b.w_photon()*b.w_phshape()*b.w_fakephoton()*b.w_nnlo()*w_years.GetScalar(b)*100;}
    return b.w_lumi()*b.w_lep()*b.w_fs_lep()*b.w_bhig()*b.w_trig()*b.w_pu()*b.w_prefire()*b.w_photon()*b.w_phshape()*b.w_fakephoton()*b.w_nnlo()*w_years.GetScalar(b);
    //return b.weight()*w_years.GetScalar(b);
  });

  vector<string> lumi_vec_str = {"61.89","7.98","7.89","26.67","26.67","17.79","17.79","9.45","9.45"};


  //This defines the plotting object used to create plots
  PlotMaker pm;

  NamedFunc sigblind = sideband();
  NamedFunc base_half = cutbitcheck(0b100000000000) && cutbitcheck(0b000100000000) && cutbitcheck(0b000010000000) && cutbitcheck(0b000001000000);//half baseline
  NamedFunc base_full = cutbitcheck(0b100000000000) && cutbitcheck(0b000100000000) && cutbitcheck(0b000010000000) &&cutbitcheck(0b000001000000) && cutbitcheck(0b000000100000) && cutbitcheck(0b000000010000) && cutbitcheck(0b000000001000) && cutbitcheck(0b000000000100) && cutbitcheck(0b000000000010);

  NamedFunc init = "pass==1 && use_event";
  NamedFunc jet2d_cut = init && cutbitcheck(0b100000000000) && cutbitcheck(0b000100000000) && cutbitcheck(0b000010000000);// && base_half;
  NamedFunc min_sel = init && base_half && cutbitcheck(0b000000100000) && sigblind;
  vector<NamedFunc> etaphimap_cuts = {gl_jet, gl_notmap, gl_jet && "!ismapvetoevt", gl_jet && "nel>0", gl_jet && "nel>2"};//, has_jet};
  vector<string> etaphimap_names = {"_no_veto","_jet_veto","_evt_veto", "test1", "test2"};//, "_all_jets"};

  NamedFunc obj_sel = init && base_half && cutbitcheck(0b000000100000) && sigblind;
  NamedFunc bdtvar_sel = init && base_full && sigblind;


  //When doing object kinematics comparisons, use these cuts to maximize stats. Compare data and MC on same plot
  vector<NamedFunc> vetoevts = {"ismapvetoevt","!ismapvetoevt"};
  vector<string> vnames = {"in_vetoevts", "in_nvetoevts"};
  

//This is the start of the plotting. the plots are seperated by "year", flavor, and control region
for(unsigned int i(0); i < procs_test.size();  i++){

  if(i>=1){//<8 2D eta phi plots. Run these on data and MC separately. Make sure to ignore the mc with "applicable evts" title. Show for uncleaned jets, with no veto, jet veto, and evt veto. For veto map
    for(unsigned int ep(0); ep < etaphimap_cuts.size(); ep++){
      pm.Push<Hist2D>(Axis(250,-4.7,4.7, "jet_eta[0]", "#eta_{j1}", {}),
                      Axis(250,-3.14,3.14, "jet_phi[0]", "#phi_{j1}", {}),
                      etaphimap_cuts[ep] && jet2d_cut, procs_test[i], ops_2D).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot35map_lead_jet_etaphi"+etaphimap_names[ep] + procs_name[i]);
    }
  }

  for(unsigned int v(0); v < vetoevts.size(); v++){

    if(i==0){//1D data MC comparisons. Run on the data+mc procs. For veto map
      //jet and met plots. Show the distribution in uncleaned and minimal cleaned jets in veto'd events. Njets is number of non-veto'd jets (fully cleaned) in the veto'd events.
/*
      pm.Push<Hist1D>(Axis(50,0,100 , "met"   ,"met [GeV]"   ,{}), vetoevts[v] && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot36map_met" + vnames[v]  + procs_name[i]);
      pm.Push<Hist1D>(Axis(50,-3.14,3.14 , "met_phi"   ,"met #phi"   ,{}), vetoevts[v] && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot36map_metphi" + vnames[v]  + procs_name[i]);

      pm.Push<Hist1D>(Axis(50,0,50 , "npv_good"   ,"N_{PV}^{good}"   ,{}), vetoevts[v] && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot__npvgood_test" + vnames[v]  + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,30,150 , "jet_pt[0]"   ,"lead jet p_{T} [GeV]"   ,{}), vetoevts[v] && has_jet && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot37map_lead_jet_pt" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,30,150 , "jet_pt[1]"   ,"sublead jet p_{T} [GeV]"   ,{}), vetoevts[v] && has_jet2 && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot37map_sublead_jet_pt" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,30,150 , getJetpT(has_clean_jet)   ,"lead jet p_{T} [GeV]"   ,{}), vetoevts[v] && has_clean_jet>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot38map_lead_cljet_pt" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,30,150 , getJetpT(has_clean_jet2)   ,"sublead jet p_{T} [GeV]"   ,{}), vetoevts[v] && has_clean_jet2>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot39map_sublead_cljet_pt" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,30,150 , getJetpT(has_clean_nvjet)   ,"lead jet p_{T} [GeV]"   ,{}), vetoevts[v] && has_clean_nvjet>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot40map_lead_cl_nv_jet_pt" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,30,150 , getJetpT(has_clean_nvjet2)   ,"sublead jet p_{T} [GeV]"   ,{}), vetoevts[v] && has_clean_nvjet2>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot41map_sublead_cl_nv_jet_pt" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,-4.5,4.5 , getJetEta(has_clean_nvjet)   ,"lead jet #eta"   ,{}), vetoevts[v] && has_clean_nvjet>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot40map_lead_cl_nv_jet_eta" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,-4.5,4.5 , getJetEta(has_clean_nvjet2)   ,"sublead jet #eta"   ,{}), vetoevts[v] && has_clean_nvjet2>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot41map_sublead_cl_nv_jet_eta" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,-3.14,3.14 , getJetPhi(has_clean_nvjet)   ,"lead jet #phi"   ,{}), vetoevts[v] && has_clean_nvjet>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot40map_lead_cl_nv_jet_phi" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(40,-3.14,3.14 , getJetPhi(has_clean_nvjet2)   ,"sublead jet #phi"   ,{}), vetoevts[v] && has_clean_nvjet2>=zero && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot41map_sublead_cl_nv_jet_phi" + vnames[v] + procs_name[i]);




      //pm.Push<Hist1D>(Axis(40,-4.5,4.5 , etaj1   ,"lead jet #eta"   ,{}), has_1min_jet && base_full, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot40map_lead_eta_forhorn" + vnames[v] + procs_name[i]);
      //pm.Push<Hist1D>(Axis(40,-4.5,4.5 , etaj2   ,"sublead jet #eta"   ,{}), has_2min_jet && base_full, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot41map_sublead_eta_forhorn" + vnames[v] + procs_name[i]);


      pm.Push<Hist1D>(Axis(6,-0.5,5.5 , "njet"   ,"njets"   ,{}), vetoevts[v] && min_sel, procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot42map_njet" + vnames[v] + procs_name[i]);

      //ll kinematics plots. These are with half the baseline selection applied, for the sake of statistics.
      pm.Push<Hist1D>(Axis(40,70,110 , "ll_m[0]"   ,"m_{ll} [GeV]"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot43map_mll" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(50,0,100 , "ll_pt[0]"   ,"pT_{ll} [GeV]"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot44map_ptll" + vnames[v] + procs_name[i]);

      pm.Push<Hist1D>(Axis(50,-4.7,4.7 , "ll_eta[0]"   ,"#eta_{ll}"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot45map_etall" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(50,-3.14,3.14 , "ll_phi[0]"   ,"#phi_{ll}"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot46map_phill" + vnames[v] + procs_name[i]);

      //photon kinematics plots. These are with half the baseline selection applied, for the sake of statistics.
      pm.Push<Hist1D>(Axis(50,10,100 , "photon_pt[0]"   ,"pT_{#gamma} [GeV]"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot47map_ptgamma" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(50,-2.5,2.5 , "photon_eta[0]"   ,"#eta_{#gamma}"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot48map_etagamma" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(50,-3.14,3.14 , "photon_phi[0]"   ,"#phi_{#gamma}"   ,{}), obj_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot49map_phigamma" + vnames[v] + procs_name[i]);

      //ggF BDT Input Variables. These are with the full baseline selection applied, for the sake of shapes. 
      pm.Push<Hist1D>(Axis(80,100,180 , "llphoton_m[0]"   ,"m_{ll#gamma} [GeV]"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot50map_llphoton_m" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(80,100,180 , "llphoton_refit_m"   ,"m_{ll#gamma} [GeV]"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot50map_llphoton_m" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,-1,1 , "llphoton_cosTheta[0]"   ,"cos #Theta"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot50map_cosTheta" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,-1,1 , "llphoton_costheta[0]"   ,"cos #theta"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot51map_costheta" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,-3.14,3.14 , "llphoton_psi[0]"   ,"#psi"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot52map_kinangle_psi" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,-3.14,3.14 , "llphoton_phi[0]"   ,"#phi"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot52map_kinangle_phi" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,3 , "photon_drmin[0]"   ,"dRmin"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot53map_dRmin" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,6 , "photon_drmax[0]"   ,"dRmax"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot54map_dRmax" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,2 , "llphoton_pt[0]/llphoton_m[0]"   ,"pT mass"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot55map_ptmass" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,0.2 ,   photon_relpterr   ,"#gamma res"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot56map_ph_res" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,-2.5,2.5 ,   lead_lepton_eta   ,"l1 #eta"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot57map_l1_eta" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,-2.5,2.5 ,   sublead_lepton_eta   ,"l2 #eta"   ,{}), bdtvar_sel && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot58map_l2_eta" + vnames[v] + procs_name[i]);

      //VBF BDT Input Variables not in ggF BDT. These are with the full baseline selection applied, for the sake of shapes. Jet variables distributions are shown with only fully cleaned jets to keep recognizable shapes, and for convenience. (n2p crafts some of these variables only from cleaned ("signal") jets. We still get information out, since it will tell us if even the good jets are messed up in these events, and in an ideal world we only veto the jets and not the events. 
      pm.Push<Hist1D>(Axis(20,0,1 ,   "photon_idmva[0]"   ,"#gamma MVA"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot59map_ph_cl_nv_mva" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,8 ,   "dijet_deta"   ,"d#eta jj"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot60map_detajj_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,3.14 ,   "dijet_dphi"   ,"d#phi jj"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot61map_dphijj_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,3.14 ,   "llphoton_dijet_dphi[0]"   ,"d#phi zgjj"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot62map_dphizgjj_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,1 ,   "llphoton_dijet_balance[0]"   ,"zgjj_bal"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot63map_zgjj_bal_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,100 ,   "llphoton_pTt[0]"   ,"pTt"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot64map_pTt_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,6 ,   "photon_zeppenfeld[0]"   ,"zepp"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot65map_ph_zep_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(20,0,4 ,   "photon_jet_mindr[0]"   ,"dRj#gamma"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot66map_phj_mindr_cl_nv" + vnames[v] + procs_name[i]);
      pm.Push<Hist1D>(Axis(50,0,2000 ,   "dijet_m"   ,"m_{jj}"   ,{}), bdtvar_sel && has_clean_nvjet2>=zero && "njet>=2" && vetoevts[v], procs_test[i], lin_datanorm).Weight(wgt_std).LuminosityTag(lumi_vec_str[i]).Tag("ShortName:plot67map_mjj_cl_nv" + vnames[v] + procs_name[i]);
*/
    }
  }

}


  //This chunk of code is needed to actually create the plots. The 1.0 can be changed to a luminosity value. 
  //Right now it is set to 2022 luminosity value but will need to be changed for the luminosity of 2022 and 2022EE
  pm.min_print_ = true;
  pm.MakePlots(1.01);
}


