#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "zgamma/zg_functions.hpp"
#include "zgamma/zg_utilities.hpp"

namespace ZgFunctions {

  //isolated dielectron triggers for run 2
  const NamedFunc HLT_pass_dielectron("dielectron triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL()||b.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ();
  });

  //isolated dimuon triggers for run 2
  const NamedFunc HLT_pass_dimuon("dimuon triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ()||b.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ();
    return b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8()||b.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8();
  });

  //isolated dilepton triggers for run 2
  const NamedFunc HLT_pass_dilepton = HLT_pass_dielectron||HLT_pass_dimuon;

  //isolated single electron triggers for run 2
  const NamedFunc HLT_pass_singleelectron("single electron triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_Ele27_WPTight_Gsf();
    if (abs(b.SampleType())==2017)
      return b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG();
    return b.HLT_Ele32_WPTight_Gsf();
  });

  //isolated single muon triggers for run 2
  const NamedFunc HLT_pass_singlemuon("single muon triggers",[](const Baby &b) -> NamedFunc::ScalarType{
    if (abs(b.SampleType())==2016)
      return b.HLT_IsoMu24();
    if (abs(b.SampleType())==2017)
      return b.HLT_IsoMu27()||b.HLT_IsoMu24();
    return b.HLT_IsoMu24();
  });

  //isolated single lepton triggers for run 2
  const NamedFunc HLT_pass_singlelepton = HLT_pass_singleelectron||HLT_pass_singlemuon;

  //year integrated lumi weights
  const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.; //data
    if (b.SampleType()==2016)
      return 36.32264; 
    else if (b.SampleType()==2017)
      return 41.52756;
    return 59.67377;
  });

  //common NamedFuncs for run 2 baseline selection
  NamedFunc zg_baseline_nolep = "(ll_m[0]>50) && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)";
  NamedFunc zg_el_cuts = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)";
  NamedFunc zg_mu_cuts = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)";
  const NamedFunc zg_baseline_el =  NamedFunc(zg_el_cuts && zg_baseline_nolep).Name("electron_baseline");
  const NamedFunc zg_baseline_mu =  NamedFunc(zg_mu_cuts && zg_baseline_nolep).Name("muon_baseline");
  const NamedFunc zg_baseline = NamedFunc((zg_el_cuts || zg_mu_cuts) && zg_baseline_nolep).Name("baseline");

  //master stitch variable, currently only has DY/ZG stitch
  const NamedFunc stitch("stitch",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.type() >= 6000 && b.type() < 7000)
      return b.stitch_dy();
    return 1.0;
  });

  //drmax of lead photon
  const NamedFunc photon_drmax("photon_drmax",[](const Baby &b) -> NamedFunc::ScalarType{
      return ZgUtilities::pdrmax(b);
  });

  //lead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  const NamedFunc lead_lepton_eta("lead_lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
      if (b.ll_lepid()->size() < 1) return 0;
      if (b.ll_lepid()->at(0) == 11) return b.el_eta()->at(b.ll_i1()->at(0));
      return b.mu_eta()->at(b.ll_i1()->at(0));
  });

  //sublead lepton eta (=lep_eta[0], but this isn't saved in slims =( )
  //only works for 2 lepton events
  const NamedFunc sublead_lepton_eta("sublead_lepton_eta",[](const Baby &b) -> NamedFunc::ScalarType{
      if (b.ll_lepid()->size() < 1) return 0;
      if (b.ll_lepid()->at(0) == 11) return b.el_eta()->at(b.ll_i2()->at(0));
      return b.mu_eta()->at(b.ll_i2()->at(0));
  });

  //and of trigger and stitch
  const NamedFunc trigger_and_stitch = stitch&&(HLT_pass_dilepton||HLT_pass_singlelepton);
}







