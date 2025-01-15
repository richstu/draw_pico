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

  //Below two functions are those that check if the triggers/trigger pT cuts are passed
  const NamedFunc pass_trigs = "trig_double_el||trig_double_mu||trig_single_el||trig_single_mu";
  const NamedFunc pass_singlelepton_trigs_and_pt("pass_singlelepton_trigs_and_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_pt_sel = false;
    if(b.trig_single_mu()){
      if(!(b.SampleTypeString().Contains("2017")) && b.nmu()>0 && b.mu_pt() -> at(0) > 25){ pass_pt_sel = true; }
      if(b.SampleTypeString().Contains("2017") && b.nmu()>0 && b.mu_pt() -> at(0) > 30){ pass_pt_sel = true; }
    }
    if(b.trig_single_el()){
      if(!(b.SampleTypeString().Contains("2016")) && b.nel()>0 && b.el_pt() -> at(0) > 35){ pass_pt_sel = true; }
      if(b.SampleTypeString().Contains("2016") && b.nel()>0 && b.el_pt() -> at(0) > 30){ pass_pt_sel = true; }
    }
    return pass_pt_sel;
  });

  const NamedFunc pass_singleelectron_trigs_and_pt("pass_singleelectron_trigs_and_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_pt_sel = false;
    if(b.trig_single_el()){
      if(!(b.SampleTypeString().Contains("2016")) && b.nel()>0 && b.el_pt() -> at(0) > 35){ pass_pt_sel = true; }
      if(b.SampleTypeString().Contains("2016") && b.nel()>0 && b.el_pt() -> at(0) > 30){ pass_pt_sel = true; }
    }
    return pass_pt_sel;
  });

  const NamedFunc pass_singlemuon_trigs_and_pt("pass_singlemuon_trigs_and_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    bool pass_pt_sel = false;
    if(b.trig_single_mu()){
      if(!(b.SampleTypeString().Contains("2017")) && b.nmu()>0 && b.mu_pt() -> at(0) > 25){ pass_pt_sel = true; }
      if(b.SampleTypeString().Contains("2017") && b.nmu()>0 && b.mu_pt() -> at(0) > 30){ pass_pt_sel = true; }
    }
    return pass_pt_sel;
  });

  const NamedFunc pass_dilepton_trigs_and_pt   = "(trig_double_mu && nmu > 1 && mu_pt[0] > 20 && mu_pt[1] > 10) || (trig_double_el && nel > 1 && el_pt[0] > 25 && el_pt[1] > 15)";
  const NamedFunc pass_dielectron_trigs_and_pt = "(trig_double_el && nel > 1 && el_pt[0] > 25 && el_pt[1] > 15)";
  const NamedFunc pass_dimuon_trigs_and_pt     = "(trig_double_mu && nmu > 1 && mu_pt[0] > 20 && mu_pt[1] > 10)";
  const NamedFunc pass_trigs_and_pt    = pass_dilepton_trigs_and_pt || pass_singlelepton_trigs_and_pt;
  const NamedFunc pass_el_trigs_and_pt = pass_dielectron_trigs_and_pt || pass_singleelectron_trigs_and_pt;
  const NamedFunc pass_mu_trigs_and_pt = pass_dimuon_trigs_and_pt || pass_singlemuon_trigs_and_pt;


  //year integrated lumi weights
  const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleTypeString().Contains("-")) 
      return 1.; //data
    if (b.SampleTypeString()=="2016APV")
      return 19.51; //from brilcalc 
    if (b.SampleTypeString()=="2016")
      return 16.80; //from brilcalc
    else if (b.SampleTypeString()=="2017")
      return 41.48;
    else if (b.SampleTypeString()=="2018")
      return 59.83;
    else if (b.SampleTypeString()=="2022")
      return 8.17;
    else if (b.SampleTypeString()=="2022EE")
      return 27.01;
    else if (b.SampleTypeString()=="2023")
      return 17.61;
    //else if (b.SampleTypeString()=="2023BPix")
    return 9.53;
  });

  //year integrated lumi weights
  //const NamedFunc w_years_noapv("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  //  if (b.SampleType()<0) return 1.; //data
  //  if (b.SampleType()==2016)
  //    return 36.32264; 
  //  else if (b.SampleType()==2017)
  //    return 41.52756;
  //  return 59.67377;
  //});

  //Run 3 weight-up
  const NamedFunc w_run3("w_run3", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.SampleType()<0) return 1.; //data
    return 2.472;
  });

  //common NamedFuncs for run 2 baseline selection
  NamedFunc zg_baseline_nolep = "nlep>=2 && nphoton>=1 && (ll_m[0]>50) && ((photon_pt[0]/llphoton_m[0])>=15.0/110.0) && ((llphoton_m[0]+ll_m[0])>=185) && (photon_drmin[0]>0.4)";
  NamedFunc zg_el_cuts = "(ll_lepid[0]==11) && (el_pt[ll_i1[0]]>25) && (el_pt[ll_i2[0]]>15)";
  NamedFunc zg_mu_cuts = "(ll_lepid[0]==13) && (mu_pt[ll_i1[0]]>20) && (mu_pt[ll_i2[0]]>10)";
  const NamedFunc zg_baseline_el = NamedFunc(zg_el_cuts && zg_baseline_nolep).Name("electron_baseline");
  const NamedFunc zg_baseline_mu = NamedFunc(zg_mu_cuts && zg_baseline_nolep).Name("muon_baseline");
  const NamedFunc zg_baseline = NamedFunc((zg_el_cuts || zg_mu_cuts) && zg_baseline_nolep).Name("baseline");

  //master stitch variable
  const NamedFunc stitch("stitch",[](const Baby &b) -> NamedFunc::ScalarType{
    //remove ZGToLLG from DYJets
    if(b.type() >= 6000 && b.type() < 7000)
      return b.stitch_dy();
    //remove DYJets from ZGToLLG
    if (b.type() >= 17000 && b.type() < 18000)
      return !b.stitch_dy();
    //remove ttG from TTJets
    if (b.type() >= 1000 && b.type() < 2000)
      //return b.stitch_photon(); currently bugged since photon doesn't exempt itself
      if (b.photon_pflavor()->size()>0)
        if (b.photon_pflavor()->at(0)==1)
          return 0;
    //remove WWG from WW
    if (b.type() >= 14000 && b.type() < 15000)
      return b.stitch_dy();
      //if (b.photon_pflavor()->size()>0)
      //  if (b.photon_pflavor()->at(0)==1)
      //    return 0;
    //remove WZG from WZ - poor man's stitch
    if (b.type() >= 15000 && b.type() < 16000)
      return b.stitch_dy();
      //if (b.photon_pflavor()->size()>0)
      //  if (b.photon_pflavor()->at(0)==1)
      //    return 0;
    //don't use ZZG because only leptonic decays
    //remove ZZG from ZZ - poor man's stitch
    //if (b.type() >= 16000 && b.type() < 17000)
    //  if (b.photon_pflavor()->size()>0)
    //    if (b.photon_pflavor()->at(0)==1)
    //      return 0;
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

  //pT/m of Higgs candidate
  const NamedFunc llphoton_rel_pt = NamedFunc("llphoton_pt[0]/llphoton_m[0]").Name("llphoton_rel_pt");

  //Below is the baseline used by HIG-19-014
  const NamedFunc hig19014_baseline     = "nllphoton>0" && pass_trigs_and_pt && "ll_m[0] > 50 && photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] + llphoton_m[0] > 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && pass";

  const NamedFunc tightened_baseline_pinnacles="(zg_cutBitMap==3582 || zg_cutBitMap==3583 || zg_cutBitMap==3070 || zg_cutBitMap==3071)";

  //Below is the tightened baseline for lassen and earlier picos
  const NamedFunc tightened_baseline= "nllphoton>0" && pass_trigs_and_pt && "photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] > 80 && ll_m[0] < 100 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && ll_m[0] + llphoton_m[0] > 185 && pass";
  const NamedFunc tightened_baseline_refit_all_sels= "nllphoton>0" && pass_trigs_and_pt && "photon_pt[0]/llphoton_refit_m > 15.0/110 && ll_refit_m > 80 && ll_refit_m < 100 && llphoton_refit_m > 100 && llphoton_refit_m < 180 && ll_refit_m + llphoton_refit_m > 185 && pass";
  const NamedFunc tightened_baseline_refit= "nllphoton>0" && pass_trigs_and_pt && "photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] > 80 && ll_m[0] < 100 && llphoton_refit_m > 100 && llphoton_refit_m < 180 && ll_refit_m + llphoton_refit_m > 185 && pass";

  const NamedFunc wgt("wgt",[](const Baby &b) -> NamedFunc::ScalarType{ 
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }

    double w_year = w_years.GetScalar(b);
    if( b.type() >= 200000 && b.type() <= 200500 ){ w_year=w_year*10;}
    if(b.SampleType() > 2020){ return b.w_lumi()*w_year;}

    return b.weight()*w_year;
  });


}






