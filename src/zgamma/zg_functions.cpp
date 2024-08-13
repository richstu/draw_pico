#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/table_row.hpp"
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


  //Turns a vector<NamedFunc> into one usable for a cutflow table
  std::vector<NamedFunc> progressive_cuts(std::vector<NamedFunc> vector_NamedFunc){
    for(unsigned int idx = 1; idx<vector_NamedFunc.size(); idx++){ vector_NamedFunc[idx] = vector_NamedFunc[idx-1] && vector_NamedFunc[idx];}
    return vector_NamedFunc;
  }

  //This function adds all selections and reverses the selection at reverse
  NamedFunc Nreverse1(std::vector<NamedFunc> vector_NamedFunc, unsigned int reverse){
    NamedFunc return_NamedFunc = "1";
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx==reverse){return_NamedFunc = return_NamedFunc && !(vector_NamedFunc[idx]); continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }

  //Returns a NamedFunc replacing one selection (marked by skip)
  NamedFunc Nreplace1(std::vector<NamedFunc> vector_NamedFunc, NamedFunc replace, unsigned int skip){
    NamedFunc return_NamedFunc = "1";//vector_NamedFunc[start];
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx==skip){return_NamedFunc = return_NamedFunc && replace; continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }

  //Returns a NamedFunc with all but one selection (marked by skip)
  NamedFunc Nminus1(std::vector<NamedFunc> vector_NamedFunc, unsigned int skip){
//    NamedFunc return_NamedFunc = "1";
//    unsigned int start = skip!=0 ?  0 : 1;
    NamedFunc return_NamedFunc = "1";//vector_NamedFunc[start];
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx==skip){continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }


  //Returns a NamedFunc without all selections in the vector skip
  NamedFunc Nminusk(std::vector<NamedFunc> vector_NamedFunc, std::vector<unsigned int> skip){
    unsigned int idx_s = 0;
    NamedFunc return_NamedFunc = "1";
    for(unsigned int idx = 0; idx<vector_NamedFunc.size(); idx++){
      if(idx_s < skip.size() && idx==skip[idx_s]){idx_s++; continue;}  
      return_NamedFunc = return_NamedFunc && vector_NamedFunc[idx];
    }
    return return_NamedFunc;
  }

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


  //This version of llgamma mass is used when nphoton can be 0. This defaults to llphoton_m[0] otherwise
  const NamedFunc mlly("mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return ZgUtilities::AssignH(b).M(); });
  const NamedFunc pTy_mlly("pTy_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.photon_pt()->at(0))/(ZgUtilities::AssignH(b).M()); });
  const NamedFunc mll_mlly("mll_mlly",[](const Baby &b) -> NamedFunc::ScalarType{ return (b.ll_m() ->at(0)) + ZgUtilities::AssignH(b).M(); });


  //This is a vector that gives the individual selections in the baseline selection
  const std::vector<NamedFunc> vector_run2_baseline      = {"nll>0", pass_dilepton_trigs_and_pt, "ll_m[0] > 50", "nphoton > 0",  pTy_mlly > 15.0/110, mll_mlly > 185, mlly > 100 && mlly < 180};
  const std::vector<NamedFunc> vector_tightened_baseline = {"nll>0", pass_trigs_and_pt, "ll_m[0] > 50", "nphoton > 0 && photon_id80[0]",  pTy_mlly > 15.0/110, "ll_m[0] > 80 && ll_m[0] < 100",
                                                            mlly > 100 && mlly < 180, mll_mlly > 185};

  //Below is the baseline used by HIG-19-014
  const NamedFunc run2_baseline     = "nllphoton>0" && pass_dilepton_trigs_and_pt && "ll_m[0] > 50 && photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] + llphoton_m[0] > 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180";

  //Below is the tightened baseline with the photon wp80 and the 80 GeV < mll < 100 GeV selection
  const NamedFunc tightened_baseline= "nllphoton>0" && pass_trigs_and_pt && "photon_id80[0] && photon_pt[0]/llphoton_m[0] > 15.0/110 && ll_m[0] > 80 && ll_m[0] < 100 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && ll_m[0] + llphoton_m[0] > 185";
  const NamedFunc tightened_baseline_refit= "nllphoton>0" && pass_trigs_and_pt && "photon_id80[0] && photon_pt[0]/llphoton_refit_m > 15.0/110 && ll_refit_m > 80 && ll_refit_m < 100 && llphoton_refit_m > 100 && llphoton_refit_m < 180 && ll_refit_m + llphoton_refit_m > 185";





  const NamedFunc wgt("wgt",[](const Baby &b) -> NamedFunc::ScalarType{ 
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }
    double w_year = w_years.GetScalar(b);
    if( b.type() >= 200000 && b.type() <= 200500 ){ w_year=w_year*10;}
    //double weight_fix = b.w_lumi()*b.w_lep()*b.w_fs_lep()*b.w_bhig_df()*b.w_isr()*b.w_pu()*b.w_prefire()*b.w_photon();
    return b.weight()*w_year;
  });

  const NamedFunc wgt_wlumi("wgt_wlumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }
    double w_year = w_years.GetScalar(b);
    if( b.type() >= 200000 && b.type() <= 200500 ){ w_year=w_year*10;}
    //double weight_fix = b.w_lumi()*b.w_lep()*b.w_fs_lep()*b.w_bhig_df()*b.w_isr()*b.w_pu()*b.w_prefire()*b.w_photon();
    return b.w_lumi()*w_year;
  });

  const NamedFunc wgt_nodiff("wgt_nodiff",[](const Baby &b) -> NamedFunc::ScalarType{ 
    if(b.SampleTypeString().Contains("-")) {
      return 1;
    }
    double w_year = w_years.GetScalar(b);
    double weight_fix = b.w_lumi()*b.w_lep()*b.w_fs_lep()*b.w_bhig_df()*b.w_isr()*b.w_pu()*b.w_prefire()*b.w_photon();
    return weight_fix*w_year;
  });


  const NamedFunc testing_sample    = "event%314159%3==0";
  const NamedFunc training_sample   = "event%314159%3==1";
  const NamedFunc validation_sample = "event%314159%3==2";

  NamedFunc add_cut(NamedFunc & current_cut, NamedFunc additional_cut) {
    current_cut = current_cut && additional_cut;
    return current_cut;
  }

  void constructCutflowTable(std::vector<TableRow> & tablerows, NamedFunc weight, int electron_or_muon, bool isMC) {
    NamedFunc current_cut = "1";
    if (isMC == true){
      current_cut = "use_event";
    }
 
    if (electron_or_muon == 0) { // el
      tablerows = std::vector<TableRow>{
        TableRow("No selection",                                  add_cut(current_cut,"1"),                                           0, 0, weight),
        TableRow("$N_e \\geq 2$",                                 add_cut(current_cut, "nel>=2 && ll_lepid[0]==11"),                  0, 0, weight),
        TableRow("e, ee trigger",                                 add_cut(current_cut, "trig_single_el||trig_double_el"),             0, 0, weight),
        TableRow("el trigger $p_{T}$ cuts",                       add_cut(current_cut, pass_el_trigs_and_pt),                         0, 0, weight),
        TableRow("$N_{\\gamma} \\geq 1$",                         add_cut(current_cut, "nphoton>=1"),                                 0, 0, weight),
        TableRow("$Charge_{ll} = 0$",                             add_cut(current_cut, "ll_charge[0] == 0"),                          0, 0, weight),
        TableRow("$80 \\text{ GeV} < m_{ee} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),                0, 0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ee\\gamma} > 15./110$",     add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),     0, 0, weight),
        TableRow("$m_{ee} + m_{ee\\gamma} > 185 \\text{ GeV}$",   add_cut(current_cut, "(llphoton_m[0] + ll_m[0]) > 185"),            0, 0, weight),
        TableRow("$100 < m_{ee\\gamma} < 180 \\text{ GeV}$",      add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"), 0, 0, weight)
      };
    } else if (electron_or_muon == 1) { // mu
      tablerows = std::vector<TableRow>{
        TableRow("No selection",                                          add_cut(current_cut, "1"),                                          0, 0, weight),
        TableRow("$N_{\\mu} \\geq 2 $",                                   add_cut(current_cut, "nmu>=2 && ll_lepid[0]==13"),                  0, 0, weight),
        TableRow("$\\mu, \\mu\\mu$ trigger",                              add_cut(current_cut, "trig_single_mu||trig_double_mu"),             0, 0, weight),
        TableRow("mu trigger $p_{T}$ cuts",                               add_cut(current_cut, pass_mu_trigs_and_pt),                         0, 0, weight),
        TableRow("$N_{\\gamma} \\geq 1$",                                 add_cut(current_cut, "nphoton>=1"),                                 0, 0, weight),
        TableRow("$Charge_{ll} = 0$",                                     add_cut(current_cut, "ll_charge[0] == 0"),                          0, 0, weight),
        TableRow("$80 \\text{ GeV} < m_{\\mu\\mu} < 100 \\text{ GeV}$",   add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),                0, 0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{\\mu\\mu\\gamma} > 15./110$",       add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),     0, 0, weight),
        TableRow("$m_{\\mu\\mu\\gamma}+m_{\\mu\\mu} > 185 \\text{ GeV}$", add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) > 185"),              0, 0, weight),
        TableRow("$100 < m_{\\mu\\mu\\gamma} < 180 \\text{ GeV}$",        add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"), 0, 0, weight)
      };
    } else { // mu + el
      tablerows = std::vector<TableRow>{
        TableRow("No selection",                                  add_cut(current_cut, "1"),                                                            0, 0, weight),
        TableRow("$N_e \\geq 2 || N_{\\mu} \\geq 2 $",            add_cut(current_cut, "(nel>=2 && ll_lepid[0]==11) || (nmu>=2 && ll_lepid[0]==13)"),   0, 0, weight),
        TableRow("e, ee trigger $||$ $\\mu, \\mu\\mu$ trigger",   add_cut(current_cut, pass_trigs),                                                     0, 0, weight),
        TableRow("mu and el trigger $p_{T}$ cuts",                add_cut(current_cut, pass_trigs_and_pt),                                              0, 0, weight),
        TableRow("$N_{\\gamma} \\geq 1$",                         add_cut(current_cut, "photon_id80[0]"),                                               0, 0, weight),
        TableRow("$Charge_{ll} = 0$",                             add_cut(current_cut, "ll_charge[0] == 0"),                                            0, 0, weight),
        TableRow("$80 \\text{ GeV} < m_{ll} < 100 \\text{ GeV}$", add_cut(current_cut, "ll_m[0] > 80 && ll_m[0]<100"),                                  0, 0, weight),
        TableRow("$p_{T}^{\\gamma}/m_{ll\\gamma} > 15./110$",     add_cut(current_cut, "(photon_pt[0]/llphoton_m[0])>(15./110)"),                       0, 0, weight),
        TableRow("$m_{ll\\gamma}+m_{ll} > 185 \\text{ GeV}$",     add_cut(current_cut, "(llphoton_m[0]+ll_m[0]) > 185"),                                0, 0, weight),
        TableRow("$100 < m_{ll\\gamma} < 180 \\text{ GeV}$",      add_cut(current_cut, "llphoton_m[0] > 100 && llphoton_m[0] < 180"),                   0, 0, weight)
      };
    }
  }



}







