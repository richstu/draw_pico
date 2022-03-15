///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  //vector<string> sigm = {}; 
  vector<string> sigm = {"175","500","950"}; 
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  //string bfolder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/";
  //if (string(getenv("LOCAL_PICO_DIR"))=="") bfolder = "/net/cms17/cms17r0/pico/NanoAODv7/";
  string bfolder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath_v3/";

  //string mc_production = "higgsino_angeles"; // higgsino_eldorado
  //string mc_production = "higgsino_eldorado";
  //string mc_production = "higgsino_humboldt";
  string mc_production = "higgsino_klamath_v3";
  set<int> years;
  HigUtilities::parseYears("run2",years);

  string mc_skim_folder("/mc/merged_higmc_higloose/");
  string sig_skim_folder("/SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/");

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  //mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mctags["other_and_single_t"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root","*_ST_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });
  
  string c_ps = "stitch";
  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, 1,
              attach_folder(bfolder, years, mc_skim_folder, mctags["other"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background, 1,
              attach_folder(bfolder, years, mc_skim_folder, mctags["single_t"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, 1,
              attach_folder(bfolder, years, mc_skim_folder, mctags["qcd"]),c_ps)); 
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, 1,
              attach_folder(bfolder, years, mc_skim_folder,mctags["wjets"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, 1,
              attach_folder(bfolder, years, mc_skim_folder,mctags["zjets"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", 
              Process::Type::background,1,
              attach_folder(bfolder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", 
              Process::Type::background,1,
              attach_folder(bfolder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  for (unsigned isig(0); isig<sigm.size(); isig++) {
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", 
        Process::Type::signal, 1, 
        attach_folder(bfolder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), 
        "1"));
  }

  string c_2bt  = "nbt>=2";
  string c_2b   = "nbt==2&&nbm==2";
  string c_3b   = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b   = "nbt>=2&&nbm>=3&&nbl>=4";
  string c_ge3b = "nbt>=2&&nbm>=3";

  string c_hig_dm = "hig_cand_dm[0]<=40";
  string c_hig_am = "hig_cand_am[0]<=200";
  string c_drmax  = "hig_cand_drmax[0]<=2.2";
  string hig = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && (hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
  string sbd = "hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200 && hig_cand_dm[0] <= 40 && !(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";

  //NamedFunc wgt = "w_lumi*w_isr";//Higfuncs::weight_higd * Higfuncs::eff_higtrig;
  //NamedFunc wgt = "weight"* Higfuncs::eff_higtrig;
  NamedFunc wgt = Higfuncs::final_weight;

  NamedFunc search_filters = Higfuncs::final_pass_filters; //pass_filters&& "met/mht<2 && met/met_calo<2&&weight<1.5"

  NamedFunc baseline = "nvlep==0 && njet>=4 && njet<=5"; //pass_muon_jet && met/met_calo<5
  NamedFunc same_cut = baseline && "met>150";

  string ncols = to_string(procs.size()+2);  
  PlotMaker pm;
  pm.Push<Table>("cutflow", vector<TableRow>{
  TableRow("1", 
    "1",0,0, wgt),
  TableRow("Filters", 
    Higfuncs::final_pass_filters,0,0, wgt),
  TableRow("MET$>$150", //included in skim
    Higfuncs::final_pass_filters && "met>150",0,0, wgt),
  TableRow("$0\\ell$", //included in skim
    Higfuncs::final_pass_filters && "nvlep==0&&met>150",0,0, wgt),
  TableRow("$N_{bT}\\geq$2", //included in skim
    Higfuncs::final_pass_filters && "nvlep==0&&met>150&&nbt>=2",0,0, wgt),
  TableRow("$\\text{4-5 jets}$", 
    Higfuncs::final_pass_filters && baseline&&"met>150&&nbt>=2",0,0, wgt),
  TableRow("DeltaPhi", 
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&nbt>=2",0,0, wgt),
  TableRow("met over mht", 
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&met/mht<2&&nbt>=2",0,0, wgt),
  TableRow("met over met calo", 
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&met/mht<2&&met/met_calo<2&&nbt>=2",0,0, wgt),
  //TableRow("nbt$>=$2", 
  //  baseline&&"met>150&&met/mht<2&&met/met_calo<2&&"+c_2bt,0,0, wgt),
  TableRow("$|\\Delta m| < 40$ GeV",     
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&met/mht<2&&met/met_calo<2&&"+c_2bt+"&&"+c_hig_dm,0,0, wgt),
  TableRow("$\\langle m_{bb} \\rangle < 200$ GeV",     
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&met/mht<2&&met/met_calo<2&&"+c_2bt+"&&"+c_hig_dm+" && "+c_hig_am,0,0, wgt),
  TableRow("$\\Delta R_{\\text{max}} < 2.2$",                    
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&met/mht<2&&met/met_calo<2&&"+c_2bt+"&&"+c_hig_dm+" && "+c_hig_am+" && "+c_drmax,0,0,wgt),
  TableRow("Track veto",                    
    Higfuncs::final_pass_filters && baseline&&"met>150&&!low_dphi_met&&met/mht<2&&met/met_calo<2&&"+c_2bt+"&&"+c_hig_dm+" && "+c_hig_am+" && "+c_drmax+" &&ntk==0",0,1,wgt),

  //TableRow("Track veto", 
  //  same_cut &&"ntk==0",0,0, wgt),
  //TableRow("pass filter", 
  //  same_cut&&Higfuncs::pass_filters,0,0, wgt),
  //TableRow("DeltaPhi",        
  //  same_cut &&"!lowDphiFix",0,1, wgt),

  // $\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$

  TableRow("\\multicolumn{"+ncols+"}{c}{HIG: $100<\\left< m \\right>\\leq140$}\\\\%", 
    "met>1e6",0,1, wgt),

  TableRow("$100<\\left< m \\right>\\leq140$ GeV", 
    baseline && "ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix &&"+hig,0,0, wgt),
  TableRow("3b + 4b", 
    baseline && "ntk==0 && met>150 &&"+c_ge3b+"&& !lowDphiFix &&"+hig,0,0,wgt),
  TableRow("4b", 
    baseline && "ntk==0 && met>150 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,0, wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>200$ GeV", 
    baseline && "ntk==0 && met>200 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,0,wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>300$ GeV", 
    baseline && "ntk==0 && met>300 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,0,wgt),
  TableRow("$p_{\\rm T}^{\\rm miss}>400$ GeV", 
    baseline && "ntk==0 && met>400 &&"+c_4b+"  && !lowDphiFix &&"+hig,0,1,wgt),

  // TableRow("\\multicolumn{"+ncols+"}{c}{SBD: $\\left< m\\right> <100$ or $140<\\left< m\\right>\\leq200$}\\\\%", 
  //   "met>1e6",0,1, wgt),

  // TableRow("CHECK:: SBD 2b, MET 150-200", 
  //   baseline + " && ntk==0 && met>150 && met<=200 &&"+c_2b+"   && !lowDphiFix &&"+sbd,0,0, wgt),
  // TableRow("$\\left< m\\right> <100$ or $140<\\left< m\\right>\\leq200$ GeV", 
  //   baseline + " && ntk==0 && met>150 &&"+c_2bt+"   && !lowDphiFix &&"+sbd,0,0, wgt),
  // TableRow("3b + 4b", 
  //   baseline +"  && ntk==0 && met>150 &&"+c_ge3b+"&& !lowDphiFix &&"+sbd,0,0,wgt),
  // TableRow("4b", 
  //   baseline + " && ntk==0 && met>150 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0, wgt),
  // TableRow("$p_{\\rm T}^{\\rm miss}>200$ GeV", 
  //   baseline + " && ntk==0 && met>200 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0,wgt),
  // TableRow("$p_{\\rm T}^{\\rm miss}>300$ GeV", 
  //   baseline + " && ntk==0 && met>300 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0,wgt),
  // TableRow("$p_{\\rm T}^{\\rm miss}>450$ GeV", 
  //   baseline + " && ntk==0 && met>450 &&"+c_4b+"  && !lowDphiFix &&"+sbd,0,0,wgt),

  },procs,0);


  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
