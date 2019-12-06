#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

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
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  bool do_allbkg = true;
  string sample_name = "search";
  float lumi = 137;
  bool unblind = true;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  PlotOpt lin_shapes_info("txt/plot_styles.txt", "CMSPaper");
  lin_shapes_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes)
    .LegendColumns(3);
  // PlotOpt log_shapes_info = lin_shapes_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_types = {lin_shapes_info};

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/mc/merged_higmc_higtight/";
  // if (sample_name=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep1/";
  // if (sample_name=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep2/";
  // if (sample_name=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higqcd/";
  string folderdata = bfolder+"";
  // if (sample_name=="ttbar") folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higlep1/";
  // if (sample_name=="zll") folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higlep2/";
  // if (sample_name=="qcd") folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higqcd/";

  set<string> alltags; 
  if (sample_name=="ttbar" || sample_name=="search") alltags = {"*TTJets_*Lept*",
                                      "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", 
                                      "*ttHTobb*.root","*_TTTT*.root"};
  // if (sample_name=="zll") alltags = {"*DYJetsToLL*.root"};
  // if (sample_name=="qcd") alltags = {//"*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
  //                               //"*QCD_HT300to500_Tune*", 
  //                                "*QCD_HT500to700_Tune*",
  //                                "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
  //                                "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"};
  if (do_allbkg) { // don't include QCD MC unless in QCD control sample_name
    if(sample_name=="qcd") alltags = {"*TTJets_*Lept*", 
                                 "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
                                 "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
                                 "*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
                                 "*QCD_HT300to500_Tune*", 
                                 "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
    else  alltags = {"*TTJets_*Lept*",
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
            "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
  }
  set<string> allfiles = attach_folder(foldermc,alltags);

  // Baseline definitions
  string baseline = "!lowDphiFix && nvlep==0 && ntk==0";
  string higtrim = "hig_cand_drmax[0]<=2.2 && hig_cand_dm[0] <= 40 && hig_cand_am[0]<=200";
  baseline += "&& njet>=4 && njet<=5 && nbt>=2 && "+higtrim;
  

  // zll skim: ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && 
  // nlep==2 && nlep>=1 && Max$(leps_pt)>30 && njets>=4&&njets<=5
  // if (sample_name=="zll") baseline = baseline+"&& nlep==2 && met<50";
  // // qcd skim - met>150 && nvlep==0 && (njets==4||njets==5)
  // if (sample_name=="qcd") baseline = baseline+"&& nvlep==0 && lowDphiFix";
  // // ttbar skim - nlep==1 && (njets==4||njets==5) && nbm>=2
  // if (sample_name=="ttbar") baseline = baseline+"&& nlep==1 && mt<100";
  // search skim - met>100 && nvlep==0 && (njets==4||njets==5) && nbm>=2
  // if (sample_name=="search") baseline = baseline+"&& nvlep==0";
    

  ////// Nb cuts
  vector<string> nbcuts;
  unsigned firstnb = 0;
  nbcuts.push_back("nbm==0");
  nbcuts.push_back("nbm==1");
  nbcuts.push_back("nbm==2");
  if (sample_name=="ttbar" || sample_name=="search") {
    firstnb = 2;
    nbcuts.push_back("nbm==3 && nbl==3");
    nbcuts.push_back("nbm>=3 && nbl>=4");
  } 

  string sname = "t#bar{t}+X";
  if (sample_name=="qcd") sname = "QCD";
  if (sample_name=="zll") sname = "Z#rightarrow ll";
  if (do_allbkg) sname = "Bkg.";

  vector<int> colors = {kGreen+3, kGreen+1, kOrange, kAzure+1, kBlue+1};

  // Cuts applied to all processes but not shown in plot title
  NamedFunc filters = HigUtilities::pass_2016 && "stitch"; 

  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    // if (sample_name=="qcd" && inb==nbcuts.size()-1) continue;
    procs.push_back(Process::MakeShared<Baby_pico>(sname+" "+RoundNumber(inb,0).Data()+"b", 
      Process::Type::background, colors[inb], allfiles, baseline+"&&"+nbcuts[inb]));
  }

  vector<int> colors_trub = {kAzure-4, kTeal-8, kOrange-4, kPink+2, kMagenta-1};
  vector<shared_ptr<Process> > procs_trub = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    procs_trub.push_back(Process::MakeShared<Baby_pico>(sname+" "+RoundNumber(inb,0).Data()+" B-hadrons", 
        Process::Type::background, colors_trub[inb], allfiles, Functions::ntrub==inb && baseline));
  }

  string lumi_s=RoundNumber(lumi,1).Data();

  vector<vector<string>> combos, combos_labels;
  // int color_data = kBlue-7;
  // if (sample_name=="zll") { // do 0b vs 1b
  //   color_data = kOrange+1;
  //   combos.push_back({"nbm==0","nbm==1"});
  //   combos.push_back({"nbm==1","nbm==2"});
  // } else if (sample_name=="qcd") { // do 0b vs 1b and 2b vs 3+b
  //   color_data = kOrange;
  //   combos.push_back({"nbm==0","nbm==1"});
  //   combos.push_back({"nbm==1","nbm==2"});
  //   combos.push_back({"nbm==2","nbm==3"});
  //   combos.push_back({"nbm==2", "nbm>=3 && nbl>=4"});
  // } else 
  if (sample_name=="search" || sample_name=="ttbar") {
    combos.push_back({"nbm==2", "nbm==3 && nbl==3"});
    combos.push_back({"nbm==2", "nbm>=3 && nbl>=4"});
  }
  for (unsigned ind(0); ind<combos.size(); ind++){
    vector<string> label;
    for(unsigned lab=0; lab<2; lab++){
      if (Contains(combos[ind][lab],"nbm==0")) label.push_back("0b");
      if (Contains(combos[ind][lab],"nbm==1")) label.push_back("1b");
      if (Contains(combos[ind][lab],"nbm==2")) label.push_back("2b");
      if (Contains(combos[ind][lab],"nbm==3")) label.push_back("3b");
      if (Contains(combos[ind][lab],"nbm>=3") && !Contains(combos[ind][lab],"nbl>=4")) label.push_back("#geq 3b");
      if (Contains(combos[ind][lab],"nbl>=4")) label.push_back("4b");
    }
    combos_labels.push_back(label);
  }

  // vector<vector<shared_ptr<Process> >> procs_data;
  // for (unsigned ind(0); ind<combos.size(); ind++){
  //   vector<string> icomb = combos[ind], ilab = combos_labels[ind];
  //   procs_data.push_back(vector<shared_ptr<Process> >());
  //   string tmpcuts = filters +"&&"+icomb[0]; //if not cast here, it crashes
  //   procs_data.back().push_back(Process::MakeShared<Baby_pico>(ilab[0]+" Data "+lumi_s+" fb^{-1}", 
  //                    Process::Type::background, kBlack, {folderdata+"*root"}, 
  //                    Higfuncs::trig_hig && tmpcuts));
  //   procs_data.back().back()->SetFillColor(color_data);
  //   procs_data.back().back()->SetLineColor(color_data);
  //   procs_data.back().back()->SetLineWidth(2);

  //   tmpcuts = filters +"&&"+icomb[1];
  //   procs_data.back().push_back(Process::MakeShared<Baby_pico>(ilab[1]+" Data "+lumi_s+" fb^{-1}", 
  //                    Process::Type::data, kBlack, 
  //     {folderdata+"*root"}, Higfuncs::trig_hig && tmpcuts));
  // }


  string metcut = "met>150";
  // if (sample_name=="zll") metcut = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>0";
  // else if (sample_name=="ttbar") metcut = "1";

  vector<string> xcuts;
  xcuts.push_back(metcut);
  xcuts.push_back("hig_cand_drmax[0]<=1.1");
  xcuts.push_back("hig_cand_drmax[0]>1.1");

  PlotMaker pm;
  NamedFunc wgt = "w_lumi*w_isr";//Higfuncs::weight_higd*Higfuncs::eff_higtrig;

  for (unsigned ic(0); ic<xcuts.size(); ic++){
    // if (sample_name!="search" || unblind) {
    //   for (unsigned i(0); i<combos.size(); i++)
    //     pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
    //       baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs_data[i], plt_types)
    // .Tag(sample_name+"_datavdata"+to_string(i)).RatioTitle(combos_labels[i][1],combos_labels[i][0]);
    // }

    pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      baseline+"&&"+xcuts[ic], procs, plt_types).Weight(wgt).Tag(sample_name+"_shape_bcats")
      .RatioTitle("Bkg. nb","Bkg. 2b");

    pm.Push<Hist1D>(Axis(10,0,200,"hig_cand_am[0]", "#LTm#GT [GeV]", {100., 140.}),
      baseline+"&&"+xcuts[ic], procs_trub, plt_types).Weight(wgt).Tag(sample_name+"_shape_trub");
  }

  pm.min_print_ = true;
  pm.MakePlots(35.9);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    
      {"unblind", no_argument, 0, 'u'},   
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "js:tau", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'u':
      unblind = true;
      break;
    case 'a':
      do_allbkg = true;
      break;
    case 's':
      sample_name = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;

    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}