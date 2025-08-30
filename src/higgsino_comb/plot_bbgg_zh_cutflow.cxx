  
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
#include "core/sample_loader.hpp"
#include "higgsino_comb/vardef.hpp"
#include "higgsino_comb/regions.hpp"
#include "higgsino_comb/weights.hpp"
#include "higgsino_comb/paths.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace regions_bbgg;
using namespace weights;
using namespace vardef;

int main() {
  gErrorIgnoreLevel = 6000;

  //trigger
  const NamedFunc diphoton_trig = triggers::diphoton_trig;
  
  //Load sample paths
  //using ProcInfo = bbgg_paths::ProcInfo;
  auto SampleInfo = bbgg_paths::SampleInfo;
  
  //Group processes
  std::set<std::string> gjetfiles;
  for (const auto& pair : SampleInfo) {
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;   
    std::string fullpath(ProcInfo.basepath + "*" + searchstring + "*");
    if (ProcInfo.proclegend.find("GJet") != ProcInfo.proclegend.npos) {
      gjetfiles.insert(fullpath);
    }
  }

  //Define processes
  // Iterate over the sample info map
  vector<shared_ptr<Process>> procs_vec;
  bool found_gjet = false;

  for (const auto& pair : SampleInfo) {
    
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
    
    string fullpath(ProcInfo.basepath + "*" + searchstring + "*");

    //Pick processes to ignore
    if ((ProcInfo.proclegend.find("-HH") != ProcInfo.proclegend.npos) || 
        (ProcInfo.proclegend.find("bbHToGG_SMH") != ProcInfo.proclegend.npos) ||
        (ProcInfo.proclegend.find("VBFHToGG_SMH") != ProcInfo.proclegend.npos) ||
        (ProcInfo.proclegend.find("QCD") != ProcInfo.proclegend.npos) ||
        //(ProcInfo.proclegend.find("_SMH") != ProcInfo.proclegend.npos) ||
        (ProcInfo.proclegend.find("data") != ProcInfo.proclegend.npos)) continue;
    
    //Grouped processes    
    if (ProcInfo.proclegend.find("GJet") != ProcInfo.proclegend.npos) {
      if (!found_gjet) {
        auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), gjetfiles, diphoton_trig && ProcInfo.special_proccuts);
        found_gjet = true;
        procs_vec.push_back(proc);
      }
      else continue;
    } 

    //Other processes
    else { 
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, diphoton_trig && ProcInfo.special_proccuts);
      procs_vec.push_back(proc);
    }

  }
  
  //Define plot options
  PlotOpt cmspaper("txt/plot_styles.txt","CMSPaper");
  cmspaper.Title(TitleType::simulation)
          .Overflow(OverflowType::none)
          .FileExtensions({"pdf"})
          .Stack(StackType::signal_overlay); //other options: shapes, data_norm
  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear), cmspaper.YAxis(YAxisType::log)};

  //Cutflow table using the cuts vector defined in the regions.cpp file
  PlotMaker pm;
  vector<TableRow> tablevec;
  stringstream cuts;
  int index = 0;

  //table for Zh signal region
  tablevec.push_back(TableRow("No cuts", "1", 0, 0, w_lumi_gmsb)); //first row without any cuts
  tablevec.push_back(TableRow("sig_decay_bbgg", sig_decay_bbgg, 0, 0, w_lumi_gmsb)); 
  // Iterate over the sample info map
  for (const auto& pair : cuts_bbgg_zh) {

    const auto& cutlabel = pair.first;
    const auto& cut = pair.second; 

    if (index == 0) { cuts << cut; }
    else { cuts << "&&" << cut; }
    index += 1;  

    tablevec.push_back(TableRow(cutlabel,cuts.str() && sig_decay_bbgg, 0, 0, w_lumi_gmsb));
  }

  //push special cuts not in the cut vector
  tablevec.push_back(TableRow("trigger", cuts.str() && sig_decay_bbgg && diphoton_trig, 0, 0, w_lumi_gmsb));  
  tablevec.push_back(TableRow("deepflavmedium_lead", cuts.str() && sig_decay_bbgg && diphoton_trig && DeepFlav_leadjet, 0, 0, w_lumi_gmsb));
  tablevec.push_back(TableRow("deepflavmedium_sublead", cuts.str() && sig_decay_bbgg && diphoton_trig && DeepFlav_leadjet && DeepFlav_subleadjet, 0, 0, w_lumi_gmsb));

  pm.Push<Table>("Cutflow_bbgg_zh", tablevec, procs_vec, false, true, false, false, false, false); //do_zbi, print_table, print_pie, print_titlepie, do_eff, do_unc
  pm.Push<Table>("Cutflow_bbgg_zh_eff", tablevec, procs_vec, false, true, false, false, true, false); //do_zbi, print_table, print_pie, print_titlepie, do_eff, do_unc

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(1);
  
}






