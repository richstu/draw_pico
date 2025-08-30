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

int main() {
  gErrorIgnoreLevel = 6000;
  
  //trigger
  const NamedFunc diphoton_trig = triggers::diphoton_trig;
  
  //Load sample paths
  //using ProcInfo = bbgg_paths::ProcInfo;
  auto SampleInfo = bbgg_paths::SampleInfo;
  
  //Define processes
  vector<shared_ptr<Process>> procs_vec;
  
  // Iterate over the sample info map
  for (const auto& pair : SampleInfo) {
    
    const auto& searchstring = pair.first;
    const auto& ProcInfo = pair.second;
    
    //Pick processes
    if (ProcInfo.proclegend.find("-HH") != ProcInfo.proclegend.npos) { //only HH signal
      string fullpath(ProcInfo.basepath + "*" + searchstring + "*");
      auto proc = Process::MakeShared<Baby_pico>(ProcInfo.proclegend, ProcInfo.type, TColor::GetColor(ProcInfo.proccolor.c_str()), {fullpath}, diphoton_trig && ProcInfo.special_proccuts);
      procs_vec.push_back(proc);
    }

  }
  
  //Define plot options
  PlotOpt cmspaper("txt/plot_styles.txt","CMSPaper");
  cmspaper.Title(TitleType::info)
          .Overflow(OverflowType::both)
          .FileExtensions({"pdf"})
          .Stack(StackType::signal_overlay); //other options: shapes, data_norm
  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear), cmspaper.YAxis(YAxisType::log)};

  //Plot
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(20,0,500,"met", "p_{T}^{miss} [GeV]", {}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");
  pm.Push<Hist1D>(Axis(8,0,8,"nvlep", "Number of veto leptons", {1}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");
  pm.Push<Hist1D>(Axis(8,0,8,"njet", "Number of jets", {2,4}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");
  pm.Push<Hist1D>(Axis(30,60,160,"bb_m", "M_{bb} [GeV]", {100,140}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");
  pm.Push<Hist1D>(Axis(30,0,4,"photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");
  pm.Push<Hist1D>(Axis(30,0,300,"photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");
  pm.Push<Hist1D>(Axis(20,100,150,"photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {}), bbgg_hh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_hh_baseline");

  pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(1);
  
}






