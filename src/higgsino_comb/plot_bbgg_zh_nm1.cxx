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
  cmspaper.Title(TitleType::info) //simulation, preliminary
          .Overflow(OverflowType::none)
          .FileExtensions({"pdf"})
          .Stack(StackType::signal_overlay); //other options: shapes, data_norm
  vector<PlotOpt> ops = {cmspaper.YAxis(YAxisType::linear), cmspaper.YAxis(YAxisType::log)};

  //Plots
  PlotMaker pm;
  
  
  pm.Push<Hist1D>(Axis(30,0,200,"met", "p_{T}^{miss} [GeV]", {}), bbgg_zh_baseline, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_met");
  pm.Push<Hist1D>(Axis(8,0,8,"nvlep", "Number of veto leptons", {1}), bbgg_zh_nm_nvl, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_nvlep");
  pm.Push<Hist1D>(Axis(8,0,8,"njet", "Number of jets", {2,4}), bbgg_zh_nm_njet, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_njet");
  pm.Push<Hist1D>(Axis(30,60,160,"bb_m", "M_{bb} [GeV]", {60,100}), bbgg_zh_nm_mbb, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_mbb");
  pm.Push<Hist1D>(Axis(30,0,4,"photonphoton_dr", "\\DeltaR_{\\gamma\\gamma}", {2.5}), bbgg_zh_nm_ggdr, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_ggdr");
  pm.Push<Hist1D>(Axis(30,0,300,"photonphoton_pt", "p_{T}^{\\gamma\\gamma} [GeV]", {75}), bbgg_zh_nm_ggpt, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_ggpt");
  pm.Push<Hist1D>(Axis(20,100,200,"photonphoton_m", "M_{\\gamma\\gamma} [GeV]", {}), bbgg_zh_nm_ggm, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_ggm");
  pm.Push<Hist1D>(Axis(20,-1,1,"photon_idmva[0]", "Photon mva id of lead photon", {}), bbgg_zh_nm_photonidlead, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_photonidlead");
  pm.Push<Hist1D>(Axis(20,-1,1,"photon_idmva[1]", "Photon mva id of sublead photon", {}), bbgg_zh_nm_photonidsublead, procs_vec, ops).Weight(w_lumi_gmsb).Tag("ShortName:bbgg_zh_nm_photonidsublead");
  

  pm.multithreaded_ = true;
  pm.min_print_ = true;
  pm.MakePlots(1);
  
}






