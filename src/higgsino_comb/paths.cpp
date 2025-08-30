  
#include <algorithm>
#include <stdlib.h>
#include <regex>
#include <string>
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/plot_opt.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino_comb/vardef.hpp"
#include "higgsino_comb/regions.hpp"
#include "higgsino_comb/paths.hpp"

using namespace std;

namespace bbgg_paths{

    /////// Define basepaths /////
    //2D locally produced fullsim signals basepath
    string bfoldersigHH("/net/cms18/cms18r0/pico/NanoAODv9/higgsino_bbgg_meru/2016/mc/unskimmed/"); 
    string bfoldersigZH("/net/cms37/data2/mlai/SMS-TChiHZ_HToGG_2D/NanoAODv9/higgsino_bbgg_meru/2016/mc/unskimmed/");

    //Background paths, March 2024 Run 2 UL production basepath
    string bfolderbkg("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/20*/mc/unskimmed/");
    string bfolderbkg_smh("/net/cms11/cms11r0/pico/NanoAODv9/higgsino_smh_rakaposhi/20*/mc/unskimmed/");
    string bfolderbkg_gjet2("/net/cms11/cms11r0/pico/NanoAODv2/bbgg_higgsino_signal_prod_11/201*/mc/unskimmed/");

    //Data basepaths
    string bfolderdata("/net/cms11/cms11r0/pico/NanoAODv9/bbgg_higgsino_signal_prod_11/");


    /////// Define map for processes /////

    //map between the search string after the basepath (basepath + "*" + searchstring + "*") and proc info
    std::map<std::string, ProcInfo> SampleInfo = { 

      { "SMS-TChiHH_mChi-150_mLSP-0_HToGG",         { bfoldersigHH, "TChiHH-150-1-HH", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-300_mLSP-0_HToGG",         { bfoldersigHH, "TChiHH-300-1-HH", "#0000FF", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-500_mLSP-0_HToGG",         { bfoldersigHH, "TChiHH-500-1-HH", "#00FF00", "1", Process::Type::signal} },

      { "SMS-TChiHZ_mChi-150_mLSP-0_HToGG",         { bfoldersigZH, "TChiHH-150-1-ZH", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHZ_mChi-300_mLSP-0_HToGG",         { bfoldersigZH, "TChiHH-300-1-ZH", "#0000FF", "1", Process::Type::signal} },
      { "SMS-TChiHZ_mChi-500_mLSP-0_HToGG",         { bfoldersigZH, "TChiHH-500-1-ZH", "#00FF00", "1", Process::Type::signal} },   

      { "DiPhotonJetsBox_",                         { bfolderbkg,  "DiPhotonJetsBox", "#BD1F01", "no_genb", Process::Type::background} },
      { "DiPhotonJetsBox1BJet",                     { bfolderbkg,  "DiPhotonJetsBox1BJet", "#FFA90E", "one_genb", Process::Type::background} },
      { "DiPhotonJetsBox2BJets",                    { bfolderbkg,  "DiPhotonJetsBox2BJets", "#FFA90E", "two_genb", Process::Type::background} },
      { "GJet_Pt-20to40",                           { bfolderbkg,  "GJet", "#92DADD", "use_event", Process::Type::background} },
      { "GJet_Pt-40toInf",                          { bfolderbkg_gjet2, "GJet", "#92DADD", "use_event", Process::Type::background} },
      { "TTGG_0Jets",                               { bfolderbkg,  "TTGG_0Jets", "#A96B59", "1", Process::Type::background} },
      { "TTGJets",                                  { bfolderbkg,  "TTGjets", "#A96B59", "1", Process::Type::background} },
      { "TTTo2L2Nu",                                { bfolderbkg,  "TTTo2L2Nu", "#A96B59", "1", Process::Type::background} },
      { "QCD_Pt-30to40",                            { bfolderbkg,  "QCD", "#0194da", "1", Process::Type::background} },
      { "QCD_Pt-40",                                { bfolderbkg,  "QCD", "#0194da", "1", Process::Type::background} },

      { "bbHToGG",                                  { bfolderbkg_smh,  "bbHToGG_SMH", "#0194da", "1", Process::Type::background} },
      { "VHToGG",                                   { bfolderbkg_smh,  "VHToGG_SMH", "#af0032", "1", Process::Type::background} },
      { "ttHJetToGG",                               { bfolderbkg_smh,  "ttHJetToGG_SMH", "#401641", "1", Process::Type::background} },
      { "GluGluHToGG",                              { bfolderbkg_smh,  "GluGluHToGG_SMH", "#99dcff", "1", Process::Type::background} },
      { "VBFHToGG",                                 { bfolderbkg_smh,  "VBFHToGG_SMH", "#ec00ff", "1", Process::Type::background} },
    
      { "2016*/data/unskimmed/*DoubleEG*",          { bfolderdata,  "data", "#000000", "1", Process::Type::data} },
      { "2017/data/unskimmed/*DoubleEG*",           { bfolderdata,  "data", "#000000", "1", Process::Type::data} },
      { "2018/data/unskimmed/*EGamma*",             { bfolderdata,  "data", "#000000", "1", Process::Type::data} },

    };


}

