  
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

namespace paths_4b{
    /////// Define map for processes /////
    string sig_path = "/net/cms18/cms18r0/pico/NanoAODv9/higgsino_4b_adelie_v3/2016/mc/unskimmed/";
    string bkg_met150 = "/net/cms18/cms18r0/pico/NanoAODv9/higgsino_4b_adelie_v4/20*/mc/skim_met150/";
    string data_met150 = "/net/cms11/cms11r0/pico/NanoAODv9/higgsino_4b_adelie_v6/20*/data/skim_met150/";
    string bkg_1lep = "/net/cms18/cms18r0/pico/NanoAODv9/higgsino_4b_adelie_v4/20*/mc/skim_higlep1T/";
    string data_1lep = "/net/cms11/cms11r0/pico/NanoAODv9/higgsino_4b_adelie_v6/20*/data/skim_higlep1T/";
    string bkg_2lep = "/net/cms18/cms18r0/pico/NanoAODv9/higgsino_4b_adelie_v4/20*/mc/skim_higlep2T/";
    string data_2lep = "/net/cms11/cms11r0/pico/NanoAODv9/higgsino_4b_adelie_v6/20*/data/skim_higlep2T/";

    //map between the search string after the basepath (basepath + "*" + searchstring + "*") and proc info

    std::map<std::string, ProcInfo> signal_info = {
      { "SMS-TChiHH_mChi-200_mLSP-0",          { sig_path, "TChiHH-200-1-HH",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-300_mLSP-0",          { sig_path, "TChiHH-300-1-HH",  "#0000FF", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-500_mLSP-0",          { sig_path, "TChiHH-500-1-HH",  "#00FF00", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-0",          { sig_path, "TChiHH-800-1-HH",  "#C849A9", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-0",         { sig_path, "TChiHH-1000-1-HH", "#FF5E02", "1", Process::Type::signal} },
    };
   
    std::map<std::vector<std::string>, paths_4b::ProcInfo> sample_info_met150 = { 
        { {"TTJets_","TTGJets","TTZToLLNuNu","TTZToQQ","TTWJetsToLNu","TTWJetsToQQ","TTTT","ttHTobb"},  { bkg_met150, "TT+X",       "#832DB6", "stitch_photon", Process::Type::background} },
        { {"_WJetsToLNu"},                       					      		{ bkg_met150, "W+Jets",     "#3F90DA", "stitch",        Process::Type::background} },
        { {"ZJetsToNuNu","DYJetsToLL"},                    				      		{ bkg_met150, "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
        { {"QCD_HT"},                          						      		{ bkg_met150, "QCD",        "#B9AC70", "stitch",        Process::Type::background} },
        { {"ST_"},                               					      		{ bkg_met150, "Single top", "#A96B59", "stitch",        Process::Type::background} },
        { {"GJets_DR-0p4_HT", "WminusH_HToBB", "WplusH_HToBB", "WWTo2L2Nu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu", "ZH_HToBB_ZToNuNu", "ZZ_TuneCP5"},                 							  	    { bkg_met150, "Others",     "#717581", "stitch",        Process::Type::background} },
        { {"*"}, 				   				              		{ data_met150,"Data", 	    "#000000", "stitch", 	      Process::Type::data} },

//        { "TTGJets",                   	       { bkg_met150,  "TT+X",	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "TTZToLLNuNu",                       { bkg_met150,  "TT+X",	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "TTZToQQ",                   	       { bkg_met150,  "TT+X",	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "TTWJetsToLNu",                      { bkg_met150,  "TT+X",	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "TTWJetsToQQ",                       { bkg_met150,  "TT+X",	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "TTTT",                              { bkg_met150,  "TT+X", 	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "ttHTobb",                           { bkg_met150,  "TT+X", 	    "#832DB6", "stitch_photon", Process::Type::background} },
//        { "DYJetsToLL",                        { bkg_met150,  "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
//        { "WminusH_HToBB",                     { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WplusH_HToBB",                      { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WWTo2L2Nu",                         { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WWTo1L1Nu2Q",                       { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WZTo1L1Nu2Q",                       { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WZTo1L3Nu",                         { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WZTo2Q2L",                          { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "WZTo3LNu",                          { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "ZH_HToBB_ZToNuNu",                  { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },
//        { "ZZ_TuneCP5",                        { bkg_met150,  "Others",     "#717581", "stitch",        Process::Type::background} },



    };
   

    std::map<std::string, paths_4b::ProcInfo> sample_info_1lep = { 
        { "TTJets_",                           { bkg_1lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTGJets",                   	       { bkg_1lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTZToLLNuNu",                       { bkg_1lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTZToQQ",                   	       { bkg_1lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTWJetsToLNu",                      { bkg_1lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTWJetsToQQ",                       { bkg_1lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTTT",                              { bkg_1lep,  "TT+X", 	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "ttHTobb",                           { bkg_1lep,  "TT+X", 	  "#832DB6", "stitch_photon", Process::Type::background} },

        { "_WJetsToLNu",                       { bkg_1lep,  "W+Jets",     "#3F90DA", "stitch",        Process::Type::background} },

        { "ZJetsToNuNu",                       { bkg_1lep,  "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
        { "DYJetsToLL",                        { bkg_1lep,  "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },

        { "QCD_HT",                            { bkg_1lep,  "QCD",        "#B9AC70", "stitch",        Process::Type::background} },

        { "ST_",                               { bkg_1lep,  "Single top", "#A96B59", "stitch",        Process::Type::background} },
  
        { "GJets_DR-0p4_HT",                   { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WminusH_HToBB",                     { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WplusH_HToBB",                      { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WWTo2L2Nu",                         { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WWTo1L1Nu2Q",                       { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo1L1Nu2Q",                       { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo1L3Nu",                         { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo2Q2L",                          { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo3LNu",                          { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "ZH_HToBB_ZToNuNu",                  { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "ZZ_TuneCP5",                        { bkg_1lep,  "Others",     "#717581", "stitch",        Process::Type::background} },

        { "*", 				       { data_1lep, "Data", 	  "#000000", "stitch",        Process::Type::data} },

    };

    std::map<std::string, paths_4b::ProcInfo> sample_info_2lep = { 
        { "TTJets_",                           { bkg_2lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTGJets",                   	       { bkg_2lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTZToLLNuNu",                       { bkg_2lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTZToQQ",                   	       { bkg_2lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTWJetsToLNu",                      { bkg_2lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTWJetsToQQ",                       { bkg_2lep,  "TT+X",	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "TTTT",                              { bkg_2lep,  "TT+X", 	  "#832DB6", "stitch_photon", Process::Type::background} },
        { "ttHTobb",                           { bkg_2lep,  "TT+X", 	  "#832DB6", "stitch_photon", Process::Type::background} },

        { "_WJetsToLNu",                       { bkg_2lep,  "W+Jets",     "#3F90DA", "stitch",        Process::Type::background} },

        { "ZJetsToNuNu",                       { bkg_2lep,  "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
        { "DYJetsToLL",                        { bkg_2lep,  "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },

        { "QCD_HT",                            { bkg_2lep,  "QCD",        "#B9AC70", "stitch",        Process::Type::background} },

        { "ST_",                               { bkg_2lep,  "Single top", "#A96B59", "stitch",        Process::Type::background} },
  
        { "GJets_DR-0p4_HT",                   { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WminusH_HToBB",                     { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WplusH_HToBB",                      { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WWTo2L2Nu",                         { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WWTo1L1Nu2Q",                       { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo1L1Nu2Q",                       { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo1L3Nu",                         { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo2Q2L",                          { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "WZTo3LNu",                          { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "ZH_HToBB_ZToNuNu",                  { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },
        { "ZZ_TuneCP5",                        { bkg_2lep,  "Others",     "#717581", "stitch",        Process::Type::background} },

        { "*", 				       { data_2lep, "Data", 	  "#000000", "stitch",        Process::Type::data} },

    };



}
