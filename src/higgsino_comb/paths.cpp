  
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
      { "SMS-TChiHH_mChi-200_mLSP-0",          { sig_path, "TChiHH-200-1",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-300_mLSP-0",          { sig_path, "TChiHH-300-1",  "#0000FF", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-500_mLSP-0",          { sig_path, "TChiHH-500-1",  "#00FF00", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-0",          { sig_path, "TChiHH-800-1",  "#C849A9", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-0",         { sig_path, "TChiHH-1000-1", "#FF5E02", "1", Process::Type::signal} },
    };
   
    std::map<std::vector<std::string>, paths_4b::ProcInfo> sample_info_met150 = { 
        { {"TTJets_","TTGJets","TTZToLLNuNu","TTZToQQ","TTWJetsToLNu","TTWJetsToQQ","TTTT","ttHTobb"},  { bkg_met150, "TT+X",       "#832DB6", "stitch_photon", Process::Type::background} },
        { {"_WJetsToLNu"},                       					      		{ bkg_met150, "W+Jets",     "#3F90DA", "stitch",        Process::Type::background} },
        { {"ZJetsToNuNu","DYJetsToLL"},                    				      		{ bkg_met150, "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
        { {"QCD_HT"},                          						      		{ bkg_met150, "QCD",        "#B9AC70", "stitch",        Process::Type::background} },
        { {"ST_"},                               					      		{ bkg_met150, "Single top", "#A96B59", "stitch",        Process::Type::background} },
        { {"GJets_DR-0p4_HT", "WminusH_HToBB", "WplusH_HToBB", "WWTo2L2Nu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu", "ZH_HToBB_ZToNuNu", "ZZ_TuneCP5"},                 							  	    { bkg_met150, "Others",     "#717581", "stitch",        Process::Type::background} },
        { {"*"}, 				   				              		{ data_met150,"Data", 	    "#000000", "stitch", 	      Process::Type::data} },
    };
   
    std::map<std::vector<std::string>, paths_4b::ProcInfo> sample_info_1lep = { 
        { {"TTJets_","TTGJets","TTZToLLNuNu","TTZToQQ","TTWJetsToLNu","TTWJetsToQQ","TTTT","ttHTobb"},  { bkg_1lep, "TT+X",       "#832DB6", "stitch_photon", Process::Type::background} },
        { {"_WJetsToLNu"},                       					      		{ bkg_1lep, "W+Jets",     "#3F90DA", "stitch",        Process::Type::background} },
        { {"ZJetsToNuNu","DYJetsToLL"},                    				      		{ bkg_1lep, "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
        { {"QCD_HT"},                          						      		{ bkg_1lep, "QCD",        "#B9AC70", "stitch",        Process::Type::background} },
        { {"ST_"},                               					      		{ bkg_1lep, "Single top", "#A96B59", "stitch",        Process::Type::background} },
        { {"GJets_DR-0p4_HT", "WminusH_HToBB", "WplusH_HToBB", "WWTo2L2Nu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu", "ZH_HToBB_ZToNuNu", "ZZ_TuneCP5"},                 							  	    { bkg_1lep, "Others",     "#717581", "stitch",        Process::Type::background} },
        { {"*"}, 				   				              		{ data_1lep,"Data", 	  "#000000", "stitch", 	      Process::Type::data} },
    };

    std::map<std::vector<std::string>, paths_4b::ProcInfo> sample_info_2lep = { 
        { {"TTJets_","TTGJets","TTZToLLNuNu","TTZToQQ","TTWJetsToLNu","TTWJetsToQQ","TTTT","ttHTobb"},  { bkg_2lep, "TT+X",       "#832DB6", "stitch_photon", Process::Type::background} },
        { {"_WJetsToLNu"},                       					      		{ bkg_2lep, "W+Jets",     "#3F90DA", "stitch",        Process::Type::background} },
        { {"ZJetsToNuNu","DYJetsToLL"},                    				      		{ bkg_2lep, "Z+Jets",     "#E76300", "stitch",        Process::Type::background} },
        { {"QCD_HT"},                          						      		{ bkg_2lep, "QCD",        "#B9AC70", "stitch",        Process::Type::background} },
        { {"ST_"},                               					      		{ bkg_2lep, "Single top", "#A96B59", "stitch",        Process::Type::background} },
        { {"GJets_DR-0p4_HT", "WminusH_HToBB", "WplusH_HToBB", "WWTo2L2Nu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2Q2L", "WZTo3LNu", "ZH_HToBB_ZToNuNu", "ZZ_TuneCP5"},                 							  	    { bkg_2lep, "Others",     "#717581", "stitch",        Process::Type::background} },
        { {"*"}, 				   				              		{ data_2lep,"Data", 	  "#000000", "stitch", 	      Process::Type::data} },
    };


    std::map<std::string, paths_4b::ProcInfo> signal_info_procs1 = {
      { "SMS-TChiHH_mChi-150_mLSP-0",            { sig_path, "TChiHH-150-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-200_mLSP-0",            { sig_path, "TChiHH-200-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-300_mLSP-0",            { sig_path, "TChiHH-300-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-300_mLSP-100",          { sig_path, "TChiHH-300-100",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-400_mLSP-0",            { sig_path, "TChiHH-400-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-400_mLSP-100",          { sig_path, "TChiHH-400-100",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-400_mLSP-200",          { sig_path, "TChiHH-400-200",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-500_mLSP-0",            { sig_path, "TChiHH-500-1",    "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs2 = {
      { "SMS-TChiHH_mChi-500_mLSP-100",          { sig_path, "TChiHH-500-100",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-500_mLSP-200",          { sig_path, "TChiHH-500-200",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-500_mLSP-300",          { sig_path, "TChiHH-500-300",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-600_mLSP-0",            { sig_path, "TChiHH-600-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-600_mLSP-100",          { sig_path, "TChiHH-600-100",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-600_mLSP-200",          { sig_path, "TChiHH-600-200",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-600_mLSP-300",          { sig_path, "TChiHH-600-300",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-600_mLSP-400",          { sig_path, "TChiHH-600-400",  "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs3 = {
      { "SMS-TChiHH_mChi-700_mLSP-0",            { sig_path, "TChiHH-700-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-700_mLSP-100",          { sig_path, "TChiHH-700-100",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-700_mLSP-200",          { sig_path, "TChiHH-700-200",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-700_mLSP-300",          { sig_path, "TChiHH-700-300",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-700_mLSP-400",          { sig_path, "TChiHH-700-400",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-700_mLSP-500",          { sig_path, "TChiHH-700-500",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-0",            { sig_path, "TChiHH-800-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-100",          { sig_path, "TChiHH-800-100",  "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs4 = {
      { "SMS-TChiHH_mChi-800_mLSP-200",          { sig_path, "TChiHH-800-200",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-300",          { sig_path, "TChiHH-800-300",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-400",          { sig_path, "TChiHH-800-400",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-500",          { sig_path, "TChiHH-800-500",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-800_mLSP-600",          { sig_path, "TChiHH-800-600",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-0",            { sig_path, "TChiHH-900-1",    "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-100",          { sig_path, "TChiHH-900-100",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-200",          { sig_path, "TChiHH-900-200",  "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs5 = {
      { "SMS-TChiHH_mChi-900_mLSP-300",          { sig_path, "TChiHH-900-300",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-400",          { sig_path, "TChiHH-900-400",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-500",          { sig_path, "TChiHH-900-500",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-600",          { sig_path, "TChiHH-900-600",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-900_mLSP-700",          { sig_path, "TChiHH-900-700",  "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-0",           { sig_path, "TChiHH-1000-1",   "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-100",         { sig_path, "TChiHH-1000-100", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-200",         { sig_path, "TChiHH-1000-200", "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs6 = {
      { "SMS-TChiHH_mChi-1000_mLSP-300",         { sig_path, "TChiHH-1000-300", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-400",         { sig_path, "TChiHH-1000-400", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-500",         { sig_path, "TChiHH-1000-500", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-600",         { sig_path, "TChiHH-1000-600", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-700",         { sig_path, "TChiHH-1000-700", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1000_mLSP-800",         { sig_path, "TChiHH-1000-800", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-0",           { sig_path, "TChiHH-1100-0",   "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-100",         { sig_path, "TChiHH-1100-100", "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs7 = {
      { "SMS-TChiHH_mChi-1100_mLSP-200",         { sig_path, "TChiHH-1100-200", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-300",         { sig_path, "TChiHH-1100-300", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-400",         { sig_path, "TChiHH-1100-400", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-500",         { sig_path, "TChiHH-1100-500", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-600",         { sig_path, "TChiHH-1100-600", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-700",         { sig_path, "TChiHH-1100-700", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-800",         { sig_path, "TChiHH-1100-800", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1100_mLSP-900",         { sig_path, "TChiHH-1100-900", "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs8 = {
      { "SMS-TChiHH_mChi-1200_mLSP-0",           { sig_path, "TChiHH-1200-1",   "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-100",         { sig_path, "TChiHH-1200-100", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-200",         { sig_path, "TChiHH-1200-200", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-300",         { sig_path, "TChiHH-1200-300", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-400",         { sig_path, "TChiHH-1200-400", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-500",         { sig_path, "TChiHH-1200-500", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-600",         { sig_path, "TChiHH-1200-600", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-700",         { sig_path, "TChiHH-1200-700", "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs9 = {
      { "SMS-TChiHH_mChi-1200_mLSP-800",         { sig_path, "TChiHH-1200-800", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-900",         { sig_path, "TChiHH-1200-900", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1200_mLSP-1000",        { sig_path, "TChiHH-1200-1000","#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-0",           { sig_path, "TChiHH-1300-1",   "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-100",         { sig_path, "TChiHH-1300-100", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-200",         { sig_path, "TChiHH-1300-200", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-300",         { sig_path, "TChiHH-1300-300", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-400",         { sig_path, "TChiHH-1300-400", "#FF0000", "1", Process::Type::signal} },
    };
    std::map<std::string, paths_4b::ProcInfo> signal_info_procs10 = {
      { "SMS-TChiHH_mChi-1300_mLSP-500",         { sig_path, "TChiHH-1300-500", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-600",         { sig_path, "TChiHH-1300-600", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-700",         { sig_path, "TChiHH-1300-700", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-800",         { sig_path, "TChiHH-1300-800", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-900",         { sig_path, "TChiHH-1300-900", "#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-1000",        { sig_path, "TChiHH-1300-1000","#FF0000", "1", Process::Type::signal} },
      { "SMS-TChiHH_mChi-1300_mLSP-1100",        { sig_path, "TChiHH-1300-1100","#FF0000", "1", Process::Type::signal} },
    };


}
