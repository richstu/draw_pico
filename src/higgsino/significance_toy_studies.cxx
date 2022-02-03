#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <random>

#include "core/gamma_params.hpp"
#include "core/utilities.hpp"

#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/DistFuncMathCore.h"
#include "RooStats/NumberCountingUtils.h"
//#include "Math/GSLRndmEngines.h"
#include "TF1.h"
#include "TPaveStats.h"

#include "RooRealVar.h"
#include "RooLognormal.h"
#include "RooDataSet.h"

#include "TStyle.h"
#include "TArrow.h"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::istringstream;
using std::setprecision;
using std::stringstream;
using std::vector;
using std::map;
using std::max;
using std::thread;
using std::ref;

void calculate_predictions(TRandom3 & random, vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > & predictions, vector<float> & pvalues, vector<TH1D*> * kappa_hists = 0, vector<TH1D*> * prediction_hists = 0, vector<TH1D*> * toy_hists = 0);
void calculate_kappas(vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > & kappas);
vector<float> get_kappa_systematics(bool fraction=true);

float to_pvalue(float significance) {return 0.5-TMath::Erf(double(significance/sqrt(2)))/2;}
float to_significance(float pvalue) {
  if (pvalue == 0) return 99;
  else if (pvalue == 1) return 0;
  else return ROOT::Math::normal_quantile_c(pvalue, 1);
}

Double_t gamma(Double_t *x, Double_t * par) {
  return ROOT::Math::gamma_pdf(x[0],par[0]+1,1);
}

// -1 sigma, 0 sigma, 1 sigma
vector<float> get_sigma_band(TH1D * hist) {
  Double_t sigma[3] = {0.5-0.3413, 0.5, 0.5+0.3413}; Double_t x_sigma[3];
  hist->GetQuantiles(3, x_sigma, sigma);
  return {static_cast<float>(x_sigma[0]), static_cast<float>(x_sigma[1]), static_cast<float>(x_sigma[2])};
}

// 0 sigma, -1 sigma diff, +1 sigma diff
vector<float> get_sigma_band(vector<float> & toy_predictions, float prediction) {
  sort(toy_predictions.begin(), toy_predictions.end());
  int prediction_index = -1;
  for (int iToy = 0; iToy < static_cast<int>(toy_predictions.size()); iToy++) {
    if (toy_predictions[iToy] > prediction) {
      prediction_index = iToy;
      break;
    }
  }
  float nSigma = 1;
  double integratedGaus = intGaus(0,1,0,nSigma);
  int prediction_low_index = prediction_index - static_cast<int>(integratedGaus * toy_predictions.size());
  if (prediction_low_index<0) prediction_low_index = 0;
  int prediction_high_index = prediction_index + static_cast<int>(integratedGaus * toy_predictions.size());
  if (prediction_high_index>static_cast<int>(toy_predictions.size())) prediction_high_index = toy_predictions.size()+1;
  float prediction_low = prediction - toy_predictions[prediction_low_index];
  float prediction_high = toy_predictions[prediction_high_index] - prediction;
  return {prediction, prediction_low, prediction_high};
}

// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(TRandom3 & rand, vector<vector<float> > const & entries, vector<vector<float> > const & weights,
    vector<float> const & powers, float &mSigma, float &pSigma, vector<float> & fKappas, 
    bool do_data, bool verbose, double syst, __attribute__((unused)) bool do_plot, int nrep, float nSigma){
  fKappas.clear();
  fKappas.reserve(nrep);
  int nbadk(0);
  double mean(0.), bignum(1e10);
  // Doing kappa variations
  for(int rep(0), irep(0); rep < nrep; rep++) {
    fKappas.push_back(1.);
    bool Denom_is0(false);
    for(unsigned obs(0); obs < powers.size(); obs++) {
      float observed(0.);
      for(unsigned sam(0); sam < entries[obs].size(); sam++) {
        // With a flat prior, the expected average of the Poisson with N observed is (Gamma(N+1,1) == Poisson(N))
        // Rounding the expected yield for data stats
        if(do_data) observed += entries[obs][sam]*weights[obs][sam];
        else observed += gsl_ran_gamma(entries[obs][sam]+1,1,rand)*weights[obs][sam];
      } // Loop over samples
      //if(do_data) observed = gsl_ran_gamma(static_cast<int>(0.5+observed)+1,1,rand);
      if(do_data) observed = gsl_ran_gamma(observed+1,1,rand);
      if(observed <= 0 && powers[obs] < 0) Denom_is0 = true;
      else fKappas[irep] *= pow(observed, powers[obs]);
    } // Loop over number of observables going into kappa

    if(syst>=0){
      double factor = exp(rand.Gaus(0,log(1+syst)));
      fKappas[irep] *= factor;
    }
    if(Denom_is0 && fKappas[irep]==0) {
      fKappas.pop_back();
      nbadk++;
    }else {
      if(Denom_is0) fKappas[irep] = bignum;
      else mean += fKappas[irep];
      irep++;
    }
  } // Loop over fluctuations of kappa (repetitions)
  int ntot(nrep-nbadk);
  mean /= static_cast<double>(ntot);

  sort(fKappas.begin(), fKappas.end());
  // integrated gaussian
  double gSigma = intGaus(0,1,0,nSigma);
  int iMedian((nrep-nbadk+1)/2-1);
  int imSigma(iMedian-static_cast<int>(gSigma*ntot)), ipSigma(iMedian+static_cast<int>(gSigma*ntot));
  float median(fKappas[iMedian]);
  mSigma = median-fKappas[imSigma]; pSigma = fKappas[ipSigma]-median;

  // Finding standard value
  float stdval(1.);
  bool infStd(false);
  for(unsigned obs(0); obs < powers.size(); obs++) {
    float stdyield(0.); // yield in bin
    if(verbose) cout<<obs<<": ";
    for(unsigned sam(0); sam < entries[obs].size(); sam++) {
      if(verbose) cout<<"Yield"<<sam<<" "<<entries[obs][sam]*weights[obs][sam]
        <<", N"<<sam<<" "<<entries[obs][sam]
          <<", avW"<<sam<<" "<<weights[obs][sam]<<". ";
      stdyield += entries[obs][sam]*weights[obs][sam];
    }
    if(verbose) cout<<"  ==> Total yield "<<stdyield<<endl;
    if(stdyield <= 0 && powers[obs] < 0) infStd = true;
    else stdval *= pow(stdyield, powers[obs]);
  } // Loop over number of observables going into kappa
  if(infStd) stdval = median;
  else {
    int istd(0);
    for(int rep(0); rep < ntot; rep++) 
      if(fKappas[rep]>stdval) {istd = rep; break;}
    imSigma = istd-static_cast<int>(gSigma*ntot);
    ipSigma = istd+static_cast<int>(gSigma*ntot);
    if(imSigma<0){ // Adjusting the length of the interval in case imSigma has less than 1sigma
      ipSigma += (-imSigma);
      imSigma = 0;
    }
    if(ipSigma>=ntot){ // Adjusting the length of the interval in case ipSigma has less than 1sigma
      imSigma -= (ipSigma-ntot+1);
      ipSigma = ntot-1;
    }
    mSigma = fabs(stdval-fKappas[imSigma]); pSigma = fKappas[ipSigma]-stdval;
  }

  if(do_plot) {
    gStyle->SetOptStat(0);              // No Stats box
    TCanvas can;
    can.SetMargin(0.15, 0.05, 0.12, 0.11);
    int nbins(100);
    double minH(stdval-3*fabs(mSigma)), maxH(stdval+3*pSigma);
    if(minH < fKappas[0]) minH = fKappas[0];
    if(maxH > fKappas[ntot-1]) maxH = fKappas[ntot-1];
    TH1D histo("h","",nbins, minH, maxH);
    TH1D herr("herr","",nbins, minH, maxH);
    for(int rep(0); rep < ntot; rep++) {
      histo.Fill(fKappas[rep]);   
      if(fKappas[rep] >= stdval - mSigma && fKappas[rep] <= stdval + pSigma)
        herr.Fill(fKappas[rep]); 
    }
    double mode(histo.GetBinLowEdge(histo.GetMaximumBin()));
    if(verbose) cout<<"Std kappa = "<<stdval<<"+"<<pSigma<<"-"<<mSigma<<".   Mean = "<<mean
      <<". Mode = "<<mode<<". Median = "<<median<<endl;
    //histo.SetBinContent(1, histo.GetBinContent(1)+nbadk);
    //histo.SetBinContent(nbins, histo.GetBinContent(nbins)+histo.GetBinContent(nbins+1));

    herr.SetLineColor(0);
    herr.SetFillColor(kGray);
    if (nSigma==2) herr.SetFillColor(kAzure+1);
    else if (nSigma==2) herr.SetFillColor(kMagenta+1);
    histo.SetTitleOffset(1.1, "X");
    histo.SetTitleOffset(1.5, "Y");
    // histo.Scale(1/histo.Integral());
    // herr.Scale(1/histo.Integral());
    histo.SetMinimum(0);
    histo.SetMaximum(histo.GetMaximum()*1.2);
    histo.SetTitleSize(0.05, "XY");
    histo.SetLineWidth(3);
    histo.Draw();
    histo.SetXTitle("Toy value");
    histo.SetYTitle("Number of toys");
    histo.Draw();
    herr.Draw("same");
    histo.Draw("same");
    histo.Draw("same axis");
    TString title;
    for(unsigned obs(0); obs < powers.size(); obs++) {
      float observed(0.);
      for(unsigned sam(0); sam < entries[obs].size(); sam++) {
        observed += entries[obs][sam]*weights[obs][sam];
      }
      if(obs>0){
        if(powers[obs]>0) title += " #times ";
        else title += " / ";
      }
      title += RoundNumber(observed,0);
    }
    TString pName = "gamma_"+title+"_"+RoundNumber(nSigma,0)+"sigma.pdf"; 
    pName.ReplaceAll("/", "_d_"); pName.ReplaceAll("#times","_");
    pName.ReplaceAll(" ","");
    title += (" #rightarrow "+RoundNumber(stdval,2)+"^{+"+RoundNumber(pSigma,2)+"}_{-"+RoundNumber(mSigma,2)+"}");
    title = "Interval for "+title;
    histo.SetTitle(title);

    int abin = (nbins * fabs(stdval-minH)/(maxH-minH))+1;
    TArrow arrow;
    arrow.SetLineColor(kRed+2); arrow.SetFillColor(kRed+2);
    arrow.SetArrowSize(0.015); arrow.SetLineWidth(4);
    arrow.DrawArrow(stdval, 0, stdval, histo.GetBinContent(abin));

    can.SaveAs("plots/"+pName);
    cout<<" open "<<"plots/"<<pName<<endl;
  } // do_plot

  return stdval;
}

// Assumes data is ordered from (plane=0,ibin=0), (plane=0,ibin=1), // ibin is split in plane
void load_counts(string const & count_in_bins_filename, vector<vector<vector<GammaParams> > > & counts) {
  ifstream count_in_bins_file(count_in_bins_filename);
  if (count_in_bins_file.is_open()) {
    string line;
    GammaParams binA_mc, binB_mc, binC_mc, binD_mc;
    GammaParams binA_data, binB_data, binC_data, binD_data;
    while ( getline(count_in_bins_file, line) ) {
      if (line.find("iplane") != std::string::npos) continue;
      else if (line.find("MC eff. entries:") != std::string::npos) {
        istringstream entries(line.substr(16));
        float binA_entries, binB_entries, binC_entries, binD_entries;
        entries >> binD_entries >> binB_entries >> binC_entries >> binA_entries;
        //cout<<line<<endl;
        //cout<<"binD_entries: "<<binD_entries<<" binB_entries: "<<binB_entries<<" binC_entries: "<<binC_entries<<" binA_entries: "<<binA_entries<<endl;
        getline(count_in_bins_file, line);
        float binA_weights, binB_weights, binC_weights, binD_weights;
        istringstream weights(line.substr(16));
        weights >> binD_weights >> binB_weights >> binC_weights >> binA_weights;
        //cout<<line<<endl;
        //cout<<"binD_weights: "<<binD_weights<<" binB_weights: "<<binB_weights<<" binC_weights: "<<binC_weights<<" binA_weights: "<<binA_weights<<endl;
        binA_mc.SetNEffectiveAndWeight(binA_entries, binA_weights);
        binB_mc.SetNEffectiveAndWeight(binB_entries, binB_weights);
        binC_mc.SetNEffectiveAndWeight(binC_entries, binC_weights);
        binD_mc.SetNEffectiveAndWeight(binD_entries, binD_weights);
      }
      else if (line.find("Data entries:") != std::string::npos) {
        istringstream entries(line.substr(13));
        float binA_entries, binB_entries, binC_entries, binD_entries;
        entries >> binD_entries >> binB_entries >> binC_entries >> binA_entries;
        //cout<<line<<endl;
        //cout<<"binD_entries: "<<binD_entries<<" binB_entries: "<<binB_entries<<" binC_entries: "<<binC_entries<<" binA_entries: "<<binA_entries<<endl;
        getline(count_in_bins_file, line);
        float binA_weights, binB_weights, binC_weights, binD_weights;
        istringstream weights(line.substr(13));
        weights >> binD_weights >> binB_weights >> binC_weights >> binA_weights;
        //cout<<line<<endl;
        //cout<<"binD_weights: "<<binD_weights<<" binB_weights: "<<binB_weights<<" binC_weights: "<<binC_weights<<" binA_weights: "<<binA_weights<<endl;
        binA_data.SetNEffectiveAndWeight(binA_entries, binA_weights);
        binB_data.SetNEffectiveAndWeight(binB_entries, binB_weights);
        binC_data.SetNEffectiveAndWeight(binC_entries, binC_weights);
        binD_data.SetNEffectiveAndWeight(binD_entries, binD_weights);

        // Save data into counts
        vector<GammaParams> mc_bins = {binA_mc, binB_mc, binC_mc, binD_mc};
        vector<GammaParams> data_bins = {binA_data, binB_data, binC_data, binD_data};
        counts.push_back({mc_bins, data_bins});
      }
    }
  } else cout<<"[Error] Unable to open "<<count_in_bins_filename<<endl;
}
void print_counts(vector<vector<vector<GammaParams> > > const & counts) {
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    cout<<"iPlane: "<<iPlane<<endl;
    cout<<"mc  (A,B,C,D): "<<counts[iPlane][0][0]<<", "<<counts[iPlane][0][1]<<", "<<counts[iPlane][0][2]<<", "<<counts[iPlane][0][3]<<", "<<endl;
    cout<<"data(A,B,C,D): "<<counts[iPlane][1][0]<<", "<<counts[iPlane][1][1]<<", "<<counts[iPlane][1][2]<<", "<<counts[iPlane][1][3]<<", "<<endl;
  }
}

void replace_counts(vector<vector<vector<GammaParams> > > const & counts_A, unsigned isData_A, unsigned isData_B, vector<vector<vector<GammaParams> > > & counts_B) {
  for (unsigned iPlane = 0; iPlane < counts_A.size(); ++iPlane) {
    counts_B[iPlane][isData_B][0] = counts_A[iPlane][isData_A][0];
    counts_B[iPlane][isData_B][1] = counts_A[iPlane][isData_A][1];
    counts_B[iPlane][isData_B][2] = counts_A[iPlane][isData_A][2];
    counts_B[iPlane][isData_B][3] = counts_A[iPlane][isData_A][3];
  }
}

void shift_kappa_counts(vector<vector<vector<GammaParams> > > const & counts, bool shiftUp, vector<vector<vector<GammaParams> > > & shifted_kappa_counts) {
  // kappa_systematics[iPlane] = fractional uncertainty
  vector<float> kappa_systematics = get_kappa_systematics();
  (void)kappa_systematics;
  // counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
  shifted_kappa_counts = counts;
  // kappas[iplane] = [0=kappa,1=kappa_down_diff,2=kappa_up_diff]
  vector<vector<float> > kappas;
  calculate_kappas(counts, kappas);
  // Shift counts by kappa 1 sigma
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    //float shiftFactor;
    //if(shiftUp) shiftFactor = kappas[iPlane][0]/(kappas[iPlane][0]+kappas[iPlane][2]);
    //else shiftFactor = kappas[iPlane][0]/(kappas[iPlane][0]-kappas[iPlane][1]);
    //shifted_kappa_counts[iPlane][0][1] = counts[iPlane][0][1]*shiftFactor;
    float shiftFactor;
    if(shiftUp) shiftFactor = kappas[iPlane][0]/(kappas[iPlane][0]*(1+kappa_systematics[iPlane]));
    else shiftFactor = kappas[iPlane][0]/(kappas[iPlane][0]*(1-kappa_systematics[iPlane]));
    shifted_kappa_counts[iPlane][0][1] = counts[iPlane][0][1]*shiftFactor;
    //cout<<iPlane<<" kappa: "<<kappas[iPlane][0]<<" +- (fractional)"<<kappa_systematics[iPlane]<<" shiftFactor: "<<shiftFactor<<endl;
    //cout<<iPlane<<" kappa_up: "<<kappas[iPlane][0]*(1+kappa_systematics[iPlane])<<" kappa_down: "<<kappas[iPlane][0]*(1-kappa_systematics[iPlane])<<endl;
    //cout<<"Before: "<<counts[iPlane][0][1].Yield()<<" After: "<<shifted_kappa_counts[iPlane][0][1].Yield()<<endl;
  }
}

// Assumes data is ordered from (plane=0,ibin=0), (plane=0,ibin=1), // ibin is split in plane
void load_fit_counts(string const & count_in_bins_filename, vector<vector<vector<GammaParams> > > & counts) {
  ifstream count_in_bins_file(count_in_bins_filename);
  if (count_in_bins_file.is_open()) {
    string line;
    GammaParams binA_data, binB_data, binC_data, binD_data;
    while ( getline(count_in_bins_file, line) ) {
      if (line.find("iplane") != std::string::npos) continue;
      else if (line.find("fit value:") != std::string::npos) {
        istringstream entries(line.substr(10));
        float binA_entries, binB_entries, binC_entries, binD_entries;
        entries >> binD_entries >> binB_entries >> binC_entries >> binA_entries;
        //cout<<line<<endl;
        //cout<<"binD_entries: "<<binD_entries<<" binB_entries: "<<binB_entries<<" binC_entries: "<<binC_entries<<" binA_entries: "<<binA_entries<<endl;
        getline(count_in_bins_file, line);
        float binA_weights, binB_weights, binC_weights, binD_weights;
        istringstream weights(line.substr(8));
        weights >> binD_weights >> binB_weights >> binC_weights >> binA_weights;
        //cout<<line<<endl;
        //cout<<"binD_weights: "<<binD_weights<<" binB_weights: "<<binB_weights<<" binC_weights: "<<binC_weights<<" binA_weights: "<<binA_weights<<endl;
        binA_data.SetYieldAndUncertainty(binA_entries, binA_weights);
        binB_data.SetYieldAndUncertainty(binB_entries, binB_weights);
        binC_data.SetYieldAndUncertainty(binC_entries, binC_weights);
        binD_data.SetYieldAndUncertainty(binD_entries, binD_weights);

        // Save data into counts
        vector<GammaParams> mc_bins = {binA_data, binB_data, binC_data, binD_data};
        vector<GammaParams> data_bins = {binA_data, binB_data, binC_data, binD_data};
        counts.push_back({mc_bins, data_bins});
      }
    }
  } else cout<<"[Error] Unable to open "<<count_in_bins_filename<<endl;
}

// counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
void make_observed_prefit_model(TRandom3 & random, vector<vector<vector<GammaParams> > > const & counts, vector<vector<vector<GammaParams> > > & observed_prefit_counts) {
  //observed_prefit_counts = counts;
  // Calculate predictions
  // predictions[iplane] = [0=prediction,1=prediction_down,2=prediction_up]
  vector<vector<float> > predictions;
  // pvalue[iplane] = pvalue for aboved observed in A
  vector<float> pvalues;
  calculate_predictions(random, counts, predictions, pvalues);

  GammaParams binA_model, binB_model, binC_model, binD_model;
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    // Get observed sideband from data
    float binB = counts[iPlane][1][1].Yield();
    float binC = counts[iPlane][1][2].Yield();
    float binD = counts[iPlane][1][3].Yield();
    float binB_unc = counts[iPlane][1][1].Uncertainty();
    float binC_unc = counts[iPlane][1][2].Uncertainty();
    float binD_unc = counts[iPlane][1][3].Uncertainty();
    // Get prefit values
    float binA = predictions[iPlane][0];
    float binA_unc = max(predictions[iPlane][1], predictions[iPlane][2]);

    binA_model.SetYieldAndUncertainty(binA, binA_unc);
    binB_model.SetYieldAndUncertainty(binB, binB_unc);
    binC_model.SetYieldAndUncertainty(binC, binC_unc);
    binD_model.SetYieldAndUncertainty(binD, binD_unc);

    //cout<<iPlane<<" prediction: "<<predictions[iPlane][0]<<" "<<predictions[iPlane][1]<<" "<<predictions[iPlane][2]<<" "<<binA_unc<<endl;
    //cout<<iPlane<<" binA unc: "<<binA_model.Uncertainty()<<endl;

    vector<GammaParams> mc_bins = {binA_model, binB_model, binC_model, binD_model};
    vector<GammaParams> data_bins = {binA_model, binB_model, binC_model, binD_model};
    observed_prefit_counts.push_back({mc_bins, data_bins});
  }
}

// counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
// kappas[iplane] = [0=kappa,1=kappa_down,2=kappa_up]
// predictions[iplane] = [0=prediction,1=prediction_down,2=prediction_up]
// pvalue[iplane] = pvalue for aboved observed in A
void printPdf_counts(vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > const & kappas, vector<vector<float> > const & predictions, vector<float> const & pvalues, vector<float> const & normalizations, string const & filename, bool print_obs_error=false) {
  vector<float> kappa_systematics = get_kappa_systematics(false);
  // Change table to vector of rows
  // | mc +- error | kappa +- error | prediction +- error | observation | Sig. |
  string table_prefix;
  string table_postfix;
  table_prefix = 
"\\documentclass[10pt,oneside]{report}\n"
"\\usepackage{graphicx,xspace,amssymb,amsmath,colordvi,colortbl,verbatim,multicol}\n"
"\\usepackage{multirow, rotating}\n"
"\\usepackage[active,tightpage]{preview}\n"
"\\begin{document}\n"
"\\begin{preview}\n";
//"\\resizebox*{!}{\\dimexpr\\textheight-2\\baselineskip\\relax}{";
//"\\resizebox*{\\textwidth}{\\dimexpr\\textheight-2\\baselineskip\\relax}{";
  table_postfix =
//"\\end{tabular}}\n"
"\\end{tabular}\n"
"\\end{preview}\n"
"\\end{document}\n";
  string table_header;
  table_header = "\\begin{tabular}{c c c c c c}\n";
  stringstream table_body;
  table_body << " & MC & $\\kappa$ & pred. & Obs. & Signi. \\\\ \n";
  table_body << "\\hline\n";
  table_body << std::fixed;
  vector<string> planeNames = {
    "$150<\\rm{MET}\\leq200, \\Delta R_{\\rm{max}} \\leq 1.1$",
    "$150<\\rm{MET}\\leq200, \\Delta R_{\\rm{max}} > 1.1$",
    "$200<\\rm{MET}\\leq300, \\Delta R_{\\rm{max}} \\leq 1.1$",
    "$200<\\rm{MET}\\leq300, \\Delta R_{\\rm{max}} > 1.1$",
    "$300<\\rm{MET}\\leq400, \\Delta R_{\\rm{max}} \\leq 1.1$",
    "$300<\\rm{MET}\\leq400, \\Delta R_{\\rm{max}} > 1.1$",
    "$400<\\rm{MET}, \\Delta R_{\\rm{max}} \\leq 1.1$",
    "$400<\\rm{MET}, \\Delta R_{\\rm{max}} > 1.1$",
  };
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    if (print_obs_error) {
      if (iPlane % 2 == 0) {
        if (normalizations.size()) table_body << "data/mc="<<setprecision(2)<<normalizations[iPlane] <<"& \\multicolumn{5}{c}{"<<planeNames[iPlane/2]<<"} \\\\ \n";
        else table_body << "& \\multicolumn{5}{c}{"<<planeNames[iPlane/2]<<"} \\\\ \n";
        table_body << "\\hline\n";
        table_body << "2b, SDB &" << setprecision(2) << counts[iPlane][0][3].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][3].Uncertainty()<<" & & & "<<counts[iPlane][1][3].Yield()<<" $\\pm$ "<<counts[iPlane][1][3].Uncertainty()<<" & \\\\ \n";
        table_body << "2b, HIG &" << setprecision(2) << counts[iPlane][0][2].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][2].Uncertainty()<<" & & & "<<counts[iPlane][1][2].Yield()<<" $\\pm$ "<<counts[iPlane][1][2].Uncertainty()<<" & \\\\ \n";
        table_body << "\\hline\n";
        table_body << "3b, SDB &" << setprecision(2) << counts[iPlane][0][1].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][1].Uncertainty()<<" & & & "<<counts[iPlane][1][1].Yield()<<" $\\pm$ "<<counts[iPlane][1][1].Uncertainty()<<" & \\\\ \n";
        table_body << "3b, HIG &" << setprecision(2) << counts[iPlane][0][0].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][0].Uncertainty()<<" & $"<<kappas[iPlane][0]<<"^{+"<<kappas[iPlane][2]<<"}_{-"<<kappas[iPlane][1]<<"}\\pm"<<kappa_systematics[iPlane]<<"$ & $"<<predictions[iPlane][0]<<"^{+"<<predictions[iPlane][2]<<"}_{-"<<predictions[iPlane][1]<<"}$ & "<<counts[iPlane][1][0].Yield()<<" $\\pm$ "<<counts[iPlane][1][0].Uncertainty()<<" & "<<setprecision(2)<<to_significance(pvalues[iPlane])<<"\\\\ \n";
        table_body << "\\hline\n";
      } else {
        table_body << "4b, SDB &" << setprecision(2) << counts[iPlane][0][1].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][1].Uncertainty()<<" & & & "<<counts[iPlane][1][1].Yield()<<" $\\pm$ "<<counts[iPlane][1][1].Uncertainty()<<" & \\\\ \n";
        table_body << "4b, HIG &" << setprecision(2) << counts[iPlane][0][0].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][0].Uncertainty()<<" & $"<<kappas[iPlane][0]<<"^{+"<<kappas[iPlane][2]<<"}_{-"<<kappas[iPlane][1]<<"}\\pm"<<kappa_systematics[iPlane]<<"$ & $"<<predictions[iPlane][0]<<"^{+"<<predictions[iPlane][2]<<"}_{-"<<predictions[iPlane][1]<<"}$ & "<<counts[iPlane][1][0].Yield()<<" $\\pm$ "<<counts[iPlane][1][0].Uncertainty()<<" & "<<setprecision(2)<<to_significance(pvalues[iPlane])<<"\\\\ \n";
        table_body << "\\hline\n";
      }
    } else {
      if (iPlane % 2 == 0) {
        if (normalizations.size()) table_body << "data/mc="<<setprecision(2)<<normalizations[iPlane] <<"& \\multicolumn{5}{c}{"<<planeNames[iPlane/2]<<"} \\\\ \n";
        else table_body << "& \\multicolumn{5}{c}{"<<planeNames[iPlane/2]<<"} \\\\ \n";
        table_body << "\\hline\n";
        table_body << "2b, SDB &" << setprecision(2) << counts[iPlane][0][3].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][3].Uncertainty()<<" & & & "<< setprecision(1) <<counts[iPlane][1][3].Yield()<<" & \\\\ \n";
        table_body << "2b, HIG &" << setprecision(2) << counts[iPlane][0][2].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][2].Uncertainty()<<" & & & "<< setprecision(1) <<counts[iPlane][1][2].Yield()<<" & \\\\ \n";
        table_body << "\\hline\n";
        table_body << "3b, SDB &" << setprecision(2) << counts[iPlane][0][1].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][1].Uncertainty()<<" & & & "<< setprecision(1) <<counts[iPlane][1][1].Yield()<<" & \\\\ \n";
        table_body << "3b, HIG &" << setprecision(2) << counts[iPlane][0][0].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][0].Uncertainty()<<" & $"<<kappas[iPlane][0]<<"^{+"<<kappas[iPlane][2]<<"}_{-"<<kappas[iPlane][1]<<"}\\pm"<<kappa_systematics[iPlane]<<"$ & $"<<predictions[iPlane][0]<<"^{+"<<predictions[iPlane][2]<<"}_{-"<<predictions[iPlane][1]<<"}$ & "<< setprecision(1) <<counts[iPlane][1][0].Yield()<<" & "<<setprecision(2)<<to_significance(pvalues[iPlane])<<"\\\\ \n";
        table_body << "\\hline\n";
      } else {
        table_body << "4b, SDB &" << setprecision(2) << counts[iPlane][0][1].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][1].Uncertainty()<<" & & & "<< setprecision(1) <<counts[iPlane][1][1].Yield()<<" & \\\\ \n";
        table_body << "4b, HIG &" << setprecision(2) << counts[iPlane][0][0].Yield() <<" $\\pm$ "<< setprecision(2) <<counts[iPlane][0][0].Uncertainty()<<" & $"<<kappas[iPlane][0]<<"^{+"<<kappas[iPlane][2]<<"}_{-"<<kappas[iPlane][1]<<"}\\pm"<<kappa_systematics[iPlane]<<"$ & $"<<predictions[iPlane][0]<<"^{+"<<predictions[iPlane][2]<<"}_{-"<<predictions[iPlane][1]<<"}$ & "<< setprecision(1) <<counts[iPlane][1][0].Yield()<<" & "<<setprecision(2)<<to_significance(pvalues[iPlane])<<"\\\\ \n";
        table_body << "\\hline\n";
      }
    }
  }

  (void)filename;
  //cout<<table_prefix+table_header+table_body.str()+table_postfix<<endl;
  std::ofstream out_file;
  out_file.open(("tables/"+filename).c_str());
  out_file <<table_prefix+table_header+table_body.str()+table_postfix<<endl;
  out_file.close();
  cout<<"pdflatex tables/"<<filename<<endl;
}

// kappas[iplane] = [0=kappa,1=kappa_down,2=kappa_up]
void calculate_kappas(vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > & kappas) {
  kappas.reserve(counts.size());
  // D, B, C, A
  vector<float> pow_kappa({ 1, -1, -1,  1});
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    vector<vector<float> > entries;
    vector<vector<float> > weights;
    entries = {{float(counts[iPlane][0][3].NEffective())}, {float(counts[iPlane][0][1].NEffective())}, {float(counts[iPlane][0][2].NEffective())}, {float(counts[iPlane][0][0].NEffective())}};
    weights = {{float(counts[iPlane][0][3].Weight())}, {float(counts[iPlane][0][1].Weight())}, {float(counts[iPlane][0][2].Weight())}, {float(counts[iPlane][0][0].Weight())}};
    //// Double uncertainty on kappa
    //if (iPlane==8) {
    //  entries = {{float(counts[iPlane][0][3].NEffective()/4)}, {float(counts[iPlane][0][1].NEffective()/4)}, {float(counts[iPlane][0][2].NEffective()/4)}, {float(counts[iPlane][0][0].NEffective()/4)}};
    //  weights = {{float(counts[iPlane][0][3].Weight()*4)}, {float(counts[iPlane][0][1].Weight()*4)}, {float(counts[iPlane][0][2].Weight()*4)}, {float(counts[iPlane][0][0].Weight()*4)}};
    //}

    //float kappa, kappa_down, kappa_up;
    //kappa = calcKappa(entries, weights, pow_kappa, kappa_down, kappa_up);
    //kappas.push_back({kappa, kappa_down, kappa_up});
    ////cout<<iPlane<<" kappa: "<<kappa<<" "<<kappa_down<<" "<<kappa_up<<endl;

    TRandom3 toy_random(1234);
    vector<float> kappa_systematics = get_kappa_systematics();
    vector<float> toy_kappas(10000);
    for (unsigned iTrial = 0; iTrial < 10000; ++iTrial) {
      float toy_mc_a = gsl_ran_gamma(entries[3][0]+1, 1, toy_random);
      float toy_mc_b = gsl_ran_gamma(entries[1][0]+1, 1, toy_random);
      float toy_mc_c = gsl_ran_gamma(entries[2][0]+1, 1, toy_random);
      float toy_mc_d = gsl_ran_gamma(entries[0][0]+1, 1, toy_random);
      float toy_kappa_mean = (toy_mc_a * weights[3][0]) / (toy_mc_b*weights[1][0]) * (toy_mc_d*weights[0][0]) / (toy_mc_c*weights[2][0]);
      float toy_kappa = toy_kappa_mean * pow(1+kappa_systematics[iPlane], toy_random.Gaus(0,1));
      toy_kappas[iTrial] = toy_kappa;
    }
    float kappa = (entries[3][0]*weights[3][0]) / (entries[1][0]*weights[1][0]) * (entries[0][0]*weights[0][0]) / (entries[2][0]*weights[2][0]);
    kappas.push_back(get_sigma_band(toy_kappas, kappa));

  }
}

float calculate_pvalue(float observed, float prediction, float prediction_up_diff, float prediction_down_diff) {
  //// Try to analytically calculate. Does not work when observed is 0.
  //(void)prediction_down_diff; (void)prediction_up_diff;
  //float test_stat = -2 * (observed * log(prediction/observed));
  //float prob = ROOT::Math::chisquared_cdf(test_stat, 1);
  //cout<<"test_stat: "<<test_stat<<" prob: "<<prob<<endl;
  //return 1-prob;

  double precision = 0.03; // Desired precision on the p-value
  double Nmin = 1/pow(precision,2), Nmax = 5e7; // Nmax controls max sigmas achievable (5e7->5.5 sigma)
  TRandom3 random(1234);

  double Nbelow=0, Nabove=0, Nequal=0;

  if (observed==0) return 1;

  while ( (std::min(Nbelow, Nabove)+Nequal)<Nmin && (Nbelow+Nabove+Nequal)<Nmax) {
    // Fluctuate mean
    double theta = random.Gaus(0,1);
    double kappa_up = 1+prediction_up_diff/prediction;
    double kappa_down = 1+prediction_down_diff/prediction;
    float fluctuated_mean;
    //// combine method
    if (prediction==0) fluctuated_mean = fabs(theta)*prediction_up_diff; //2-sided Gaussian uncertainty
    else if (theta>=0) fluctuated_mean = prediction * pow(kappa_up, theta);
    //else fluctuated_mean = prediction * pow(kappa_down, theta);
    else if(prediction_down_diff<0.8*prediction) fluctuated_mean = prediction * pow(kappa_down, theta);
    else { // Make fluctuation be a Gaussian fluctuation(larger) which saves time.
      fluctuated_mean = max(0., prediction + theta*prediction_down_diff);
    }

    //// Util's method
    //if(theta>=0) fluctuated_mean = prediction * exp(theta*log(kappa_up));
    //else if(prediction_down_diff<0.8*prediction) fluctuated_mean = prediction * exp(-theta*log(2-kappa_down));
    //else fluctuated_mean = max(0., prediction + theta*prediction_down_diff);
    // Make toy
    //float observed_toy = gsl_ran_gamma(fluctuated_mean+1, 1, random); // Mean is fluctuated_mean+1
    //observed_toy = observed_toy - 1
    //if (fluctuated_mean==0) fluctuated_mean = fabs(theta)*prediction_up_diff;
    //float observed_toy = gsl_ran_gamma(fluctuated_mean, 1, random); // Mean is fluctuated_mean
    float observed_toy = random.Poisson(fluctuated_mean); // shifts mean to high

    // Calculate test statistic

    // Method 1
    if (observed_toy>=observed) Nabove++;
    else Nbelow++;

    //// Method 2
    //if (observed_toy>observed) Nabove++;
    //else if (observed_toy==observed) Nequal++;
    //else Nbelow++;

    //// Method 3
    //float prediction_toy;
    //double theta_prediction = random.Gaus(0,1);
    //// Fluctuate prediction
    //if (prediction==0) prediction_toy = fabs(theta_prediction)*prediction_up_diff; //2-sided Gaussian uncertainty
    //if (theta>=0) prediction_toy = prediction * pow(kappa_up, theta_prediction);
    //else if(prediction_down_diff<0.8*prediction) prediction_toy = prediction * pow(kappa_down, theta_prediction);
    //else prediction_toy = max(0.,prediction + theta_prediction*prediction_down_diff); // Make fluctuation be a Gaussian fluctuation(larger) which saves time.
    //// Count
    //float test_stat = observed - prediction;
    //float test_stat_toy = observed_toy - prediction_toy;
    //if (test_stat_toy >= test_stat) Nabove++;
    ////else if (test_stat_toy == test_stat) Nequal++;
    //else Nbelow++;
  }
  //cout<<"observed: "<<observed<<" prediction: "<<prediction<<" prediction_up_diff: "<<prediction_up_diff<<" prediction_down_diff: "<<prediction_down_diff;
  //cout<<"  Nabove: "<<Nabove<<" Nequal: "<<Nequal<<" Nbelow: "<<Nbelow<<endl;

  if(Nabove+Nequal==0) return 1./Nmax;
  else if (Nbelow+Nequal==0) return 1-1./Nmax;
  return (Nabove+Nequal*1./2)/(Nabove+Nbelow+Nequal);
  //return (Nabove+Nequal)/(Nabove+Nbelow+Nequal);
}

//// Count how many toys there are above observed
//void calculate_predictions(TRandom3 & random, vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > & predictions, vector<float> & pvalues) {
//  predictions.reserve(counts.size());
//  pvalues.reserve(counts.size());
//  // D, B, C, D_mc, B_mc, C_mc, A_mc
//  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});
//  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
//    float prediction, prediction_down, prediction_up;
//    vector<vector<float> > entries;
//    vector<vector<float> > weights;
//    // Don't use odd plane D and C to calculate predictions. This correlates A and A'.
//    if (iPlane%2==0) {
//      entries = {{float(counts[iPlane][1][3].Yield())}, {float(counts[iPlane][1][1].Yield())}, {float(counts[iPlane][1][2].Yield())},
//        {float(counts[iPlane][0][3].NEffective())}, {float(counts[iPlane][0][1].NEffective())}, {float(counts[iPlane][0][2].NEffective())}, {float(counts[iPlane][0][0].NEffective())}};
//      weights = {{1}, {1}, {1},
//        {float(counts[iPlane][0][3].Weight())}, {float(counts[iPlane][0][1].Weight())}, {float(counts[iPlane][0][2].Weight())}, {float(counts[iPlane][0][0].Weight())}};
//    // Don't use odd plane D and C to calculate predictions. This correlates A and A'.
//    } else {
//      entries = {{float(counts[iPlane-1][1][3].Yield())}, {float(counts[iPlane][1][1].Yield())}, {float(counts[iPlane-1][1][2].Yield())},
//        {float(counts[iPlane-1][0][3].NEffective())}, {float(counts[iPlane][0][1].NEffective())}, {float(counts[iPlane-1][0][2].NEffective())}, {float(counts[iPlane][0][0].NEffective())}};
//      weights = {{1}, {1}, {1},
//        {float(counts[iPlane-1][0][3].Weight())}, {float(counts[iPlane][0][1].Weight())}, {float(counts[iPlane-1][0][2].Weight())}, {float(counts[iPlane][0][0].Weight())}};
//    }
//    //// Double uncertainty on kappa
//    //if (iPlane == 8) {
//    //  entries = {{float(counts[iPlane][1][3].Yield())}, {float(counts[iPlane][1][1].Yield())}, {float(counts[iPlane][1][2].Yield())},
//    //    {float(counts[iPlane][0][3].NEffective()/4)}, {float(counts[iPlane][0][1].NEffective()/4)}, {float(counts[iPlane][0][2].NEffective()/4)}, {float(counts[iPlane][0][0].NEffective()/4)}};
//    //  weights = {{1}, {1}, {1},
//    //    {float(counts[iPlane][0][3].Weight()*4)}, {float(counts[iPlane][0][1].Weight()*4)}, {float(counts[iPlane][0][2].Weight()*4)}, {float(counts[iPlane][0][0].Weight()*4)}};
//    //}
//    //// Double uncertainty on prediction
//    //if (iPlane == 8) {
//    //  prediction_down = prediction_down*2;
//    //  prediction_up = prediction_up*2;
//    //}
//    int nrep = 10000;
//    vector<float> fKappas;
//    prediction = calcKappa(random, entries, weights, pow_totpred, prediction_down, prediction_up, fKappas, /*do_data=false*/ false, /*verbose=false*/ false, /*syst=-1*/ -1, /*do_plot=false*/ false, /*nrep=100000*/ nrep, /*float nSigma=1*/ 1);
//    predictions.push_back({prediction, prediction_down, prediction_up});
//    //cout<<iPlane<<" prediction: "<<prediction<<" down: "<<prediction_down<<" up: "<<prediction_up<<endl;
//
//    //// calculate pvlaue strangly
//    //float observed_a = counts[iPlane][1][0].NEffective();
//    //// Redo kappa if not observed
//    //int nAboveObserved;
//    //int nrep_redo = nrep;
//    //for (unsigned iRedo = 0; iRedo < 4; ++iRedo) {
//    //  nAboveObserved = 0; // initialize
//    //  for(unsigned iKappa = 0; iKappa < fKappas.size(); ++iKappa) if(fKappas[iKappa]>=observed_a) nAboveObserved++;
//    //  if (nAboveObserved >=20) break; // Can make a p-value
//    //  // clear memory
//    //  // redo kappa calcualtion
//    //  nrep_redo = nrep * pow(10,iRedo+1);
//    //  //cout<<"nAboveObserved: "<<nAboveObserved<<". Redoing calcKappa with "<<nrep_redo<<endl;
//    //  calcKappa(random, entries, weights, pow_totpred, prediction_down, prediction_up, fKappas, /*do_data=false*/ false, /*verbose=false*/ false, /*syst=-1*/ -1, /*do_plot=false*/ false, /*nrep=100000*/ nrep_redo, /*float nSigma=1*/ 1);
//    //}
//    //float pvalue = nAboveObserved*1./nrep_redo;
//
//    // Use Poisson toy with mean that has log-normal fluctuation
//    float observed_a = counts[iPlane][1][0].Yield();
//    if (prediction_down>prediction) prediction_down = prediction; // Prediction can be 0 if, B or C is 0. Ignore prediction_down
//    //float significance = Significance(observed_a, prediction, prediction_up, prediction_down);
//    //float pvalue = 0.5-TMath::Erf(double(significance/sqrt(2)))/2;
//    float pvalue = calculate_pvalue(observed_a, prediction, prediction_up, prediction_down);
//    float significance = to_significance(pvalue); (void)significance;
//
//    //cout<<iPlane<<" significance: "<<significance<<" pvalue: "<<pvalue<<endl;
//
//    //cout<<iPlane<<" predicted: "<<prediction<<" -"<<prediction_down<<" +"<<prediction_up<<" observed_a: "<<observed_a<<" pvalue: "<<pvalue<<" sigma: "<<ROOT::Math::normal_quantile_c(pvalue,1)<<endl;
//    //if(prediction == 0) {
//    //  cout<<iPlane<<" predicted: "<<prediction<<" -"<<prediction_down<<" +"<<prediction_up<<" observed_a: "<<observed_a<<" pvalue: "<<pvalue<<" sigma: "<<ROOT::Math::normal_quantile_c(pvalue,1)<<endl;
//    //}
//
//    //// Use Poisson distribution to calculate significance, Can't handle case when 0 observed, 0 predicted
//    //TF1 poisson("gamma", gamma, 0,10000, 1);
//    //float observed_a = counts[iPlane][1][0].NEffective();
//    //poisson.SetParameter(0,observed_a);
//    //float pvalue = poisson.Integral(0, prediction);
//
//    //cout<<iPlane<<" predicted: "<<prediction<<" observed_a: "<<observed_a<<" pvalue: "<<pvalue<<" sigma: "<<ROOT::Math::normal_quantile_c(pvalue,1)<<endl;
//    pvalues.push_back(pvalue);
//  }
//}

// Count how many toys there are above observed
void calculate_predictions(TRandom3 & random, vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > & predictions, vector<float> & pvalues, 
  vector<TH1D*> * kappa_hists, vector<TH1D*> * prediction_hists, vector<TH1D*> * toy_hists) {
  if (prediction_hists != 0) {
    for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
      string planeName = "plane"+std::to_string(iPlane);
      kappa_hists->push_back(new TH1D(("kappa_"+planeName).c_str(), ("kappa_"+planeName).c_str(), 100, 0, 0));
      prediction_hists->push_back(new TH1D(("prediction_"+planeName).c_str(), ("prediction_"+planeName).c_str(), 100, 0, 0));
      if (iPlane == 8) {
        toy_hists->push_back(new TH1D(("toy_"+planeName).c_str(), ("toy_"+planeName).c_str(), 6, -0.5, 5.5));
      } else toy_hists->push_back(new TH1D(("toy_"+planeName).c_str(), ("toy_"+planeName).c_str(), 100, 5, 5));
    }
  }
  (void)random;
  predictions.reserve(counts.size());
  pvalues.reserve(counts.size());
  // D, B, C, D_mc, B_mc, C_mc, A_mc
  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    // entries[D,B,C,D_mc,B_mc,C_mc,A_mc] = [yield]
    // weights[D,B,C,D_mc,B_mc,C_mc,A_mc] = [yield]
    vector<vector<float> > entries;
    vector<vector<float> > weights;
    // Don't use odd plane D and C to calculate predictions. This correlates A and A'.
    if (iPlane%2==0) {
      entries = {{float(counts[iPlane][1][3].Yield())}, {float(counts[iPlane][1][1].Yield())}, {float(counts[iPlane][1][2].Yield())},
        {float(counts[iPlane][0][3].NEffective())}, {float(counts[iPlane][0][1].NEffective())}, {float(counts[iPlane][0][2].NEffective())}, {float(counts[iPlane][0][0].NEffective())}};
      weights = {{1}, {1}, {1},
        {float(counts[iPlane][0][3].Weight())}, {float(counts[iPlane][0][1].Weight())}, {float(counts[iPlane][0][2].Weight())}, {float(counts[iPlane][0][0].Weight())}};
    // Don't use odd plane D and C to calculate predictions. This correlates A and A'.
    } else {
      entries = {{float(counts[iPlane-1][1][3].Yield())}, {float(counts[iPlane][1][1].Yield())}, {float(counts[iPlane-1][1][2].Yield())},
        {float(counts[iPlane-1][0][3].NEffective())}, {float(counts[iPlane][0][1].NEffective())}, {float(counts[iPlane-1][0][2].NEffective())}, {float(counts[iPlane][0][0].NEffective())}};
      weights = {{1}, {1}, {1},
        {float(counts[iPlane-1][0][3].Weight())}, {float(counts[iPlane][0][1].Weight())}, {float(counts[iPlane-1][0][2].Weight())}, {float(counts[iPlane][0][0].Weight())}};
    }

    float observed_a = counts[iPlane][1][0].Yield();

    // Use entries and weights to get prediction (and uncertainty) and pvalue
    TRandom3 toy_random(1234);
    double precision = 0.03; // Desired precision on the p-value
    double Nmin = 1/pow(precision,2), Nmax = 5e7; // Nmax controls max sigmas achievable (5e7->5.5 sigma)
    float pvalue = -1;
    // prediction_info = [prediction, prediction_down_diff, prediction_up_diff]
    vector<float> prediction_info = {-1, -1, -1};

    if (observed_a==0) Nmax=Nmin; // For calculating prediction distribution

    double Nbelow=0, Nabove=0;
    vector<float> toy_predictions; toy_predictions.reserve(Nmax);

    // kappa = A_mc / B_mc * D_mc / C_mc
    float kappa = (entries[6][0]*weights[6][0]) / (entries[4][0]*weights[4][0]) * (entries[3][0]*weights[3][0]) / (entries[5][0]*weights[5][0]);
    // kappa_systematics[iPlane] = fractional uncertainty
    vector<float> kappa_systematics = get_kappa_systematics();

    //// Print input
    //cout<<"----------------"<<endl;
    //cout<<"Input to significance for plane: "<<iPlane<<endl;
    //cout<<"A_obs: "<<observed_a<<" B_obs: "<<entries[1][0]<<" C_obs: "<<entries[2][0]<<" D_obs: "<<entries[0][0]<<endl;
    //cout<<"A_mc: "<<entries[6][0]<<" B_obs: "<<entries[4][0]<<" C_obs: "<<entries[5][0]<<" D_obs: "<<entries[3][0]<<endl;
    //cout<<"A_weight: "<<weights[6][0]<<" B_obs: "<<weights[4][0]<<" C_obs: "<<weights[5][0]<<" D_obs: "<<weights[3][0]<<endl;
    //cout<<"kappa_systematics: "<<kappa_systematics[iPlane]<<endl;
    //cout<<"----------------"<<endl;

    // Generate toys
    while ( (std::min(Nbelow, Nabove))<Nmin && (Nbelow+Nabove)<Nmax) {
      // Generate B, C, D toy yields
      float toy_b = gsl_ran_gamma(entries[1][0]+1, 1, toy_random);
      float toy_c = gsl_ran_gamma(entries[2][0]+1, 1, toy_random);
      float toy_d = gsl_ran_gamma(entries[0][0]+1, 1, toy_random);
      // Generate kappa with systematics from statistics. Covers MET systematics.
      float toy_mc_a = gsl_ran_gamma(entries[6][0]+1, 1, toy_random);
      float toy_mc_b = gsl_ran_gamma(entries[4][0]+1, 1, toy_random);
      float toy_mc_c = gsl_ran_gamma(entries[5][0]+1, 1, toy_random);
      float toy_mc_d = gsl_ran_gamma(entries[3][0]+1, 1, toy_random);
      float toy_kappa_mean = (toy_mc_a * weights[6][0]) / (toy_mc_b*weights[4][0]) * (toy_mc_d*weights[3][0]) / (toy_mc_c*weights[5][0]);
      //float toy_kappa = toy_kappa_mean;
      // Add in other systematics to kappa with log-normal
      float toy_kappa = toy_kappa_mean * pow(1+kappa_systematics[iPlane], toy_random.Gaus(0,1));
      //float toy_kappa = kappa; (void)toy_kappa_mean;

      // Calculate toy prediction
      float toy_prediction = (toy_b*weights[1][0]) * (toy_c*weights[2][0]) / (toy_d*weights[0][0]) * toy_kappa;
      toy_predictions.push_back(toy_prediction);
      // Generate toy_observed using toy prediction
      float toy_observed = toy_random.Poisson(toy_prediction);

      if (prediction_hists != 0) {
        kappa_hists->at(iPlane)->Fill(toy_kappa);
        prediction_hists->at(iPlane)->Fill(toy_prediction);
        toy_hists->at(iPlane)->Fill(toy_observed);
      }

      if (toy_observed >= observed_a) Nabove++;
      else Nbelow++;
    }

    //cout<<"b: "<<entries[1][0]<<" c: "<<entries[2][0]<<" d: "<<entries[0][0]<<" kappa: "<<kappa<<" kappa uncertainty: "<<kappa_systematics[iPlane]<<" obs: "<<observed_a<<" Nabove: "<<Nabove<<" Nbelow: "<<Nbelow<<endl;
    //if (prediction_hists != 0) {
    //  cout<<toy_hists->at(iPlane)->GetEntries()<<endl;
    //}
    
    // Calculate pvalue
    {
      if (observed_a==0) pvalue = 1;
      else if (Nabove==0) pvalue = 1./Nmax;
      else if (Nbelow==0) pvalue = 1-1./Nmax;
      else pvalue = Nabove/(Nabove+Nbelow);
    }

    // Calculate prediction uncertainty
    {
      float prediction = kappa * entries[1][0] * entries[2][0] / entries[0][0];
      prediction_info = get_sigma_band(toy_predictions, prediction);
      //cout<<"iPlane "<<iPlane<<" "<<prediction_info[0]<<" "<<prediction_info[1]<<" "<<prediction_info[2]<<endl;

      //sort(toy_predictions.begin(), toy_predictions.end());
      //int prediction_index = -1;
      //for (int iToy = 0; iToy < static_cast<int>(toy_predictions.size()); iToy++) {
      //  if (toy_predictions[iToy] > prediction) {
      //    prediction_index = iToy;
      //    break;
      //  }
      //}
      //float nSigma = 1;
      //double integratedGaus = intGaus(0,1,0,nSigma);
      //int prediction_low_index = prediction_index - static_cast<int>(integratedGaus * toy_predictions.size());
      //if (prediction_low_index<0) prediction_low_index = 0;
      //int prediction_high_index = prediction_index + static_cast<int>(integratedGaus * toy_predictions.size());
      //if (prediction_high_index>static_cast<int>(toy_predictions.size())) prediction_high_index = toy_predictions.size()+1;
      //float prediction_low = prediction - toy_predictions[prediction_low_index];
      //float prediction_high = toy_predictions[prediction_high_index] - prediction;
      //prediction_info = {prediction, prediction_low, prediction_high};
      ////cout<<"prediction: "<<prediction<<" toy_predictions[prediction_index]: "<<toy_predictions[prediction_index]<<endl;
      ////cout<<"prediction: "<<prediction<<" low: "<<toy_predictions[prediction_low_index]<<" high: "<<toy_predictions[prediction_high_index]<<endl;
      ////cout<<"prediction index pred: "<<toy_predictions[prediction_index]<<endl;
      ////cout<<"prediction: "<<prediction<<" pvalue: "<<pvalue<<" significance: "<<to_significance(pvalue)<<endl;
    }

    pvalues.push_back(pvalue);
    predictions.push_back(prediction_info);
  }
}

void calculate_test_stats(vector<vector<vector<GammaParams> > > const & counts, vector<vector<float> > const & predictions, vector<float> & test_stats) {
  test_stats.reserve(predictions.size());
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    float prediction = predictions[iPlane][0];
    float prediction_uncertainty = max(predictions[iPlane][1], predictions[iPlane][2]);
    float observed = counts[iPlane][1][0].Yield();
    //float test_stat = fabs(observed-prediction)/sqrt(prediction_uncertainty*prediction_uncertainty+observed);
    float test_stat = fabs(observed-prediction)/sqrt(prediction_uncertainty*prediction_uncertainty);
    // Above Prediction distribution is not a gaussian so test_stat doesn't work that well.
    // When calculating prediction, give observed value to also calculate p-value.
    //cout<<iPlane<<" prediction: "<<prediction<<" "<<prediction_uncertainty<<" observed: "<<observed<<" test_stat: "<<test_stat<<endl;
    test_stats.push_back(test_stat);
  }
}

float get_max_test_stat(vector<float> const & test_stats) {
  float max_test_stat = 0;
  for (unsigned iPlane = 0; iPlane < test_stats.size(); ++iPlane) if (max_test_stat < test_stats[iPlane]) max_test_stat = test_stats[iPlane];
  return max_test_stat;
}

float get_max_significance(vector<float> const & pvalues) {
  float max_significance = 0;
  for (unsigned iPlane = 0; iPlane < pvalues.size(); ++iPlane) {
    float significance;
    if (pvalues[iPlane] == 0) significance = 99;
    else if (pvalues[iPlane] == 1) significance = 0;
    else significance = ROOT::Math::normal_quantile_c(pvalues[iPlane],1);

    if (max_significance < significance) max_significance = significance;
    //cout<<iPlane<<" pvalue: "<<pvalues[iPlane]<<" significance: "<<significance<<endl;
  }
  return max_significance;
}

void normalize_mc_to_data(vector<vector<vector<GammaParams> > > const & counts, vector<vector<vector<GammaParams> > > & normalized_counts, vector<float> & normalizations) {
  normalized_counts = counts;
  normalizations.resize(counts.size());
  // Calculate data/mc factor
  GammaParams total_mc;
  GammaParams total_data;
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    total_mc += counts[iPlane][0][0] + counts[iPlane][0][1] +counts[iPlane][0][2] +counts[iPlane][0][3];
    total_data += counts[iPlane][1][0] + counts[iPlane][1][1] +counts[iPlane][1][2] +counts[iPlane][1][3];
  }
  float data_over_mc = total_data.Yield()/total_mc.Yield();
  //float data_over_mc_unc = sqrt(pow(sqrt(total_data.Yield())/total_mc.Yield(),2) + pow(total_data.Yield()*total_mc.Uncertainty()/total_mc.Yield()/total_mc.Yield(),2));
  //cout<<"total data: "<<total_data.Yield()<<" total_mc: "<<total_mc<<" data_over_mc: "<<data_over_mc<<" "<<data_over_mc_unc<<endl;
  // Apply normalization on mc
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    normalized_counts[iPlane][0][0] *= data_over_mc;
    normalized_counts[iPlane][0][1] *= data_over_mc;
    normalized_counts[iPlane][0][2] *= data_over_mc;
    normalized_counts[iPlane][0][3] *= data_over_mc;
    normalizations[iPlane] = data_over_mc;
  }
}

void normalize_mc_to_data_by_met(vector<vector<vector<GammaParams> > > const & counts, vector<vector<vector<GammaParams> > > & normalized_counts, vector<float> & normalizations) {
  normalized_counts = counts;
  normalizations.resize(counts.size());
  unsigned factor = 4; // 4 is by met, 2 is by plane
  // Calculate data/mc factor
  vector<float> data_over_mc(counts.size()/factor);
  vector<GammaParams> total_mc(counts.size()/factor);
  vector<GammaParams> total_data(counts.size()/factor);
  //GammaParams total_mc;
  //GammaParams total_data;
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    if (iPlane%2 == 0) {
      total_mc[iPlane/factor] += counts[iPlane][0][0] + counts[iPlane][0][1] +counts[iPlane][0][2] +counts[iPlane][0][3];
      total_data[iPlane/factor] += counts[iPlane][1][0] + counts[iPlane][1][1] +counts[iPlane][1][2] +counts[iPlane][1][3];
    } else { 
      total_mc[iPlane/factor] += counts[iPlane][0][0] + counts[iPlane][0][1];
      total_data[iPlane/factor] += counts[iPlane][1][0] + counts[iPlane][1][1];
    }
    if (iPlane%factor == (factor-1)) {
      data_over_mc[iPlane/factor] = total_data[iPlane/factor].Yield()/total_mc[iPlane/factor].Yield();
      //cout<<total_mc[iPlane/factor]<<" "<<total_data[iPlane/factor]<<endl;
    }
  }
  //float data_over_mc_unc = sqrt(pow(sqrt(total_data.Yield())/total_mc.Yield(),2) + pow(total_data.Yield()*total_mc.Uncertainty()/total_mc.Yield()/total_mc.Yield(),2));
  //cout<<"total data: "<<total_data.Yield()<<" total_mc: "<<total_mc<<" data_over_mc: "<<data_over_mc<<" "<<data_over_mc_unc<<endl;
  // Apply normalization on mc
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    normalized_counts[iPlane][0][0] *= data_over_mc[iPlane/factor];
    normalized_counts[iPlane][0][1] *= data_over_mc[iPlane/factor];
    normalized_counts[iPlane][0][2] *= data_over_mc[iPlane/factor];
    normalized_counts[iPlane][0][3] *= data_over_mc[iPlane/factor];
    normalizations[iPlane] = data_over_mc[iPlane/factor];
  }
}

float generateYield(TRandom3 & random, GammaParams binB, GammaParams binC, GammaParams binD, 
  GammaParams binA_mc, GammaParams binB_mc, GammaParams binC_mc, GammaParams binD_mc,
  float kappa_systematics,
  float & true_value
) {
  float toy_b = gsl_ran_gamma(binB.Yield()+1, 1, random);
  float toy_c = gsl_ran_gamma(binC.Yield()+1, 1, random);
  float toy_d = gsl_ran_gamma(binD.Yield()+1, 1, random);
  // Generate kappa with systematics from statistics
  float toy_mc_a = gsl_ran_gamma(binA_mc.NEffective()+1, 1, random);
  float toy_mc_b = gsl_ran_gamma(binB_mc.NEffective()+1, 1, random);
  float toy_mc_c = gsl_ran_gamma(binC_mc.NEffective()+1, 1, random);
  float toy_mc_d = gsl_ran_gamma(binD_mc.NEffective()+1, 1, random);
  float toy_kappa_mean = (toy_mc_a * binA_mc.Weight()) / (toy_mc_b*binB_mc.Weight()) * (toy_mc_d*binD_mc.Weight()) / (toy_mc_c*binC_mc.Weight());
  float toy_kappa = toy_kappa_mean * pow(1+kappa_systematics, random.Gaus(0,1));
  //float toy_kappa = binA_mc.Yield() / binB_mc.Yield() * binD_mc.Yield() / binC_mc.Yield(); (void)kappa_systematics;
  // Calculate toy prediction
  float toy_prediction = (toy_b) * (toy_c) / (toy_d) * toy_kappa;
  true_value = toy_prediction;
  // Generate toy_observed using toy prediction
  float toy_observed = random.Poisson(toy_prediction);
  return toy_observed;
}

// mode: 0 poisson
// mode: 1 poisson x Gamma distribution
// mode: 2 poisson x gaussian with unc. from poisson (Model post-fit)
// mode: 3 poisson x gaussian with unc. multiplied by factor
// mode: 4 poisson x log-normal 
// uncorrelated_uncertinaty and correlated_uncertainty work only with mode 4
float generateYield(TRandom3 & random, float mean, float uncertainty, int mode, 
  float & true_value,
  float uncorrelated_uncertinaty=0, float correlated_uncertainty=0, float correlated_theta=0,
  float weight=1 // used for gamma function
) {
  if (mode == 0) {
    true_value = mean;
    return random.Poisson(mean);
  }
  // Fluctuate mean with Gamma
  else if (mode == 1) {
    //float fluctuated_mean=-1;
    //fluctuated_mean= gsl_ran_gamma(mean+1, 1, random);

    //// Additional uncertinaty following log-normal
    //double theta = random.Gaus(0,1);
    //double kappa = 1+uncertainty/mean/weight; 
    //double uncorrelated_theta = random.Gaus(0,1);

    //true_value = fluctuated_mean*weight*pow(kappa,theta)*pow(1+correlated_uncertainty,correlated_theta)*pow(1+uncorrelated_uncertinaty, uncorrelated_theta);
    ////cout<<"mean: "<<mean<<" uncertainty: "<<uncertainty<<" theta: "<<theta<<" uncertainty/mean: "<<uncertainty/mean/weight<<" theta: "<<theta<<" fluctuated_mean: "<<true_value<<endl;
    //return random.Poisson(true_value);
    (void) weight;
    true_value = mean;
    return gsl_ran_gamma(mean+1, 1, random);
  }
  // Fluctuate mean with gaussian, where gaussian uncertainty is sqrt(mean) [Try to model post-fit]
  else if (mode == 2) {
    float gaussian_unc = sqrt(mean);
    float fluctuated_mean=-1;
    while (fluctuated_mean<0) fluctuated_mean= random.Gaus(mean, gaussian_unc);
    true_value = fluctuated_mean;
    return random.Poisson(fluctuated_mean);
  }
  // Fluctuate mean with gaussian
  else if (mode == 3) {
    // Include systematics
    float fluctuated_mean=-1;
    double theta = random.Gaus(0,1);
    double uncorrelated_theta = random.Gaus(0,1);
    fluctuated_mean = mean + mean*(theta*(uncertainty/mean)) + mean*(uncorrelated_theta*uncorrelated_uncertinaty) + mean*(correlated_theta*correlated_uncertainty);
    if(fluctuated_mean<=0) fluctuated_mean = 0.001;
    true_value = fluctuated_mean;
    //cout<<"mean: "<<mean<<" theta: "<<theta<<" uncertainty/mean: "<<uncertainty/mean<<" uncorrelated_theta: "<<uncorrelated_theta<<" uncorrelated_uncertinaty: "<<uncorrelated_uncertinaty<<" correlated_theta: "<<correlated_theta<<" correlated_uncertainty: "<<correlated_uncertainty<<" fluctuated_mean: "<<fluctuated_mean<<endl;

    //cout<<"mean: "<<mean<<" uncertainty: "<<uncertainty<<" fluctuated_mean: "<<fluctuated_mean<<endl;

    return random.Poisson(fluctuated_mean);
  }
  // Use log normal to generate fluctuated_mean
  else if (mode == 4) {
    // https://cds.cern.ch/record/1379837/files/NOTE2011_005.pdf
    double theta = random.Gaus(0,1);

    // First method
    double kappa = 1+uncertainty/mean; 
    if (mean ==0) kappa = 1;
    double uncorrelated_theta = random.Gaus(0,1);
    //// Original 
    //double fluctuated_mean = mean*pow(kappa,theta);
    //// Example to show multiplication works
    //double fluctuated_mean = mean*pow(kappa,theta)*pow(1+correlated_uncertainty,correlated_theta)*pow(kappa, uncorrelated_theta);
    double fluctuated_mean = mean*pow(kappa,theta)*pow(1+correlated_uncertainty,correlated_theta)*pow(1+uncorrelated_uncertinaty, uncorrelated_theta);
    //// Second method
    //double kappa = 1+(sqrt(pow(uncertainty,2)+pow(uncorrelated_uncertinaty*mean,2)))/mean; 
    //double fluctuated_mean = mean*pow(kappa,theta)*pow(1+correlated_uncertainty,correlated_theta);

    true_value = fluctuated_mean;

    //cout<<"mean: "<<mean<<" uncertainty: "<<uncertainty<<" kappa: "<<kappa<<" theta: "<<theta<<" fluctuated_mean: "<<fluctuated_mean<<endl;

    return random.Poisson(fluctuated_mean);
  }
  else return 0;
}

//float generateYield(TRandom3 & random, float mean, float uncertainty, int mode, 
//  float & true_value,
//  float uncorrelated_uncertinaty=0, float correlated_uncertainty=0, float correlated_theta=0,
//  float weight=1 // used for gamma function
//) {
//}

// correlated_uncertainty is in %.
// counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
void generateToy (vector<vector<vector<GammaParams> > > const & counts, int seed, 
    vector<vector<vector<GammaParams> > > & counts_toy, vector<vector<float> > & predictions_toy, vector<float> & pvalues_toy, float & max_test_stat_toy,
    vector<vector<vector<float> > > & counts_mean_toy, 
    int model_mode = 0, float factor_uncerainty = 0,
    float uncorrelated_uncertinaty = 0, float correlated_uncertainty = 0) 
{
  (void)uncorrelated_uncertinaty; (void)correlated_uncertainty;
  TRandom3 random_toy(seed);
  // counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
  counts_toy = counts; // Copy structure
  counts_mean_toy = vector<vector<vector<float> > > (counts.size(), vector<vector<float> > (counts[0].size(), vector<float>(counts[0][0].size()) ) );
  float true_A, true_B, true_C, true_D;
  float generation_mode = 4; // 1 is gamma, 3 is gaussian, 4 is log-normal
  (void)generation_mode;
  float correlated_theta = random_toy.Gaus(0,1);
  (void)correlated_theta;
  // Replace observed with Poisson(MC)
  for (unsigned iPlane = 0; iPlane < counts.size(); ++iPlane) {
    float obs_A = counts[iPlane][1][0].Yield();
    float obs_B = counts[iPlane][1][1].Yield();
    float obs_C = counts[iPlane][1][2].Yield();
    float obs_D = counts[iPlane][1][3].Yield();
    (void)obs_A;(void)obs_B;(void)obs_C;(void)obs_D;
    float obs_A_unc = counts[iPlane][1][0].Uncertainty();
    float obs_B_unc = counts[iPlane][1][1].Uncertainty();
    float obs_C_unc = counts[iPlane][1][2].Uncertainty();
    float obs_D_unc = counts[iPlane][1][3].Uncertainty();
    (void)obs_A_unc;(void)obs_B_unc;(void)obs_C_unc;(void)obs_D_unc;
    float obs_A_neffective = counts[iPlane][1][0].NEffective();
    float obs_B_neffective = counts[iPlane][1][1].NEffective();
    float obs_C_neffective = counts[iPlane][1][2].NEffective();
    float obs_D_neffective = counts[iPlane][1][3].NEffective();
    (void)obs_A_neffective;(void)obs_B_neffective;(void)obs_C_neffective;(void)obs_D_neffective;
    float obs_A_weight = counts[iPlane][1][0].Weight();
    float obs_B_weight = counts[iPlane][1][1].Weight();
    float obs_C_weight = counts[iPlane][1][2].Weight();
    float obs_D_weight = counts[iPlane][1][3].Weight();
    (void)obs_A_weight;(void)obs_B_weight;(void)obs_C_weight;(void)obs_D_weight;
    //float factor = sqrt(1); // sqrt(3) is to include systematics
    //(void)factor;

    if (model_mode==0) {
      // Use Poisson x gamma (MC stats) to generate toy yields. 
      counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, obs_A_neffective, obs_A_unc*factor_uncerainty, /*mode*/ 1, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_A_weight), 1);
      counts_toy[iPlane][1][1].SetNEffectiveAndWeight(generateYield(random_toy, obs_B_neffective, obs_B_unc*factor_uncerainty, /*mode*/ 1, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_B_weight), 1);
      counts_toy[iPlane][1][2].SetNEffectiveAndWeight(generateYield(random_toy, obs_C_neffective, obs_C_unc*factor_uncerainty, /*mode*/ 1, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_C_weight), 1);
      counts_toy[iPlane][1][3].SetNEffectiveAndWeight(generateYield(random_toy, obs_D_neffective, obs_D_unc*factor_uncerainty, /*mode*/ 1, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_D_weight), 1);
      //// Use Poisson to generate toy yields
      //counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, obs_A, obs_A_unc*factor_uncerainty, /*mode*/ 0, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_A_weight), 1);
      //counts_toy[iPlane][1][1].SetNEffectiveAndWeight(generateYield(random_toy, obs_B, obs_B_unc*factor_uncerainty, /*mode*/ 0, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_B_weight), 1);
      //counts_toy[iPlane][1][2].SetNEffectiveAndWeight(generateYield(random_toy, obs_C, obs_C_unc*factor_uncerainty, /*mode*/ 0, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_C_weight), 1);
      //counts_toy[iPlane][1][3].SetNEffectiveAndWeight(generateYield(random_toy, obs_D, obs_D_unc*factor_uncerainty, /*mode*/ 0, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_D_weight), 1);
    } else if (model_mode==1) {
      // Use Poisson x log-normal(Uncertainty) to generate toy yields
      counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, obs_A, obs_A_unc*factor_uncerainty, /*mode*/ 4, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_A_weight), 1);
      counts_toy[iPlane][1][1].SetNEffectiveAndWeight(generateYield(random_toy, obs_B, obs_B_unc*factor_uncerainty, /*mode*/ 4, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_B_weight), 1);
      counts_toy[iPlane][1][2].SetNEffectiveAndWeight(generateYield(random_toy, obs_C, obs_C_unc*factor_uncerainty, /*mode*/ 4, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_C_weight), 1);
      counts_toy[iPlane][1][3].SetNEffectiveAndWeight(generateYield(random_toy, obs_D, obs_D_unc*factor_uncerainty, /*mode*/ 4, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_D_weight), 1);
    } else if (model_mode==2) {
      // Use Poisson x log-normal(Uncertainty) to generate toy yields
      vector<float> kappa_systematics = get_kappa_systematics();
      //counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, obs_A, obs_A_unc*factor_uncerainty, /*mode*/ 4, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_A_weight), 1);
      counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, counts[iPlane][1][1], counts[iPlane][1][2], counts[iPlane][1][3], counts[iPlane][0][0], counts[iPlane][0][1], counts[iPlane][0][2], counts[iPlane][0][3], kappa_systematics[iPlane], true_A),1);
      counts_toy[iPlane][1][1].SetNEffectiveAndWeight(generateYield(random_toy, obs_B, obs_B_unc*factor_uncerainty, /*mode*/ 1, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_B_weight), 1);
      counts_toy[iPlane][1][2].SetNEffectiveAndWeight(generateYield(random_toy, obs_C, obs_C_unc*factor_uncerainty, /*mode*/ 1, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_C_weight), 1);
      counts_toy[iPlane][1][3].SetNEffectiveAndWeight(generateYield(random_toy, obs_D, obs_D_unc*factor_uncerainty, /*mode*/ 1, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_D_weight), 1);
      //cout<<"toy A: "<<counts_toy[iPlane][1][0].Yield()<<" "<<counts_toy[iPlane][1][0].NEffective()<<endl;
      //// Use Poisson to generate toy yields
      //counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, obs_A, obs_A_unc*factor_uncerainty, /*mode*/ 0, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_A_weight), 1);
      //counts_toy[iPlane][1][1].SetNEffectiveAndWeight(generateYield(random_toy, obs_B, obs_B_unc*factor_uncerainty, /*mode*/ 0, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_B_weight), 1);
      //counts_toy[iPlane][1][2].SetNEffectiveAndWeight(generateYield(random_toy, obs_C, obs_C_unc*factor_uncerainty, /*mode*/ 0, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_C_weight), 1);
      //counts_toy[iPlane][1][3].SetNEffectiveAndWeight(generateYield(random_toy, obs_D, obs_D_unc*factor_uncerainty, /*mode*/ 0, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta, obs_D_weight), 1);
    }

    //// Use Poisson x log-normal(statistical uncertianty+systematics uncertainty) to generate toy yields
    //counts_toy[iPlane][1][0].SetNEffectiveAndWeight(generateYield(random_toy, obs_A, obs_A_unc*factor_uncerainty, /*mode*/ generation_mode, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), 1);
    //counts_toy[iPlane][1][1].SetNEffectiveAndWeight(generateYield(random_toy, obs_B, obs_B_unc*factor_uncerainty, /*mode*/ generation_mode, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), 1);
    //counts_toy[iPlane][1][2].SetNEffectiveAndWeight(generateYield(random_toy, obs_C, obs_C_unc*factor_uncerainty, /*mode*/ generation_mode, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), 1);
    //counts_toy[iPlane][1][3].SetNEffectiveAndWeight(generateYield(random_toy, obs_D, obs_D_unc*factor_uncerainty, /*mode*/ generation_mode, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), 1);

    // Fill true value
    counts_mean_toy[iPlane][1][0] = true_A;
    counts_mean_toy[iPlane][1][1] = true_B;
    counts_mean_toy[iPlane][1][2] = true_C;
    counts_mean_toy[iPlane][1][3] = true_D;

    //cout<<iPlane<<" obs_A:     "<<obs_A<<    " obs_B:     "<<obs_B    <<" obs_C:     "<<obs_C    <<" obs_D: "    <<obs_D<<endl;
    //cout<<iPlane<<" obs_A_unc: "<<obs_A_unc<<" obs_B_unc: "<<obs_B_unc<<" obs_C_unc: "<<obs_C_unc<<" obs_D_unc: "<<obs_D_unc<<endl;


    //// Also Change MC.
    float mc_A = counts[iPlane][0][0].Yield();
    float mc_B = counts[iPlane][0][1].Yield();
    float mc_C = counts[iPlane][0][2].Yield();
    float mc_D = counts[iPlane][0][3].Yield();
    (void)mc_A;(void)mc_B;(void)mc_C;(void)mc_D;
    float mc_A_unc = counts[iPlane][0][0].Uncertainty();
    float mc_B_unc = counts[iPlane][0][1].Uncertainty();
    float mc_C_unc = counts[iPlane][0][2].Uncertainty();
    float mc_D_unc = counts[iPlane][0][3].Uncertainty();
    (void)mc_A_unc;(void)mc_B_unc;(void)mc_C_unc;(void)mc_D_unc;
    float mc_A_neffective = counts[iPlane][0][0].NEffective();
    float mc_B_neffective = counts[iPlane][0][1].NEffective();
    float mc_C_neffective = counts[iPlane][0][2].NEffective();
    float mc_D_neffective = counts[iPlane][0][3].NEffective();
    (void)mc_A_neffective;(void)mc_B_neffective;(void)mc_C_neffective;(void)mc_D_neffective;
    float mc_A_weight = counts[iPlane][0][0].Weight();
    float mc_B_weight = counts[iPlane][0][1].Weight();
    float mc_C_weight = counts[iPlane][0][2].Weight();
    float mc_D_weight = counts[iPlane][0][3].Weight();
    (void)mc_A_weight;(void)mc_B_weight;(void)mc_C_weight;(void)mc_D_weight;

    // Case when not changing MC
    // Fill true value
    counts_mean_toy[iPlane][0][0] = mc_A_neffective*mc_A_weight;
    counts_mean_toy[iPlane][0][1] = mc_B_neffective*mc_B_weight;
    counts_mean_toy[iPlane][0][2] = mc_C_neffective*mc_C_weight;
    counts_mean_toy[iPlane][0][3] = mc_D_neffective*mc_D_weight;
    //// Case when model is from MC
    //counts_toy[iPlane][0][0].SetNEffectiveAndWeight(generateYield(random_toy, mc_A_neffective, mc_A_unc/mc_A_weight*factor_uncerainty, /*mode*/ generation_mode, true_A, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), mc_A_weight);
    //counts_toy[iPlane][0][1].SetNEffectiveAndWeight(generateYield(random_toy, mc_B_neffective, mc_B_unc/mc_B_weight*factor_uncerainty, /*mode*/ generation_mode, true_B, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), mc_B_weight);
    //counts_toy[iPlane][0][2].SetNEffectiveAndWeight(generateYield(random_toy, mc_C_neffective, mc_C_unc/mc_C_weight*factor_uncerainty, /*mode*/ generation_mode, true_C, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), mc_C_weight);
    //counts_toy[iPlane][0][3].SetNEffectiveAndWeight(generateYield(random_toy, mc_D_neffective, mc_D_unc/mc_D_weight*factor_uncerainty, /*mode*/ generation_mode, true_D, uncorrelated_uncertinaty, correlated_uncertainty, correlated_theta), mc_D_weight);
    //// Fill true value
    //counts_mean_toy[iPlane][0][0] = true_A*mc_A_weight;
    //counts_mean_toy[iPlane][0][1] = true_B*mc_B_weight;
    //counts_mean_toy[iPlane][0][2] = true_C*mc_C_weight;
    //counts_mean_toy[iPlane][0][3] = true_D*mc_D_weight;
    //// Case when model is from fitted toys
    //counts_toy[iPlane][0][0].SetNEffectiveAndWeight(generateYield(random_toy, mc_A, mc_A_unc*factor_uncerainty, /*mode*/ 4, true_A), 1);
    //counts_toy[iPlane][0][1].SetNEffectiveAndWeight(generateYield(random_toy, mc_B, mc_B_unc*factor_uncerainty, /*mode*/ 4, true_B), 1);
    //counts_toy[iPlane][0][2].SetNEffectiveAndWeight(generateYield(random_toy, mc_C, mc_C_unc*factor_uncerainty, /*mode*/ 4, true_C), 1);
    //counts_toy[iPlane][0][3].SetNEffectiveAndWeight(generateYield(random_toy, mc_D, mc_D_unc*factor_uncerainty, /*mode*/ 4, true_D), 1);
    //// Fill true value
    //counts_mean_toy[iPlane][0][0] = true_A;
    //counts_mean_toy[iPlane][0][1] = true_B;
    //counts_mean_toy[iPlane][0][2] = true_C;
    //counts_mean_toy[iPlane][0][3] = true_D;

  }

  // Calculate max_test_stat
  calculate_predictions(random_toy, counts_toy, predictions_toy, pvalues_toy);
  //// test_stats[iplane] = test_stat
  //vector<float> test_stats_toy;
  //calculate_test_stats(counts_toy, predictions_toy, test_stats_toy);
  //// max_test_stat
  //float max_test_stat_toy = get_max_test_stat(test_stats_toy);
  max_test_stat_toy = get_max_significance(pvalues_toy);

}

void exampleToy(int seed, float & max_test_stat_toy) {
  TRandom3 random(seed);
  max_test_stat_toy = random.Poisson(1);
}

void makeTransparent(TCanvas & canvas) {
  canvas.SetFillStyle(4000);
  //canvas.SetFrameFillStyle(4000);
  //canvas.SetFrameFillColor(0);
  //canvas.SetFrameBorderMode(0);
}

// Returns fractional systematics
vector<float> get_kappa_systematics(bool fractional) {
  vector<float> kappa(16);
  kappa[0] = 1.04805;
  kappa[1] = 1.26201;
  kappa[2] = 1.08626;
  kappa[3] = 1.13387;
  kappa[4] = 1.04333;
  kappa[5] = 1.21014;
  kappa[6] = 0.922012;
  kappa[7] = 0.955641;
  kappa[8] = 0.715363;
  kappa[9] = 2.35375;
  kappa[10] = 0.944108;
  kappa[11] = 0.891454;
  kappa[12] = 1.23673;
  kappa[13] = 0.938487;
  kappa[14] = 0.975312;
  kappa[15] = 0.918603;



  // controlSystematics[iPlane][control] = fractional systematic
  vector<map<string, float> > controlSystematics(16);
  controlSystematics[0]["ttbar"] = 1.11; 
  controlSystematics[1]["ttbar"] = 1.18;
  controlSystematics[2]["ttbar"] = 1.02;
  controlSystematics[3]["ttbar"] = 1.07;
  controlSystematics[4]["ttbar"] = 1.10;
  controlSystematics[5]["ttbar"] = 1.17;
  controlSystematics[6]["ttbar"] = 1.02;
  controlSystematics[7]["ttbar"] = 1.07;
  controlSystematics[8]["ttbar"] = 1.06;
  controlSystematics[9]["ttbar"] = 1.13;
  controlSystematics[10]["ttbar"] = 1.01;
  controlSystematics[11]["ttbar"] = 1.06;
  controlSystematics[12]["ttbar"] = 1.06;
  controlSystematics[13]["ttbar"] = 1.09;
  controlSystematics[14]["ttbar"] = 1.01;
  controlSystematics[15]["ttbar"] = 1.07;

  controlSystematics[0]["vjets"] = 1.01; 
  controlSystematics[1]["vjets"] = 1.01;
  controlSystematics[2]["vjets"] = 1.00;
  controlSystematics[3]["vjets"] = 1.00;
  controlSystematics[4]["vjets"] = 1.04;
  controlSystematics[5]["vjets"] = 1.05;
  controlSystematics[6]["vjets"] = 1.01;
  controlSystematics[7]["vjets"] = 1.00;
  controlSystematics[8]["vjets"] = 1.10;
  controlSystematics[9]["vjets"] = 1.06;
  controlSystematics[10]["vjets"] = 1.01;
  controlSystematics[11]["vjets"] = 1.01;
  controlSystematics[12]["vjets"] = 1.06;
  controlSystematics[13]["vjets"] = 1.10;
  controlSystematics[14]["vjets"] = 1.02;
  controlSystematics[15]["vjets"] = 1.02;

  controlSystematics[0]["qcd"] = 1.01; 
  controlSystematics[1]["qcd"] = 1.00;
  controlSystematics[2]["qcd"] = 1.00;
  controlSystematics[3]["qcd"] = 1.00;
  controlSystematics[4]["qcd"] = 1.00;
  controlSystematics[5]["qcd"] = 1.02;
  controlSystematics[6]["qcd"] = 1.00;
  controlSystematics[7]["qcd"] = 1.00;
  controlSystematics[8]["qcd"] = 1.00;
  controlSystematics[9]["qcd"] = 1.00;
  controlSystematics[10]["qcd"] = 1.00;
  controlSystematics[11]["qcd"] = 1.00;
  controlSystematics[12]["qcd"] = 1.00;
  controlSystematics[13]["qcd"] = 1.00;
  controlSystematics[14]["qcd"] = 1.00;
  controlSystematics[15]["qcd"] = 1.01;

  vector<float> kappa_systematics(controlSystematics.size()); 
  for (unsigned iPlane = 0; iPlane < controlSystematics.size(); ++iPlane) {
    float kappa_ttbar_systematic = kappa[iPlane] * (controlSystematics[iPlane]["ttbar"]-1);
    float kappa_zjets_systematic = kappa[iPlane] * (controlSystematics[iPlane]["vjets"]-1);
    float kappa_qcd_systematic = kappa[iPlane] * (controlSystematics[iPlane]["qcd"]-1);
    kappa_systematics[iPlane] = sqrt(pow(kappa_ttbar_systematic,2)+pow(kappa_zjets_systematic,2)+pow(kappa_qcd_systematic,2));
    if(fractional) kappa_systematics[iPlane] = kappa_systematics[iPlane]/kappa[iPlane];
    //cout<<iPlane<<" ttbar: "<<kappa_ttbar_systematic<<" zjets: "<<kappa_zjets_systematic<<" qcd: "<<kappa_qcd_systematic<<" kappa systematic: "<<kappa_systematics[iPlane]<<" kappa: "<<kappa[iPlane]<<endl;
    //cout<<iPlane<<" ttbar: "<<controlSystematics[iPlane]["ttbar"]<<" vjets: "<<controlSystematics[iPlane]["vjets"]<<" qcd: "<<controlSystematics[iPlane]["qcd"]<<" kappa_systematics: "<<kappa_systematics[iPlane]<<" kappa: "<<kappa[iPlane]<<endl;
  }
  return kappa_systematics;
}

int main(/* int argc, char *argv[] */){
  std::chrono::high_resolution_clock::time_point begTime;
  begTime = std::chrono::high_resolution_clock::now();

  // Load data
  // counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
  vector<vector<vector<GammaParams> > > counts;
  // Can generate counts_in_bins.txt by two methods
  // 1. scripts/datacard_to_counts.py
  // 2. ./run/higgsino/plot_kappas.exe --sample search --year run2 --unblind --scen data --string_options save_entries_weights_to_file 
  //string count_in_bins_filename = "./txt/counts_in_bins/counts_in_bins.txt";
  //string count_in_bins_filename = "./txt/counts_in_bins/counts_in_bins_mc.txt";
  string count_in_bins_filename = "./txt/counts_in_bins/counts_in_bins_correctBTagSF.txt";
  load_counts(count_in_bins_filename, counts);
  //print_counts(counts);

  // Calculate observed test-stat
  // kappas[iplane] = [0=kappa,1=kappa_down_diff,2=kappa_up_diff]
  vector<vector<float> > kappas;
  calculate_kappas(counts, kappas);
  // predictions[iplane] = [0=prediction,1=prediction_down_diff,2=prediction_up_diff]
  vector<vector<float> > predictions;
  // pvalue[iplane] = pvalue for aboved observed in A
  vector<float> pvalues;
  TRandom3 random;
  vector<TH1D*> kappa_hists;
  vector<TH1D*> prediction_hists;
  vector<TH1D*> toy_hists;
  calculate_predictions(random, counts, predictions, pvalues, &kappa_hists, &prediction_hists, &toy_hists);

  // Draw kappa and prediction distributions
  for (unsigned iPlane=0; iPlane < counts.size(); ++iPlane) {
    TCanvas * canvas_plane = new TCanvas(("plane"+std::to_string(iPlane)).c_str(), ("plane"+std::to_string(iPlane)).c_str(), 1500, 500);
    canvas_plane->Divide(3,1);
    canvas_plane->cd(1);
    kappa_hists[iPlane]->SetTitle(("#kappa for bin "+std::to_string(iPlane)+";#kappa;# toys").c_str());
    kappa_hists[iPlane]->Draw();
    vector<float> kappa_sigma_band = get_sigma_band(kappa_hists[iPlane]);
    //cout<<"kappa["<<iPlane<<"]: "<<kappa_sigma_band[1]-kappa_sigma_band[0]<<" "<<kappa_sigma_band[2]-kappa_sigma_band[1]<<endl;
    TLine * kappa_m1 = new TLine(kappa_sigma_band[0], 0, kappa_sigma_band[0], kappa_hists[iPlane]->GetMaximum()); kappa_m1->Draw();
    TLine * kappa_p1 = new TLine(kappa_sigma_band[2], 0, kappa_sigma_band[2], kappa_hists[iPlane]->GetMaximum()); kappa_p1->Draw();
    canvas_plane->cd(2);
    prediction_hists[iPlane]->SetTitle(("Pre-fit prediction for bin "+std::to_string(iPlane)+";Pre-fit prediction;# toys").c_str());
    prediction_hists[iPlane]->Draw();
    vector<float> prediction_sigma_band = get_sigma_band(prediction_hists[iPlane]);
    TLine * prediction_m1 = new TLine(prediction_sigma_band[0], 0, prediction_sigma_band[0], prediction_hists[iPlane]->GetMaximum()); prediction_m1->Draw();
    TLine * prediction_p1 = new TLine(prediction_sigma_band[2], 0, prediction_sigma_band[2], prediction_hists[iPlane]->GetMaximum()); prediction_p1->Draw();
    //cout<<"prediction["<<iPlane<<"]: "<<prediction_sigma_band[1]-prediction_sigma_band[0]<<" "<<prediction_sigma_band[2]-prediction_sigma_band[1]<<endl;
    canvas_plane->cd(3);
    gPad->SetLogy();
    toy_hists[iPlane]->SetTitle(("Toy yields for bin A in bin "+std::to_string(iPlane)+";Toy yield in bin A;# toys").c_str());
    toy_hists[iPlane]->Draw();
    float observed_yield = counts[iPlane][1][0].Yield()-0.5;
    //vector<float> toy_sigma_band = get_sigma_band(toy_hists[iPlane]);
    TLine * observed_line = new TLine(observed_yield, 0, observed_yield, toy_hists[iPlane]->GetMaximum()); observed_line->Draw();
    //TLine * toy_m1 = new TLine(toy_sigma_band[0], 0, toy_sigma_band[0], toy_hists[iPlane]->GetMaximum()); toy_m1->Draw();
    //TLine * toy_p1 = new TLine(toy_sigma_band[2], 0, toy_sigma_band[2], toy_hists[iPlane]->GetMaximum()); toy_p1->Draw();
    canvas_plane->SaveAs(("plots/kappa_prediction_toy_plane_"+std::to_string(iPlane)+".pdf").c_str());
    cout<<"open plots/kappa_prediction_toy_plane_"+std::to_string(iPlane)+".pdf"<<endl;
  }

  //// test_stats[iplane] = test_stat
  //vector<float> test_stats;
  //calculate_test_stats(counts, predictions, test_stats);
  //// max_test_stat
  //float max_test_stat = get_max_test_stat(test_stats);
  float max_test_stat = get_max_significance(pvalues);
  cout<<"[Observed] max_test_stat: "<<max_test_stat<<" p-value: "<<to_pvalue(max_test_stat)<<endl;
  int number_planes=16;
  cout<<"  expected global significance: "<<to_significance(to_pvalue(max_test_stat)*number_planes)<< " expected p-value: "<<to_pvalue(max_test_stat)*number_planes<<endl;

  printPdf_counts(counts, kappas, predictions, pvalues, vector<float>(), "counts_base.tex");

  TH1F * hist_significance = new TH1F("hist_significance","hist_significance", 10, -5, 5);
  TH1F * hist_significance_3b = new TH1F("hist_significance_3b","hist_significance_3b", 10, -5, 5);
  TH1F * hist_significance_4b = new TH1F("hist_significance_4b","hist_significance_4b", 10, -5, 5);
  for (unsigned iPlane = 0; iPlane<pvalues.size(); ++iPlane) {
    float significance;
    //cout<<iPlane<<" pvalue: "<<pvalues[iPlane]<<endl;
    significance = to_significance(pvalues[iPlane]);
    hist_significance->Fill(significance);
    if (iPlane%2==0) hist_significance_3b->Fill(significance);
    if (iPlane%2==1) hist_significance_4b->Fill(significance);
  }

  // Generate toys
  bool generate_toys = true;
  bool use_threads = true;
  unsigned nToys = 10000;
  float histNBin = 10000;
  float histMax = 100;

  // 0: MC, 1: background-only fit, 2: b,c,d observed+prefit (default)
  int model_mode = 2;
  float correlated_uncertainty = 0.0; //fraction ex) 0.1 is 10%
  float uncorrelated_uncertinaty = 0.00; //fraction: ex) 0.05 is 5%
  float factor_uncerainty = 0; // Uncertainty factor of log-normal with statistical uncertainty ex) 0 is no log-normal uncertianty. Used for Gamma model_mode.
  if (model_mode != 0) factor_uncerainty = 1; // Use log-normal uncertainty
  int shift_kappa_by_systematics = -1; // -1: no shift, 0: shift down, 1: shift up

  int probe_plane = 8;
  int probe_max_a = 10;
  int probe_max_b = 30;
  int probe_max_c = 30;
  int probe_max_d = 100;

  // Model related
  // Normalize mc to data
  // counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
  vector<vector<vector<GammaParams> > > model_counts;
  // normalization[iplane] = normalization
  model_counts = counts;
  vector<float > normalizations;
  //normalize_mc_to_data(counts, model_counts, normalizations);
  //normalize_mc_to_data_by_met(counts, model_counts, normalizations);
  //print_counts(model_counts);
  if (model_mode == 0) {
    // Replace data with MC for model
    replace_counts(model_counts, 0, 1, model_counts);
  } else if (model_mode == 1) {
    // Load fit data
    vector<vector<vector<GammaParams> > > fit_counts;
    //load_fit_counts("counts_in_bins_fit.txt", fit_counts); // background only
    //load_fit_counts("counts_in_bins_fit_signal500.txt", fit_counts);
    //load_fit_counts("counts_in_bins_fit_signal450.txt", fit_counts); // r is 0.15
    load_fit_counts("counts_in_bins_fit_signal450_r1.txt", fit_counts); // Takes too much resources
    //load_fit_counts("counts_in_bins_fit_signal450_r0p5.txt", fit_counts);
    //load_fit_counts("counts_in_bins_fit_signal450_r0p3.txt", fit_counts);
    replace_counts(counts, 0, 0, fit_counts);
    model_counts = fit_counts;
  } else if (model_mode == 2) {
    // Made model with observed+prefit. Can use counts to do it.
    vector<vector<vector<GammaParams> > > observed_prefit_counts;
    make_observed_prefit_model(random, counts, observed_prefit_counts);
    replace_counts(counts, 0, 0, observed_prefit_counts); // Replace MC
    model_counts = observed_prefit_counts;
  }

  if (shift_kappa_by_systematics!=-1) {
    // Shift kappa by one sigma
    vector<vector<vector<GammaParams> > > shifted_kappa_counts;
    shift_kappa_counts(model_counts, shift_kappa_by_systematics, shifted_kappa_counts);
    model_counts = shifted_kappa_counts;
  }

  vector<vector<float> > model_kappas;
  calculate_kappas(model_counts, model_kappas);
  vector<vector<float> > model_predictions;
  vector<float> model_pvalues;
  calculate_predictions(random, model_counts, model_predictions, model_pvalues);
  printPdf_counts(model_counts, model_kappas, model_predictions, model_pvalues, normalizations, "model_counts_base.tex", /*print_obs_error*/ true);

  vector<thread> threads;
  // Common storage for threads
  vector<vector<vector<vector<GammaParams> > > > counts_toy_threads(nToys);
  vector<vector<vector<vector<float> > > > counts_true_mean_toy_threads(nToys);
  vector<vector<vector<float> > > predictions_toy_threads(nToys);
  vector<vector<float> > pvalues_toy_threads(nToys);
  vector<float> max_test_stat_toy_threads(nToys);

  if (generate_toys) {
    TH1F * hist_probe_a = new TH1F("hist_probe_a","hist_probe_a", 100, 0, probe_max_a);
    TH1F * hist_probe_b = new TH1F("hist_probe_b","hist_probe_b", 100, 0, probe_max_b);
    TH1F * hist_probe_c = new TH1F("hist_probe_c","hist_probe_c", 100, 0, probe_max_c);
    TH1F * hist_probe_d = new TH1F("hist_probe_d","hist_probe_d", 100, 0, probe_max_d);
    TH1F * hist_probe_a_mc = new TH1F("hist_probe_a_mc","hist_probe_a_mc", 100, 0, probe_max_a);
    TH1F * hist_probe_b_mc = new TH1F("hist_probe_b_mc","hist_probe_b_mc", 100, 0, probe_max_b);
    TH1F * hist_probe_c_mc = new TH1F("hist_probe_c_mc","hist_probe_c_mc", 100, 0, probe_max_c);
    TH1F * hist_probe_d_mc = new TH1F("hist_probe_d_mc","hist_probe_d_mc", 100, 0, probe_max_d);
    TH1F * hist_probe_true_a = new TH1F("hist_probe_true_a","hist_probe_true_a", 100, 0, probe_max_a);
    TH1F * hist_probe_true_b = new TH1F("hist_probe_true_b","hist_probe_true_b", 100, 0, probe_max_b);
    TH1F * hist_probe_true_c = new TH1F("hist_probe_true_c","hist_probe_true_c", 100, 0, probe_max_c);
    TH1F * hist_probe_true_d = new TH1F("hist_probe_true_d","hist_probe_true_d", 100, 0, probe_max_d);
    TH1F * hist_probe_true_a_mc = new TH1F("hist_probe_true_a_mc","hist_probe_true_a_mc", 100, 0, probe_max_a);
    TH1F * hist_probe_true_b_mc = new TH1F("hist_probe_true_b_mc","hist_probe_true_b_mc", 100, 0, probe_max_b);
    TH1F * hist_probe_true_c_mc = new TH1F("hist_probe_true_c_mc","hist_probe_true_c_mc", 100, 0, probe_max_c);
    TH1F * hist_probe_true_d_mc = new TH1F("hist_probe_true_d_mc","hist_probe_true_d_mc", 100, 0, probe_max_d);
    TH1F * hist_probe_prediction = new TH1F("hist_probe_prediction", "hist_probe_prediction", 100, 0, probe_max_a);
    TH1F * hist_probe_significance = new TH1F("hist_probe_significance", "hist_probe_significance", 100, -5, 5);

    TH1F * hist_max_test_stat = new TH1F("hist_max_test_stat", "hist_max_test_stat", histNBin, 0, histMax);
    TH1F * hist_significance_toy = new TH1F("hist_significance_toy","hist_significance_toy", 10, -5, 5);
    TH1F * hist_significance_toy_nozeroobs = new TH1F("hist_significance_toy_nozeroobs","hist_significance_toy_nozeroobs", 10, -5, 5);
    TH1F * hist_significance_toy_3b = new TH1F("hist_significance_toy_3b","hist_significance_toy_3b", 10, -5, 5);
    TH1F * hist_significance_toy_4b = new TH1F("hist_significance_toy_4b","hist_significance_toy_4b", 10, -5, 5);
    TH1F * hist_delta_toy = new TH1F("hist_delta_toy","hist_delta_toy", 20, -1, 1);


    if (use_threads) {
      // Use multiple threads to run toys
      for (unsigned iToy = 0; iToy < nToys; ++iToy) {
        // input: model_counts, iToy
        // output: counts_toy, predictions_toy, pvalues_toy, max_test_stat_toy
        threads.push_back( thread(generateToy, model_counts, iToy+1, 
              ref(counts_toy_threads[iToy]), ref(predictions_toy_threads[iToy]), ref(pvalues_toy_threads[iToy]), ref(max_test_stat_toy_threads[iToy]), ref(counts_true_mean_toy_threads[iToy]), 
              model_mode, factor_uncerainty, uncorrelated_uncertinaty, correlated_uncertainty));
      };
      for (auto & thread : threads) thread.join();
      // Collect toy result
      for (unsigned iToy = 0; iToy < nToys; ++iToy) {
        //cout<<"iToy: "<<iToy<<" max_test_stat_toy: "<<max_test_stat_toy_threads[iToy]<<endl;
        hist_max_test_stat->Fill(max_test_stat_toy_threads[iToy]);
        hist_probe_a->Fill(counts_toy_threads[iToy][probe_plane][1][0].Yield());
        hist_probe_b->Fill(counts_toy_threads[iToy][probe_plane][1][1].Yield());
        hist_probe_c->Fill(counts_toy_threads[iToy][probe_plane][1][2].Yield());
        hist_probe_d->Fill(counts_toy_threads[iToy][probe_plane][1][3].Yield());
        hist_probe_a_mc->Fill(counts_toy_threads[iToy][probe_plane][0][0].Yield());
        hist_probe_b_mc->Fill(counts_toy_threads[iToy][probe_plane][0][1].Yield());
        hist_probe_c_mc->Fill(counts_toy_threads[iToy][probe_plane][0][2].Yield());
        hist_probe_d_mc->Fill(counts_toy_threads[iToy][probe_plane][0][3].Yield());
        hist_probe_true_a->Fill(counts_true_mean_toy_threads[iToy][probe_plane][1][0]);
        hist_probe_true_b->Fill(counts_true_mean_toy_threads[iToy][probe_plane][1][1]);
        hist_probe_true_c->Fill(counts_true_mean_toy_threads[iToy][probe_plane][1][2]);
        hist_probe_true_d->Fill(counts_true_mean_toy_threads[iToy][probe_plane][1][3]);
        hist_probe_true_a_mc->Fill(counts_true_mean_toy_threads[iToy][probe_plane][0][0]);
        hist_probe_true_b_mc->Fill(counts_true_mean_toy_threads[iToy][probe_plane][0][1]);
        hist_probe_true_c_mc->Fill(counts_true_mean_toy_threads[iToy][probe_plane][0][2]);
        hist_probe_true_d_mc->Fill(counts_true_mean_toy_threads[iToy][probe_plane][0][3]);
        hist_probe_prediction->Fill(predictions_toy_threads[iToy][probe_plane][0]);
        hist_probe_significance->Fill(ROOT::Math::normal_quantile_c(pvalues_toy_threads[iToy][probe_plane],1));
        // counts[iplane][0=mc/1=data][0=A,1=B,2=C,3=D] = GammaParams
        for (unsigned iPlane = 0; iPlane<pvalues_toy_threads[iToy].size(); ++iPlane) {
          float significance;
          significance = ROOT::Math::normal_quantile_c(pvalues_toy_threads[iToy][iPlane],1);
          hist_significance_toy->Fill(significance);
          if (iPlane%2==0) hist_significance_toy_3b->Fill(significance);
          if (iPlane%2==1) hist_significance_toy_4b->Fill(significance);
          // 0 yields makes histogram skew to left(less observed). Because 0 effect.
          // Above 1 makes histogram skew to right(more observed). Probably because observation is Poisson.
          if ( counts_toy_threads[iToy][iPlane][1][0].Yield()>=1 && iPlane==3) {
            hist_significance_toy_nozeroobs->Fill(significance);
          }
          if (counts_toy_threads[iToy][iPlane][1][0].Yield()>=1) {
            hist_delta_toy->Fill((counts_toy_threads[iToy][iPlane][1][0].Yield()-predictions_toy_threads[iToy][iPlane][0])/counts_toy_threads[iToy][iPlane][1][0].Yield());
          }
        }
      }
      vector<vector<float> > kappas_toy;
      calculate_kappas(counts_toy_threads[0], kappas_toy);
      printPdf_counts(counts_toy_threads[0], kappas_toy, predictions_toy_threads[0], pvalues_toy_threads[0], vector<float>(), "toy_counts.tex");
    } else {
      // No thread version
      for (unsigned iToy = 0; iToy < nToys; ++iToy) {
        // Generate toy
        vector<vector<vector<GammaParams> > >counts_toy;
        vector<vector<vector<float> > >counts_mean_toy;
        vector<vector<float> > predictions_toy;
        vector<float> pvalues_toy;
        vector<vector<float> > kappas_toy;
        float max_test_stat_toy;
        generateToy(model_counts, iToy, counts_toy, predictions_toy, pvalues_toy, max_test_stat_toy, counts_mean_toy, 
          model_mode, factor_uncerainty, uncorrelated_uncertinaty, correlated_uncertainty);
        calculate_kappas(counts_toy, kappas_toy);
        cout<<"iToy: "<<iToy<<" max_test_stat_toy: "<<max_test_stat_toy<<endl;
        //print_counts(counts_toy);
        //printPdf_counts(counts_toy, kappas_toy, predictions_toy, pvalues_toy, vector<float>(), "toy_counts.tex");
        hist_max_test_stat->Fill(max_test_stat_toy);
      }
    }

    // Save results
    TFile test_stat("test_stat.root", "RECREATE");
    hist_max_test_stat->Write();
    hist_probe_a->Write();
    hist_probe_b->Write();
    hist_probe_c->Write();
    hist_probe_d->Write();
    hist_probe_a_mc->Write();
    hist_probe_b_mc->Write();
    hist_probe_c_mc->Write();
    hist_probe_d_mc->Write();
    hist_probe_true_a->Write();
    hist_probe_true_b->Write();
    hist_probe_true_c->Write();
    hist_probe_true_d->Write();
    hist_probe_true_a_mc->Write();
    hist_probe_true_b_mc->Write();
    hist_probe_true_c_mc->Write();
    hist_probe_true_d_mc->Write();
    hist_probe_prediction->Write();
    hist_probe_significance->Write();
    hist_significance_toy->Write();
    hist_significance_toy_3b->Write();
    hist_significance_toy_4b->Write();
    hist_significance_toy_nozeroobs->Write();
    hist_significance->Write();
    hist_delta_toy->Write();
    // TODO: Delete above histograms
  }

  TFile test_stat("test_stat.root");
  TH1F * hist_max_test_stat = static_cast<TH1F*> (test_stat.Get("hist_max_test_stat"));
  TH1F * hist_significance_toy = static_cast<TH1F*> (test_stat.Get("hist_significance_toy"));
  TH1F * hist_significance_toy_3b = static_cast<TH1F*> (test_stat.Get("hist_significance_toy_3b"));
  TH1F * hist_significance_toy_4b = static_cast<TH1F*> (test_stat.Get("hist_significance_toy_4b"));
  TH1F * hist_significance_toy_nozeroobs = static_cast<TH1F*> (test_stat.Get("hist_significance_toy_nozeroobs"));
  TH1F * hist_delta_toy = static_cast<TH1F*> (test_stat.Get("hist_delta_toy"));

  TH1F * hist_probe_a = static_cast<TH1F*> (test_stat.Get("hist_probe_a"));
  TH1F * hist_probe_b = static_cast<TH1F*> (test_stat.Get("hist_probe_b"));
  TH1F * hist_probe_c = static_cast<TH1F*> (test_stat.Get("hist_probe_c"));
  TH1F * hist_probe_d = static_cast<TH1F*> (test_stat.Get("hist_probe_d"));
  TH1F * hist_probe_a_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_a_mc"));
  TH1F * hist_probe_b_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_b_mc"));
  TH1F * hist_probe_c_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_c_mc"));
  TH1F * hist_probe_d_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_d_mc"));
  TH1F * hist_probe_true_a = static_cast<TH1F*> (test_stat.Get("hist_probe_true_a"));
  TH1F * hist_probe_true_b = static_cast<TH1F*> (test_stat.Get("hist_probe_true_b"));
  TH1F * hist_probe_true_c = static_cast<TH1F*> (test_stat.Get("hist_probe_true_c"));
  TH1F * hist_probe_true_d = static_cast<TH1F*> (test_stat.Get("hist_probe_true_d"));
  //TH1F * hist_probe_true_a_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_true_a_mc"));
  //TH1F * hist_probe_true_b_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_true_b_mc"));
  //TH1F * hist_probe_true_c_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_true_c_mc"));
  //TH1F * hist_probe_true_d_mc = static_cast<TH1F*> (test_stat.Get("hist_probe_true_d_mc"));
  TH1F * hist_probe_prediction = static_cast<TH1F*> (test_stat.Get("hist_probe_prediction"));
  TH1F * hist_probe_significance = static_cast<TH1F*> (test_stat.Get("hist_probe_significance"));

  // Calculate significance
  float number_above_toys = hist_max_test_stat->Integral(int(histNBin/histMax*max_test_stat),histNBin);
  cout<<"Number of toys: "<<hist_max_test_stat->Integral(0,histNBin)<<endl;
  cout<<"Number of toy above "<<int(histNBin/histMax*max_test_stat)*histMax/histNBin<<": "<<number_above_toys<<endl;
  float pvalue = number_above_toys/hist_max_test_stat->Integral(0,histNBin);
  float pvalue_up = (sqrt(number_above_toys)+number_above_toys)/hist_max_test_stat->Integral(0,histNBin);
  float pvalue_down = (-sqrt(number_above_toys)+number_above_toys)/hist_max_test_stat->Integral(0,histNBin);
  float sigma = ROOT::Math::normal_quantile_c(pvalue,1);
  float sigma_up = ROOT::Math::normal_quantile_c(pvalue_up,1);
  float sigma_down = ROOT::Math::normal_quantile_c(pvalue_down,1);
  cout<<"Sigma: "<<sigma<<" +"<<sigma_down-sigma<<" "<<sigma_up-sigma<<endl;
  cout<<"  pvalue: "<<pvalue<<" +"<<pvalue_up-pvalue<<" "<<pvalue_down-pvalue<<endl;

  TCanvas c1("c1", "c1", 500, 500);
  c1.SetLogy();
  hist_max_test_stat->SetTitle("Max. significance between 16 bins;Max. significance (#sigma);# toys/0.1 #sigma");
  hist_max_test_stat->Rebin(10);
  hist_max_test_stat->GetXaxis()->SetRangeUser(0, 10);
  hist_max_test_stat->Draw();
  TLine * observed_test_stat = new TLine(max_test_stat,0,max_test_stat,hist_max_test_stat->GetMaximum()); 
  observed_test_stat->SetLineColor(kRed);observed_test_stat->Draw();
  c1.SaveAs("plots/test_stat.pdf");
  cout<<"open plots/test_stat.pdf"<<endl;

  gStyle->SetOptStat(111111);              // No Stats box
  TCanvas c2("c2", "c2", 500, 500);
  c2.Divide(2,2);
  c2.cd(1);
  hist_probe_a->SetStats(1);
  hist_probe_a->Draw();
  c2.cd(2);
  hist_probe_b->SetStats(1);
  hist_probe_b->Draw();
  c2.cd(3);
  hist_probe_c->SetStats(1);
  hist_probe_c->Draw();
  c2.cd(4);
  hist_probe_d->SetStats(1);
  hist_probe_d->Draw();
  c2.SaveAs("plots/test_stat_toy_yields.pdf");
  cout<<"open plots/test_stat_toy_yields.pdf"<<endl;

  TCanvas c3("c3", "c3", 500, 500);
  c3.Divide(2,2);
  c3.cd(1);
  hist_probe_a_mc->SetStats(1);
  hist_probe_a_mc->Draw();
  c3.cd(2);
  hist_probe_b_mc->SetStats(1);
  hist_probe_b_mc->Draw();
  c3.cd(3);
  hist_probe_c_mc->SetStats(1);
  hist_probe_c_mc->Draw();
  c3.cd(4);
  hist_probe_d_mc->SetStats(1);
  hist_probe_d_mc->Draw();
  c3.SaveAs("plots/test_stat_toy_yields_mc.pdf");
  cout<<"open plots/test_stat_toy_yields_mc.pdf"<<endl;

  TCanvas c4("c4", "c4", 500, 500);
  c4.Divide(2,2);
  c4.cd(1);
  hist_probe_true_a->SetStats(1);
  hist_probe_true_a->Draw();
  c4.cd(2);
  hist_probe_true_b->SetStats(1);
  hist_probe_true_b->Draw();
  c4.cd(3);
  hist_probe_true_c->SetStats(1);
  hist_probe_true_c->Draw();
  c4.cd(4);
  hist_probe_true_d->SetStats(1);
  hist_probe_true_d->Draw();
  c4.SaveAs("plots/test_stat_toy_true_yields.pdf");
  cout<<"open plots/test_stat_toy_true_yields.pdf"<<endl;

  TCanvas c5("c5", "c5", 500, 500);
  gStyle->SetOptStat(220001111);              // No Stats box
  hist_significance_toy->SetStats(1);
  hist_significance_toy->SetLineColor(kRed);
  TH1F* hist_significance_toy_norm = static_cast<TH1F*>(hist_significance_toy->DrawNormalized());
  c5.Update();
  TPaveStats * stat_significance_toy_norm =  static_cast<TPaveStats*>(hist_significance_toy_norm->FindObject("stats"));
  stat_significance_toy_norm->SetY1NDC(0.3);
  stat_significance_toy_norm->SetY2NDC(0.6);
  hist_significance_toy_norm->Scale(16);
  hist_significance_toy_norm->SetMaximum(6.5);
  //hist_significance_toy_nozeroobs->SetLineColor(kGreen);
  //hist_significance_toy_nozeroobs->DrawNormalized("same");
  (void)hist_significance_toy_nozeroobs;
  hist_significance->SetStats(1);
  hist_significance->Draw("sames hist");
  //hist_significance->DrawNormalized("sames hist");
  //hist_significance->DrawNormalized("same")->Scale(0.5);
  c5.Modified(); c5.Update();
  c5.SaveAs("plots/hist_significances.pdf");
  cout<<"open plots/hist_significances.pdf"<<endl;

  TCanvas c9("c9", "c9", 500, 500);
  gStyle->SetOptStat(220001111);              // No Stats box
  hist_significance_toy_3b->SetStats(1);
  hist_significance_toy_3b->SetLineColor(kRed);
  TH1F* hist_significance_toy_norm_3b = static_cast<TH1F*>(hist_significance_toy_3b->DrawNormalized());
  c9.Update();
  TPaveStats * stat_significance_toy_norm_3b =  static_cast<TPaveStats*>(hist_significance_toy_norm_3b->FindObject("stats"));
  stat_significance_toy_norm_3b->SetY1NDC(0.3);
  stat_significance_toy_norm_3b->SetY2NDC(0.6);
  hist_significance_toy_norm_3b->Scale(8);
  hist_significance_toy_norm_3b->SetMaximum(3.5);
  hist_significance_3b->SetStats(1);
  hist_significance_3b->Draw("sames hist");
  //hist_significance_3b->DrawNormalized("sames hist");
  //hist_significance->DrawNormalized("same")->Scale(0.5);
  c9.Modified(); c9.Update();
  c9.SaveAs("plots/hist_significances_3b.pdf");
  cout<<"open plots/hist_significances_3b.pdf"<<endl;

  TCanvas c10("c10", "c10", 500, 500);
  gStyle->SetOptStat(220001111);              // No Stats box
  hist_significance_toy_4b->SetStats(1);
  hist_significance_toy_4b->SetLineColor(kRed);
  TH1F* hist_significance_toy_norm_4b = static_cast<TH1F*>(hist_significance_toy_4b->DrawNormalized());
  c10.Update();
  TPaveStats * stat_significance_toy_norm_4b =  static_cast<TPaveStats*>(hist_significance_toy_norm_4b->FindObject("stats"));
  stat_significance_toy_norm_4b->SetY1NDC(0.3);
  stat_significance_toy_norm_4b->SetY2NDC(0.6);
  hist_significance_toy_norm_4b->Scale(8);
  hist_significance_toy_norm_4b->SetMaximum(3.5);
  hist_significance_4b->SetStats(1);
  hist_significance_4b->Draw("sames hist");
  //hist_significance->DrawNormalized("same")->Scale(0.5);
  c10.Modified(); c10.Update();
  c10.SaveAs("plots/hist_significances_4b.pdf");
  cout<<"open plots/hist_significances_4b.pdf"<<endl;

  gStyle->SetOptStat(220001111);              // No Stats box
  TCanvas c6("c6", "c6", 500, 500);
  hist_delta_toy->SetStats(1);
  hist_delta_toy->Draw();
  c6.SaveAs("plots/hist_delta_toy.pdf");
  cout<<"open plots/hist_delta_toy.pdf"<<endl;

  TCanvas c7("c7", "c7", 500, 500);
  makeTransparent(c7);
  hist_probe_a->Draw();
  hist_probe_prediction->SetLineColor(kRed);
  hist_probe_prediction->Scale(5);
  hist_probe_prediction->Draw("same hist");
  c7.SaveAs("plots/hist_probe_prediction.pdf");
  cout<<"open plots/hist_probe_prediction.pdf"<<endl;

  TCanvas c8("c8", "c8", 500, 500);
  makeTransparent(c8);
  gStyle->SetOptStat(220001111);              // No Stats box
  hist_probe_significance->SetStats(1);
  hist_probe_significance->Draw();
  c8.SaveAs("plots/hist_probe_significance.pdf");
  cout<<"open plots/hist_probe_significance.pdf"<<endl;
  float probe_all = hist_probe_significance->Integral(0,100);
  float probe_above = hist_probe_significance->Integral(int(100./10*to_significance(pvalues[probe_plane]))+50, 100);
  float probe_pvalue = probe_above / probe_all;
  cout<<"Number of significance: "<<probe_all<<endl;
  cout<<"Number of significance above significance: "<<probe_above<<endl;
  cout<<"Calculated significance: "<<to_significance(pvalues[probe_plane])<<" toys significance: "<<to_significance(probe_pvalue)<<endl;

  double seconds = (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;

  return 0;
}
