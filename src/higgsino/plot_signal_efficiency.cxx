#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TStyle.h"
#include <iostream>
#include <set>
#include <string>

using std::string;
using std::to_string;
using std::set;
using std::cout;
using std::endl;


void setMaximum(float max)
{
  TList * list = gPad->GetListOfPrimitives();
  TIter next(list);
  while (TObject * obj = next())
  {
    std::string className = obj->ClassName();
    if (className.find("TH1") != std::string::npos)
    {
      TH1 * th1 = static_cast<TH1*>(obj);
      th1->SetMaximum(max);
    }
    if (className.find("THStack") != std::string::npos)
    {
      THStack * thstack = static_cast<THStack*>(obj);
      thstack->SetMaximum(max);
    }
  }
  gPad->Modified();
  gPad->Update();

  //TH1 * obj = (TH1*)(gPad->GetListOfPrimitives()->First());
  //obj->SetMaximum(max);
  //gPad->Update();
}

double getMaximumTH1()
{
  TList * list = gPad->GetListOfPrimitives();
  TIter next(list);
  int index = 0;
  double max = 0;
  while (TObject * obj = next())
  {
    std::string className = obj->ClassName();
    if (className.find("TH1") != std::string::npos)
    {
      TH1 * th1 = static_cast<TH1*>(obj);
      double t_max = th1->GetMaximum();
      if (t_max>max || index==0) max = t_max;
    }
    index++;
  }
  return max;
}

double setMaximumTH1(double maxFraction = 1.05)
{
  double max = getMaximumTH1() * maxFraction;
  setMaximum(max);
  return 1;
}

TCanvas * newCanvas(string const & name = "", int size = 500)
{
  TSeqCollection * canvases = gROOT->GetListOfCanvases();
  double iCanvas = canvases->GetEntries();
  string canvasName;
  if (name == "") canvasName = "c_g_" + to_string(iCanvas++);
  else canvasName = name;
  return new TCanvas(canvasName.c_str(), canvasName.c_str(), size, size);
}

void draw_1D_compare_eff(TChain * tree, int dist_nbins, float dist_min, float dist_max, float lumi, string x_axis, string analysis1_cut, string analysis2_cut, string weight, string folder_name, string prefix)
{
  TH1D unskimmed((prefix+"_unskimmed").c_str(), (prefix+"_unskimmed").c_str(), dist_nbins, dist_min, dist_max);
  tree->Draw((x_axis+">>"+prefix+"_unskimmed").c_str(),weight.c_str(),"goff");
  TH1D analysis1((prefix+"_analysis1").c_str(), (prefix+"_analysis1").c_str(), dist_nbins, dist_min, dist_max);
  tree->Draw((x_axis+">>"+prefix+"_analysis1").c_str(), ("("+analysis1_cut+")*("+weight+")").c_str(), "goff");
  TH1D analysis2((prefix+"_analysis2").c_str(), (prefix+"_analysis2").c_str(), dist_nbins, dist_min, dist_max);
  tree->Draw((x_axis+">>"+prefix+"_analysis2").c_str(), ("("+analysis2_cut+")*("+weight+")").c_str(), "goff");

  TH1D analysis1_eff((prefix+"_analysis1_eff").c_str(),(prefix+"_analysis1_eff").c_str(), dist_nbins, dist_min, dist_max);
  analysis1_eff.Divide(&analysis1,&unskimmed);
  TH1D analysis2_eff((prefix+"_analysis2_eff").c_str(),(prefix+"_analysis2_eff").c_str(), dist_nbins, dist_min, dist_max);
  analysis2_eff.Divide(&analysis2,&unskimmed);

  cout<<"number unskimmed: "<<unskimmed.Integral()*lumi<<endl;
  cout<<"number analysis1:  "<<analysis1.Integral()*lumi<<endl;
  cout<<"number analysis2:   "<<analysis2.Integral()*lumi<<endl;

  gStyle->SetOptStat(0);

  TCanvas * t_c  = 0;
  t_c= newCanvas();
  unskimmed.Draw("hist");
  t_c->SaveAs((folder_name+prefix+"_comp_unskimmed.pdf").c_str());
  t_c = newCanvas();
  analysis1.Draw("hist");
  t_c->SaveAs((folder_name+prefix+"_comp_analysis1.pdf").c_str());
  t_c = newCanvas();
  analysis2.Draw("hist");
  t_c->SaveAs((folder_name+prefix+"_comp_analysis2.pdf").c_str());

  t_c = newCanvas();
  analysis1_eff.SetTitle("Efficiency of analysis");
  analysis1_eff.Draw();
  analysis2_eff.Draw("same");
  analysis1_eff.SetLineColor(kRed);

  setMaximumTH1();

  t_c->SaveAs((folder_name+prefix+"_comp_eff.pdf").c_str());
}

void draw_1D_eff(TChain * tree, int dist_nbins, float dist_min, float dist_max, float lumi, string x_axis, string analysis_cut, string weight, string folder_name, string prefix)
{
  TH1D unskimmed((prefix+"_unskimmed").c_str(), (prefix+"_unskimmed").c_str(), dist_nbins, dist_min, dist_max);
  tree->Draw((x_axis+">>"+prefix+"_unskimmed").c_str(),weight.c_str(),"goff");
  TH1D analysis((prefix+"_analysis").c_str(), (prefix+"_analysis").c_str(), dist_nbins, dist_min, dist_max);
  tree->Draw((x_axis+">>"+prefix+"_analysis").c_str(), ("("+analysis_cut+")*("+weight+")").c_str(), "goff");

  TH1D analysis_eff((prefix+"_analysis_eff").c_str(),(prefix+"_analysis_eff").c_str(), dist_nbins, dist_min, dist_max);
  analysis_eff.Divide(&analysis,&unskimmed);

  cout<<"number unskimmed: "<<unskimmed.Integral()*lumi<<endl;
  cout<<"number analysis:  "<<analysis.Integral()*lumi<<endl;

  gStyle->SetOptStat(0);

  TCanvas * t_c  = 0;
  t_c= newCanvas();
  unskimmed.Draw("hist");
  t_c->SaveAs((folder_name+prefix+"_unskimmed.pdf").c_str());
  t_c = newCanvas();
  analysis.Draw("hist");
  t_c->SaveAs((folder_name+prefix+"_analysis.pdf").c_str());

  t_c = newCanvas();
  analysis_eff.SetTitle("Efficiency of analysis");
  analysis_eff.Draw();

  setMaximumTH1();

  t_c->SaveAs((folder_name+prefix+"_eff.pdf").c_str());
}


void draw_2D_eff(TChain * tree, int dist_nbins, float dist_min, float dist_max, float lumi, string x_axis, string y_axis, string analysis_cut, string weight, string folder_name, string prefix)
{
  TH2D unskimmed((prefix+"_unskimmed").c_str(), (prefix+"_unskimmed").c_str(), dist_nbins, dist_min, dist_max, dist_nbins, dist_min, dist_max);
  tree->Draw((y_axis+":"+x_axis+">>"+prefix+"_unskimmed").c_str(),weight.c_str(),"goff");
  TH2D analysis((prefix+"_analysis").c_str(), (prefix+"_analysis").c_str(), dist_nbins, dist_min, dist_max, dist_nbins, dist_min, dist_max);
  tree->Draw((y_axis+":"+x_axis+">>"+prefix+"_analysis").c_str(), ("("+analysis_cut+")*("+weight+")").c_str(), "goff");

  TH2D analysis_eff((prefix+"_analysis_eff").c_str(),(prefix+"_analysis_eff").c_str(), dist_nbins, dist_min, dist_max, dist_nbins, dist_min, dist_max);
  analysis_eff.Divide(&analysis,&unskimmed);

  cout<<"number unskimmed: "<<unskimmed.Integral()*lumi<<endl;
  cout<<"number analysis:  "<<analysis.Integral()*lumi<<endl;

  gStyle->SetOptStat(0);

  TCanvas * t_c  = 0;
  t_c= newCanvas();
  unskimmed.Draw("colz");
  t_c->SaveAs((folder_name+prefix+"_unskimmed.pdf").c_str());
  t_c = newCanvas();
  analysis.Draw("colz");
  t_c->SaveAs((folder_name+prefix+"_analysis.pdf").c_str());

  t_c = newCanvas();
  analysis_eff.SetTitle("Efficiency of analysis");
  analysis_eff.Draw("colz");
  t_c->SaveAs((folder_name+prefix+"_eff.pdf").c_str());
}

int main(){
  TChain * tree = new TChain("tree");

  //tree->Add("/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/merged_higmc_higloose/*mChi-700_mLSP-1_*.root");
  //tree->Add("/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/unskimmed/*mChi-700_mLSP-1_*.root");
  //tree->Add("/net/cms29/cms29r0/pico/NanoAODv5/higgsino_angeles/2016/TChiHH/unskimmed/*mChi-*_mLSP-1_*.root");
  tree->Add("/net/cms29/cms29r0/pico/NanoAODv5/higgsino_humboldt/2016/SMS-TChiHH_2D/unskimmed/*mChi-*_mLSP-0_*.root");

  int dist_nbins = 60;
  float dist_min = 0;
  float dist_max = 1500;
  float lumi = 35.9;
  string folder_name = "plots/";

  string weight = "1";

  string resolved_baseline = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150";
  string resolved_bcuts = "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)||(nbt==2&&nbm==2))";
  string resolved_baseline_fourjet = "njet>=4&&njet<=5&&Alt$(hig_cand_drmax[0]<2.2,0)&&Alt$(hig_cand_am[0]<200,0)&&Alt$(hig_cand_dm[0]<40,0)"; 
  string resolved = resolved_baseline+"&&"+resolved_bcuts+"&&"+resolved_baseline_fourjet;

  string boosted_baseline = "met>150 && ht>300 && !low_dphi_met && nvlep==0 && ntk==0";
  string boosted_baseline_ak8 = "Alt$(fjet_pt[0]>300,0) && Alt$(fjet_pt[1]>300,0) && Alt$(fjet_msoftdrop[0]>85,0) && Alt$(fjet_msoftdrop[0]<135,0) && Alt$(fjet_msoftdrop[1]>85,0) && Alt$(fjet_msoftdrop[1]<135,0)";
  //string boosted_bcuts = "(fjet_mva_hbb_btv[0]>0.3||fjet_mva_hbb_btv[1]>0.3)";
  string boosted_doubleh = "((met<300&&Alt$(fjet_mva_hbb_btv[0]>0.6,0)&&Alt$(fjet_mva_hbb_btv[1]>0.6,0))||(met>300&&Alt$(fjet_mva_hbb_btv[0]>0.3,0)&&Alt$(fjet_mva_hbb_btv[1]>0.3,0)))";
  string boosted_singleh = "((Alt$(fjet_mva_hbb_btv[0]>0.3,0)&&Alt$(fjet_mva_hbb_btv[1]<=0.3,0))||(Alt$(fjet_mva_hbb_btv[0]<=0.3,0)&&Alt$(fjet_mva_hbb_btv[1]>0.3,0)))";
  string boosted_bcuts = "("+boosted_doubleh + "||" + boosted_singleh+")";
  string boosted = boosted_baseline+"&&"+boosted_bcuts+"&&"+boosted_baseline_ak8;

  //string x_axis = "MaxIf$(mc_pt,mc_id==25)";
  //string x_axis = "MinIf$(mc_pt,mc_id==25)";

  // Depending on mass of higgsino
  draw_1D_compare_eff(tree, dist_nbins, dist_min, dist_max, lumi, "mprod", resolved, boosted, weight, folder_name, "mprod");
  draw_1D_compare_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MaxIf$(mc_pt,mc_id==25)", resolved, boosted, weight, folder_name, "higgs_leading");
  draw_1D_compare_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MinIf$(mc_pt,mc_id==25)", resolved, boosted, weight, folder_name, "higgs_subleading");

  // Depending on leading higgs pt
  draw_1D_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MaxIf$(mc_pt,mc_id==25)", resolved, weight, folder_name, "resolved_higgs_leading");
  draw_1D_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MaxIf$(mc_pt,mc_id==25)", boosted, weight, folder_name, "boosted_higgs_leading");
  draw_1D_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MaxIf$(mc_pt,mc_id==25)", "("+boosted+")||("+resolved+")", weight, folder_name, "combined_higgs_leading");

  draw_2D_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MinIf$(mc_pt,mc_id==25)", "MaxIf$(mc_pt,mc_id==25)", resolved, weight, folder_name, "resolved_2d");
  draw_2D_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MinIf$(mc_pt,mc_id==25)", "MaxIf$(mc_pt,mc_id==25)", boosted, weight, folder_name, "boosted_2d");
  draw_2D_eff(tree, dist_nbins, dist_min, dist_max, lumi, "MinIf$(mc_pt,mc_id==25)", "MaxIf$(mc_pt,mc_id==25)", "("+boosted+")||("+resolved+")", weight, folder_name, "combined_2d");

}
