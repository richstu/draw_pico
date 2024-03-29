#include "higgsino/limit_scan.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"

#include "core/utilities.hpp"
#include "core/styles.hpp"
#include "core/cross_sections.hpp"

using namespace std;

namespace{
  int num_smooth_ = 0; // Number of times to smooth TH2D
  string filename_ = "txt/t1tttt_limit_scan.txt";
  string tag = "2D";
  //string model_ = "T1tttt";
  string model_ = "CN";
  bool unblind = false;
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);
  // styles style("2Dscan"); style.setDefaultStyle();

  if (Contains(filename_,"5tttt")) model_ = "T5tttt";

  if(filename_ == "") ERROR("No input file provided");

  vector<double> vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown, vsigobs, vsigexp;
  ReadPoints(vmx, vmy, vxsec, vobs, vobsup, vobsdown, vexp, vup, vdown, vsigobs, vsigexp);

  vector<double> vlim(vxsec.size());
  for(size_t i = 0; i < vxsec.size(); ++i){
    vlim.at(i) = vxsec.at(i) * vobs.at(i);
  }

  TH2D hsigobs = MakeObservedSignificancePlot(vmx, vmy, vsigobs);
  TH2D hsigexp = MakeExpectedSignificancePlot(vmx, vmy, vsigexp);

  MakeLimitPlot(vmx, vmy, vlim,
                vobs, vobsup, vobsdown,
                vexp, vup, vdown,
                hsigobs, hsigexp);
}

void ReadPoints(vector<double> &vmx,
                vector<double> &vmy,
                vector<double> &vxsec,
                vector<double> &vobs,
                vector<double> &vobsup,
                vector<double> &vobsdown,
                vector<double> &vexp,
                vector<double> &vup,
                vector<double> &vdown,
                vector<double> &vsigobs,
                vector<double> &vsigexp){
  ifstream infile(filename_);
  string line;

  while(getline(infile, line)){
    istringstream iss(line);
    //double pmx, pmy, pxsec, pxsecunc, pobs, pobsup, pobsdown, pexp, pup, pdown, pup2, pdown2, sigobs, sigexp;
    //iss >> pmx >> pmy >> pxsec >> pxsecunc >> pobs >> pobsup >> pobsdown >> pexp >> pup >> pdown >> pup2 >> pdown2 >> sigobs >> sigexp;
    double pmx, pmy, pxsec, pxsecunc, pobs, pexp, pup, pdown, pup2, pdown2, sigobs, sigexp;
    iss >> pmx >> pmy >> pxsec >> pxsecunc >> pobs >> pexp >> pup >> pdown >> pup2 >> pdown2 >> sigobs >> sigexp;
    int mglu(pmx);
    double xsec, exsec;
    //xsec::higgsinoCrossSection(mglu, xsec, exsec);
    xsec::higgsino2DCrossSection(mglu, xsec, exsec);
    if(Contains(model_, "T5HH")) xsec::gluinoCrossSection(mglu, xsec, exsec);
    // int factor(50), mlsp(pmy);
    // if((mglu%factor!=0 || mlsp%factor!=0) && mglu-mlsp!=225 && mlsp!=1450) continue;
    // if(mglu-mlsp==225 && mglu%factor!=0) continue;
    vmx.push_back(pmx);
    vmy.push_back(pmy);
    if(Contains(model_, "CN") || Contains(model_, "N1N2")) vxsec.push_back(xsec);
    else if(Contains(model_, "T5HH")) vxsec.push_back(xsec);
    else vxsec.push_back(pxsec);
    vobs.push_back(pobs);
    vobsup.push_back(pobs/(1+pxsecunc));
    vobsdown.push_back(pobs/(1-pxsecunc));
    vexp.push_back(pexp);
    vup.push_back(pup);
    vdown.push_back(pdown);
    vsigobs.push_back(sigobs);
    vsigexp.push_back(sigexp);
  }
  infile.close();

  if(vmx.size() <= 2) ERROR("Need at least 3 model_s to draw scan");
  if(vmx.size() != vmy.size()
     || vmx.size() != vxsec.size()
     || vmx.size() != vobs.size()
     || vmx.size() != vobsup.size()
     || vmx.size() != vobsdown.size()
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size()
     || vmx.size() != vsigobs.size()
     || vmx.size() != vsigexp.size()) ERROR("Error parsing text file. Model_ point not fully specified");
}

TH2D MakeObservedSignificancePlot(vector<double> vmx,
                                  vector<double> vmy,
                                  vector<double> vobs){
  SetupSignedColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV];Observed Significance";

  TGraph2D g("", title.c_str(), vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));

  double the_max = 3.;
  // for(int i = 0; i < g.GetN(); ++i){
  //   double z = fabs(g.GetZ()[i]);
  //   if(z>the_max) the_max = z;
  // }
  g.SetMinimum(-the_max);
  g.SetMaximum(the_max);

  g.SetNpx((2600.-800.)/12.5);
  g.SetNpy(1600/12.5);

  g.GetHistogram()->SetTitle(title.c_str());
  g.GetHistogram()->SetTickLength(0., "Z");

  TCanvas c("can","can", 800, 900);
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
              "#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
               "#scale[0.8]{137 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);
  TLatex model = GetModelLabel(c.GetLeftMargin()+0.03, 1.-c.GetTopMargin()-0.03);

  g.Draw("colz");

  vector<TLine> lines;
  for(double z = 0.; z < the_max; z+=0.5){
    int style = min(5, static_cast<int>(2.*fabs(z))+1);
    double width = max(1., 3.-fabs(z));
    DrawContours(g, 1, style, width, 0, z);
    double x1 = 1.-c.GetRightMargin()+0.0047;
    double x2 = 1.-c.GetRightMargin()+0.05;
    double ybot = c.GetBottomMargin();
    double ytop = 1.-c.GetTopMargin();
    double zpos = ybot+(ytop-ybot)*(the_max+z)/(2.*the_max);
    lines.emplace_back(x1, zpos, x2, zpos);
    lines.back().SetLineColor(1);
    lines.back().SetLineStyle(style);
    lines.back().SetLineWidth(width);
    lines.back().SetNDC(true);
    if(z != 0.){
      DrawContours(g, 1, style, width, 0, -z);
      double zneg = ybot+(ytop-ybot)*(the_max-z)/(2.*the_max);
      lines.emplace_back(x1, zneg, x2, zneg);
      lines.back().SetLineColor(1);
      lines.back().SetLineStyle(style);
      lines.back().SetLineWidth(width);
      lines.back().SetNDC(true);
    }
  }
  for(auto &l: lines){
    l.Draw("same");
  }
  
  model.Draw("same");
  ltitle.Draw("same");
  rtitle.Draw("same");

  c.Print(("plots/"+model_+"_sigobs_"+tag+".pdf").c_str());
  cout<<"open "<<"plots/"+model_+"_sigobs_"+tag+".pdf"<<endl;
  c.Print(("plots/"+model_+"_sigobs_"+tag+".root").c_str());

  TH2D h = *g.GetHistogram();
  h.SetTitle("Observed Significance");
  return h;
}

TH2D MakeExpectedSignificancePlot(vector<double> vmx,
                                  vector<double> vmy,
                                  vector<double> vobs){
  SetupColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV]; Expected Significance";

  TGraph2D g("", title.c_str(), vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));

  double the_max = 0.;
  for(int i = 0; i < g.GetN(); ++i){
    double z = g.GetZ()[i];
    if(z>the_max) the_max = z;
  }
  if(the_max > 6.) the_max = 6.;
  g.SetMinimum(0.);
  g.SetMaximum(the_max);

  g.SetNpx((2600.-800.)/12.5);
  g.SetNpy(1600/12.5);

  g.GetHistogram()->SetTitle(title.c_str());

  TCanvas c;
  c.cd();

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
              "#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
               "#scale[0.8]{137 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);
  TLatex model = GetModelLabel(c.GetLeftMargin()+0.03, 1.-c.GetTopMargin()-0.03);

  g.Draw("colz");

  for(double z = 0.; z < the_max; z+=1.){
    int style = (z == 5. ? 1 : 2);
    double width = (z == 5. ? 2. : 1.);
    DrawContours(g, 1, style, width, 0, z);
  }
  
  model.Draw("same");
  ltitle.Draw("same");
  rtitle.Draw("same");
  
  c.Print(("plots/"+model_+"_sigexp_"+tag+".pdf").c_str());
  cout<<"open "<<"plots/"+model_+"_sigexp_"+tag+".pdf"<<endl;
  c.Print(("plots/"+model_+"_sigexp_"+tag+".root").c_str());

  TH2D h = *g.GetHistogram();
  h.SetTitle("Expected Significance");
  return h;
}

void MakeLimitPlot(vector<double> vmx,
                   vector<double> vmy,
                   vector<double> vlim,
                   vector<double> vobs,
                   vector<double> vobsup,
                   vector<double> vobsdown,
                   vector<double> vexp,
                   vector<double> vup,
                   vector<double> vdown,
                   const TH2D &hsigobs,
                   const TH2D &hsigexp){
  SetupColors();

  string xparticle, yparticle;
  GetParticleNames(xparticle, yparticle);
  string title = ";m_{"+xparticle+"} [GeV];m_{"+yparticle+"} [GeV];95% CL upper limit on cross section [pb]";
  
  TGraph2D glim("", title.c_str(), vlim.size(), &vmx.at(0), &vmy.at(0), &vlim.at(0));
  TGraph2D gobs("", "Observed Limit", vobs.size(), &vmx.at(0), &vmy.at(0), &vobs.at(0));
  TGraph2D gobsup("", "Observed +1#sigma Limit", vobsup.size(), &vmx.at(0), &vmy.at(0), &vobsup.at(0));
  TGraph2D gobsdown("", "Observed -1#sigma Limit", vobsdown.size(), &vmx.at(0), &vmy.at(0), &vobsdown.at(0));
  TGraph2D gexp("", "Expected Limit", vexp.size(), &vmx.at(0), &vmy.at(0), &vexp.at(0));
  TGraph2D gup("", "Expected +1#sigma Limit", vup.size(), &vmx.at(0), &vmy.at(0), &vup.at(0));
  TGraph2D gdown("", "Expected -1#sigma Limit", vdown.size(), &vmx.at(0), &vmy.at(0), &vdown.at(0));

  glim.SetMinimum(0.00001);
  if (model_ == "T5HH") glim.SetMinimum(0.001);
  glim.SetMaximum(2);

  glim.SetNpx((2600.-800.)/12.5);
  glim.SetNpy(1600/12.5);
  //glim.SetNpx(800/12.5);
  //glim.SetNpy(700/12.5);
  //glim.SetNpx(800/25);
  //glim.SetNpy(700/25);

  glim.SetTitle(title.c_str());

  TLegend l(gStyle->GetPadLeftMargin(), 1.-2.*gStyle->GetPadTopMargin(),
            1.-gStyle->GetPadRightMargin(), 1.-gStyle->GetPadTopMargin());
  l.SetNColumns(2);
  l.SetTextSize(0.05);
  l.SetBorderSize(0);
  
  TCanvas c("","",1200, 950);
  c.cd();
  c.SetRightMargin(0.15);

  TLatex ltitle(c.GetLeftMargin(), 1.-0.5*c.GetTopMargin(),
                "#font[62]{CMS}#scale[0.76]{#font[52]{ Supplementary}}");
  TLatex rtitle(1.-c.GetRightMargin(), 1.-0.5*c.GetTopMargin(),
                "#scale[0.8]{137 fb^{-1} (13 TeV)}");
  ltitle.SetNDC();
  rtitle.SetNDC();
  ltitle.SetTextAlign(12);
  rtitle.SetTextAlign(32);

  c.SetTopMargin(2.*c.GetTopMargin());
  c.SetLeftMargin(0.12);
  c.SetLogz();
  glim.Draw("colz");

  float limit_width = 2.0; //was 5.0
  TGraph cup = DrawContours(gup, 2, 2, limit_width, num_smooth_);
  TGraph cdown = DrawContours(gdown, 2, 2, limit_width, num_smooth_);
  TGraph cexp = DrawContours(gexp, 2, 1, limit_width, num_smooth_, 1.);
  TGraph cobsup, cobsdown, cobs;
  if (unblind) {
    std::cout << "DEBUG: unblinding" << std::endl;
    cobsup = DrawContours(gobsup, 1, 2, limit_width, num_smooth_);
    cobsdown = DrawContours(gobsdown, 1, 2, limit_width, num_smooth_);
    cobs = DrawContours(gobs, 1, 1, limit_width, num_smooth_, 1.);
  }

  l.AddEntry(&cexp, "Expected", "l");
  if (unblind) l.AddEntry(&cobs, "Observed", "l");

  l.Draw("same");

  TLatex model = GetModelLabel(c.GetLeftMargin()+0.03, 1.-c.GetTopMargin()-0.03);
  model.Draw("same");

  ltitle.Draw("same");
  rtitle.Draw("same");
  
  string filebase = model_+"_limit_scan";
  if(num_smooth_>0){
    filebase += "_smooth";
    filebase += to_string(num_smooth_);
  }

  gPad->Update();
  glim.GetZaxis()->SetTitleOffset(1.3);
  glim.GetXaxis()->SetRangeUser(0,700);
  glim.GetYaxis()->SetRangeUser(0,800);
  if (model_=="T5HH") {
    glim.GetXaxis()->SetRangeUser(1000,2600);
    glim.GetYaxis()->SetRangeUser(0,1600);
  }

  c.Print(("plots/"+filebase+"_"+tag+".pdf").c_str());
  cout<<"open "<<"plots/"+filebase+"_"+tag+".pdf"<<endl;
  
  TFile file(("plots/"+filebase+"_"+tag+".root").c_str(), "recreate");
  glim.GetHistogram()->Write((model_+"ObservedExcludedXsec").c_str());
  if (unblind) {
    cobs.Write((model_+"ObservedLimit").c_str());
    cobsup.Write((model_+"ObservedLimitUp").c_str());
    cobsdown.Write((model_+"ObservedLimitDown").c_str());
  }
  cexp.Write((model_+"ExpectedLimit").c_str());
  cup.Write((model_+"ExpectedLimitUp").c_str());
  cdown.Write((model_+"ExpectedLimitDown").c_str());
  (void) hsigobs;
  (void) hsigexp;
  //if(false){ // Significances not saved together with limits as per recommendations
  //  hsigobs.Write("ObservedSignificance");
  //  hsigexp.Write("ExpectedSignificance");
  //}
  file.Close();
  cout << "\nSaved limit curves in " << filebase << "_"<<tag<<".root\n" << endl;
}

int GetNumBins(const vector<double> &pts, double width){
  double pmin = *min_element(pts.cbegin(), pts.cend());
  double pmax = *max_element(pts.cbegin(), pts.cend());
  return max(1, min(500, static_cast<int>(ceil((pmax-pmin)/width))));
}

void GetParticleNames(string &xparticle, string &yparticle){
  if(model_=="N1N2"){
    xparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{2}}}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  } else if(model_=="T1tttt"){
  //if(model_=="T1tttt"){
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }else if(model_=="T5tttt"){
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }else if(model_=="T2tt"){
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }else if(model_=="T6ttWW"){
    xparticle = "sbottom";
    yparticle = "chargino";
  } else if(model_=="T5HH"){
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }else{
    DBG(("Unknown model: "+model_));
    xparticle = "#tilde{g}";
    yparticle = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  }
}

TLatex GetModelLabel(double x, double y){
  string lsp = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  string label = "";
  if(model_=="T1tttt"){
    label = "pp #rightarrow #tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #rightarrow t#kern[0.4]{#bar{t}}#kern[0.4]{"+lsp+"}";
  }else if(model_=="T5tttt"){
    label = "#splitline{pp #rightarrow #tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #rightarrow #tilde{t}_{1}t, #tilde{t}_{1} #rightarrow #bar{t}#kern[0.4]{"
      +lsp+"}}{(m#kern[0.3]{_{#lower[-0.12]{#tilde{t}_{1}}}} - m#kern[0.12]{_{"+lsp+"}} = 175 GeV)}";
  }else if(model_=="N1N2") {
    TString chii= "#lower[-0.12]{#tilde{#chi}}#kern[+0.1]{#lower[0.2]{#scale[0.99]{^{0,#kern[+0.2]{#lower[-0.15]{#pm}}}}}}#kern[-4]{#scale[0.99]{_{i}}}";
    TString chij= "#lower[-0.12]{#tilde{#chi}}#kern[+0.1]{#lower[0.2]{#scale[0.99]{^{0,#kern[+0.2]{#lower[0.15]{#mp}}}}}}#kern[-4]{#scale[0.99]{_{j}}}";
    TString chi10= "#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{1}}}";
    TString chi20= "#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{2}}}";
    TString chi30= "#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{3}}}";
    TString xsoft= "X#lower[-0.2]{#scale[0.85]{_{soft}}}";
    //label = "pp #rightarrow "+chii+"#kern[0.7]{"+chij+"}  #rightarrow "+chi10+"#kern[0.3]{"+chi10+"} + "
    label = "pp #rightarrow "+chi30+"#kern[0.3]{"+chi20+"} #rightarrow HH#kern[0.3]{"+chi10+"}#kern[0.3]{"+chi10+"}";
  }else if(model_=="T5HH") {
    TString chi10= "#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{1}}}";
    TString chi20= "#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{2}}}";
    TString xsoft= "X#lower[-0.2]{#scale[0.85]{_{soft}}}";
    //label = "pp #rightarrow "+chii+"#kern[0.7]{"+chij+"}  #rightarrow "+chi10+"#kern[0.3]{"+chi10+"} + "
    label = "pp #rightarrow #tilde{g}#kern[0.3]{#tilde{g}}, #tilde{g} #rightarrow"+chi20+"#kern[0.3]{q}#kern[0.3]{#bar{q}}, "+chi20+" #rightarrow H#kern[0.3]{"+chi10+"}";
  }else{
    DBG(("Unknown model: "+model_));
    label = "";
  }
  TLatex l(x, y, label.c_str());
  l.SetNDC();
  l.SetTextAlign(13);
  l.SetTextSize(0.032);
  return l;
}

void Style(TGraph *g, int color, int style, float width){
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}

TGraph DrawContours(TGraph2D &g2, int color, int style, double width,
                    int n_smooth, double val){
  TGraph graph;

  TList *l;
  //// Finding the TH2D, smoothing it, and creating a TGraph2D to get a new Delauny interpolation
  if(n_smooth>0){
    TH2D *histo2d = g2.GetHistogram();
    TH2D htemp("", "",
               100, histo2d->GetXaxis()->GetXmin(), histo2d->GetXaxis()->GetXmax(),
               100, histo2d->GetYaxis()->GetXmin(), histo2d->GetYaxis()->GetXmax());
    for(int binx=1; binx<=htemp.GetNbinsX(); ++binx){
      double x = htemp.GetXaxis()->GetBinCenter(binx);
      for(int biny=1; biny<=htemp.GetNbinsY(); ++biny){
        double y = htemp.GetYaxis()->GetBinCenter(biny);
        double z = g2.Interpolate(x,y);
        if(z!=0.){
          htemp.SetBinContent(htemp.GetBin(binx, biny), z);
        }
      }
    }
    
    for(int ind=0; ind<n_smooth; ++ind){
      htemp.Smooth(1,"k5b");
    }
    
    vector<double> vx, vy, vz;
    double glu_lsp = 225;
    for(int binx=1; binx<=htemp.GetNbinsX(); ++binx){
      double x = htemp.GetXaxis()->GetBinCenter(binx);
      for(int biny=1; biny<=htemp.GetNbinsY(); ++biny){
        double y = htemp.GetYaxis()->GetBinCenter(biny);
        double z = htemp.GetBinContent(htemp.GetBin(binx,biny));
        
        vx.push_back(x);
        vy.push_back(y);
        int thresh = glu_lsp+30;
        if (model_=="T5tttt") thresh = glu_lsp+50;
        if(x-y>thresh){
          vz.push_back(z);
        }else{
          vz.push_back(g2.Interpolate(x,y));
        }
      }
    }
    
    TGraph2D gsmooth("gsmooth", "Cross-Section Limit", vx.size(), &vx.at(0), &vy.at(0), &vz.at(0));
    gsmooth.GetHistogram();
    l = gsmooth.GetContourList(val);
  } else {
    g2.GetHistogram();
    l = g2.GetContourList(val);
  }
  if(l == nullptr) return graph;
  int max_points = -1;
  vector<TGraph*> g;
  for(int i = 0; i < l->GetSize(); ++i){
    g.push_back(static_cast<TGraph*>(l->At(i)));
    Style(g[i], color, style, width);
    if(g[i] == nullptr) continue;
    int n_points = g[i]->GetN();
    if(n_points > max_points){
      if(n_smooth>0) FixGraph(*(g[i]));
      TGraph* old_graph = static_cast<TGraph*>(graph.Clone());
      if (i>0) ReverseGraph(*(g[i]));
      graph = *(joinGraphs(old_graph, g[i]));
      graph.RemovePoint(graph.GetN()-1);
      max_points = n_points;
    }
    g[i]->Draw("L same");
  }

  graph.SetTitle(g2.GetTitle());
  graph.SetLineColor(color);
  graph.SetLineWidth(width);
  graph.SetLineStyle(style);
  return graph;
}

TGraph* joinGraphs(TGraph *graph1, TGraph *graph2){
  TGraph *graph = new TGraph;
  double mglu, mlsp;
  for(int point(0); point < graph1->GetN(); point++) {
    graph1->GetPoint(point, mglu, mlsp);
    graph->SetPoint(graph->GetN(), mglu, mlsp);
  } // Points in graph1
  for(int point(0); point < graph2->GetN(); point++) {
    graph2->GetPoint(point, mglu, mlsp);
    graph->SetPoint(graph->GetN(), mglu, mlsp);
  } // Points in graph1
  graph1->GetPoint(0, mglu, mlsp);
  graph->SetPoint(graph->GetN(), mglu, mlsp);
  TString gname = graph1->GetName(); gname += graph2->GetName();
  graph->SetName(gname);

  return graph;
}

void FixGraph(TGraph &graph){
  double glu_lsp = 225;
  if(model_=="T5tttt") glu_lsp = 265.;
  int np(graph.GetN());
  double iniglu, endglu, inilsp, endlsp;

  graph.GetPoint(0, iniglu, inilsp);
  graph.GetPoint(np-1, endglu, endlsp);

  // Reversing graph if printed towards decreasing mgluino
  if(inilsp < endlsp) {
    ReverseGraph(graph);
    endglu = iniglu;
    endlsp = inilsp;
  }

  // Adding a point so that it goes down to mLSP = 0, but not for WZ,SOS
  if(endlsp<30){
    graph.SetPoint(graph.GetN(), endglu, 0);
    np++;
  }

  ReverseGraph(graph);
  // Adding a point at mLSP = 0, and removing points beyond the diagonal
  for(int point(0); point < np; point++){
    double mglu, mlsp;
    graph.GetPoint(point, mglu, mlsp);
    if(mlsp > mglu-glu_lsp-5){
      while(point <= graph.GetN() && point!=0) {
        graph.RemovePoint(graph.GetN()-1);
        np--;
      }
      break;
    }
  }
}

void ReverseGraph(TGraph &graph){
  int np(graph.GetN());
  double mglu, mlsp;
  vector<double> mglus, mlsps;
  for(int point(np-1); point >= 0; point--){
    graph.GetPoint(point, mglu, mlsp);
    mglus.push_back(mglu);
    mlsps.push_back(mlsp);
  }
  for(int point(0); point < np; point++)
    graph.SetPoint(point, mglus[point], mlsps[point]);
}

void SetupColors(){
  const unsigned num = 5;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[num] = {0.50, 0.50, 1.00, 1.00, 1.00};
  double green[num] = {0.50, 1.00, 1.00, 0.60, 0.50};
  double blue[num] = {1.00, 1.00, 0.50, 0.40, 0.50};
  /*const unsigned num = 6;
    double red[num] =   {1.,0.,0.,0.,1.,1.};
    double green[num] = {0.,0.,1.,1.,1.,0.};
    double blue[num] =  {1.,1.,1.,0.,0.,0.};
    double stops[num] = {0.,0.2,0.4,0.6,0.8,1.};*/
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    colors[i] = fi+i;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);
}

void SetupSignedColors(){
  const unsigned num = 3;
  const int bands = 255;
  int colors[bands];
  double stops[num] = {0.0, 0.5, 1.0};
  double red[num]   = {0.0, 1.0, 1.0};
  double green[num] = {0.0, 1.0, 0.0};
  double blue[num]  = {1.0, 1.0, 0.0};
  int fi = TColor::CreateGradientColorTable(num,stops,red,green,blue,bands);
  for(int i = 0; i < bands; ++i){
    colors[i] = fi+i;
  }
  gStyle->SetNumberContours(bands);
  gStyle->SetPalette(bands, colors);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"num_smooth", required_argument, 0, 's'},
      {"model", required_argument, 0, 'm'},
      {"tag", required_argument, 0, 't'},
      {"unblind", no_argument, 0, 0},
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:m:s:t:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      num_smooth_ = atoi(optarg);
      break;
    case 'm':
      model_ = optarg;
      break;
    case 't':
      tag = optarg;
      break;
    case 'f':
      filename_ = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
        printf("Bad option! Found option name %s\n", optname.c_str());
      } else if (optname == "unblind") {
        unblind = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
