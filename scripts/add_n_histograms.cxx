#include <unordered_map>

//content rfrom tdrstyle.C
#include "TStyle.h"
TStyle *tdrStyle;

// tdrGrid: Turns the grid lines on (true) or off (false)

void tdrGrid(bool gridOn) {
  tdrStyle->SetPadGridX(gridOn);
  tdrStyle->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.0125);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.75);//1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

}
// /tdrstyle.C

//content from CMS_lumi.C/CMS_lumi.h
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText_sim      = "Simulation";//"Preliminary";
TString extraText_prelim   = "Preliminary";
TString extraText_simprelim= "#splitline{Simulation}{Preliminary}";
TString extraText_simsupp  = "#splitline{Simulation}{Supplementary}";
TString extraText_supp     = "Supplementary";

float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.5;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.76;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.7;

TString lumi_13TeV = "138 fb^{-1}";
TString lumi_13p6TeV = "62 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_sqrtS = "";

bool drawLogo      = false;

#include <iostream>

void 
CMS_lumi( TPad* pad, int iPeriod, int iPosX, int extraTextSelect)
{            
  TString extraText = extraText_sim;
  if(extraTextSelect==1){         extraText = extraText_prelim;
  } else if(extraTextSelect==2){  extraText = extraText_simprelim;
  } else if(extraTextSelect==3){  extraText = extraText_supp;
  } else if(extraTextSelect==4){  extraText = extraText_simsupp;
  }

  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 ) 
    {
      lumiText = lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
    }
  else if ( iPeriod==7 )
    { 
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV; 
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
 else if ( iPeriod==8 )
    { 
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV; 
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_13p6TeV; 
      lumiText += " (13.6 TeV)";
      if( outOfFrame) lumiText += "}";
    }
 else if ( iPeriod==9 )
    { 
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV; 
      lumiText += " (13 TeV)";
      if( outOfFrame) lumiText += "}";
    }

  else if ( iPeriod==12 )
    {
      lumiText += "8 TeV";
    }
  else if ( iPeriod==0 )
    {
      lumiText += lumi_sqrtS;
    }
   
  std::cout << lumiText << endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }
  
  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
	{
	  posX_ =   l + 0.045*(1-l-r)*W/H;
	  posY_ = 1-t - 0.045*(1-t-b);
	  float xl_0 = posX_;
	  float yl_0 = posY_ - 0.15;
	  float xl_1 = posX_ + 0.15*H/W;
	  float yl_1 = posY_;
	  //TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
	  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
	  pad_logo->Draw();
	  pad_logo->cd();
	  //CMS_logo->Draw("X");
	  pad_logo->Modified();
	  pad->cd();
	}
      else
	{
	  latex.SetTextFont(cmsTextFont);
	  latex.SetTextSize(cmsTextSize*t);
	  latex.SetTextAlign(align_);
	  latex.DrawLatex(posX_, posY_, cmsText);
	  if( writeExtraText ) 
	    {
	      latex.SetTextFont(extraTextFont);
	      latex.SetTextAlign(align_);
	      latex.SetTextSize(extraTextSize*t);
	      latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
	    }
	}
    }
  else if( writeExtraText )
    {
      if( iPosX==0) 
	{
	  posX_ =   l +  relPosX*(1-l-r);
	  posY_ =   1-t+lumiTextOffset*t;
	}
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);      
    }
  return;
}
// /CMS_lumi

void overlay_histograms(std::vector<TString> filenames, std::vector<std::string> labels, TString x_axis_lab = "",TString outfile = "outfile", bool normalize = false, bool rerange=false, string proc_name = "H#rightarrowZ#gamma"){
  TFile   *new_file      = nullptr;
  std::vector<TH1D*> new_histogram = {};
  TCanvas *new_canvas    = nullptr;
  TLegend *legend = nullptr;
  bool first = true;
  double scale = 1;
  for(unsigned int idx_file = 0; idx_file < filenames.size(); idx_file++){
//  for(TString filename :filenames){
    new_histogram.push_back(nullptr);
    new_file = new TFile(filenames[idx_file]);
    new_canvas = static_cast<TCanvas*>(new_file -> Get("canvas"));


    for(int idx_h=0; idx_h<50000; idx_h++){
      if(new_histogram[idx_file] != nullptr){break;}
      new_histogram[idx_file] = static_cast<TH1D*>( new_canvas->FindObject( ("sig_" + proc_name + "_" + std::to_string(idx_h)).c_str() ) );
    }
    scale = 0.1;

    if(normalize){
      scale = 1.0/(new_histogram[idx_file] -> Integral());
    } 

    new_histogram[idx_file] -> Scale(scale); 
  }

  int Ncolors = 7;
  std::vector<Color_t> color_wheel = {(kBlue),(kRed),(kGreen+2),(kOrange+9),(kViolet-1),(kAzure+1),(kCyan-6)};

  TString pdf= outfile + ".pdf";
  TCanvas *canvas = new TCanvas("c1", "outplot", 1920, 1440);

  canvas -> Update();
  gStyle -> SetOptFit();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);
  gPad   -> SetTopMargin(0.1); 
  gPad   -> SetRightMargin(0.03);


  double max = 0;
  for(unsigned int idx_file = 0; idx_file < new_histogram.size(); idx_file++){
    double hist_max = new_histogram[idx_file] -> GetMaximum();
    max = hist_max > max ? hist_max : max;
  }

  for(unsigned int idx_file = 0; idx_file < new_histogram.size(); idx_file++){
    new_histogram[idx_file] -> SetLineColor(color_wheel[idx_file%Ncolors]);
    new_histogram[idx_file] -> SetLineWidth(2);
    
    //Set range around the mass of the higgs
    if(rerange){
      new_histogram[idx_file] -> GetXaxis() -> SetRange(0,26*4*5);
    }

    if(!first){
      new_histogram[idx_file] -> Draw("HIST SAME");
    } else {

      
      if(filenames[idx_file].Contains("lly_m")){
        new_histogram[idx_file] -> SetAxisRange(110,140);
      }
      
      new_histogram[idx_file] -> Draw("HIST");
      new_histogram[idx_file] -> GetXaxis() -> SetTitle(x_axis_lab);
      new_histogram[idx_file] -> GetYaxis() -> SetTitle("Events/0.5 GeV");
      new_histogram[idx_file] -> GetYaxis() -> SetTitleOffset(0.75);
      new_histogram[idx_file] -> GetXaxis() -> SetTitleOffset(0.65);
      new_histogram[idx_file] -> SetMaximum(1.1*max);
      new_histogram[idx_file] -> GetXaxis() -> SetLabelSize(0.05);
      new_histogram[idx_file] -> GetYaxis() -> SetLabelSize(0.05);
      new_histogram[idx_file] -> GetXaxis() -> SetTitleSize(0.07);
      new_histogram[idx_file] -> GetYaxis() -> SetTitleSize(0.07);  
      legend = new TLegend(0.63,0.67,0.93,0.87);
      legend -> SetTextSize(0.04); 
      first=false;
    }

    legend->AddEntry(new_histogram[idx_file],labels[idx_file].c_str(),"l");
  }
  CMS_lumi( canvas, 8, 10, 4);
  setTDRStyle();


  legend->Draw();
  canvas -> Print(pdf);
  canvas -> Print(pdf+"]");

  delete new_canvas;
  delete canvas;
}

std::vector<TString> filename_vector(TString prefix, std::vector<TString> &filenames){
  std::vector<TString> ret_vector = {};
  for(unsigned int idx = 0; idx < filenames.size(); idx++){
    ret_vector.push_back(prefix + filenames[idx]);
  }
  return ret_vector;
}


void overlay_bkg_histograms(std::vector<TString> filenames, std::vector<std::string> labels, TString x_axis_lab = "",TString outfile = "outfile", bool normalize = false, bool rerange=false){
  TFile   *new_file      = nullptr;
  std::vector<TH1D*> new_histogram = {};
  TCanvas *new_canvas    = nullptr;
  TLegend *legend = nullptr;
  bool first = true;
  double scale = 1;
  for(unsigned int idx_file = 0; idx_file < filenames.size(); idx_file++){
//  for(TString filename :filenames){
    new_histogram.push_back(nullptr);
    new_file = new TFile(filenames[idx_file]);
    new_canvas = static_cast<TCanvas*>(new_file -> Get("canvas"));


    for(int idx_h=0; idx_h<50000; idx_h++){
      if(new_histogram[idx_file] != nullptr){break;}
      new_histogram[idx_file] = static_cast<TH1D*>( new_canvas->FindObject( ("bkg_Z+FakePhoton_"+std::to_string(idx_h)).c_str() ) );
    }
    scale = 0.1;

    if(normalize){
      scale = 1.0/(new_histogram[idx_file] -> Integral());
    } 

    new_histogram[idx_file] -> Scale(scale); 
  }

  int Ncolors = 7;
  std::vector<Color_t> color_wheel = {(kBlue),(kRed),(kGreen+2),(kOrange+9),(kViolet-1),(kAzure+1),(kCyan-6)};

  TString pdf= outfile + ".pdf";
  TCanvas *canvas = new TCanvas("c1", "outplot", 1920, 1440);

  canvas -> Update();
  gStyle -> SetOptFit();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);
  gPad   -> SetTopMargin(0.1); 
  gPad   -> SetRightMargin(0.03);

  double max = 0;
  for(unsigned int idx_file = 0; idx_file < new_histogram.size(); idx_file++){
    double hist_max = new_histogram[idx_file] -> GetMaximum();
    max = hist_max > max ? hist_max : max;
  }

  for(unsigned int idx_file = 0; idx_file < new_histogram.size(); idx_file++){
    new_histogram[idx_file] -> SetLineColor(color_wheel[idx_file%Ncolors]);
    new_histogram[idx_file] -> SetLineWidth(2);
    new_histogram[idx_file]->SetFillColorAlpha(kBlue,0);//SetFillStyle(4000);

  
    //Set range around the mass of the higgs
    //if(rerange){
    //  new_histogram[idx_file] -> GetXaxis() -> SetRange(0,26);
    //}

    if(!first){
      new_histogram[idx_file] -> Draw("HIST SAME");
    } else {

      
      if(filenames[idx_file].Contains("lly_m")){
        new_histogram[idx_file] -> SetAxisRange(100,180);
      }
      
      new_histogram[idx_file] -> Draw("HIST");
      new_histogram[idx_file] -> GetXaxis() -> SetTitle(x_axis_lab);
      new_histogram[idx_file] -> GetYaxis() -> SetTitle("Events (%)");
      new_histogram[idx_file] -> GetYaxis() -> SetTitleOffset(0.75);
      new_histogram[idx_file] -> GetXaxis() -> SetTitleOffset(0.65);
      new_histogram[idx_file] -> SetMaximum(1.4*max);
      new_histogram[idx_file] -> GetXaxis() -> SetLabelSize(0.05);
      new_histogram[idx_file] -> GetYaxis() -> SetLabelSize(0.04);
      new_histogram[idx_file] -> GetXaxis() -> SetTitleSize(0.07);
      new_histogram[idx_file] -> GetYaxis() -> SetTitleSize(0.07);  
      legend = new TLegend(0.63,0.67,0.93,0.87);
      legend -> SetTextSize(0.04); 
      first=false;
    }

    legend->AddEntry(new_histogram[idx_file],labels[idx_file].c_str(),"l");
  }
  CMS_lumi( canvas, 8, 10,4);
  setTDRStyle();


  legend->Draw();
  canvas -> Print(pdf);
  canvas -> Print(pdf+"]");

  delete new_canvas;
  delete canvas;
}




void add_n_histograms(){
  //TString plot_folder = "./plots/an_redwoodsv1_kinematicrefit/";// "plots/an_kinematic_refit/";
  TString plot_folder = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/PushControlRegions/draw_pico/plots/an_redwoodsv1_kinematicrefit/";// "plots/an_kinematic_refit/";
  
  //--------------------------------------------------Constrained Fit Comparisons-----------------------------------------------------------//
  vector<string> histogram_labels = {"w/out Z-mass fit","w/ Z-mass fit"};

  string outfile_mee_overlay = "plots/constrainedzfit-mee-overlay";
  string outfile_mmumu_overlay = "plots/constrainedzfit-mmumu-overlay";
  vector<TString> mee_histograms = {plot_folder + "an_kinematicfit_ee_m_nomll__ll_m0__wgt__lumi_nonorm_lin.root",
                                    plot_folder + "an_kinematicfit_ee_refit_m_nomll__ll_refit_m__wgt__lumi_nonorm_lin.root"};
  vector<TString> mmumu_histograms = {plot_folder + "an_kinematicfit_mumu_m_nomll__ll_m0__wgt__lumi_nonorm_lin.root",
                                      plot_folder + "an_kinematicfit_mumu_refit_m_nomll__ll_refit_m__wgt__lumi_nonorm_lin.root"};


  string outfile_meegamma_overlay = "plots/constrainedzfit-meegamma-overlay";
  string outfile_mmumugamma_overlay = "plots/constrainedzfit-mmumugamma-overlay";
  vector<TString> meegamma_histograms = {plot_folder + "an_kinematicfit_eephoton_m_nomll__llphoton_m0__wgt__lumi_nonorm_lin.root",
                                        plot_folder + "an_kinematicfit_refit_eephoton_m_nmll__llphoton_refit_m__wgt__lumi_nonorm_lin.root"};
  vector<TString> mmumugamma_histograms = {plot_folder + "an_kinematicfit_mumuphoton_m_nomll__llphoton_m0__wgt__lumi_nonorm_lin.root",
                                        plot_folder + "an_kinematicfit_refit_mumuphoton_m_nmll__llphoton_refit_m__wgt__lumi_nonorm_lin.root"};



  setTDRStyle();
  overlay_histograms(mee_histograms,   histogram_labels, "m_{ee}",     outfile_mee_overlay);
  overlay_histograms(mmumu_histograms, histogram_labels, "m_{#mu#mu}", outfile_mmumu_overlay);

  overlay_histograms(meegamma_histograms,   histogram_labels, "m_{ee#gamma}",     outfile_meegamma_overlay, false, true);
  overlay_histograms(mmumugamma_histograms, histogram_labels, "m_{#mu#mu#gamma}", outfile_mmumugamma_overlay, false, true);
  //--------------------------------------------END OF Constrained Fit Comparisons-----------------------------------------------------------//


  /*
  //----------------------------------------OTHER COMPARISONS COMMENTED OUT FOR CLARITY---------------------------------------------------//
  vector<TString> flavor_comparison = {plot_folder + "an_kinematicfit_refit_eephoton_m_nmll__llphoton_refit_m__wgt__lumi_nonorm_lin.root",
                                       plot_folder + "an_kinematicfit_refit_mumuphoton_m_nmll__llphoton_refit_m__wgt__lumi_nonorm_lin.root"};
  vector<string> histogram_labels_2 = {"Z #rightarrow e^{+}e^{-}","Z #rightarrow #mu^{+}#mu^{-}"};
  string outfile_flav_comp = "plots/constrainedzfit-flav-comp";

  vector<TString> bkg_ee_comparison = {plot_folder + "an_kinematicfit_ee_llphoton_bkg__llphoton_m0__wgt__lumi_nonorm_lin.root",
                                       plot_folder + "an_kinematicfit_ee_refit_llphoton_m_bkg__llphoton_refit_m__wgt__lumi_nonorm_lin.root"};
  vector<TString> bkg_mumu_comparison = {plot_folder + "an_kinematicfit_mumu_llphoton_bkg__llphoton_m0__wgt__lumi_nonorm_lin.root",
                                       plot_folder + "an_kinematicfit_mumu_refit_llphoton_m_bkg__llphoton_refit_m__wgt__lumi_nonorm_lin.root"};
 

  string outfile_meegamma_bkg_overlay = "plots/constrainedzfit-meegamma-bkg-overlay";
  string outfile_mmumugamma_bkg_overlay = "plots/constrainedzfit-mmumugamma-bkg-overlay";
 
  //overlay_histograms(flavor_comparison, histogram_labels_2, "m_{ll#gamma}", outfile_flav_comp, true, true);

  overlay_bkg_histograms(bkg_ee_comparison, histogram_labels, "m_{ee#gamma}",  outfile_meegamma_bkg_overlay, true, true);
  overlay_bkg_histograms(bkg_mumu_comparison, histogram_labels, "m_{#mu#mu#gamma}",  outfile_mmumugamma_bkg_overlay, true, true);


  //-------------------------------------------------------------------------------------------------------------------------------------//
  vector<TString> mll_histograms = {plot_folder + "an_kinematicfit_ll_m_nomll__ll_m0__wgt__lumi_nonorm_lin.root",
                                    plot_folder + "an_kinematicfit_ll_refit_m_nomll__ll_refit_m__wgt__lumi_nonorm_lin.root"};

  string outfile_mll_overlay = "plots/constrainedzfit-mll-overlay";
  overlay_histograms(mll_histograms, histogram_labels, "m_{ll}",  outfile_mll_overlay);//, true, true); 

 
  //Make plots testing pt uncertainty changes for Z -> ee
  plot_folder = "./plots/refit_test_diff_unc/";
  string plot_ee_part1 = "kf_test_unc_mlly_test_";
  string plot_ee_part2 = "__mll_refit_test__wgt__lumi_nonorm_lin.root";
  string plot_eeg_part1 = "kf_test_unc_mlly_test_";
  string plot_eeg_part2 = "__mlly_refit_test__wgt__lumi_nonorm_lin.root";
 
  vector<TString> rft_mee_hists = {plot_folder + plot_ee_part1 + to_string(0) + plot_ee_part2,
                                   plot_folder + plot_ee_part1 + to_string(4) + plot_ee_part2,
                                   plot_folder + plot_ee_part1 + to_string(8) + plot_ee_part2,
                                   plot_folder + plot_ee_part1 + to_string(12) + plot_ee_part2,
                                   plot_folder + plot_ee_part1 + to_string(16) + plot_ee_part2};
  
  vector<TString> rft_meeg_hists = {plot_folder + plot_eeg_part1 + to_string(0) + plot_eeg_part2,
                                    plot_folder + plot_eeg_part1 + to_string(4) + plot_eeg_part2,
                                    plot_folder + plot_eeg_part1 + to_string(8) + plot_eeg_part2,
                                    plot_folder + plot_eeg_part1 + to_string(12) + plot_eeg_part2,
                                    plot_folder + plot_eeg_part1 + to_string(16) + plot_eeg_part2};
 
  string outfile_ee_test = "plots/constrainedzfit-test-ee-unc"; 
  string outfile_eeg_test = "plots/constrainedzfit-test-eeg-unc";

  vector<string> histogram_rf_labels = {"1x unc","0.8x unc", "0.6x unc", "0.4x unc", "0.2x unc"};

  //overlay_histograms(rft_mee_hists,  histogram_rf_labels, "m_{ll}",       outfile_ee_test,  true); 
  //overlay_histograms(rft_meeg_hists, histogram_rf_labels, "m_{ll#gamma}", outfile_eeg_test, true);
  //--------------------------------------------------------------------------------------------------//

  rft_mee_hists = {plot_folder + plot_ee_part1 + to_string(0) + plot_ee_part2,
                   plot_folder + plot_ee_part1 + to_string(1) + plot_ee_part2,
                   plot_folder + plot_ee_part1 + to_string(2) + plot_ee_part2,
                   plot_folder + plot_ee_part1 + to_string(3) + plot_ee_part2,
                   plot_folder + plot_ee_part1 + to_string(4) + plot_ee_part2};
  
  rft_meeg_hists = {plot_folder + plot_eeg_part1 + to_string(0) + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + to_string(1) + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + to_string(2) + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + to_string(3) + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + to_string(4) + plot_eeg_part2};


  outfile_ee_test = "plots/constrainedzfit-test-ee-unc-closeup"; 
  outfile_eeg_test = "plots/constrainedzfit-test-eeg-unc-closeup";

  histogram_rf_labels = {"1x unc","0.95x unc", "0.9x unc", "0.85x unc", "0.8x unc"};

  //overlay_histograms(rft_mee_hists,  histogram_rf_labels, "m_{ll}",       outfile_ee_test,  true); 
  //overlay_histograms(rft_meeg_hists, histogram_rf_labels, "m_{ll#gamma}", outfile_eeg_test, true);
  
  //Test larger uncertainties
  rft_mee_hists = {plot_folder + plot_ee_part1 + "0" + plot_ee_part2,
                   plot_folder + plot_ee_part1 + "m3" + plot_ee_part2,
                   plot_folder + plot_ee_part1 + "m5" + plot_ee_part2,
                   plot_folder + plot_ee_part1 + "m7" + plot_ee_part2,
                   plot_folder + plot_ee_part1 + "m9" + plot_ee_part2,
                   "plots/an_kinematic_refit/an_kinematicfit_mZtruth__mZtruth__wgt__lumi_nonorm_lin.root"};
  
  rft_meeg_hists = {plot_folder + plot_eeg_part1 + "0" + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + "m3" + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + "m5" + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + "m7" + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + "m9" + plot_eeg_part2};
 


  outfile_ee_test = "plots/constrainedzfit-test-ee-unc-larger"; 
  outfile_eeg_test = "plots/constrainedzfit-test-eeg-unc-larger";

  histogram_rf_labels = {"1x unc","1.15x unc", "1.25x unc", "1.35x unc", 
                         "1.45x unc",
                         "truth m_{ll}"};

  //overlay_histograms(rft_mee_hists,  histogram_rf_labels, "m_{ll}",       outfile_ee_test,  true); 

  histogram_rf_labels = {"1x unc","1.15x unc", "1.25x unc", "1.35x unc", 
                         "1.45x unc"};


  overlay_histograms(rft_meeg_hists, histogram_rf_labels, "m_{ll#gamma}", outfile_eeg_test, true);

  //Test particular scalign plot
  rft_mee_hists = {plot_folder + plot_ee_part1 + "0" + plot_ee_part2,
                   //plot_folder + plot_ee_part1 + "5" + plot_ee_part2,
                   "plots/an_kinematic_refit/an_kinematicfit_mZtruth__mZtruth__wgt__lumi_nonorm_lin.root"};
  
  rft_meeg_hists = {plot_folder + plot_eeg_part1 + "0" + plot_eeg_part2,
                    plot_folder + plot_eeg_part1 + "5" + plot_eeg_part2};
 


  outfile_ee_test = "plots/constrainedzfit-test-ee-unc-0p75"; 
  outfile_eeg_test = "plots/constrainedzfit-test-eeg-unc-0p75";

  histogram_rf_labels = {//"1x unc",//"1.15x unc", "1.25x unc", "1.35x unc", 
                         "Z-mass fit",
                         "truth m_{ll}"};

  //overlay_histograms(rft_mee_hists,  histogram_rf_labels, "m_{ll}",       outfile_ee_test,  true); 
  overlay_histograms(rft_mee_hists,  histogram_rf_labels, "m_{ll}",       "plots/constrainedzfit-truth-comp",  true); 

  histogram_rf_labels = {"1x unc",//"1.15x unc", "1.25x unc", "1.35x unc", 
                         "0.75x unc"};


  //overlay_histograms(rft_meeg_hists, histogram_rf_labels, "m_{ll#gamma}", outfile_eeg_test, true);
 


  //Run 2 vs. Run 3 Z --> mumu test
  plot_folder = "./plots/refit_test_diff_unc/";
 
  vector<TString> rft_mmumu_hists = {plot_folder + "kf_test_unc_mlly_test_mumu_r2__mll_refit_test__wgt__lumi_nonorm_lin.root",
                                    plot_folder + "kf_test_unc_mlly_test_mumu_r3__mll_refit_test__wgt__lumi_nonorm_lin.root"};
  
  vector<TString> rft_mmumug_hists = {plot_folder + "kf_test_unc_mlly_test_mumug_r2__mlly_refit_test__wgt__lumi_nonorm_lin.root",
                                     plot_folder + "kf_test_unc_mlly_test_mumug_r3__mlly_refit_test__wgt__lumi_nonorm_lin.root"};


  string outfile_mumu_test = "plots/constrainedzfit-r2vr3-mumu-refit"; 
  string outfile_mumug_test = "plots/constrainedzfit-r2vr3-mumug-refit";

  histogram_rf_labels = {"Run 2", "Run 3"};

  //overlay_histograms(rft_mmumu_hists,  histogram_rf_labels, "m_{ll}",       outfile_mumu_test,  true); 
  //overlay_histograms(rft_mmumug_hists, histogram_rf_labels, "m_{ll#gamma}", outfile_mumug_test, true);
 

  //ee vs mumu test
  plot_folder = "./plots/refit_test_diff_unc/";
 
  vector<TString> rft_mee_mmumu_hists = {plot_folder + "kf_test_unc_mlly_test_0__mll_refit_test__wgt__lumi_nonorm_lin.root",
                                     plot_folder + "kf_test_unc_mlly_test_mumu_r2__mll_refit_test__wgt__lumi_nonorm_lin.root"};
  
  vector<TString> rft_meeg_mmumug_hists = {plot_folder + "kf_test_unc_mlly_test_0__mlly_refit_test__wgt__lumi_nonorm_lin.root",
                                      plot_folder + "kf_test_unc_mlly_test_mumug_r2__mlly_refit_test__wgt__lumi_nonorm_lin.root"};


  string outfile_ll_test = "plots/constrainedzfit-ee-vs-mumu-refit"; 
  string outfile_llg_test = "plots/constrainedzfit-eeg-vs-mumug-refit";

  histogram_rf_labels = {"Z #rightarrow e^{+}e^{-}", "Z #rightarrow #mu^{+}#mu^{-}"};

  //overlay_histograms(rft_mee_mmumu_hists,  histogram_rf_labels, "m_{ll}",       outfile_ll_test,  true); 
  //overlay_histograms(rft_meeg_mmumug_hists, histogram_rf_labels, "m_{ll#gamma}", outfile_llg_test, true);
 

  //ttH lep refit vs ggF test
  plot_folder = "./plots/convener_comments_120125/";
 
  vector<TString> rft_catcomp_hists = {plot_folder + "an_conv_response_mllgamma_refit_ggF_noimini__llphoton_refit_m__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_refit_tthlep__llphoton_refit_m__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_refit_tthlep_noimini__llphoton_refit_m__wgt__shapes_lin.root"
                                       };

  vector<TString> catcomp_hists     = {plot_folder + "an_conv_response_mllgamma_ggF_noimini__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_noimini__llphoton_m0__wgt__shapes_lin.root"
                                       };

  vector<TString> catcomp_ptg30_hists = {plot_folder + "an_conv_response_mllgamma_ggF_ptg30__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_ptg30__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_noimini_ptg30__llphoton_m0__wgt__shapes_lin.root"
                                       };


  vector<TString> catcomp_ptg50_hists = {plot_folder + "an_conv_response_mllgamma_ggF_ptg50__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_ptg50__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_noimini_ptg50__llphoton_m0__wgt__shapes_lin.root"
                                       };


  vector<TString> tth_noimini_hists = {
                                       plot_folder + "an_conv_response_mllgamma_tthlep_noimini__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_noimini_ptg30__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_noimini_ptg50__llphoton_m0__wgt__shapes_lin.root"
                                       };

  vector<TString> tth_hists = {
                                       plot_folder + "an_conv_response_mllgamma_tthlep__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_ptg30__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_tthlep_ptg50__llphoton_m0__wgt__shapes_lin.root"
                                       };


  vector<TString> ggF_hists = {
                                       plot_folder + "an_conv_response_mllgamma_ggF_noimini__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_ggF_ptg30__llphoton_m0__wgt__shapes_lin.root",
                                       plot_folder + "an_conv_response_mllgamma_ggF_ptg50__llphoton_m0__wgt__shapes_lin.root"
                                       };

  vector<TString> mll_ohdm_hists     = {plot_folder + "an_conv_response_ohdm_mll__ll_m0__wgt__shapes_lin.root",
                                        plot_folder + "an_conv_response_ohdm_mll_refit__ll_refit_m__wgt__shapes_lin.root"
                                       };


  vector<TString> mllgamma_ohdm_hists = {plot_folder + "an_conv_response_ohdm_mllgamma__llphoton_m0__wgt__shapes_lin.root",
                                         plot_folder + "an_conv_response_ohdm_mllgamma_refit__llphoton_refit_m__wgt__shapes_lin.root"
                                        };
 
  string outfile_refit_llg_catcomp = "plots/constrainedzfit-catcomp-refit";
  string outfile_llg_catcomp       = "plots/constrainedzfit-catcomp";
  string outfile_llg_catcomp_ptg30 = "plots/constrainedzfit-catcomp-ptg30";
  string outfile_llg_catcomp_ptg50 = "plots/constrainedzfit-catcomp-ptg50";
  string outfile_llg_ggf           = "plots/constrainedzfit-ggf"; 
  string outfile_llg_tthlep        = "plots/constrainedzfit-tthlep";
  string outfile_llg_tthlep_noimini= "plots/constrainedzfit-tthlep-noimini"; 
  string outfile_ll_ohdm           = "plots/constrainedzfit-ll-ohdm";
  string outfile_llg_ohdm          = "plots/constrainedzfit-llg-ohdm";

  vector<string> histogram_mll_rf_labels  = {"m_{ll}","m_{ll,refit}"}; 
  vector<string> histogram_mllg_rf_labels = {"m_{ll#gamma}","m_{ll#gamma,refit}"};
  vector<string> histogram_catcomp_labels = {"ggF","ttH lep.", "ttH lep. w/ I_{mini} < 0.1"};
  vector<string> histogram_pt_labels = {"p_{T}(#gamma) > 15","p_{T}(#gamma) > 30","p_{T}(#gamma) > 50"};


  overlay_histograms(rft_catcomp_hists,   histogram_catcomp_labels, "m_{ll#gamma}", outfile_refit_llg_catcomp, false, false); 
  overlay_histograms(catcomp_hists,       histogram_catcomp_labels, "m_{ll#gamma}", outfile_llg_catcomp,       false, false);
  overlay_histograms(catcomp_ptg30_hists, histogram_catcomp_labels, "m_{ll#gamma}", outfile_llg_catcomp_ptg30, false, false);
  overlay_histograms(catcomp_ptg50_hists, histogram_catcomp_labels, "m_{ll#gamma}", outfile_llg_catcomp_ptg50, false, false);
 
  overlay_histograms(mll_ohdm_hists,       histogram_mll_rf_labels,  "m_{ll}",       outfile_ll_ohdm,  true, true, "H#rightarrow#mu#mu");
  overlay_histograms(mllgamma_ohdm_hists,  histogram_mllg_rf_labels, "m_{ll#gamma}", outfile_llg_ohdm, true, false, "H#rightarrow#mu#mu");


  overlay_histograms(ggF_hists,         histogram_pt_labels, "m_{ll#gamma}", outfile_llg_ggf,            false, false);
  overlay_histograms(tth_hists,         histogram_pt_labels, "m_{ll#gamma}", outfile_llg_tthlep,         false, false);
  overlay_histograms(tth_noimini_hists, histogram_pt_labels, "m_{ll#gamma}", outfile_llg_tthlep_noimini, false, false);
 */

} //END OF SCRIPT




/*
  //Loops to cover
  std::vector<TString> flavor          = {"_ee","_mumu"};
  std::vector<TString> topology        = {"ggF_", "VBF_", "WH_3l_", "ZH_MET_", "ttH_had_", "ttH_lep_"};
  std::vector<TString> photon_location = {"_barrel", "_crack", "_endcap"};
  TString file_start="plots/kinematicfit_0327/KinFit_cat_";
//  std::vector< std::tuple<std::vector<std::string>, std::string, std::string, std::string, bool> > plot_vector = {};
  
//  plot_vector.push_back(std::make_tuple({"_lly_m_Refit_nCut80__kinFitH__lumi_nonorm_lin.root","_80_lly_m_nCut__llphoton_m0__lumi_nonorm_lin.root"},) );
 

  //filenames to loop
  std::vector<std::vector<TString>> file_lists = {{"_lly_m_Refit_llmCut80__kinFitH__lumi_nonorm_lin.root","_80_lly_m_llmCut__llphoton_m0__lumi_nonorm_lin.root"},
                                                  {"_delta_mll_truth_refit__mll_diff_refit__lumi_nonorm_lin.root","_delta_mll__mll_diff__lumi_nonorm_lin.root"},
                                                  {"_deltaPt_lminus_refit__l1_refit_ptdiff__lumi_nonorm_lin.root","_deltaPt_lminus__l1_ptdiff__lumi_nonorm_lin.root"},
                                                  {"_deltaPt_lplus_refit__l2_refit_ptdiff__lumi_nonorm_lin.root","_deltaPt_lplus__l2_ptdiff__lumi_nonorm_lin.root"},
                                                  {"_ll_m_Refit_nCut__kinFit__lumi_nonorm_lin.root","_ll_m_nCut__ll_m0__lumi_nonorm_lin.root"},
                                                  {"_ll_m_Refit_nCut__kinFit__lumi_nonorm_lin.root","_llm_truth__Z_truth__lumi_nonorm_lin.root"},

                                                  {"_ll_m_Refit_nCut__kinFit__lumi_nonorm_lin.root","_ll_m_nCut__ll_m0__lumi_nonorm_lin.root","_llm_truth__Z_truth__lumi_nonorm_lin.root"},

                                                  {"_ll_m_Refit_nCut_barrel__kinFit__lumi_nonorm_lin.root","_ll_m_Refit_nCut_crack__kinFit__lumi_nonorm_lin.root",
                                                   "_ll_m_Refit_nCut_endcap__kinFit__lumi_nonorm_lin.root"},

                                                  {"_lly_m_Refit_nCut_barrel_80__kinFitH__lumi_nonorm_lin.root","_lly_m_Refit_nCut_crack_80__kinFitH__lumi_nonorm_lin.root",
                                                   "_lly_m_Refit_nCut_endcap_80__kinFitH__lumi_nonorm_lin.root"},

                                                  {"_lly_m_Refit_nCut80__kinFitH__lumi_nonorm_lin.root","_80_lly_m_nCut__llphoton_m0__lumi_nonorm_lin.root","_lly_m_Refit_llmCut80__kinFitH__lumi_nonorm_lin.root"}

                                                 };

  std::vector<std::vector<std::string>> labels = {{"m_{ll#gamma,refit}","m_{ll#gamma}"},
                                                  {"#Delta m(ll_{refit},ll_{truth})","#Delta m(ll,ll_{truth})"},
                                                  {"#Delta p_{T}(l^{-}_{refit},l^{-}_{truth})","#Delta p_{T}(l^{-},l^{-}_{truth})"},
                                                  {"#Delta p_{T}(l^{+}_{refit},l^{+}_{truth})","#Delta p_{T}(l^{+},l^{+}_{truth})"},
                                                  {"m_{ll,refit}","m_{ll}"},
                                                  {"m_{ll,refit}","m_{ll,truth}"},

                                                  {"m_{ll,refit}","m_{ll}","m_{ll,truth}"},

                                                  {"m_{ll,refit} - barrel","m_{ll,refit} - crack","m_{ll,refit} - endcap"},

                                                  {"m_{ll#gamma,refit} - barrel","m_{ll#gamma,refit} - crack","m_{ll#gamma,refit} - endcap"},

                                                  {"m_{ll#gamma,refit}","m_{ll#gamma}","m_{ll#gamma,refit} - w/ m_{ll} selection"}
                                                 };


  std::vector<std::string>  x_axis_labels       = {"m_{ll#gamma}",
                                                  "#Delta m(ll,ll_{truth})",
                                                  "#Delta p_{T}(l^{-},l^{-}_{truth})",
                                                  "#Delta p_{T}(l^{+},l^{+}_{truth})",
                                                  "m_{ll}",
                                                  "m_{ll}",

                                                  "m_{ll}",

                                                  "m_{ll,refit}",

                                                  "m_{ll#gamma,refit}",
                          
                                                  "m_{ll#gamma}"
                                                  };


  std::vector<std::string> outfile             = {"_mlly_wmllsel",
                                                  "_mll_diff",
                                                  "_delta_ptm",
                                                  "_delta_ptp",
                                                  "_mll_refit_reco",
                                                  "_mll_refit_tru",

                                                  "_mll_refit_reco_tru",
                                                   
                                                  "_mll_bce_comp",
                                                  "_mlly_bce_comp",
                                                  
                                                  "_mlly_mllcomp"

                                                 };


  //Loop making plot overlays
  std::vector<TString> filenames = {};
  TString prefix = "";
  for(TString flav : flavor){
    for(TString top : topology){
      for(unsigned int idx_files = 0; idx_files < file_lists.size(); idx_files++){
        prefix = file_start + top + flav;

        filenames = filename_vector(prefix,file_lists[idx_files]);

        //make overlays
        if(filenames[0].Contains("barrel")){
         //overlay_histograms(filenames, labels[idx_files], x_axis_labels[idx_files], "plots/overlay_" + outfile[idx_files] + "_" + top + flav,true );
         continue; 
        }
        //overlay_histograms(filenames, labels[idx_files], x_axis_labels[idx_files], "plots/overlay_" + outfile[idx_files] + "_" + top + flav );
      }
    }
  }


}

*/
