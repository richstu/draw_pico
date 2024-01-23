#include "core/utilities.hpp"

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>

#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <glob.h>
#include <libgen.h>
#include <sys/stat.h>

#include "TMath.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TArrow.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;

mutex Multithreading::root_mutex;

set<string> Glob(const string &pattern){
  glob_t glob_result;
  glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result);
  set<string> ret;
  for(size_t i=0; i<glob_result.gl_pathc; ++i){
    ret.emplace(realpath(glob_result.gl_pathv[i], nullptr));
  }
  globfree(&glob_result);
  return ret;
}

string Basename(const string &filename){
  vector<char> c(filename.cbegin(), filename.cend());
  c.push_back(0);
  return string(basename(&c.at(0)));
}

bool Contains(const string &str, const string &pat){
  return str.find(pat) != string::npos;
}

bool StartsWith(const string &str, const string &pat){
  return str.find(pat) == 0;
}

void ReplaceAll(string &str, const string &orig, const string &rep){
  size_t loc = 0;
  while ((loc = str.find(orig, loc)) != string::npos) {
    str.replace(loc, orig.length(), rep);
    loc += rep.length();
  }
}

string CopyReplaceAll(string str, const string &orig, const string &rep){
  ReplaceAll(str, orig, rep);
  return str;
}

string LeftStrip(string s){
  s.erase(s.begin(), find_if(s.begin(), s.end(), [](char c){return isspace(c)==0;}));
  return s;
}

string RightStrip(string s){
  s.erase(find_if(s.rbegin(), s.rend(), [](char c){return isspace(c)==0;}).base(), s.end());
  return s;
}

string Strip(string s){
  return LeftStrip(RightStrip(s));
}

bool FileExists(const string &path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

string CodeToPlainText(string code){
  ReplaceAll(code, " ", "");
  ReplaceAll(code, "stitch &&", "");
  ReplaceAll(code, "stitch_htmet &&", "");
  ReplaceAll(code, "stitch_ht &&", "");
  ReplaceAll(code, "nbm==0","0b");
  ReplaceAll(code, "nbm==1","1b");
  ReplaceAll(code, "nbt==2&&nbm==2","2b");
  ReplaceAll(code, "nbt>=2&&nbm==3&&nbl==3","3b");
  ReplaceAll(code, "nbt>=2&&nbm>=3&&nbl>=4","4b");
  ReplaceAll(code, "nbt>=2&&nbm>=3","ge3b");
  ReplaceAll(code, "nbm==2","2b");
  ReplaceAll(code, "nbm==3&&nbl==3","3b");
  ReplaceAll(code, "nbm>=3&&nbl>=4","4b");
  ReplaceAll(code, "nbm>=3","ge3b");
  ReplaceAll(code, "nbt>=2","ge2b");
  ReplaceAll(code, "!(hig_cand_am[0]>100&&hig_cand_am[0]<=140)","SBD");
  ReplaceAll(code, "hig_cand_am[0]>100&&hig_cand_am[0]<=140","HIG");
  ReplaceAll(code, "boostedRegionIdx==0","0H, SBD");
  ReplaceAll(code, "boostedRegionIdx==1","1H, SBD");
  ReplaceAll(code, "boostedRegionIdx==2","2H, SBD");
  ReplaceAll(code, "boostedRegionIdx==3","0H, HIG");
  ReplaceAll(code, "boostedRegionIdx==4","1H, HIG");
  ReplaceAll(code, "boostedRegionIdx==5","2H, HIG");
  
  ReplaceAll(code, ".", "p");
  ReplaceAll(code, "(", "");
  ReplaceAll(code, ")", "");
  ReplaceAll(code, "[", "");
  ReplaceAll(code, "]", "");
  ReplaceAll(code, "{", "");
  ReplaceAll(code, "}", "");
  ReplaceAll(code, "+", "p");
  ReplaceAll(code, "-", "m");
  ReplaceAll(code, "*", "x");
  ReplaceAll(code, "/", "d");
  ReplaceAll(code, "%", "_");
  ReplaceAll(code, "!", "n");
  ReplaceAll(code, "&&", "_");
  ReplaceAll(code, "||", "_");
  ReplaceAll(code, "==", "");
  ReplaceAll(code, "<=", "le");
  ReplaceAll(code, ">=", "ge");
  ReplaceAll(code, ">", "g");
  ReplaceAll(code, "<", "l");
  ReplaceAll(code, "=", "");
  ReplaceAll(code, "&", "_");
  ReplaceAll(code, "|", "_");
  ReplaceAll(code, "^", "_");
  ReplaceAll(code, "~", "_");
  ReplaceAll(code, "___", "__");
  for(size_t i = 0; i < code.size(); ++i){
    if(isalnum(code.at(i)) || code.at(i) == '.' || code.at(i) == '_') continue;
    code = code.substr(0,i)+code.substr(i+1);
  }

  return code;
}

string CodeToRootTex(string code){
  ReplaceAll(code, " ", "");
  ReplaceAll(code, "stitch&&", "");
  ReplaceAll(code, "pass&&stitch", "ps");
  ReplaceAll(code, "&&1", "");
  ReplaceAll(code, "weight", "w");
  ReplaceAll(code, "nbm==0","0b");
  ReplaceAll(code, "nbm==1","1b");
  ReplaceAll(code, "nbt==2&&nbm==2","2b");
  ReplaceAll(code, "nbt>=2&&nbm==3&&nbl==3","3b");
  ReplaceAll(code, "nbt>=2&&nbm>=3&&nbl>=4","4b");
  ReplaceAll(code, "nbt>=2&&nbm>=3","#geq3b");
  ReplaceAll(code, "nbt>=2","#geq2b");
  ReplaceAll(code, "nbm==2","2b");
  ReplaceAll(code, "nbm==3&&nbl==3","3b");
  ReplaceAll(code, "nbm>=3&&nbl>=4","4b");
  ReplaceAll(code, "nbm>=3","#geq3b");
  ReplaceAll(code, "!(hig_cand_am[0]>100&&hig_cand_am[0]<=140)","SBD");
  ReplaceAll(code, "hig_cand_am[0]>100&&hig_cand_am[0]<=140","HIG");
  ReplaceAll(code, "hig_cand_am[0]", "#LTm#GT");
  ReplaceAll(code, "hig_cand_drmax[0]<=1.1", "#DeltaR_{max} < 1.1");
  ReplaceAll(code, "hig_cand_drmax[0]>1.1", "1.1 < #DeltaR_{max} < 2.2");
  ReplaceAll(code, "hig_cand_drmax[0]", "#DeltaR_{max}");
  ReplaceAll(code, "hig_cand_dm[0]", "#Delta m");

  ReplaceAll(code, "met>0&&met<=75", "0<p_{T}^{miss}#leq 75");
  ReplaceAll(code, "met>75&&met<=150", "75<p_{T}^{miss}#leq 150");
  ReplaceAll(code, "met>100&&met<=150", "100<p_{T}^{miss}#leq 150");
  ReplaceAll(code, "met>150&&met<=200", "150<p_{T}^{miss}#leq 200");
  ReplaceAll(code, "met>200&&met<=300", "200<p_{T}^{miss}#leq 300");
  ReplaceAll(code, "met>300&&met<=400", "300<p_{T}^{miss}#leq 400");
  ReplaceAll(code, "met>300&&met<=500", "300<p_{T}^{miss}#leq 500");
  ReplaceAll(code, "met>450&&met<=700", "450<p_{T}^{miss}#leq 700");
  ReplaceAll(code, "met>500&&met<=700", "500<p_{T}^{miss}#leq 700");
  ReplaceAll(code, "met>200&&met<=500", "200<p_{T}^{miss}#leq 500");

  ReplaceAll(code, "ll_pt[0]>0&&ll_pt[0]<=75","0<p_{T}^{Z}#leq 75");
  ReplaceAll(code, "ll_pt[0]>75&&ll_pt[0]<=150","75<p_{T}^{Z}#leq 150");
  ReplaceAll(code, "ll_pt[0]>150&&ll_pt[0]<=200","150<p_{T}^{Z}#leq 200");
  ReplaceAll(code, "ll_pt[0]>200&&ll_pt[0]<=300","200<p_{T}^{Z}#leq 300");
  ReplaceAll(code, "ll_pt[0]>300&&ll_pt[0]<=400","300<p_{T}^{Z}#leq 400");
  ReplaceAll(code, "ll_pt[0]","p_{T}^{Z}");
  
  ReplaceAll(code, "fjet_pt[0]>300&&fjet_pt[1]>300", "p_{T,J1},p_{T,J2}> 300");
  ReplaceAll(code, "((met<=300&&fjet_deep_md_hbb_btv[0]>0.6&&fjet_deep_md_hbb_btv[1]>0.6)||(met>300&&fjet_deep_md_hbb_btv[0]>0.3&&fjet_deep_md_hbb_btv[1]>0.3))", "DDB_{J1},DDB_{J2}> 0.3 (0.6)");
  ReplaceAll(code, "fjet_msoftdrop[0]>50&&fjet_msoftdrop[0]<=250&&fjet_msoftdrop[1]>50&&fjet_msoftdrop[1]<=250", "50<m_{J1},m_{J2}#leq 250");
  ReplaceAll(code, "fjet_msoftdrop[0]>85&&fjet_msoftdrop[0]<=135&&fjet_msoftdrop[1]>85&&fjet_msoftdrop[1]<=135", "85<m_{J1},m_{J2}#leq 135");
  ReplaceAll(code, "fjet_deep_md_hbb_btv[0]>0.7&&fjet_deep_md_hbb_btv[1]>0.7", "DDB_{J1},DDB_{J2}> 0.7");
  ReplaceAll(code, "!lowDphiFix", "hi-#Delta#phi");
  ReplaceAll(code, "!low_dphi_met", "hi-#Delta#phi");
  ReplaceAll(code, "low_dphi_met", "low #Delta#phi");
  ReplaceAll(code, "hig_drmax", "#DeltaR^{max}_{bb}");
  ReplaceAll(code, "leadingSignalLeptonPt", "lep. p_{T}");
  ReplaceAll(code, "leadingSignalMuonPt", "#mu p_{T}");
  ReplaceAll(code, "leadingSignalElectronPt", "e p_{T}");
  ReplaceAll(code, "ntk==0", "0 trk");
  ReplaceAll(code, "ntk", "N_{tks}");
  ReplaceAll(code, "nlep", "N_{lep}");
  ReplaceAll(code, "nvlep==0", "0L");
  ReplaceAll(code, "nvlep", "N_{lep}");
  ReplaceAll(code, "nmu", "N_{#mu}");
  ReplaceAll(code, "nel", "N_{e}");
  ReplaceAll(code, "nphoton", "N_{#gamma}");
  ReplaceAll(code, "nvmu", "N^{veto}_{#mu}");
  ReplaceAll(code, "nvel", "N^{veto}_{e}");
  ReplaceAll(code, "ntrulep", "N^{true}_{lep}");
  ReplaceAll(code, "npv", "N_{PV}");

  ReplaceAll(code, "resolved_baseline", "Resolved topology baseline selection");

  ReplaceAll(code, "njet>=4&&njet<=5","4-5j");
  ReplaceAll(code, "njet","N_{jet}");
  ReplaceAll(code, "nfjet","N_{J}");
  ReplaceAll(code, "abs(lep_id)==13&&","");
  ReplaceAll(code, ">=", " #geq ");
  ReplaceAll(code, ">", " > ");
  ReplaceAll(code, "<=", " #leq ");
  ReplaceAll(code, "<", " < ");
  ReplaceAll(code, "&&", ", ");
  ReplaceAll(code, "==", " = ");
  ReplaceAll(code, "met_calo", "E_{T,calo}^{miss}");
  ReplaceAll(code, "met", "p_{T}^{miss}");
  ReplaceAll(code, "mht", "H_{T}^{miss}");
  ReplaceAll(code, "ht", "H_{T}");
  ReplaceAll(code, "mt", "m_{T}");
  ReplaceAll(code, "nbm","N_{b,M}");
  ReplaceAll(code, "nbt","N_{b,T}");
  ReplaceAll(code, "nbl","N_{b,L}");

  ReplaceAll(code, "jetys_nob_pt", "ISR p_{T}");
  ReplaceAll(code, "(", "");
  ReplaceAll(code, ")", "");

  return code;
}

string CodeToLatex(string code){
  code = CodeToRootTex(code);
  ReplaceAll(code, "#", "\\");
  ReplaceAll(code, "\\DeltaR", "\\Delta R");
  return code;
}

vector<string> Tokenize(const string& input,
                        const string& tokens){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

string MakeDir(string prefix){
  prefix += "XXXXXX";
  char *dir_name = new char[prefix.size()];
  if(dir_name == nullptr) ERROR("Could not allocate directory name");
  strcpy(dir_name, prefix.c_str());
  (void)! mkdtemp(dir_name);
  prefix = dir_name;
  delete[] dir_name;
  return prefix;
}

string MakeTemp(string prefix){
  prefix += "XXXXXX";
  char *file_name = new char[prefix.size()];
  if(file_name == nullptr) ERROR("Could not allocate file name");
  strcpy(file_name, prefix.c_str());
  (void)! mkstemp(file_name);
  prefix = file_name;
  delete[] file_name;
  return prefix;
}

void AdjustDensityForBinWidth(TH1D &h){
  double entries = h.GetEntries();
  int nbins = h.GetNbinsX();
  double low = h.GetBinLowEdge(1);
  double high = h.GetBinLowEdge(nbins+1);
  double width = (high-low)/nbins;
  for(int bin = 1; bin <= nbins; ++bin){
    double content = h.GetBinContent(bin);
    double error = h.GetBinError(bin);
    double this_width = h.GetBinWidth(bin);
    double scale = width/this_width;
    h.SetBinContent(bin, content*scale);
    h.SetBinError(bin, error*scale);
  }
  h.SetEntries(entries);
}

string ChangeExtension(string path, const string &new_ext){
  auto pos = path.rfind(".");
  if(pos == string::npos){
    path += new_ext;
  }else{
    auto count = path.size() - pos;
    path = path.replace(pos, count, new_ext);
  }
  return path;
}


void getLegendBoxes(TLegend &leg, vector<vector<float> > &boxes){
  int nRows = leg.GetNRows();
  boxes = vector<vector<float> > (nRows, vector<float>(4,0));
  float x1 = leg.GetX1NDC(), y1 = leg.GetY1NDC(), x2 = leg.GetX2NDC(), y2 = leg.GetY2NDC();
  float rowH = (y2-y1)/nRows;
  for(int row=0; row<nRows; row++){
    float bx1 = x1+0.038*(x2-x1);
    float bx2 = x1+0.21*(x2-x1);
    float by1 = y2-row*rowH-0.2*rowH;
    float by2 = y2-row*rowH-0.8*rowH;
    boxes[row] = vector<float>({bx1, by1, bx2, by2});
  } // Loop over rows
}

void Normalize(TH1D &h, double normalization, bool norm_per_avg_width){
  int nbins = h.GetNbinsX();
  double low = h.GetBinLowEdge(1);
  double high = h.GetBinLowEdge(nbins+1);
  double width = (high-low)/nbins;
  if(norm_per_avg_width) normalization *= width;
  double integral = h.Integral("width");
  if (integral != 0)
    h.Scale(normalization/integral);
}

void MergeOverflow(TH1D &h, bool merge_underflow, bool merge_overflow){
  if(merge_underflow){
    h.SetBinContent(1, h.GetBinContent(0)+h.GetBinContent(1));
    h.SetBinContent(0, 0.);
    h.SetBinError(1, hypot(h.GetBinError(0), h.GetBinError(1)));
    h.SetBinError(0, 0.);
  }
  int nbins = h.GetNbinsX();
  if(merge_overflow){
    h.SetBinContent(nbins, h.GetBinContent(nbins)+h.GetBinContent(nbins+1));
    h.SetBinContent(nbins+1, 0.);
    h.SetBinError(nbins, hypot(h.GetBinError(nbins), h.GetBinError(nbins+1)));
    h.SetBinError(nbins+1, 0.);
  }
}

string FixedDigits(double x, int n_digits){
  int digits_left = max(floor(log10(x))+1., 0.);
  int digits_right = max(n_digits-digits_left, 0);

  double multiplier = pow(10., digits_right);

  ostringstream oss;
  oss << setprecision(numeric_limits<double>::digits10) << round(x*multiplier)/multiplier << flush;
  string out = oss.str();
  if(out.substr(0,2) == "0."){
    out = out.substr(1);
  }
  return out;
}

string FullTitle(const TH1 &h){
  return string(h.GetTitle())
    +";"+h.GetXaxis()->GetTitle()
    +";"+h.GetYaxis()->GetTitle()
    +";"+h.GetZaxis()->GetTitle();
}

std::vector<double> ScaleBins(const TAxis &axis, double scale){
  vector<double> bins(axis.GetNbins()+1);
  for(size_t i = 0; i < bins.size(); ++i){
    bins.at(i) = scale * axis.GetBinLowEdge(i+1);
  }
  return bins;
}

TH1D ScaleAxes(const TH1 &h, double scale, const std::string &axes){
  double xscale = (Contains(axes,"x") || Contains(axes,"X")) ? scale : 1.;
  double yscale = (Contains(axes,"y") || Contains(axes,"Y")) ? scale : 1.;
  //double zscale = (Contains(axes,"z") || Contains(axes,"Z")) ? scale : 1.;
  
  vector<double> xbins = ScaleBins(*h.GetXaxis(), xscale);
  
  TH1D hout("",FullTitle(h).c_str(),
            xbins.size()-1, &xbins.front());
  
  int nbins = h.GetNcells();
  for(int bin = 0; bin < nbins; ++bin){
    hout.SetBinContent(bin, yscale*h.GetBinContent(bin));
    hout.SetBinError(bin, yscale*h.GetBinError(bin));
  }

  CopyStyle(h, hout);

  return hout;
}

TH2D ScaleAxes(const TH2 &h, double scale, const std::string &axes){
  double xscale = (Contains(axes,"x") || Contains(axes,"X")) ? scale : 1.;
  double yscale = (Contains(axes,"y") || Contains(axes,"Y")) ? scale : 1.;
  double zscale = (Contains(axes,"z") || Contains(axes,"Z")) ? scale : 1.;
  
  vector<double> xbins = ScaleBins(*h.GetXaxis(), xscale);
  vector<double> ybins = ScaleBins(*h.GetYaxis(), yscale);
  
  TH2D hout("",FullTitle(h).c_str(),
            xbins.size()-1, &xbins.front(),
            ybins.size()-1, &ybins.front());
  
  int nbins = h.GetNcells();
  for(int bin = 0; bin < nbins; ++bin){
    hout.SetBinContent(bin, zscale*h.GetBinContent(bin));
    hout.SetBinError(bin, zscale*h.GetBinError(bin));
  }

  CopyStyle(h, hout);

  return hout;
}

TH3D ScaleAxes(const TH3 &h, double scale, const std::string &axes){
  double xscale = (Contains(axes,"x") || Contains(axes,"X")) ? scale : 1.;
  double yscale = (Contains(axes,"y") || Contains(axes,"Y")) ? scale : 1.;
  double zscale = (Contains(axes,"z") || Contains(axes,"Z")) ? scale : 1.;
  
  vector<double> xbins = ScaleBins(*h.GetXaxis(), xscale);
  vector<double> ybins = ScaleBins(*h.GetYaxis(), yscale);
  vector<double> zbins = ScaleBins(*h.GetZaxis(), zscale);
  
  TH3D hout("",FullTitle(h).c_str(),
            xbins.size()-1, &xbins.front(),
            ybins.size()-1, &ybins.front(),
            zbins.size()-1, &zbins.front());
  
  int nbins = h.GetNcells();
  for(int bin = 0; bin < nbins; ++bin){
    hout.SetBinContent(bin, h.GetBinContent(bin));
    hout.SetBinError(bin, h.GetBinError(bin));
  }

  CopyStyle(h, hout);

  return hout;
}

void CopyStyle(const TH1 &hin, TH1 &hout){
  hout.SetLineStyle(hin.GetLineStyle());
  hout.SetLineWidth(hin.GetLineWidth());
  hout.SetLineColor(hin.GetLineColor());

  hout.SetFillColor(hin.GetFillColor());
  hout.SetFillStyle(hin.GetFillStyle());

  hout.SetMarkerStyle(hin.GetMarkerStyle());
  hout.SetMarkerSize(hin.GetMarkerSize());
  hout.SetMarkerColor(hin.GetMarkerColor());
}

TString HoursMinSec(float fseconds){
  long seconds = round(fseconds);
  int minutes((seconds/60)%60), hours(seconds/3600);
  TString hhmmss("");
  if(hours<10) hhmmss += "0";
  hhmmss += hours; hhmmss += ":";
  if(minutes<10) hhmmss += "0";
  hhmmss += minutes; hhmmss += ":";
  if((seconds%60)<10) hhmmss += "0";
  hhmmss += seconds%60; 

  return hhmmss;
}

TString AddCommas(double num){
  TString result(""); result += num;
  int posdot(result.First('.'));
  if(posdot==-1) posdot = result.Length();
  for(int ind(posdot-3); ind > 0; ind -= 3)
    result.Insert(ind, ",");
  return result;
}


TString RoundNumber(double num, int decimals, double denom){
  if(denom==0 || !isfinite(num) || !isfinite(denom)) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  if(abs(num) > 1e16) return "-";
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  if(result.Length()>15) cout<<"num "<<num<<", denom "<<denom<<"  ---->  "<<result<<endl;
  return result;
}

// Code from http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
// Parameter b could be theta...
double gsl_ran_gamma(const double a, const double b, TRandom3 &rand){
  if (a < 1){
    double u = rand.Uniform(1);
    return gsl_ran_gamma(1.0 + a, b, rand) * pow (u, 1.0 / a);
  }

  double x, v, u;
  double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt (d);
  
  while (1) {
    do {
      x = rand.Gaus(0, 1.0);
      v = 1.0 + c * x;
    }
    while (v <= 0);
      
    v = v * v * v;
    u = rand.Uniform(1);

    if (u < 1 - 0.0331 * x * x * x * x) 
      break;

    if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
      break;
  }
    
  return b * d * v;
}

double intGaus(double mean, double sigma, double minX, double maxX){
  return (TMath::Erf((maxX-mean)/sigma/sqrt(2.))-TMath::Erf((minX-mean)/sigma/sqrt(2.)))/2.;
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
  return hypot(TVector2::Phi_mpi_pi(phi2-phi1), eta2-eta1);
}

double deltaPhi(double phi1, double phi2){
  const double PI = acos(-1.);
  double dphi = fmod(fabs(phi2-phi1), 2.*PI);
  return dphi>PI ? 2.*PI-dphi : dphi;
}

float utilities::calculate_pvalue(float observed, float prediction, float prediction_up_diff, float prediction_down_diff) {
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

    float observed_toy = random.Poisson(fluctuated_mean); // shifts mean to high

    // Calculate test statistic
    if (observed_toy>=observed) Nabove++;
    else Nbelow++;

  }

  if(Nabove==0) return 1./Nmax;
  else if (Nbelow==0) return 1-1./Nmax;
  return (Nabove)/(Nabove+Nbelow);
}
float utilities::to_pvalue(float significance) {return 0.5-TMath::Erf(double(significance/sqrt(2)))/2;}
float utilities::to_significance(float pvalue) {return ROOT::Math::normal_quantile_c(pvalue, 1);}

// Below is an old method, use utilities::to_significance(utilities::calculate_pvalue()) instead.
// Finds significance of observation Nbkg, for an expected background Nbkg+Eup_bkg-Edown_bkg
// The mean of the Poisson is Nbkg convolved with the asymmetric lognormal(Eup_bkg, Edown_bkg)
double Significance(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg){
  double precision = 0.03; // Desired precision on the p-value
  double Nmin = 1/pow(precision,2), Nmax = 5e7; // Nmax controls max sigmas achievable (5e7->5.5 sigma)
  double Nbelow=0, Nabove=0, Nequal=0;
  if(Edown_bkg<0) Edown_bkg = Eup_bkg; // If down uncertainty not specified, symmetric lognormal
  if (Nbkg==0) Edown_bkg = 0;
  if(Edown_bkg>Nbkg) {
    cout<<"Down uncertainty ("<<Edown_bkg<<") has to be smaller than Nbkg ("<<Nbkg<<")"<<endl;
    return -999.;
  }
  //if(Nobs==Nbkg) return 0.;
  TRandom3 rand(1234);
  double mu, valG;
  while( (min(Nbelow,Nabove)+Nequal)<Nmin && (Nbelow+Nequal+Nabove)<Nmax){
    // Convolving expected bkg with log-normal
    mu = Nbkg;
    valG = rand.Gaus(0,1);
    if(mu==0) mu = fabs(valG)*Eup_bkg; // Apply 2-sided Gaussian uncertainty
    else if(valG>=0) mu *= exp(valG*log(1+Eup_bkg/Nbkg));
    else if(Edown_bkg<0.8*Nbkg) mu *= exp(-valG*log(1-Edown_bkg/Nbkg));
    else mu = max(0., mu + valG*Edown_bkg);
    // Finding if toy above the observed yield
    double valPois = rand.PoissonD(mu);
    if(valPois>Nobs) Nabove++;
    else if(valPois==Nobs) Nequal++;
    else Nbelow++;
  }

  //cout<<"Nobs: "<<Nobs<<" Nbkg: "<<Nbkg<<" Eup_bkg: "<<Eup_bkg<<" Edown_bkg: "<<Edown_bkg<<" Nabove: "<<Nabove<<" Nequal: "<<Nequal<<" Nbelow: "<<Nbelow<<endl;

  if(Nabove+Nequal==0){
    cout<<"No toys above or at Nobs="<<Nobs<<" for Nbkg "<<Nbkg<<"+"<<Eup_bkg<<"-"<<Edown_bkg<<". Returning "
	<<RooStats::PValueToSignificance(1/Nmax)<<endl;
    return RooStats::PValueToSignificance(1/Nmax);
  }
  if(Nbelow+Nequal==0){
    cout<<"No toys below or at Nobs="<<Nobs<<" for Nbkg "<<Nbkg<<"+"<<Eup_bkg<<"-"<<Edown_bkg<<". Returning "
	<<RooStats::PValueToSignificance(1-1/Nmax)<<endl;
    return RooStats::PValueToSignificance(1-1/Nmax);
  }
  return RooStats::PValueToSignificance((Nabove+Nequal/2.)/(Nbelow+Nequal+Nabove));
}

// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(vector<vector<float> > &entries, vector<vector<float> > &weights,
                 vector<float> &powers, float &mSigma, float &pSigma, bool do_data,
                 bool verbose, double syst, bool do_plot, int nrep, float nSigma){
  TRandom3 rand(1234);
  int nbadk(0);
  vector<float> fKappas;
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
    float stdyield(0.);
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
  if(do_plot) {
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

set<string> attach_folder(string folder, set<string> &fileset) {
  set<string> fset = set<string>();
  for (auto &ifile: fileset) {
    fset.insert(folder+ifile);
    cout<<folder+ifile<<endl;
  }
  return fset; 
}

set<string> attach_folder(string base_folder, set<int> years, string sample_folder, set<string> fileset) {
  set<string> fset = set<string>();
  for (auto & year: years) {
    for (auto &ifile: fileset) {
      fset.insert(base_folder+"/"+to_string(year)+"/"+sample_folder+"/"+ifile);
       cout<<base_folder+"/"+to_string(year)+"/"+sample_folder+"/"+ifile<<endl;
    }
  }
  return fset; 
}

set<string> attach_folder(string base_folder, set<string> years, string sample_folder, set<string> fileset) {
  set<string> fset = set<string>();
  for (auto & year: years) {
    for (auto &ifile: fileset) {
      fset.insert(base_folder+"/"+year+"/"+sample_folder+"/"+ifile);
       cout<<base_folder+"/"+year+"/"+sample_folder+"/"+ifile<<endl;
    }
  }
  return fset; 
}


void parseMasses(const string &prs, int &mglu, int &mlsp){
  mglu = stoi(prs.substr(prs.find("Chi-")+4,prs.find("_mLSP")-prs.find("ino-")-4));
  mlsp = stoi(prs.substr(prs.find("LSP-")+4,prs.find("_Tune")-prs.find("LSP-")-4));
}


void parseMassesGluino(const string &prs, int &mglu, int &mlsp){
  std::cout << prs.substr(prs.find("ino-")+4,prs.find("_mLSP")-prs.find("ino-")-4) << std::endl;
  std::cout << prs.substr(prs.find("LSP-")+4,prs.find("_Tune")-prs.find("LSP-")-4) << std::endl;
  mglu = stoi(prs.substr(prs.find("ino-")+4,prs.find("_mLSP")-prs.find("ino-")-4));
  mlsp = stoi(prs.substr(prs.find("LSP-")+4,prs.find("_Tune")-prs.find("LSP-")-4));
}
