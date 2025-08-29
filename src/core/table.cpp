#include "core/table.hpp"

#include <fstream>
#include <iomanip>

#include <sys/stat.h>

#include "RooStats/NumberCountingUtils.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPie.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

#include "core/utilities.hpp"

using namespace std;

namespace{
  std::string ToLatex(std::string x){
    if (x.find("!!") != string::npos) {
      ReplaceAll(x,"!!","");
      ReplaceAll(x,"_","\\_");
      return x;
    }
    ReplaceAll(x, "#", "\\");
    auto pos_1 = x.find("\\");
    auto pos_2 = x.find("_");
    auto pos_3 = x.find("^");
    auto pos_4 = x.find("{");
    auto pos_5 = x.find("}");
    if(pos_1 != string::npos
        || pos_2 != string::npos
        || pos_3 != string::npos
        || pos_4 != string::npos
        || pos_5 != string::npos){
      x = "$"+x+"$";
    }
    return x;
  }
}

Table::TableColumn::TableColumn(const Table &table,
    const shared_ptr<Process> &process):
  FigureComponent(table, process),
  sumw_(table.rows_.size(), 0.),
  sumw2_(table.rows_.size(), 0.),
  proc_and_table_cut_(table.rows_.size(), process->cut_),
  cut_vector_(),
  wgt_vector_(),
  val_vector_(){
    for(size_t irow = 0; irow < table.rows_.size(); ++irow){
      proc_and_table_cut_.at(irow) = table.rows_.at(irow).cut_ && process->cut_;
    }
  }

void Table::TableColumn::RecordEvent(const Baby &baby){
  const Table& table = static_cast<const Table&>(figure_);

  bool have_vector;
  size_t min_vec_size;
  for(size_t irow = 0; irow < table.rows_.size(); ++irow){
    have_vector = false;
    min_vec_size = 0;

    const TableRow& row = table.rows_.at(irow);
    if(!row.is_data_row_) continue;
    const NamedFunc &cut = proc_and_table_cut_.at(irow);
    const NamedFunc &wgt = row.weight_;

    if(cut.IsScalar()){
      if(!cut.GetScalar(baby)) continue;

    }else{
      cut_vector_ = cut.GetVector(baby);
      if(!have_vector || cut_vector_.size() < min_vec_size){
        have_vector = true;
        min_vec_size = cut_vector_.size();
      }
    }

    NamedFunc::ScalarType wgt_scalar = 0.;
    if(wgt.IsScalar()){
      wgt_scalar = wgt.GetScalar(baby);
    }else{
      wgt_vector_ = wgt.GetVector(baby);
      if(!have_vector || wgt_vector_.size() < min_vec_size){
        have_vector = true;
        min_vec_size = wgt_vector_.size();
      }
    }

    if(!have_vector){
      sumw_.at(irow) += wgt_scalar;
      sumw2_.at(irow) += wgt_scalar*wgt_scalar;
    }else{
      for(size_t iobject = 0; iobject < min_vec_size; ++iobject){
        NamedFunc::ScalarType this_cut = cut.IsScalar() ? true : cut_vector_.at(iobject);
        if(!this_cut) continue;
        NamedFunc::ScalarType this_wgt = wgt.IsScalar() ? wgt_scalar : wgt_vector_.at(iobject);
        sumw_.at(irow) += this_wgt;
        sumw2_.at(irow) += this_wgt*this_wgt;
      }
    }
  }
}

Table::Table(const string &name,
    const vector<TableRow> &rows,
    const vector<shared_ptr<Process> > &processes,
    bool do_zbi,
    bool print_table,
    bool print_pie,
    bool print_titlepie, 
    bool do_eff,
    bool do_unc,
    bool rel_prev_cut):
  Figure(),
  name_(name),
  rows_(rows),
  do_zbi_(do_zbi),
  print_table_(print_table),
  print_pie_(print_pie),
  print_titlepie_(print_titlepie),
  do_eff_(do_eff),
  do_unc_(do_unc),
  rel_prev_cut_(rel_prev_cut),
  plot_options_({PlotOpt("txt/plot_styles.txt", "Pie")}),
  backgrounds_(),
  signals_(),
  datas_(){
    precision_ = -1;
    for(const auto &process: processes){
      switch(process->type_){
        case Process::Type::data:
          datas_.emplace_back(new TableColumn(*this, process));
          break;
        case Process::Type::background:
          backgrounds_.emplace_back(new TableColumn(*this, process));
          break;
        case Process::Type::signal:
          signals_.emplace_back(new TableColumn(*this, process));
          break;
        default:
          break;
      }
    }


  }

void Table::Print(double luminosity,
    const string &subdir){
  if(!print_table_) return;
  if(subdir != "") mkdir(("tables/"+subdir).c_str(), 0777);
  string fmt_lumi = CopyReplaceAll(RoundNumber(luminosity,1).Data(),".","p");
  if (luminosity_tag_ != "") fmt_lumi = CopyReplaceAll(luminosity_tag_,".","p");
  string file_name = subdir != ""
    ? "tables/"+subdir+"/"+name_+"_lumi_"+fmt_lumi+".tex"
    : "tables/"+name_+"_lumi_"+fmt_lumi+".tex";
  if (Contains(name_, "FixName:")) {
    string tagName=name_;
    ReplaceAll(tagName, "FixName:", "");
    file_name = subdir != ""
      ? "tables/"+subdir+"/"+tagName+".tex"
      : "tables/"+tagName+".tex";
  }
  std::ofstream file(file_name);
  if (print_pie_) file << fixed << setprecision(4);
  else file << fixed << setprecision(4);
  PrintHeader(file, luminosity);
  for(size_t i = 0; i < rows_.size(); ++i){
    PrintRow(file, i, luminosity);
  }
  PrintFooter(file);
  file << flush;
  file.close();
  cout << "pdflatex " << file_name << " &> /dev/null; #./python/texify.py tables" << endl;
}

vector<GammaParams> Table::Yield(const Process *process, double luminosity) const{
  const auto &component_list = GetComponentList(process);
  const TableColumn *col = nullptr;
  for(const auto &component: component_list){
    if(component->process_.get() == process){
      col = static_cast<const TableColumn *>(component.get());
    }
  }
  if(col == nullptr) return vector<GammaParams>();
  vector<GammaParams> yields(rows_.size());
  for(size_t i = 0; i < yields.size(); ++i){
    yields.at(i).SetYieldAndUncertainty(luminosity*col->sumw_.at(i), luminosity*sqrt(col->sumw2_.at(i)));
  }
  return yields;
}

vector<GammaParams> Table::BackgroundYield(double luminosity) const{
  vector<GammaParams> yields(rows_.size());  
  auto procs = GetProcesses();
  for(const auto &proc: procs){
    if(proc->type_ != Process::Type::background) continue;
    vector<GammaParams> proc_yields = Yield(proc, luminosity);
    for(size_t i = 0; i < proc_yields.size(); ++i){
      yields.at(i) += proc_yields.at(i);
    }
  }
  return yields;
}

vector<GammaParams> Table::DataYield() const{
  vector<GammaParams> yields(rows_.size());  
  auto procs = GetProcesses();
  for(const auto &proc: procs){
    if(proc->type_ != Process::Type::data) continue;
    vector<GammaParams> proc_yields = Yield(proc, 1.);
    for(size_t i = 0; i < proc_yields.size(); ++i){
      yields.at(i) += proc_yields.at(i);
    }
  }
  return yields;
}

set<const Process*> Table::GetProcesses() const{
  set<const Process*> processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: signals_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: datas_){
    processes.insert(proc->process_.get());
  }
  return processes;
}

Figure::FigureComponent * Table::GetComponent(const Process *process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_.get() == process){
      return component.get();
    }
  }
  DBG("Could not find histogram for process "+process->name_+".");
  return nullptr;
}

const vector<unique_ptr<Table::TableColumn> >& Table::GetComponentList(const Process *process) const{
  switch(process->type_){
    case Process::Type::data:
      return datas_;
    case Process::Type::background:
      return backgrounds_;
    case Process::Type::signal:
      return signals_;
    default:
      ERROR("Did not understand process type "+to_string(static_cast<long>(process->type_))+".");
      return backgrounds_;
  }
}
void Table::PrintHeader(ofstream &file, double luminosity) const{
  file << "\\documentclass[10pt,oneside]{report}\n";
  file << "\\usepackage{graphicx,xspace,amssymb,amsmath,colordvi,colortbl,verbatim,multicol}\n";
  file << "\\usepackage{multirow, rotating}\n\n";
  file << "\\usepackage[active,tightpage]{preview}\n\n";
  file << "\\usepackage{siunitx}\n";
  file << "\\sisetup{round-mode = figures, round-precision=2}\n\n";
  file << "\\renewcommand{\\arraystretch}{1.1}\n\n";

  file << "\\begin{document}\n";
  file << "\\begin{preview}\n";
  file << "  \\begin{tabular}{ l";

  if(backgrounds_.size() > 1){
    file << " | ";
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      file << 'r';
    }
    file << " | r";
  }else if(backgrounds_.size() == 1){
    file << " | r";
  }

  if(datas_.size() > 1){
    file << " | ";
    for(size_t i = 0; i < datas_.size(); ++i){
      file << 'r';
    }
    file << " | r";
  }else if(datas_.size() == 1){
    file << " | r";
  }

  for(size_t i = 0; i < signals_.size(); ++i){
    if(do_zbi_) file << " | rrr"; // ***EDITED LINE CHANGED rr --> rrr***
    else file << " | r";
  }

  file << " }\n";
  file << "    \\hline\\hline\n";
  if (luminosity_tag_ != "") file <<" \\multicolumn{1}{c|}{${\\cal L} = "<<setprecision(1)<<luminosity_tag_<<"$ fb$^{-1}$} ";
  else file <<" \\multicolumn{1}{c|}{${\\cal L} = "<<setprecision(1)<<luminosity<<"$ fb$^{-1}$} ";
  if (precision_<0) {
    if(do_unc_)
      file <<setprecision(2);
    else
      file <<setprecision(1);
  } else {
    file <<setprecision(precision_);
  }

  if(backgrounds_.size() > 1){
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      file << " & " << ToLatex(backgrounds_.at(i)->process_->name_);
    }
    file << " & SM Bkg.";
  }else if(backgrounds_.size() == 1){
    file << " & " <<ToLatex(backgrounds_.front()->process_->name_);
  }

  if(datas_.size() > 1){ 
    file << " & ";
    for(size_t i = 0; i < datas_.size(); ++i){
      file << " & " << ToLatex(datas_.at(i)->process_->name_);
    }
    file << " & Data Tot.";
  }else if(datas_.size() == 1){
    file << " & " << ToLatex(datas_.front()->process_->name_);
  }

  for(size_t i = 0; i < signals_.size(); ++i){
    file << " & " << ToLatex(signals_.at(i)->process_->name_);
    if(do_zbi_){
      file << " & $Z_{\\text{Bi}}$";
      file << " & $\\text{Cowan}$"; // ***NEW LINE***
    }
  }

  file << "\\\\\n";
  file << "\\hline\n";
}

void Table::PrintRow(ofstream &file, size_t irow, double luminosity) const{
  const TableRow& row = rows_.at(irow);
  if(row.lines_before_ > 0){
    file << "    ";
    for(size_t i = 0; i < row.lines_before_; ++i){
      file << "\\hline";
    }
    file << "\n";
  }

  size_t cb = 1; //countback
  if (rel_prev_cut_ == 0) cb = irow;
  double totyield(0), totyield_err(0);
  double eff(0), eff_relUnc(0);

  if(row.is_data_row_){
    file << "    " << setw(52) << row.label_;



    if(backgrounds_.size() >= 1){//backgrounds
      totyield = luminosity*GetYield(backgrounds_, irow);
      totyield_err = luminosity*GetError(backgrounds_, irow);
      if(irow!=0){
        eff = 100*totyield / (luminosity*GetYield(backgrounds_, irow-cb));
        if(do_eff_) eff_relUnc = hypot(luminosity*GetError(backgrounds_,irow)/totyield, GetError(backgrounds_,irow-cb)/GetYield(backgrounds_,irow-cb));
      }

      if((backgrounds_.size() == 1 && irow == 0) || (backgrounds_.size() == 1 && do_eff_ != 1)){//Single background table entries for eff(y/n) and unc(y/n). Single background meaningless for pie chart
        WriteSeg(file, totyield, totyield_err, do_unc_,0);
      }else if(backgrounds_.size() == 1 && irow != 0 && do_eff_ == 1) {
        WriteSeg(file, eff, eff*eff_relUnc, do_unc_,0);
      }else if(backgrounds_.size() > 1){//multi background tables, w/ pie chart details

        for(size_t i = 0; i < backgrounds_.size(); ++i){//backgrounds loop
          if(print_pie_){
            eff = luminosity*backgrounds_.at(i)->sumw_.at(irow)/totyield;
            eff_relUnc = luminosity*sqrt(backgrounds_.at(i)->sumw2_.at(irow))/totyield;
            WriteSeg(file, eff, eff_relUnc, do_unc_,0);//Fraction of total yield
          }
          if(do_eff_ && irow == 0){
            eff = luminosity*backgrounds_.at(i)->sumw_.at(irow);
            WriteSeg(file, eff, eff_relUnc, do_unc_,0);
          } else if(do_eff_ && irow !=0){
            eff = 100*backgrounds_.at(i)->sumw_.at(irow)/backgrounds_.at(i)->sumw_.at(irow-cb);
            eff_relUnc = hypot(sqrt(backgrounds_.at(i)->sumw2_.at(irow))/backgrounds_.at(i)->sumw_.at(irow),sqrt(backgrounds_.at(i)->sumw2_.at(irow-1))/backgrounds_.at(i)->sumw_.at(irow-1));
            WriteSeg(file, eff, eff_relUnc, do_unc_,0);
          } else if(do_eff_ != 1){
            WriteSeg(file, luminosity*backgrounds_.at(i)->sumw_.at(irow), luminosity*sqrt(backgrounds_.at(i)->sumw2_.at(irow)), do_unc_, 1);
          }
        }//backgrounds loop

        if(do_eff_ && irow == 0){//Does the last column of table with multiple bkg processes, or only row of table with 1 process
          WriteSeg(file, totyield, totyield_err, do_unc_, 0);
        } else if(do_eff_){
          WriteSeg(file, eff, eff_relUnc, do_unc_, 0);
        } else {
          WriteSeg(file, totyield, totyield_err, do_unc_, 1);
        }

      }//multi background tables 
    }//backgrounds

    if(signals_.size() >= 1){//signals
      totyield = luminosity*GetYield(signals_, irow);
      totyield_err = luminosity*GetError(signals_, irow);
      if(irow!=0){
        eff = 100*totyield / (luminosity*GetYield(signals_, irow-cb));
        if(do_eff_) eff_relUnc = hypot(luminosity*GetError(signals_,irow)/totyield, GetError(signals_,irow-cb)/GetYield(signals_,irow-cb));
      }

      if((signals_.size() == 1 && irow == 0) || (signals_.size() == 1 && do_eff_ != 1)){//Single signal table entries for eff(y/n) and unc(y/n). Single signal meaningless for pie chart
        WriteSeg(file, totyield, totyield_err, do_unc_,0);
      }else if(signals_.size() == 1 && irow != 0 && do_eff_ == 1) {
        WriteSeg(file, eff, eff*eff_relUnc, do_unc_,0);
      }else if(signals_.size() > 1){//multi signal tables, w/ pie chart details

        for(size_t i = 0; i < signals_.size(); ++i){//signals loop
          if(print_pie_){
            eff = luminosity*signals_.at(i)->sumw_.at(irow)/totyield;
            eff_relUnc = luminosity*sqrt(signals_.at(i)->sumw2_.at(irow))/totyield;
            WriteSeg(file, eff, eff_relUnc, do_unc_,0);//Fraction of total yield
          }
          if(do_eff_ && irow == 0){
            eff = luminosity*signals_.at(i)->sumw_.at(irow);
            WriteSeg(file, eff, eff_relUnc, do_unc_,0);
          } else if(do_eff_ && irow !=0){
            eff = 100*signals_.at(i)->sumw_.at(irow)/signals_.at(i)->sumw_.at(irow-cb);
            eff_relUnc = hypot(sqrt(signals_.at(i)->sumw2_.at(irow))/signals_.at(i)->sumw_.at(irow),sqrt(signals_.at(i)->sumw2_.at(irow-1))/signals_.at(i)->sumw_.at(irow-1));
            WriteSeg(file, eff, eff_relUnc, do_unc_,0);
          } else if(do_eff_ != 1){
            WriteSeg(file, luminosity*signals_.at(i)->sumw_.at(irow), luminosity*sqrt(signals_.at(i)->sumw2_.at(irow)), do_unc_, 1);
          }
          if(do_zbi_){//haven't updated this
            file << " & " << RooStats::NumberCountingUtils::BinomialExpZ(luminosity*signals_.at(i)->sumw_.at(irow),
                luminosity*GetYield(backgrounds_, irow),
                GetError(backgrounds_, irow)/GetYield(backgrounds_, irow));
            //hypot(GetError(backgrounds_, irow)/GetYield(backgrounds_, irow), 0.3));
            //double sigma_b = hypot(luminosity*GetError(backgrounds_, irow), 0.3*luminosity*GetYield(backgrounds_, irow));
            double sigma_b = luminosity*GetError(backgrounds_, irow);
            double signalYield = luminosity*signals_.at(i)->sumw_.at(irow);
            double bkgdYield = luminosity*GetYield(backgrounds_, irow);
            double cowan1 = log((signalYield + bkgdYield)*(bkgdYield + sigma_b*sigma_b)/(bkgdYield*bkgdYield + (signalYield + bkgdYield)*sigma_b*sigma_b));
            cowan1 = cowan1 * (signalYield + bkgdYield);
            double cowan2 = log(1 + sigma_b*sigma_b*signalYield/(bkgdYield*(bkgdYield + sigma_b*sigma_b)));
            cowan2 = cowan2 * bkgdYield*bkgdYield / (sigma_b*sigma_b);
            double cowan = sqrt(2 * (cowan1 - cowan2));
            file << " & " <<  cowan; // ***NEW LINE***
        }
 
        }//signals loop

        if(do_eff_ && irow == 0){//Does the last column of table with multiple sig processes, or only row of table with 1 process
          WriteSeg(file, totyield, totyield_err, do_unc_, 0);
        } else if(do_eff_){
          WriteSeg(file, eff, eff_relUnc, do_unc_, 0);
        } else {
          WriteSeg(file, totyield, totyield_err, do_unc_, 1);
        }
      }//multi signal tables
    }//signals



    if(datas_.size() >= 1){
      bool pad = false;
      totyield = GetYield(datas_, irow);
      if(datas_.size() > 1) pad = 1;
      for(size_t i = 0; i < datas_.size(); ++i){
        WriteSeg(file, datas_.at(i)->sumw_.at(irow), totyield_err, 0, pad);
      }
      if(datas_.size() > 1) WriteSeg(file, totyield, totyield_err, 0, pad);
    }
    
  }else{
    file << "    \\multicolumn{" << NumColumns() << "}{c}{" << row.label_ << "}";
  }

  file << "\\\\\n";

  if(row.lines_after_ > 0){
    file << "    ";
    for(size_t i = 0; i < row.lines_after_; ++i){
      file << "\\hline";
    }
    file << "\n";
  }
  if(print_pie_) PrintPie(irow, luminosity);
} // PrintRow

void Table::WriteSeg(ofstream &file, double val, double valUnc = 0.0, bool unc = false, bool pad = false) const{
  if(pad == 0){
    if(unc == 1){
      file << " & " << val << " $\\pm$ " << valUnc;
    } else{
      file << " & " << val;
    }
  } else{
    if(unc == 1){
      file << " & " << setw(10) << val << " $\\pm$ " << valUnc;
    } else{
      file << " & " << setw(10) << val;
    }
  }
}

void Table::PrintPie(std::size_t irow, double luminosity) const{
  size_t Nbkg = backgrounds_.size();
  float Yield_tt = 0, Yield_tot = luminosity*GetYield(backgrounds_, irow);
  vector<double> counts(Nbkg);
  vector<int> colors(Nbkg);
  vector<const char*> labels(Nbkg);
  vector<TH1D> histos(Nbkg, TH1D("","",1,-1.,1.));
  TLegend leg(0., 0., 1., 1.); leg.SetFillStyle(0); leg.SetBorderSize(0);
  bool print_ttbar = false;
  for(size_t ind = 0; ind < Nbkg; ++ind){
    counts[ind] = luminosity*backgrounds_.at(ind)->sumw_.at(irow);
    histos[ind].SetFillColor(backgrounds_.at(ind)->process_->GetFillColor());
    colors.at(ind) = backgrounds_.at(ind)->process_->GetFillColor();
    string label = backgrounds_.at(ind)->process_->name_;
    leg.AddEntry(&histos[ind], label.c_str(), "f");
    labels.at(ind) =  label.c_str();
    if(Contains(label,"t#bar{t}") && !Contains(label,"t#bar{t}V")){
      Yield_tt += counts[ind];
      print_ttbar = true;
    }
  } // Loop over backgrounds

  // For now, use only the first PlotOpt in the vector
  gStyle->SetTitleW(0.95);
  TCanvas can("", "", plot_options_[0].CanvasWidth(), plot_options_[0].CanvasHeight());
  can.SetFillColorAlpha(0, 0.);
  can.SetFillStyle(4000);
  string plot_name;

  // Printing legend
  if(irow==0){
    leg.Draw();
    plot_name = "plots/pie_"+name_+"_legend_lumi"+RoundNumber(luminosity,0)+".pdf";
    if (Contains(name_,"ShortName:")) {
      string tagName=name_;
      ReplaceAll(tagName, "ShortName:", "");
      plot_name = "plots/pie_"+tagName+"_legend_lumi"+RoundNumber(luminosity,0)+".pdf";
    } else if (Contains(name_,"FixName:")) {
      string tagName=name_;
      ReplaceAll(tagName, "FixName:", "");
      plot_name = "plots/"+tagName+"_legend.pdf";
    }
    can.SaveAs(plot_name.c_str());
    cout<<" open "<<plot_name<<endl;
  }

  // Define piechart
  TPie pie("", "", Nbkg, &counts.at(0), &colors.at(0), &labels.at(0));
  pie.SetCircle(0.5, 0.48, 0.35);
  if(print_titlepie_){
    TString title = CodeToRootTex(rows_.at(irow).cut_.Name())+" (N="+RoundNumber(Yield_tot,1);
    if(print_ttbar) title += ", t#bar{t}="+RoundNumber(Yield_tt*100,1,Yield_tot)+"%";
    title += ")";
    pie.SetTitle(title);
  }

  // Printing pie chart with percentages
  pie.SetLabelFormat("%perc");
  pie.Draw();
  TLatex total(0.68,0.5,RoundNumber(Yield_tot,1));
  if(print_titlepie_) total.Draw();
  plot_name = "plots/pie_"+name_+"_"+CodeToPlainText(rows_.at(irow).cut_.Name())+"_perc_lumi"+RoundNumber(luminosity,0).Data()+".pdf";
  if (Contains(name_,"ShortName:")) {
    string tagName=name_;
    ReplaceAll(tagName, "ShortName:", "");
    plot_name = "plots/pie_"+tagName+"_perc_lumi"+RoundNumber(luminosity,0).Data()+".pdf";
  } else if (Contains(name_,"FixName:")) {
    string tagName=name_;
    ReplaceAll(tagName, "FixName:", "");
    plot_name = "plots/"+tagName+".pdf";
  }
  can.SaveAs(plot_name.c_str());
  cout<<" open "<<plot_name<<endl;

} // PrintPie


void Table::PrintFooter(ofstream &file) const{
  file << "    \\hline\n";
  file << "    ";

  if(backgrounds_.size() > 1){
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      file << " & " << ToLatex(backgrounds_.at(i)->process_->name_);
    }
    file << " & SM Bkg.";
  }else if(backgrounds_.size() == 1){
    file << " & " << ToLatex(backgrounds_.front()->process_->name_);
  }

  if(datas_.size() > 1){ 
    file << " & ";
    for(size_t i = 0; i < datas_.size(); ++i){
      file << " & " << ToLatex(datas_.at(i)->process_->name_);
    }
    file << " & Data Tot.";
  }else if(datas_.size() == 1){
    file << " & " << ToLatex(datas_.front()->process_->name_);
  }

  for(size_t i = 0; i < signals_.size(); ++i){
    file << " & " << ToLatex(signals_.at(i)->process_->name_);
    if(do_zbi_){
      file << " & $Z_{\\text{Bi}}$";
      file << " & $\\text{Cowan}$"; // ***NEW LINE***
    }
  }

  file << "\\\\\n";
  file << "    \\hline\\hline\n";
  file << "  \\end{tabular}\n";
  file << "\\end{preview}\n";
  file << "\\end{document}\n";
}

size_t Table::NumColumns() const{
  return 1
    + (backgrounds_.size() <= 1 ? backgrounds_.size() : backgrounds_.size()+1)
    + (datas_.size() <= 1 ? datas_.size() : datas_.size()+1)
    + (do_zbi_ ? 3 : 1)*signals_.size(); // ***EDITED LINE CHANGED 2-->3***
}

double Table::GetYield(const vector<unique_ptr<TableColumn> > &columns,
    size_t irow){
  double yield = 0.;
  for(const auto &column: columns){
    yield += column->sumw_.at(irow);
  }
  return yield;
}

double Table::GetError(const vector<unique_ptr<TableColumn> > &columns,
    size_t irow){
  double error = 0.;
  for(const auto &column: columns){
    error += column->sumw2_.at(irow);
  }
  return sqrt(error);
}

void Table::SetLuminosityTag(const string &tag) {
  this->LuminosityTag(tag);
}
