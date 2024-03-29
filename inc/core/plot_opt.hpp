#ifndef H_PLOT_OPT
#define H_PLOT_OPT

#include <cstddef>

#include <set>
#include <string>

namespace PlotOptTypes{
  enum class BottomType{off, ratio, diff, sorb, sorb_cut_upper};
  enum class YAxisType{linear, log};
  enum class TitleType{info, preliminary, simulation, simulation_preliminary, simulation_supplementary, supplementary, data};
  enum class StackType{signal_overlay, signal_on_top, data_norm, lumi_shapes, shapes};
  enum class OverflowType{none, underflow, overflow, both};
}

class PlotOpt{
public:
  PlotOpt();
  PlotOpt(const std::string &file_name,
          const std::string &config_name);
  PlotOpt(const PlotOpt &) = default;
  PlotOpt& operator=(const PlotOpt &) = default;
  PlotOpt(PlotOpt &&) = default;
  PlotOpt& operator=(PlotOpt &&) = default;
  ~PlotOpt() = default;

  PlotOpt operator()() const;

  PlotOpt & LoadOptions(const std::string &file_name,
                        const std::string &config_name);

  PlotOpt & Bottom(PlotOptTypes::BottomType bottom_type);
  PlotOptTypes::BottomType Bottom() const;

  PlotOpt & YAxis(PlotOptTypes::YAxisType y_axis_type);
  PlotOptTypes::YAxisType YAxis() const;

  PlotOpt & XAxis(PlotOptTypes::YAxisType x_axis_type);
  PlotOptTypes::YAxisType XAxis() const;

  PlotOpt & Title(PlotOptTypes::TitleType title_type);
  PlotOptTypes::TitleType Title() const;

  PlotOpt & Stack(PlotOptTypes::StackType stack_type);
  PlotOptTypes::StackType Stack() const;

  PlotOpt & Overflow(PlotOptTypes::OverflowType overflow_type);
  PlotOptTypes::OverflowType Overflow() const;

  PlotOpt & FileExtensions(const std::set<std::string> &file_extensions);
  const std::set<std::string> & FileExtensions() const;

  PlotOpt & LabelSize(double label_size);
  double LabelSize() const;

  PlotOpt & TitleSize(double title_size);
  double TitleSize() const;

  PlotOpt & ExtraLabelSize(double extra_label_size);
  double ExtraLabelSize() const;

  PlotOpt & XTitleOffset(double x_title_offset);
  double XTitleOffset() const;

  PlotOpt & YTitleOffset(double y_title_offset);
  double YTitleOffset() const;

  PlotOpt & ZTitleOffset(double z_title_offset);
  double ZTitleOffset() const;

  PlotOpt & AutoYAxis(bool auto_y_axis);
  bool AutoYAxis() const;

  PlotOpt & ErrorOnZeroData(bool error_on_zero_data);
  bool ErrorOnZeroData() const;

  PlotOpt & TitleInFrame(bool title_in_frame);
  bool TitleInFrame() const;

  PlotOpt & CanvasSize(int width, int height);
  PlotOpt & CanvasWidth(int width);
  int CanvasWidth() const;
  PlotOpt & CanvasHeight(int height);
  int CanvasHeight() const;

  PlotOpt & Margin(double left, double right, double bottom, double top);
  PlotOpt & LeftMargin(double left);
  double LeftMargin() const;
  PlotOpt & RightMargin(double right);
  double RightMargin() const;
  PlotOpt & BottomMargin(double bottom);
  double BottomMargin() const;
  PlotOpt & TopMargin(double top);
  double TopMargin() const;

  PlotOpt & BottomHeight(double bottom_height);
  double BottomHeight() const;

  PlotOpt & LegendColumns(int columns);
  int LegendColumns() const;
  PlotOpt & LegendEntryHeight(double height);
  double LegendEntryHeight() const;
  PlotOpt & LegendMaxHeight(double height);
  double LegendMaxHeight() const;
  PlotOpt & LegendMarkerWidth(double width);
  double LegendMarkerWidth() const;
  PlotOpt & LegendPad(double pad);
  double LegendPad() const;
  PlotOpt & LegendLeftPad(double left_pad);
  double LegendLeftPad() const;
  PlotOpt & LegendLeftColumnOffset(double left_column_offset);
  double LegendLeftColumnOffset() const;
  PlotOpt & LegendDensity(double density);
  double LegendDensity() const;

  PlotOpt & LogMinimum(double log_minimum);
  double LogMinimum() const;

  PlotOpt & RatioMinimum(double ratio_minimum);
  double RatioMinimum() const;
  PlotOpt & RatioMaximum(double ratio_maximum);
  double RatioMaximum() const;

  PlotOpt & NDivisions(int n_divisions);
  int NDivisions() const;
  PlotOpt & NDivisionsBottom(int n_divisions);
  int NDivisionsBottom() const;

  PlotOpt & Font(int font);
  int Font() const;

  PlotOpt & ShowBackgroundError(bool show_background_error);
  bool ShowBackgroundError() const;

  PlotOpt & UseCMYK(bool use_cmyk);
  bool UseCMYK() const;

  PlotOpt & PrintVals(bool print_vals);
  bool PrintVals() const;

  double TopToGlobalYNDC(double top_y) const;
  double GlobalToTopYNDC(double global_y) const;
  double BottomToGlobalYNDC(double bottom_y) const;
  double GlobalToBottomYNDC(double global_y) const;

  double TrueLegendHeight(std::size_t num_entries) const;
  double TrueLegendEntryHeight(std::size_t num_entries) const;
  double TrueLegendWidth(std::size_t num_entries) const;

  bool BackgroundsStacked() const;
  bool DisplayLumiEntry() const;

  std::string TypeString() const;

  void MakeSane();

private:
  PlotOptTypes::BottomType bottom_type_;
  PlotOptTypes::YAxisType y_axis_type_;
  PlotOptTypes::YAxisType x_axis_type_;
  PlotOptTypes::TitleType title_type_;
  PlotOptTypes::StackType stack_type_;
  PlotOptTypes::OverflowType overflow_type_;
  std::set<std::string> file_extensions_;
  double title_size_, extra_label_size_, label_size_;
  double x_title_offset_, y_title_offset_, z_title_offset_;
  bool auto_y_axis_;
  bool error_on_zero_data_;
  int canvas_width_, canvas_height_;
  double left_margin_, right_margin_, bottom_margin_, top_margin_;
  double bottom_height_;
  int legend_columns_;
  double legend_entry_height_, legend_max_height_;
  double legend_marker_width_, legend_pad_, legend_density_;
  double legend_left_pad_, legend_left_column_offset_;
  double log_minimum_;
  double ratio_minimum_, ratio_maximum_;
  int n_divisions_, n_divisions_bottom_;
  int font_;
  bool show_background_error_;
  bool use_cmyk_;
  bool print_vals_;
  bool title_in_frame_;

  void SetProperty(const std::string &property_name,
                   const std::string &value_string);
};

#endif
