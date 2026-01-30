#/usr/bin/env python
'''
script for generating scale factor NamedFuncs from draw_pico plots
'''

from array import array
from math import sqrt
from subprocess import check_call
import ROOT

ROOT.gInterpreter.Declare("""

TH1D* cast_tobject_to_th1d(TObject* object)
{ return static_cast<TH1D*>(object); }

""")

def get_sfs(filename, den_histname, num_histname, var_name):
  hist_file = ROOT.TFile(filename,'read')
  canvas = hist_file.canvas
  den_hist = ROOT.cast_tobject_to_th1d(canvas.FindObject(den_histname)).Clone()
  num_hist = ROOT.cast_tobject_to_th1d(canvas.FindObject(num_histname)).Clone()
  den_hist.Scale(1.0/den_hist.Integral())
  num_hist.Scale(1.0/num_hist.Integral())
  num_hist.Divide(den_hist)
  axis = num_hist.GetXaxis()
  for ibin in range(axis.GetNbins()):
    if (ibin != 0):
      print('else '),
    print('if ('+var_name+'>'+str(axis.GetBinLowEdge(ibin+1))+' && '+var_name+'<'+str(axis.GetBinUpEdge(ibin+1))+') '),
    print(' return '+str(num_hist.GetBinContent(ibin+1))+';')
  hist_file.Close()

def get_sfs_twodim(num_filename, num_histname, den_filename, den_histname, 
                   xvar_name, yvar_name):
  num_hist_file = ROOT.TFile(num_filename,'read')
  den_hist_file = ROOT.TFile(den_filename,'read')
  num_canvas = num_hist_file.canvas
  den_canvas = den_hist_file.canvas
  den_hist = ROOT.cast_tobject_to_th1d(
      den_canvas.FindObject(den_histname)).Clone()
  num_hist = ROOT.cast_tobject_to_th1d(
      num_canvas.FindObject(num_histname)).Clone()
  num_hist.SetDirectory(ROOT.nullptr)
  den_hist.SetDirectory(ROOT.nullptr)
  num_hist.Scale(1.0/num_hist.Integral())
  den_hist.Scale(1.0/den_hist.Integral())
  xaxis = num_hist.GetXaxis()
  yaxis = num_hist.GetYaxis()
  first = True
  for ibinx in range(xaxis.GetNbins()):
    for ibiny in range(yaxis.GetNbins()):
      xlo = xaxis.GetBinLowEdge(ibinx+1)
      xhi = xaxis.GetBinUpEdge(ibinx+1)
      ylo = yaxis.GetBinLowEdge(ibiny+1)
      yhi = yaxis.GetBinUpEdge(ibiny+1)
      num = num_hist.GetBinContent(ibinx+1, ibiny+1)
      den = den_hist.GetBinContent(ibinx+1, ibiny+1)
      ratio = 1
      if den > 0:
        ratio = num/den
      if first:
        first = False
      else:
        print('else ',end='')
      print(f'if ({xvar_name}>={xlo} && {xvar_name}<{xhi} && {yvar_name}>={ylo}'
           +f'&& {yvar_name}<{yhi}) return {ratio:.3f};')

if __name__ == '__main__':
  #get_sfs('plots/validation/zgvalidate_kingscanyonv1_2018__ll_pt0__Ztoee_CR__weightxw_years__lumi_lin.root','bkg_Z+FakePhoton_108','dat_Data','b.ll_pt()->at(0)')
  get_sfs_twodim('plots/zgcurrent_photon_corr_data__lead_photon_abseta__photon_pt0__ZtoMuMuGammaCR__weightxw_years__lumi_nonorm_lin.root', '', 'plots/zgcurrent_photon_corr_simu__lead_photon_abseta__photon_pt0__ZtoMuMuGammaCR__weightxw_years__lumi_nonorm_lin.root', '', 'photon_pt', 'photon_abseta')

