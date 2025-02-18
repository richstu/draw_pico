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

if __name__ == '__main__':
  get_sfs('plots/validation/zgvalidate_kingscanyonv1_2018__ll_pt0__Ztoee_CR__weightxw_years__lumi_lin.root','bkg_Z+FakePhoton_108','dat_Data','b.ll_pt()->at(0)')

