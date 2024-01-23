#/usr/bin/env python

from array import array
from math import sqrt
from subprocess import check_call
import ROOT

ROOT.gInterpreter.Declare("""

RooHist* cast_tobject_to_roohist(TObject* object)
{ return static_cast<RooHist*>(object); }

RooCurve* cast_tobject_to_roocurve(TObject* object)
{ return static_cast<RooCurve*>(object); }

""")

def get_value_tgraph(graph, x):
  '''Estimates the value of a TGraph at a point x via linear interpolation 
  between the nearest points

  graph  TGraph or similar (ex. RooCurve)
  x      point to get value at
  '''
  npoints = graph.GetN();
  xi = ROOT.Double(0.0)
  yi = ROOT.Double(0.0)
  graph.GetPoint(0,xi,yi)
  if (x < xi):
    raise ValueError('Requested x '+str(x)+' is less than graph minimum '+str(xi))
  graph.GetPoint(npoints-1,xi,yi)
  if (x > xi):
    raise ValueError('Requested x '+str(x)+' is greater than graph maximum '+str(xi))
  if (x == xi):
    return (yi+0.0)
  idx = 0
  graph.GetPoint(idx,xi,yi)
  while (xi <= x):
    graph.GetPoint(idx,xi,yi)
    idx += 1
  xi_prev = ROOT.Double(0.0)
  yi_prev = ROOT.Double(0.0)
  graph.GetPoint(idx-2,xi_prev,yi_prev)
  return (x-xi_prev)/(xi-xi_prev)*(yi-yi_prev)+yi_prev

def get_tgraph_maximum(graph):
  npoints = graph.GetN();
  xi = ROOT.Double(0.0)
  yi = ROOT.Double(0.0)
  maximum = 0.0
  for i in range(npoints):
    graph.GetPoint(i,xi,yi)
    if yi > maximum:
      maximum = yi+0.0
  return maximum

def get_stddev(lst):
  '''Get std deviation from Python list'''
  lst_sum = 0
  lst_sumsq = 0
  for i in lst:
    lst_sum += i
    lst_sumsq += i*i
  lst_len = float(len(lst))
  return sqrt(lst_sumsq/lst_len-(lst_sum/lst_len)**2)

if __name__ == '__main__':
  check_call('combine -M MultiDimFit datacards/test_datacard.txt --saveWorkspace -n .bestfit'.split())
  
  f = ROOT.TFile('higgsCombine.bestfit.MultiDimFit.mH120.root')
  w = f.Get('w')
  w.Print('v')
  
  n_bins = 55
  lower_bound = 105
  upper_bound = 160
  binning = ROOT.RooFit.Binning(n_bins,lower_bound,upper_bound)

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetErrorX(0)
  
  for channel in [('cat_ggh4',0),('cat_ggh3',1),('cat_ggh2',2),('cat_ggh1',3),
                  ('cat_vbf3',4),('cat_vbf2',5),('cat_vbf1',6),('cat_vh3l',7),
                  ('cat_vhmet',8),('cat_tthhad',9),('cat_tthlep',10)]:
    #Convert esoteric RooFit formats into regular ROOT
    can_rf = ROOT.TCanvas()
    plot = w.var('mllg_'+channel[0]).frame()
    w.data('data_obs').reduce('CMS_channel=='+str(channel[1])).plotOn( plot, binning )
    sb_model = w.pdf('model_s').getPdf(channel[0])
    w.loadSnapshot('MultiDimFit')
    sb_model.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name('postfit_sb') )
    r_bestfit = w.var('r').getVal()
    w.var('r').setVal(0)
    sb_model.plotOn( plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name('postfit_b') )
    plot.Draw()
    can_rf.Update()
    can_rf.SaveAs('plots/zgamma_roofit_'+channel[0]+'.pdf')
    data_obs_roohist = ROOT.cast_tobject_to_roohist(can_rf.FindObject('h_data_obs'))
    postfitsb_roocurve = ROOT.cast_tobject_to_roocurve(can_rf.FindObject('postfit_sb'))
    postfitb_roocurve = ROOT.cast_tobject_to_roocurve(can_rf.FindObject('postfit_b'))

    #generate histograms and curves
    difference_hist = ROOT.TH1D('','',n_bins,lower_bound,upper_bound)
    data_hist = ROOT.TH1D('','',n_bins,lower_bound,upper_bound)
    diffs = []
    for i in range(n_bins):
      mllg = i*1.0+(lower_bound+0.5)
      obs_yield = get_value_tgraph(data_obs_roohist,mllg)
      if (obs_yield < 0):
        obs_yield = 0
      diff = obs_yield-get_value_tgraph(postfitb_roocurve,mllg)
      diffs.append(diff)
      difference_hist.SetBinContent(i+1,diff)
      difference_hist.SetBinError(i+1,sqrt(obs_yield))
      data_hist.SetBinContent(i+1,obs_yield)
      data_hist.SetBinError(i+1,sqrt(obs_yield))

    difference_x = array('d')
    difference_b = array('d')
    difference_sb = array('d')
    for i in range(postfitsb_roocurve.GetN()):
      xi = ROOT.Double(0.0)
      yi = ROOT.Double(0.0)
      postfitsb_roocurve.GetPoint(i,xi,yi)
      difference_x.append(xi+0.0)
      difference_b.append(0.0)
      difference_sb.append(yi-get_value_tgraph(postfitb_roocurve,xi+0.0))
    difference_b_plot = ROOT.TGraph(len(difference_x),difference_x,difference_b)
    difference_sb_plot = ROOT.TGraph(len(difference_x),difference_x,difference_sb)
  
    #make things nice
    data_hist.SetMarkerStyle(ROOT.kFullCircle)
    data_hist.SetMarkerColor(ROOT.kBlack)
    data_hist.SetLineColor(ROOT.kBlack)
    data_hist.SetLineWidth(2)
    postfitsb_roocurve.SetLineColor(ROOT.kRed)
    postfitb_roocurve.SetLineColor(ROOT.kRed)
    postfitb_roocurve.SetLineStyle(ROOT.kDashed)
    postfitsb_roocurve.SetLineWidth(3)
    postfitb_roocurve.SetLineWidth(3)
    difference_hist.GetYaxis().SetTitle('Data-Fit')
    difference_hist.GetYaxis().SetNdivisions(606)
    difference_hist.GetXaxis().SetTitle('m_{ll#gamma} [GeV]')
    difference_hist.SetLabelSize(0.035,'y')
    difference_hist.SetMinimum(-3.0*get_stddev(diffs))
    difference_hist.SetMaximum(3.0*get_stddev(diffs))
    difference_hist.SetLineColor(ROOT.kBlack)
    difference_hist.SetLineWidth(2)
    difference_hist.SetMarkerColor(ROOT.kBlack)
    difference_hist.SetMarkerStyle(ROOT.kFullCircle)
    difference_b_plot.SetLineColor(ROOT.kRed)
    difference_b_plot.SetLineWidth(3)
    difference_b_plot.SetLineStyle(ROOT.kDashed)
    difference_sb_plot.SetLineColor(ROOT.kRed)
    difference_sb_plot.SetLineWidth(3)
    can = ROOT.TCanvas()
    #hist formatting
    #top plot settings
    data_hist.SetLabelSize(0,'x')
    data_hist.SetLabelSize(0.035,'y')
    data_hist.GetYaxis().SetTitle('Events/(1 GeV)')
    data_hist.GetYaxis().SetNdivisions(606)
    data_hist.SetMinimum(0.0)
    data_hist.SetMaximum(get_tgraph_maximum(data_obs_roohist)*1.2)
    #pad stuff
    top_pad = ROOT.TPad('top_pad_'+channel[0],'',0.0,0.0,1.0,1.0)
    bot_pad = ROOT.TPad('bot_pad_'+channel[0],'',0.0,0.0,1.0,1.0)
    top_pad.SetTicks(1,1)
    bot_pad.SetTicks(1,1)
    top_pad.SetMargin(0.19,0.06,0.3,0.05)
    bot_pad.SetMargin(0.19,0.06,0.1,0.7)
    top_pad.SetFillStyle(4000)
    bot_pad.SetFillStyle(4000)
    top_pad.Draw()
    top_pad.cd()
    #draw stuff here
    data_hist.Draw()
    postfitb_roocurve.Draw('same')
    postfitsb_roocurve.Draw('same')
    data_hist.Draw('same')
    leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
    leg.AddEntry('postfit_sb', 'Postfit S+B model (#mu=%.2f)'%r_bestfit, 'L')
    leg.AddEntry('postfit_b', 'Postfit B model (#mu=0)', 'L')
    leg.SetBorderSize(0)
    leg.Draw('same')
    label = ROOT.TLatex()
    label.SetTextSize(0.035)
    label.SetNDC(ROOT.kTRUE)
    label.SetTextAlign(11)
    label.DrawLatex(0.20,0.96,"#font[62]{CMS} #scale[0.8]{#font[52]{Work in Progress}}")
    label.SetTextAlign(31)
    label.SetTextSize(0.03)
    label.DrawLatex(0.93,0.96,"#font[42]{60 fb^{-1} (13 TeV)}")
    top_pad.Modified()
    can.cd()
    bot_pad.Draw('same')
    bot_pad.cd()
    difference_hist.Draw()
    difference_b_plot.Draw('same')
    difference_sb_plot.Draw('same')
    difference_hist.Draw('same')
    bot_pad.Modified()
    can.Draw()
    can.SaveAs('plots/zgamma_fit_'+channel[0]+'.pdf')

  f.Close()
  #check_call('rm higgsCombine.bestfit.MultiDimFit.mH120.root'.split())

