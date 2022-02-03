#! /usr/bin/env python
import os
import math
from ROOT import *
import sys
import gc

def read_datacard_results(bkg_pre, ebkg_pre_up, ebkg_pre_dn, bkg_post, ebkg_post_up, ebkg_post_dn, data, signal_500_raw):
  #fit_bin_ordering = [50,52,54,56,58,60,62,64,49,51,53,55,57,59,61,63]
  fit_bin_ordering = [51,53,55,57,59,61,63,65,50,52,54,56,58,60,62,64]
  #data_bin_ordering = [11,23,35,47,12,24,36,48,5,17,29,41,6,18,30,42]
  data_bin_ordering = [24,36,48,60,25,37,49,61,18,30,42,54,19,31,43,55,4,5,6,1,2,3]
  input_dir = 'input'
  prefit_filename = 'multidimfit_prefit.root'
  postfit_filename = 'multidimfit_postfit.root'
  boosted_filename = 'boosted_results.root'
  datacard_filename = '1DTChiHH450_LSP1_Data_Combo.txt'
  prefit_file = TFile.Open(input_dir+'/'+prefit_filename,'READ')
  postfit_file = TFile.Open(input_dir+'/'+postfit_filename,'READ')
  boosted_file = TFile.Open(input_dir+'/'+boosted_filename,'READ')
  prefit_args = []
  postfit_args = []
  data_args = []
  prefit_result = prefit_file.fit_mdf.floatParsFinal()
  for bin_idx in fit_bin_ordering:
    prefit_args.append(prefit_result.at(bin_idx))
    bkg_pre.append(prefit_result.at(bin_idx).getValV())
    ebkg_pre_up.append(prefit_result.at(bin_idx).getErrorHi())
    ebkg_pre_dn.append(-1.0*prefit_result.at(bin_idx).getErrorLo())
  postfit_result = postfit_file.fit_mdf.floatParsFinal()
  for bin_idx in fit_bin_ordering:
    postfit_args.append(postfit_result.at(bin_idx))
    bkg_post.append(postfit_result.at(bin_idx).getValV())
    ebkg_post_up.append(postfit_result.at(bin_idx).getErrorHi())
    ebkg_post_dn.append(-1.0*postfit_result.at(bin_idx).getErrorLo())
  for bin_idx in [0,1,2,3,4,5]:
    prefit_args.append(boosted_file.prefit.at(bin_idx))
    bkg_pre.append(boosted_file.prefit.at(bin_idx).getValV())
    ebkg_pre_up.append(boosted_file.prefit.at(bin_idx).getErrorHi())
    ebkg_pre_dn.append(-1.0*boosted_file.prefit.at(bin_idx).getErrorLo())
    postfit_args.append(boosted_file.postfit.at(bin_idx))
    bkg_post.append(boosted_file.postfit.at(bin_idx).getValV())
    ebkg_post_up.append(boosted_file.postfit.at(bin_idx).getErrorHi())
    ebkg_post_dn.append(-1.0*boosted_file.postfit.at(bin_idx).getErrorLo())
  datacard = open(input_dir+'/'+datacard_filename,'r')
  datacard_lines = datacard.read().split('\n')
  datacard_binnames = datacard_lines[67].split()
  datacard_yields = datacard_lines[68].split()
  for bin_idx in data_bin_ordering:
    data_args.append(RooRealVar(datacard_binnames[bin_idx],datacard_binnames[bin_idx],float(datacard_yields[bin_idx])))
    data.append(float(datacard_yields[bin_idx]))
  datacard.close()
  #prefit_values = RooArgList('prefit_values')
  #postfit_values = RooArgList('postfit_values')
  #data_values = RooArgList('observed_values')
  #for prefit_arg in prefit_args:
  #  prefit_values.add(prefit_arg)
  #for postfit_arg in postfit_args:
  #  postfit_values.add(postfit_arg)
  #for data_arg in data_args:
  #  data_values.add(data_arg)
  #output_file = TFile.Open('SUS-20-007_fitresults.root','RECREATE')
  #prefit_values.Write()
  #postfit_values.Write()
  #data_values.Write()
  #output_file.Close()
  prefit_file.Close()
  postfit_file.Close()
  boosted_file.Close()
  #try to fix naming issues with manual garbage collection
  del prefit_file
  del postfit_file
  del boosted_file
  gc.collect()

plot_pulls = False
plot_sig = False
preliminary = False
bin_labels = False

for arg in sys.argv:
  if (arg == '--preliminary'):
    preliminary = True
  if (arg == '--pulls'):
    pulls = True
  if (arg == '--signal'):
    plot_sig = True
  if (arg == '--binlabels'):
    bin_labels = True

#tag1_lbl = 'Pre-fit'
tag1_lbl = 'Pred'
tag1_color = kPink+2

#tag2_lbl = 'Post-fit'
tag2_lbl = 'Fit'
tag2_color = kAzure+1
tag2_color_2 = kAzure+2

bkg_pre = []
ebkg_pre_up = []
ebkg_pre_dn = []
bkg_post = []
ebkg_post_up = []
ebkg_post_dn = []
data = []
signal_500_raw = []
read_datacard_results(bkg_pre, ebkg_pre_up, ebkg_pre_dn, bkg_post, ebkg_post_up, ebkg_post_dn, data, signal_500_raw)

#14 - D met0 low drmax      20 - D met0 hidrmax
#15 - B nb3 met0 low drmax  21 - B nb3 met0 hidrmax
#16 - B nb4 met0 low drmax  22 - B nb4 met0 hidrmax
#17 - C met0 low drmax      23 - C met0 hidrmax
#18 - A nb3 met0 low drmax  24 - A nb3 met0 hidrmax
#19 - A nb4 met0 low drmax  25 - A nb4 met0 hidrmax

#signal yields by hand for now? from datacard?
#signal_500_raw = [1.36, 5.06, 6.49, 5.17, 2.53, 9.76, 12.12, 10.76, 1.40, 4.18, 3.77, 2.43, 3.45, 10.18, 8.94, 6.15] + [3.21, 0.61, 0.13, 9.76, 1.46, 0.24]
#despite name is 450,1 TChiHH-G signal
signal_500_raw = [2.61, 9.63, 9.34, 4.25, 4.92, 18.40, 19.14, 8.45, 2.07, 6.11, 4.18, 1.35, 5.38, 15.11, 11.13, 3.12] + [3.17, 0.45, 0.09, 9.14, 1.06, 0.17]

print(ebkg_pre_up)
print(ebkg_pre_dn)
print(ebkg_post_up)
print(ebkg_post_dn)

signal_500 = []
edata_up = []
edata_dn = []
esignal_up = []
esignal_dn = []
pull_pre = []
pull_sig = []
pull_post = []
#pulls as likelihood ratio significances
#pull_pre = [-0.9, 0.0, 0.6, 0.1,   0.0, 1.1, 0.6, -1.9,   0.9, 0.0, 3.5, -1.1,   -1.0, 1.5, 0.6, 0.3,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
for ibin in range(len(data)):
  #calculate expected signal+bkg
  signal_500.append(signal_500_raw[ibin] + bkg_pre[ibin])
  #make Poisson errors
  alpha = 1.0 - 0.6827
  l = 0
  if (data[ibin] != 0):
    l = ROOT.Math.gamma_quantile(alpha/2,data[ibin],1.0)
  u =  ROOT.Math.gamma_quantile_c(alpha/2,data[ibin]+1,1.0)
  edata_up.append(u-data[ibin])
  edata_dn.append(data[ibin]-l)
  #
  if (data[ibin] != 0):
    l = ROOT.Math.gamma_quantile(alpha/2,signal_500[ibin],1.0)
  u =  ROOT.Math.gamma_quantile_c(alpha/2,signal_500[ibin]+1,1.0)
  esignal_up.append(u-signal_500[ibin])
  esignal_dn.append(signal_500[ibin]-l)
  #make pulls
  w_pre = 1.0/ebkg_pre_up[ibin]**2
  w_pos = 1.0/ebkg_post_up[ibin]**2
  e_dat = edata_dn[ibin]
  if data[ibin] < bkg_pre[ibin]:
    e_dat = edata_up[ibin]
    w_pre = 1.0/ebkg_pre_dn[ibin]**2
    w_pos = 1.0/ebkg_post_dn[ibin]**2
  w_dat = 1.0/e_dat**2
  e_sig = esignal_dn[ibin]
  w_sig = 1.0/e_sig**2
  w_pre_sig = 1.0/ebkg_pre_up[ibin]**2
  e_com_pre = 1.0/math.sqrt(w_pre+w_dat)
  e_com_pos = 1.0/math.sqrt(w_pos+w_dat)
  e_com_sig = 1.0/math.sqrt(w_pre_sig+w_sig)
  l_com_pre = (bkg_pre[ibin]*w_pre+data[ibin]*w_dat)/(w_pre+w_dat)
  l_com_pos = (bkg_post[ibin]*w_pos+data[ibin]*w_dat)/(w_pos+w_dat)
  l_com_sig = (bkg_pre[ibin]*w_pre_sig+signal_500[ibin]*w_sig)/(w_pre_sig+w_sig)
  this_pull_pre = (data[ibin]-l_com_pre)/math.sqrt(e_dat**2-e_com_pre**2)
  this_pull_pos = (data[ibin]-l_com_pos)/math.sqrt(e_dat**2-e_com_pos**2)
  this_pull_sig = (signal_500[ibin]-l_com_sig)/math.sqrt(e_sig**2-e_com_sig**2)
  pull_pre.append(this_pull_pre)
  pull_post.append(this_pull_pos)
  pull_sig.append(this_pull_sig)

  #naive-pulls
  #pull_pre.append((data[ibin]-bkg_pre[ibin])/math.sqrt(bkg_pre[ibin]+ebkg_pre[ibin]*ebkg_pre[ibin]))
  #pull_post.append((data[ibin]-bkg_post[ibin])/math.sqrt(bkg_post[ibin]+ebkg_post[ibin]*ebkg_post[ibin]))
  #poisson significance
  #if (ibin > 15):
  #  pull_pre.append(math.sqrt(2)*TMath.ErfInverse(-1.0+2.0*ROOT.Math.poisson_cdf(int(data[ibin]),bkg_pre[ibin])))
  #  pull_post.append(math.sqrt(2)*TMath.ErfInverse(-1.0+2.0*ROOT.Math.poisson_cdf(int(data[ibin]),bkg_post[ibin])))

#print(pull_pre)
#print(pull_post)

nbins = len(data)
nhbins = len(data)

#         Plotting the data
#-------------------------------------
miny, maxy = 0.015, 7000
if bin_labels:
  maxy = 1000
htopdummy = TH1D("","",nhbins,0,nhbins+1)
grbkg_pre = TGraphAsymmErrors(nhbins)
grbkg_post = TGraphAsymmErrors(nhbins)
grbkg_horiz_pre = TGraphAsymmErrors(nhbins)
grbkg_horiz_post = TGraphAsymmErrors(nhbins)
grdata = TGraphAsymmErrors(nhbins)
grpull_pre = TGraph(nhbins)
grpull_post = TGraph(nhbins)
grsignal = TGraphAsymmErrors(nhbins)
grpull_sig = TGraph(nhbins)
histdata = TH1D('histdata','data',22,0.5,22.5)
i = 0
print(nbins)
for ibin in range(nbins):
  #if 'r4' not in bins[ibin]: continue
  i +=1
  grbkg_pre.SetPoint(i-1, i,bkg_pre[ibin])
  grbkg_pre.SetPointEYhigh(i-1, ebkg_pre_up[ibin])    
  grbkg_pre.SetPointEYlow(i-1, ebkg_pre_dn[ibin])    
  grbkg_pre.SetPointEXhigh(i-1, 0.5)
  grbkg_pre.SetPointEXlow(i-1, 0.5)
  grbkg_horiz_pre.SetPoint(i-1, i,bkg_pre[ibin])
  grbkg_horiz_pre.SetPointEXhigh(i-1, 0.5)
  grbkg_horiz_pre.SetPointEXlow(i-1, 0.5)
  grpull_pre.SetPoint(i-1, i, pull_pre[ibin])

  grbkg_post.SetPoint(i-1, i, bkg_post[ibin])
  grbkg_post.SetPointEYhigh(i-1, ebkg_post_up[ibin])
  grbkg_post.SetPointEYlow(i-1, ebkg_post_dn[ibin])
  grbkg_post.SetPointEXhigh(i-1, 0.5)
  grbkg_post.SetPointEXlow(i-1, 0.5)
  grbkg_horiz_post.SetPoint(i-1, i, bkg_post[ibin])
  grbkg_horiz_post.SetPointEXhigh(i-1, 0.5)
  grbkg_horiz_post.SetPointEXlow(i-1, 0.5)
  grpull_post.SetPoint(i-1, i, pull_post[ibin])

  #if ((not plot_pulls) and (data[ibin] == 0)):
  #  grdata.SetPoint(i-1, i, miny)
  #else:
  grdata.SetPoint(i-1, i, data[ibin])
  grdata.SetPointEYhigh(i-1, edata_up[ibin])
  grdata.SetPointEYlow(i-1, edata_dn[ibin])
  for idata in range(int(data[ibin])):
    histdata.Fill(i)
  grsignal.SetPoint(i-1, i, signal_500[ibin])
  grsignal.SetPointEYhigh(i-1, 0)
  grsignal.SetPointEYlow(i-1, 0)
  grsignal.SetPointEXhigh(i-1, 0.5)
  grsignal.SetPointEXlow(i-1, 0.5)
  grpull_sig.SetPoint(i-1, i, pull_sig[ibin])
  #grsignal.SetPointEYhigh(i-1, esignal_up[ibin])
  #grsignal.SetPointEYlow(i-1, esignal_dn[ibin])
  #grdata.SetPointEXhigh(i, 0.)
  #grdata.SetPointEXlow(i, 0.)

# for ibin in range(nbins):
#     print bkg_pre[ibin], ebkg_pre[ibin]

can = TCanvas('c','c',1000,500)
can.cd()

left_margin = 0.06
right_margin = 0.02

an_region_size = 0.05
nb_title_size = 0.045
met_title_size = 0.04
met_range_size = 0.03
tpad_bottom = 0.3
tpad_bottom_margin = 0.0
top_axis_label_size = 0.06
top_axis_title_size = 0.075
top_axis_title_offset = 0.5
if not plot_pulls:
  an_region_size = 0.041
  nb_title_size = 0.038
  met_title_size = 0.037
  met_range_size = 0.0325
  tpad_bottom = 0.0
  tpad_bottom_margin = 0.1
  if bin_labels:
    tpad_bottom_margin = 0.3
  top_axis_label_size = 0.04
  top_axis_title_size = 0.05
  top_axis_title_offset = 0.565

gStyle.SetOptStat(0)
top = TPad("top_pad", "top_pad", 0., tpad_bottom, 1., 1.)
can.SetMargin(0., 0., 0., 0.);
can.SetFillStyle(4000);

top.SetTopMargin(0.1)
top.SetBottomMargin(tpad_bottom_margin)
top.SetLeftMargin(left_margin)
top.SetRightMargin(right_margin)
top.SetFillStyle(4000);
top.Draw()

#     Top pad
# -----------------------------

top.cd()
top.SetLogy()
htopdummy.GetYaxis().SetLabelSize(top_axis_label_size)
htopdummy.GetYaxis().SetRangeUser(miny,maxy)
htopdummy.GetYaxis().SetTitle("Events / bin")
htopdummy.GetYaxis().CenterTitle()
htopdummy.GetYaxis().SetTitleSize(top_axis_title_size)
htopdummy.GetYaxis().SetTitleOffset(top_axis_title_offset)
#htopdummy.GetXaxis().SetRangeUser(0.01, nhbins+0.99)
htopdummy.GetXaxis().SetLimits(0.01, nhbins+0.99)
htopdummy.Draw()

ptmiss = "#font[52]{p}#font[42]{#lower[-0.1]{_{T}}#kern[-0.25]{#scale[1.0]{#lower[0.2]{^{miss}}}}}";

if not plot_pulls:
  htopdummy.GetXaxis().SetLabelSize(top_axis_label_size)
  #htopdummy.GetXaxis().SetLabelOffset(0.02)
  htopdummy.GetXaxis().SetNdivisions(24);
  if not bin_labels:
    htopdummy.GetXaxis().SetTitle("Bin number")
  htopdummy.GetXaxis().CenterTitle()
  htopdummy.GetXaxis().SetTitleSize(top_axis_title_size)
  htopdummy.GetXaxis().SetTickLength(0.02)
  htopdummy.GetXaxis().SetTitleOffset(0.85)
  htopdummy.GetXaxis().SetLimits(0.01,nhbins+0.99)
  
  if bin_labels:
    #htopdummy.GetXaxis().SetBinLabel(1,'150<'+ptmiss+'<200 GeV, 3b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(2,'200<'+ptmiss+'<300 GeV, 3b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(3,'300<'+ptmiss+'<400 GeV, 3b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(4,ptmiss+'>400 GeV, 3b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(5,'150<'+ptmiss+'<200 GeV, 4b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(6,'200<'+ptmiss+'<300 GeV, 4b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(7,'300<'+ptmiss+'<400 GeV, 4b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(8,ptmiss+'>400 GeV, 4b, Resolved High #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(9,'150<'+ptmiss+'<200 GeV, 3b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(10,'200<'+ptmiss+'<300 GeV, 3b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(11,'300<'+ptmiss+'<400 GeV, 3b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(12,ptmiss+'>400 GeV, 3b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(13,'150<'+ptmiss+'<200 GeV, 4b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(14,'200<'+ptmiss+'<300 GeV, 4b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(15,'300<'+ptmiss+'<400 GeV, 4b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(16,ptmiss+'>400 GeV, 4b, Resolved Low #DeltaR#lower[-0.1]{_{max}}')
    #htopdummy.GetXaxis().SetBinLabel(17,'300<'+ptmiss+'<500 GeV, 1H, Boosted')
    #htopdummy.GetXaxis().SetBinLabel(18,'500<'+ptmiss+'<700 GeV, 1H, Boosted')
    #htopdummy.GetXaxis().SetBinLabel(19,ptmiss+'>700 GeV, 1H, Boosted')
    #htopdummy.GetXaxis().SetBinLabel(20,'300<'+ptmiss+'<500 GeV, 2H, Boosted')
    #htopdummy.GetXaxis().SetBinLabel(21,'500<'+ptmiss+'<700 GeV, 2H, Boosted')
    #htopdummy.GetXaxis().SetBinLabel(22,ptmiss+'>700 GeV, 2H, Boosted')
    htopdummy.GetXaxis().SetBinLabel(1,'150<'+ptmiss+'<200 GeV')
    htopdummy.GetXaxis().SetBinLabel(2,'200<'+ptmiss+'<300 GeV')
    htopdummy.GetXaxis().SetBinLabel(3,'300<'+ptmiss+'<400 GeV')
    htopdummy.GetXaxis().SetBinLabel(4,ptmiss+'>400 GeV')
    htopdummy.GetXaxis().SetBinLabel(5,'150<'+ptmiss+'<200 GeV')
    htopdummy.GetXaxis().SetBinLabel(6,'200<'+ptmiss+'<300 GeV')
    htopdummy.GetXaxis().SetBinLabel(7,'300<'+ptmiss+'<400 GeV')
    htopdummy.GetXaxis().SetBinLabel(8,ptmiss+'>400 GeV')
    htopdummy.GetXaxis().SetBinLabel(9,'150<'+ptmiss+'<200 GeV')
    htopdummy.GetXaxis().SetBinLabel(10,'200<'+ptmiss+'<300 GeV')
    htopdummy.GetXaxis().SetBinLabel(11,'300<'+ptmiss+'<400 GeV')
    htopdummy.GetXaxis().SetBinLabel(12,ptmiss+'>400 GeV')
    htopdummy.GetXaxis().SetBinLabel(13,'150<'+ptmiss+'<200 GeV')
    htopdummy.GetXaxis().SetBinLabel(14,'200<'+ptmiss+'<300 GeV')
    htopdummy.GetXaxis().SetBinLabel(15,'300<'+ptmiss+'<400 GeV')
    htopdummy.GetXaxis().SetBinLabel(16,ptmiss+'>400 GeV')
    htopdummy.GetXaxis().SetBinLabel(17,'300<'+ptmiss+'<500 GeV')
    htopdummy.GetXaxis().SetBinLabel(18,'500<'+ptmiss+'<700 GeV')
    htopdummy.GetXaxis().SetBinLabel(19,ptmiss+'>700 GeV')
    htopdummy.GetXaxis().SetBinLabel(20,'300<'+ptmiss+'<500 GeV')
    htopdummy.GetXaxis().SetBinLabel(21,'500<'+ptmiss+'<700 GeV')
    htopdummy.GetXaxis().SetBinLabel(22,ptmiss+'>700 GeV')
    hgopdummy.LabelsOption('v')

leg = TLegend(0.35, 0.9, 0.75, 0.98)
#leg = TLegend(0.3, 0.9, 0.7, 0.98) #w/ signal
# leg.SetTextSize(30)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetNColumns(3)
if plot_sig:
  leg.SetNColumns(4)

grbkg_pre.SetLineColor(tag1_color)
grbkg_pre.SetLineWidth(1)
grbkg_pre.SetFillColor(0)
grbkg_pre.SetFillStyle(0)
grbkg_pre.Draw('2')
grbkg_horiz_pre.SetLineColor(tag1_color)
grbkg_horiz_pre.SetLineWidth(1)
grbkg_horiz_pre.SetFillColor(0)
grbkg_horiz_pre.SetFillStyle(0)
grbkg_horiz_pre.SetMarkerColorAlpha(tag1_color, 1.0)
grbkg_horiz_pre.Draw('P')
leg.AddEntry(grbkg_pre, tag1_lbl.replace("--","-"), "F")

grbkg_post.SetFillColor(tag2_color)
grbkg_post.SetFillStyle(3144)
grbkg_post.SetLineWidth(0)
grbkg_post.Draw('2')
grbkg_horiz_post.SetLineColor(tag2_color_2)
grbkg_horiz_post.SetFillStyle(3144)
grbkg_horiz_post.SetLineWidth(1)
grbkg_horiz_post.SetMarkerColorAlpha(tag2_color, 1.0)
grbkg_horiz_post.Draw('P')
#grbkg_post.SetLineWidth(1)
#grbkg_post.SetLineColor(tag2_color_2)
leg.AddEntry(grbkg_post, tag2_lbl.replace("--","-"), "F")

if plot_sig:
  grsignal.SetMarkerStyle(20)
  grsignal.SetMarkerColorAlpha(kGreen+3, 1.0)
  grsignal.SetLineColor(kGreen+3)
  grsignal.Draw('P')
  leg.AddEntry(grsignal, "TChiHH-G(450,1)+bkg", "PL")

grdata.SetMarkerStyle(20)
grdata.Draw('P')
leg.AddEntry(grdata, "Data", "eP")

leg.Draw()

gPad.RedrawAxis()

cmslabel = TLatex()
cmslabel.SetTextSize(0.09)
cmslabel.SetNDC(kTRUE)
cmslabel.SetTextAlign(11)
if not plot_pulls:
  cmslabel.SetTextSize(0.07)
if preliminary:
  cmslabel.DrawLatex(top.GetLeftMargin()+0.005, 0.92,"#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}")
elif plot_sig:
  cmslabel.DrawLatex(top.GetLeftMargin()+0.005, 0.92,"#font[62]{CMS} #scale[0.8]{#font[52]{Supplementary}}")
else:
  cmslabel.DrawLatex(top.GetLeftMargin()+0.005, 0.92,"#font[62]{CMS}")
cmslabel.SetTextAlign(31)
cmslabel.SetTextSize(0.05)
cmslabel.DrawLatex(1-top.GetRightMargin()-0.005, 0.92,"#font[42]{137 fb^{-1} (13 TeV)}")

merge_category_y = 3500
if bin_labels:
  merge_category_y = 500

binlabel = TLatex()
binlabel.SetTextSize(an_region_size)
# binlabel.SetNDC(kTRUE)
binlabel.SetTextAlign(21)
binlabel.DrawLatex(4.5, merge_category_y,"Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2")
binlabel.DrawLatex(12.5, merge_category_y,"Resolved, #DeltaR#lower[-0.1]{_{max}} < 1.1")
binlabel.DrawLatex(19.5, merge_category_y,"Boosted")
binlabel.SetTextSize(nb_title_size)

binlabel.DrawLatex(2.5, merge_category_y/1.944,"#font[42]{N_{b}=3}")
binlabel.DrawLatex(6.5, merge_category_y/1.944,"#font[42]{N_{b}=4}")
binlabel.DrawLatex(10.5, merge_category_y/1.944,"#font[42]{N_{b}=3}")
binlabel.DrawLatex(14.5, merge_category_y/1.944,"#font[42]{N_{b}=4}")
binlabel.DrawLatex(18, merge_category_y/1.944,"#font[42]{N_{H}=1}")
binlabel.DrawLatex(21, merge_category_y/1.944,"#font[42]{N_{H}=2}")

if not bin_labels:
  for i in range(4):
      binlabel.SetTextSize(met_title_size)
      binlabel.SetTextAngle(0.0)
      binlabel.DrawLatex(2.5+i*4, 900,ptmiss+" #font[42]{[GeV]}")
      binlabel.SetTextSize(met_range_size)
      binlabel.SetTextAngle(-40.0)
      binlabel.DrawLatex(1+i*4, 350,"#font[42]{150-200}")
      binlabel.DrawLatex(2+i*4, 350,"#font[42]{200-300}")
      binlabel.DrawLatex(3+i*4, 350,"#font[42]{300-400}")
      binlabel.DrawLatex(4+i*4, 350,"#font[42]{>400}")
  
  for i in range(2):
      binlabel.SetTextSize(met_title_size)
      binlabel.SetTextAngle(0.0)
      binlabel.DrawLatex(17+1+i*3, 900,ptmiss+" #font[42]{[GeV]}")
      binlabel.SetTextSize(met_range_size)
      binlabel.SetTextAngle(-40.0)
      binlabel.DrawLatex(17+0+i*3, 350,"#font[42]{300-500}")
      binlabel.DrawLatex(17+1+i*3, 350,"#font[42]{500-700}")
      binlabel.DrawLatex(17+2+i*3, 350,"#font[42]{>700}")

#for i in range(2):
#    binlabel.DrawLatex(4+i*18, 350,"#font[52]{200 < "+ptmiss+"#leq 350 GeV}")
#    binlabel.DrawLatex(10+i*18, 350,"#font[52]{350 < "+ptmiss+"#leq 500 GeV}")
#    binlabel.DrawLatex(16+i*18, 350,"#font[52]{"+ptmiss+"#geq 500 GeV}")
#
#binlabel.SetTextSize(0.045)
#for i in range(6):
#    binlabel.DrawLatex(2+i*6, 150,"#font[52]{1b}")
#    binlabel.DrawLatex(4+i*6, 150,"#font[52]{2b}")
#    binlabel.DrawLatex(6+i*6, 150,"#font[52]{#geq3b}")

#binlabel.SetTextAlign(11)
#binlabel.SetTextSize(0.045)
#binlabel.DrawLatex(31.3, 50,"#font[52]{Low N#lower[-0.1]{_{jets}} (odd bin #)}")
#binlabel.DrawLatex(31.3, 25,"#font[52]{High N#lower[-0.1]{_{jets}} (even bin #)}")

a = TLine()
a.SetLineWidth(1)
a.SetLineStyle(3)
a.SetLineColor(kBlack)
a.DrawLine(8.5,miny,8.5,maxy)
a.DrawLine(16.5,miny,16.5,maxy)
a.DrawLine(4.5,miny,4.5,0.4*maxy)
a.DrawLine(12.5,miny,12.5,0.4*maxy)
a.DrawLine(19.5,miny,19.5,0.4*maxy)

#     Bottom pad
# -----------------------------
if plot_pulls:
  can.cd()
  bottom = TPad("bottom_pad", "bottom_pad", 0., 0., 1., 0.3)
  bottom.SetTopMargin(0.)
  bottom.SetBottomMargin(0.3)
  bottom.SetLeftMargin(left_margin)
  bottom.SetRightMargin(right_margin)
  bottom.SetFillStyle(4000);
  bottom.Draw()

  bottom.cd()
  hbotdummy = TH1D("","",nhbins,0,nhbins+1.0)
  hbotdummy.GetYaxis().SetRangeUser(-2.9,2.9)
  hbotdummy.GetYaxis().SetLabelSize(0.12)
  hbotdummy.GetYaxis().SetTitle("Pull")
  hbotdummy.GetYaxis().SetNdivisions(206)
  hbotdummy.GetYaxis().CenterTitle()
  hbotdummy.GetYaxis().SetTitleSize(0.15)
  hbotdummy.GetYaxis().SetTitleOffset(0.2)
  
  hbotdummy.GetXaxis().SetLabelSize(0.12)
  hbotdummy.GetXaxis().SetLabelOffset(0.02)
  hbotdummy.GetXaxis().SetTitle("Bin number")
  hbotdummy.GetXaxis().SetNdivisions(24);
  hbotdummy.GetXaxis().CenterTitle()
  hbotdummy.GetXaxis().SetTitleSize(0.15)
  hbotdummy.GetXaxis().SetTitleOffset(0.9)
  #hbotdummy.GetXaxis().SetRangeUser(0.01,nhbins+0.99)
  hbotdummy.GetXaxis().SetLimits(0.01,nhbins+0.99)

  hbotdummy.Draw()
  
  grpull_pre.SetLineColor(tag1_color)
  grpull_pre.SetLineWidth(1)
  grpull_pre.SetFillColor(0)
  grpull_pre.SetFillStyle(0)
  grpull_pre.Draw('B')

  grpull_sig.SetLineColor(kRed+2)
  grpull_sig.SetLineWidth(1)
  grpull_sig.SetFillColor(0)
  grpull_sig.SetFillStyle(0)
  #grpull_sig.Draw('B')
  
  grpull_post.SetFillColor(tag2_color)
  grpull_post.SetFillStyle(3144)
  #grpull_post.Draw('B')
  #not well defined for pulls or likelihood ratio significance
  
  a = TLine()
  a.SetLineWidth(1)
  a.SetLineStyle(3)
  a.DrawLine(8.5, -2.9, 8.5, 2.9)
  a.DrawLine(16.5, -2.9, 16.5, 2.9)
  a.SetLineColor(kGray+1)
  a.DrawLine(4.5, -2.9, 4.5, 2.9)
  a.DrawLine(12.5, -2.9, 12.5, 2.9)
  a.DrawLine(19.5, -2.9, 19.5, 2.9)
  
  b = TLine()
  b.SetLineWidth(1)
  b.SetLineColor(kGray+1)
  b.SetLineStyle(1)
  b.DrawLine(0,0, nhbins+1,0)
  b.SetLineStyle(2)
  for i in [-1.,1.]:
      b.DrawLine(0,i, nhbins+1,i)
  b.SetLineStyle(3)
  for i in [-2.,2.]:
      b.DrawLine(0,i, nhbins+1,i)
  #b.SetLineStyle(3)
  #for i in [-3.,3.]:
  #    b.DrawLine(0,i, nhbins+1,i)

pname = 'plots/results_plot.pdf'
if preliminary:
  pname = 'plots/results_plot_preliminary.pdf'
if plot_sig:
  pname = 'plots/results_plot_with_signal.pdf'
can.Print(pname)

print('open '+pname)

output_file = TFile.Open('tables/CMS-SUS-20-004_Figure_010.root','RECREATE')
grbkg_pre.Write('background_prefit')
grbkg_post.Write('background_postfit')
histdata.GetXaxis().SetBinLabel(1,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 3b, 150<'+ptmiss+'<200 GeV')
histdata.GetXaxis().SetBinLabel(2,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 3b, 200<'+ptmiss+'<300 GeV')
histdata.GetXaxis().SetBinLabel(3,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 3b, 300<'+ptmiss+'<400 GeV')
histdata.GetXaxis().SetBinLabel(4,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 3b, '+ptmiss+'>400 GeV')
histdata.GetXaxis().SetBinLabel(5,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 4b, 150<'+ptmiss+'<200 GeV')
histdata.GetXaxis().SetBinLabel(6,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 4b, 200<'+ptmiss+'<300 GeV')
histdata.GetXaxis().SetBinLabel(7,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 4b, 300<'+ptmiss+'<400 GeV')
histdata.GetXaxis().SetBinLabel(8,'Resolved, 1.1 < #DeltaR#lower[-0.1]{_{max}} < 2.2, 4b, '+ptmiss+'>400 GeV')
histdata.GetXaxis().SetBinLabel(9,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 3b, 150<'+ptmiss+'<200 GeV')
histdata.GetXaxis().SetBinLabel(10,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 3b, 200<'+ptmiss+'<300 GeV')
histdata.GetXaxis().SetBinLabel(11,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 3b, 300<'+ptmiss+'<400 GeV')
histdata.GetXaxis().SetBinLabel(12,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 3b, '+ptmiss+'>400 GeV')
histdata.GetXaxis().SetBinLabel(13,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 4b, 150<'+ptmiss+'<200 GeV')
histdata.GetXaxis().SetBinLabel(14,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 4b, 200<'+ptmiss+'<300 GeV')
histdata.GetXaxis().SetBinLabel(15,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 4b, 300<'+ptmiss+'<400 GeV')
histdata.GetXaxis().SetBinLabel(16,'Resolved, Low #DeltaR#lower[-0.1]{_{max}} < 1.1, 4b, '+ptmiss+'>400 GeV')
histdata.GetXaxis().SetBinLabel(17,'Boosted, 1H, 300<'+ptmiss+'<500 GeV')
histdata.GetXaxis().SetBinLabel(18,'Boosted, 1H, 500<'+ptmiss+'<700 GeV')
histdata.GetXaxis().SetBinLabel(19,'Boosted, 1H, '+ptmiss+'>700 GeV')
histdata.GetXaxis().SetBinLabel(20,'Boosted, 2H, 300<'+ptmiss+'<500 GeV')
histdata.GetXaxis().SetBinLabel(21,'Boosted, 2H, 500<'+ptmiss+'<700 GeV')
histdata.GetXaxis().SetBinLabel(22,'Boosted, 2H, '+ptmiss+'>700 GeV')
histdata.Write('data')
output_file.Close()
print('open tables/CMS-SUS-20-004_Figure_010.root')
