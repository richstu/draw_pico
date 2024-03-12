import ROOT
import os
import sys
import argparse

#sys.path.append(os.path.abspath("scripts/zg_fitutils"))
from bkg_functions_fit import *
from bkg_functions_class import *
from Xc_Minimizer import *
from plot_utility import *

def singleBernFTest(x, gauss_mu, histogram, cat = "", method = "Chi2", e_type = "Poisson", eps = 0.1, offset = False, strategy = 0):
    bern2_model = Bern2Class(x, gauss_mu, cat, 10, 0.3, 10, 7., 105.)
    bern3_model = Bern3Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern4_model = Bern4Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    bern5_model = Bern5Class(x, gauss_mu, cat, 10, 0.3, 50, 7., 105.)
    if e_type == "Poisson": error = ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson)
    elif e_type == "SumW2": error = ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2)

    if method == "Chi2": 
        stat1 = ROOT.RooChi2Var("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram, error)
        stat2 = ROOT.RooChi2Var("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram, error)
        stat3 = ROOT.RooChi2Var("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram, error)
        stat4 = ROOT.RooChi2Var("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram, error)
    elif method == "NLL":
        stat1 = ROOT.RooNLLVar("stat_bern2_" + cat, "stat bern2 " + cat, bern2_model.pdf,  histogram)
        stat2 = ROOT.RooNLLVar("stat_bern3_" + cat, "stat bern3 " + cat, bern3_model.pdf,  histogram)
        stat3 = ROOT.RooNLLVar("stat_bern4_" + cat, "stat bern4 " + cat, bern4_model.pdf,  histogram)
        stat4 = ROOT.RooNLLVar("stat_bern5_" + cat, "stat bern5 " + cat, bern5_model.pdf,  histogram)

    stats = [stat1, stat2, stat3, stat4]
    if method == "Chi2": 
        for entry in stats:
            print(entry.GetTitle())
            Minimizer_Chi2(entry, -1, 100, False, strategy)
            r = Minimizer_Chi2(entry, -1, eps, offset, strategy)
            r.Print("V")
            # print("Cov q = ", r.covQual(), " status = ", r.status(), end="\n\n")
    elif method == "NLL": 
        for entry in stats:
            print(entry.GetTitle())
            Minimizer_NLL(entry, -1, 100, False, strategy)
            r = Minimizer_NLL(entry, -1, eps, offset, strategy).Print("V")
            r.Print("V")
    output = [stat1.getVal(), stat2.getVal(), stat3.getVal(), stat4.getVal()]
    fs = []
    if method == "Chi2":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.fdistribution_cdf_c((output[i] - output[i+1])*(260-5-i)/output[i+1], 1, 260-5-i))
    elif method == "NLL":
        for i in range(len(output) - 1):
            fs.append(ROOT.Math.chisquared_cdf_c(2*(output[i] - output[i+1]), 1))

    print(method, " = ", output)
    print("P-value = ", fs)
    x_argset = ROOT.RooArgSet(x)
    #plotClass(x_argset, histogram, bern2_model)

if __name__=='__main__':
  parser = argparse.ArgumentParser(description = "F test method (NLL or Chi2)")
  parser.add_argument("-i","--input_file",default='datacards/test_datacard.root')
  parser.add_argument("-m","--method",default="Chi2")
  args = parser.parse_args()
  if not(args.method == "Chi2" or args.method == "NLL") :
      print("Please use the correct method.")
      sys.exit(1)

  datacard_file = ROOT.TFile(args.input_file,'READ')
  ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

  #Definition of categories
  categories = ['ggh4','ggh3','ggh2','ggh1','vbf3','vbf2','vbf1','vh3l','vhmet','tthhad','tthlep']

  for category in categories:
    ws = datacard_file.Get('WS_data_obs_cat_'+category)
    data_hist = ws.data('data_obs_cat_'+category)
    mllg = ws.var('mllg_cat_'+category)
    mu_gauss = ROOT.RooRealVar("mu_gauss","always 0"       ,0.)
    singleBernFTest(mllg, mu_gauss, data_hist, "u1", args.method, "Poisson", False)

  datacard_file.Close()



  
  
