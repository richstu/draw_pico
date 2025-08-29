"""Quickly generates a ROOT file with templates for quick testing of datacards
with combine without having to run the full HtoZg_fitting machinery

MO 26.02.2025
"""

import ROOT
from argparse import ArgumentParser

#def get_datacard_names(datacard_filename: str) -> list[tuple[str,str]]:
def get_datacard_names(datacard_filename):
  """Extracts workspace and datset/PDF names from combine datacard

  Args:
    datacard_filename: name of datacard txt file

  Returns:
    list of tuples (workspace name, DataSet/PDF name)
  """
  ret = []
  with open(datacard_filename) as datacard:
    line = datacard.readline()
    while line != '':
      if (line[:6]=='shapes'):
        ws_name, pdf_name = line.split()[4].split(':')
        ret.append((ws_name, pdf_name))
      line = datacard.readline()
  return ret

#def get_mllg_varname(dataset_name: str) -> str:
def get_mllg_varname(dataset_name):
  """Extracts name of fitting variable from PDF name

  Args:
    dataset_name: name of PDF/dataset

  Returns:
    name of associated mllg fit variable
  """
  name_split = dataset_name.split('_')
  cat_name =  ''
  for ifragment in range(len(name_split)):
    if name_split[ifragment]=='cat':
      cat_name = name_split[ifragment+1]
      break
  if cat_name=='':
    raise ValueError('string "cat" not present in dataset_name')
  return 'mllg_cat_'+cat_name

def round_bins(hist):
  """Rounds bins of a TH1D to the nearest integer, also removing negative 
  yields

  Args:
    hist: histogram to round (modified in place)
  """
  for ibin in range(hist.GetNbinsX()+2):
    rounded_value = float(round(hist.GetBinContent(ibin)))
    if (rounded_value < 0.0):
      rounded_value= 0.0
    hist.SetBinContent(ibin, rounded_value)

def remove_negative_bins(hist):
  """Removes negative bins in histogram

  Args:
    hist: histogram to round (modified in place)
  """
  for ibin in range(hist.GetNbinsX()+2):
    normed_value = hist.GetBinContent(ibin)
    if (normed_value < 0.0):
      normed_value= 0.0
    hist.SetBinContent(ibin, normed_value)

if __name__=='__main__':
  #parse args and set up
  argument_parser = ArgumentParser(prog='datacard_quickconvert',
      description='generates simple workspaces and ROOT file for datacard')
  argument_parser.add_argument('-i','--input_datacard')
  args = argument_parser.parse_args()
  if (args.input_datacard[-4:] != '.txt'):
    raise ValueError('Input should be datacard txt file.')
  input_filename = args.input_datacard[:-4]+'_rawdata.root'
  output_filename = args.input_datacard[:-4]+'.root'
  datacard_names = get_datacard_names(args.input_datacard)

  #generate output
  input_file = ROOT.TFile(input_filename, 'READ')
  for ws_name, pdf_name in datacard_names:
    if ('data_obs' in ws_name):
      #save RooDataSet/RooDataHist
      data = input_file.Get(ws_name).data(pdf_name)
      mllg = input_file.Get(ws_name).var(get_mllg_varname(pdf_name))
      data_hist = ROOT.TH1D('data_hist_'+pdf_name, '', 80, 100, 180)
      data.fillHistogram(data_hist, ROOT.RooArgList(mllg))
      round_bins(data_hist)
      datahist = ROOT.RooDataHist(pdf_name, pdf_name, ROOT.RooArgList(mllg), 
                                  data_hist)
      ws = ROOT.RooWorkspace(ws_name)
      getattr(ws,'import')(datahist)
      ws.writeToFile(output_filename,False)
    else:
      #save RooHistPdf
      ws_name_clean = ws_name.replace('background','mc')
      pdf_name_clean = pdf_name.replace('background','mc').replace('pdf_','')
      mcdata = input_file.Get(ws_name_clean).data(
          'mcdata_'+pdf_name_clean+'_nominal')
      mllg = input_file.Get(ws_name_clean).var(get_mllg_varname(pdf_name))
      hist_name = pdf_name+'_datahist'
      pdf_hist = ROOT.TH1D('data_hist_'+pdf_name, '', 80, 100, 180)
      data.fillHistogram(pdf_hist, ROOT.RooArgList(mllg))
      remove_negative_bins(pdf_hist)
      datahist = ROOT.RooDataHist(pdf_name, pdf_name, ROOT.RooArgList(mllg), 
                                  pdf_hist)
      pdf = ROOT.RooHistPdf(pdf_name, pdf_name, ROOT.RooArgSet(mllg), datahist)
      ws = ROOT.RooWorkspace(ws_name)
      getattr(ws,'import')(pdf)
      if 'background' in pdf_name:
        nentries = mcdata.sumEntries()
        norm = ROOT.RooRealVar(pdf_name+'_norm', pdf_name, nentries, 0, 
                               3*nentries)
        getattr(ws,'import')(norm)
      ws.writeToFile(output_filename,False)
  input_file.Close()
