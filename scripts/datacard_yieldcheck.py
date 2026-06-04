#script for quickly checking data sideband yields

from argparse import ArgumentParser
import ROOT

CATS = ['ggf4','ggf3','ggf2','ggf1','vbf4','vbf3','vbf2','vbf1','vh3l','vhmet','tthlep','tthhad','untagged']

ROOT.gInterpreter.Declare("""

RooAbsData* cast_rooabsdata_to_roodatahist(RooAbsData* object)
{ return static_cast<RooDataHist*>(object); }

""")

if __name__=='__main__':
  argument_parser = ArgumentParser()
  argument_parser.add_argument('-i','--input_data')
  args = argument_parser.parse_args()
  input_file = ROOT.TFile(args.input_data, 'READ')
  for cat in CATS:
    ws = getattr(input_file, f'WS_data_obs_cat_{cat}')
    hist = ROOT.cast_rooabsdata_to_roodatahist(ws.data(f'data_obs_cat_{cat}'))
    total_yield = 0.0
    for ibin in range(340):
      #skip bins in blinded window
      if (ibin >= 100 and ibin < 140):
        continue
      #unclear what volume actually is, we'll assume it is yield*width
      total_yield += hist.weight(ibin)
    print(f'Cat {cat} yield {total_yield}') 

  input_file.Close()


