import argparse
import glob
import ROOT

if __name__=='__main__':
  #parse arguments
  argument_parser = argparse.ArgumentParser(prog='get_nano_sumw',
      description='Gets the total sum of weights for a set of NanoAOD files')
  argument_parser.add_argument('-i','--input_filenames')
  args = argument_parser.parse_args()

  total_sumw = 0.0
  input_filenames = glob.glob(args.input_filenames)
  for input_filename in input_filenames:
    input_file = ROOT.TFile(input_filename,'READ')
    for event in input_file.Runs:
      total_sumw += event.genEventSumw

  print 'Total sum of weights is',
  print(total_sumw)

