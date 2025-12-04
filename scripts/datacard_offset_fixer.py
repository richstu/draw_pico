#Fix asymmetric systematics in datacard

from statistics import mean, stdev
#from argparse import ArgumentParser

#def line_item(str_list: list[str], idx: int) -> str:
def line_item(str_list, idx):
  if idx >= len(str_list):
    return ''
  return str_list[idx]

if __name__ == '__main__':
  #todo add arguments
  input_filename = 'datacards/hzg_datacard_v1p2.txt'
  output_filename = 'datacards/hzg_datacard_v1p2_cleaned.txt'
  output_content = ''
  with open(input_filename, 'r') as datacard:
    input_content = datacard.read().split('\n')
    for line in input_content:
      line_split = line.split()
      #ignore non-systematic lines
      if line_item(line_split, 1) != 'lnN':
        output_content += line
        output_content += '\n'
      else:
        output_content += line_split[0].ljust(35)
        output_content += line_split[1].ljust(9)
        for icat in range(2, len(line_split)):
          if not ('/' in line_split[icat]):
            output_content += line_split[icat].ljust(19)
          else:
            up_var = float(line_split[icat].split('/')[0])
            dn_var = float(line_split[icat].split('/')[1])
            if ((up_var > 1.0 and dn_var < 1.0) 
                or (up_var < 1.0 and dn_var > 1.0)):
              output_content += line_split[icat].ljust(19)
            elif abs(up_var-dn_var) < 0.00001:
              output_content += '-'.ljust(19)
            else:
              avr = ((up_var-1.0) + (dn_var-1.0))/2.0
              up_var -= avr
              dn_var -= avr
              print(f'Shifting {line_split[0]} cat {icat}')
              syst_string = '{:.5f}'.format(up_var)+'/'+'{:.5f}'.format(dn_var)
              output_content += syst_string.ljust(19)
        output_content += '\n'
  with open(output_filename, 'w') as datacard:
    datacard.write(output_content)
    print(f'Wrote output to {output_filename}')

