

#signalStrengths = {'ggF 1' : , 'ggF 2' : , 'ggF 3' : , 'ggF 4' : ,
#                   'VBF 1' : , 'VBF 1' : , 'VBF 1' : , 'VBF 1' : ,
#                   'ggF 1' : , 'ggF 1' : , 'ggF 1' : , 'ggF 1' : , 
#                   'ggF 1' : , }

#def createErrorBar(color,):

def yAxisDistance(yAxisSpacing = 3.3, offset=0):
  ydist = [10, 10, 10, 10,  7,  7,  7,  7,  4,  4, 4, 4, 4, 2]
  ydist.reverse()
  y = [yAxisSpacing*idx + ydist[idx] for idx in range(len(ydist))]
  y.reverse()
  y = [3*entry for entry in y]
  y = [entry - offset for entry in y]
  return y




categories = ['ggF 1', 'ggF 2', 'ggF 3', 'ggF 4', 'VBF 1', 'VBF 2', 'VBF 3', 'VBF 4', r'VH $N_{lep} \geq 3$', 'VH $p_{T}^{miss}$', r'$\text{t}\overline{\text{t}}\text{H}$ leptonic', r'$\text{t}\overline{\text{t}}\text{H}$ hadronic', 'Untagged', 'Simultaneous fit']
x = [1.451, 0.23, 0.07, 3.26, 0.62, 0.42, 1.99, 8.95, 7.22, 0.351, 1.06, -2.22, -1.36, 1.10]

xTotErrHigh = [2.06, 2.49, 1.05, 1.95, 1.22, 1.08, 1.46, 4.43, 6.98, 13.47, 4.04, 6.18, 5.28, 0.52]
xTotErrLow  = [1.58, 1.36, 0.99, 1.92, 0.89, 0.97, 1.36, 4.22, 5.75, 13.70, 3.06, 5.45, 5.35, 0.61]

xSystErrHigh = [0.499, 0.399, 0.109, 0.470, 0.229, 0.194, 0.396, 1.357, 3.221, 2.429, 0.572, 0.609, 0.357, 0.16]
xSystErrLow  = [0.066, 0.123, 0.008, 0.112, 0.118, 0.079, 0.204, 0.398, 0.327, 1.633, 0.289, 1.194, 0.778, 0.10]

xStatErrHigh = [2.003, 2.459, 1.042, 1.901, 1.202, 1.064, 1.407, 4.221, 6.190, 13.25, 4.001, 6.146, 5.268, 0.50]
xStatErrLow  = [1.583, 1.356, 0.986, 1.903, 0.880, 0.964, 1.343, 4.197, 5.739, 13.60, 3.042, 5.313, 5.293, 0.60]
