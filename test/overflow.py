from petbox import dca
import numpy as np
np.seterr(all='print')

print(dca.MH(qi=1, Di=0.594, bi=1e-03, Dterm=.00001))
dca.MH(qi=1, Di=0.594, bi=1e-03, Dterm=.01)

print(dca.MH(qi=1, Di=0.594, bi=1e-04, Dterm=.01))
dca.MH(qi=1, Di=0.594, bi=1e-04, Dterm=.01)

print(dca.MH(qi=1, Di=0.594, bi=1e-05, Dterm=.01))
dca.MH(qi=1, Di=0.594, bi=1e-05, Dterm=.01)

print(dca.MH(qi=1, Di=0.594, bi=1e-06, Dterm=.01))
dca.MH(qi=1, Di=0.594, bi=1e-06, Dterm=.01)

print(dca.MH(qi=1, Di=0.594, bi=1e-010, Dterm=.0000001))
dca.MH(qi=1, Di=0.594, bi=1e-07, Dterm=.01)
