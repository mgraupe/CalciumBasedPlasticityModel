from collections import *

# fixed parameter and parameter ranges
thetaD = 1.
thetaPH = 1.3884 # 1.24133883 # :4 1.3884 :1 # 1.2885 # :3
thetaPJ = 1.63069609#1.3
nonlinear =2. #2.
w0 = 0.5

# short-term plasticity parameters
# first value is from Loebel et al. 2009, second is from fit of  Sjoestrom 2003
who = 'henry' # 'henry' or 'jesper'
Npresentations = [10,15]  # henry, jesper
presentationInterval = [4,10] # in s henry, jesper
tauRec = [0.525,0.1489192] # henry, jesper
#U = [0.,0.]
U = [0.46,0.383753]  # henry, jesper
Nvesicles = 10
Npairs = 5

# fitting loops 
Nruns = 1000000
threshold = 1. # lower threshold to include parameter set

# smoothnessWeight = 1.

# parameterLimits
limits = OrderedDict([
    ('tauCa',[0.015,0.1]),
    ('Cpre',[0.1,4.]),
    ('Cpost',[0.3,4.]),
    ('thetaP',[1.2,4.1]),
    ('gammaD',[20,1000.]),
    ('gammaP',[100,1000.]),
    ('tau',[1.,50000.]),
    ('D',[0.0,0.015]),
])

