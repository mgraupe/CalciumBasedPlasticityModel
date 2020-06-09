from numpy import *

###########################################################################
# full simulation

sJFullSim0 = (array([3.83492083e-02, 3.99132241e+00, 1.12940834e+00, 1.63069609e+00,
        1.11320539e+02, 5.64392975e+02, 2.99877800e+02, 9.23545841e-03]),
 0.08000169726924594,
 1469,
 2273,
 0)

sJFullNoSTDSim0 = (array([3.21900754e-02, 1.60681037e+00, 1.12436420e+00, 3.19759883e+01,
        1.61987985e+02, 7.99756573e+01, 5.75272377e-03]),
 0.031440976841944995,
 1592,
 2452,
 0)


sJFullNonlinSim0 = (array([4.16066039e-02, 3.00753343e+00, 1.20342240e+00, 2.10000000e+00,
        2.00766828e+02, 9.99990054e+02,  7.19895333e+02, 9.99429160e-03]),
 0.07806673255109622,
 1512,
 2387,
 0)

############################################################################
sHFullSim1 = (array([4.89774484e-02, 2.41618557e+00, 1.38836494e+00, 1.38843434e+00,
        1.76541097e+02, 5.79578738e+02, 1.43096290e+02, 1.00700540e-02]),
 0.008389786038669585,
 2331,
 3456,
 0)

sHFullSim2 = (array([6.66711823e-02, 1.98364185e+00, 9.31190036e-01, 1.18809086e+00,
        2.85953581e+02, 8.71155123e+02, 4.77417713e+02, 5.00000000e-03]),
 0.011051309270841778,
 2932,
 4243,
 0)

sHFullSim3 = (array([4.85950258e-02, 2.33853752e+00, 1.28841747e+00, 1.28848216e+00,
        3.41941772e+02, 9.99999694e+02,  2.43391252e+02, 1.04677682e-02]),
 0.011106812148214289,
 1778,
 2771,
 0)



sHFullSim4 = (array([6.29867685e-02, 1.89727782e+00, 1.03131758e+00, 1.24133883e+00,
        2.77229124e+02, 8.73563447e+02, 3.74056796e+02, 5.00000000e-03]),
 0.011459122657145741,
 1976,
 2964,
 0)


sHFullNoSTDSim1 = (array([3.45093185e-02, 5.07766907e-01, 1.39498534e+00, 1.35988823e+00,
        1.89868323e+02, 7.04119738e+02, 7.27739718e+01, 8.18018412e-03]),
 1.9563232946153946e-05,
 2114,
 3110,
 0)

sHFullNoSTDSim3 = (array([4.33362310e-02, 1.97759135e+00, 4.79345689e-01, 1.60761785e+02,
        5.00510565e+02, 1.79574188e+02, 6.53234809e-03]),
 0.003620827154466764,
 1526,
 2353,
 0)
sHFullNoSTDSim3b =  (array([3.28590422e-02, 4.35222738e-01, 1.29706333e+00, 1.56233834e+02,
        5.72058546e+02, 4.65163078e+01, 6.94446366e-03]),
 0.00012036275965056429,
 1696,
 2566,
 0)

 # new parameter set of choice with nonlinear calicum dynamics
sHFullNonlinSim1 = (array([4.27662771e-02, 2.51325173e+00, 1.40326199e+00, 1.40328478e+00,
        2.55848921e+02, 8.27284991e+02, 2.11114257e+02, 9.13482505e-03]),
 0.010169857354502439,
 2937,
 4261,
 0)

############################################################################
# solTP11n0c , with nl=1

linearCaModel = (array([  2.22721171e-02,   8.44990637e-01,   1.62137900e+00,
           2.00928899e+00,   1.37758631e+02,   5.97089216e+02,
           5.20761286e+02,   9.53708736e-03]),
  20.117666357585907,
  1757,
  2661,
  0)



############################################################################
# solNL10n10,  with nl=2

nonlinearCaModel = (array([  1.89304409e-02,   8.64671725e-01,   2.30814844e+00,
          4.99779542e+00,   1.11825154e+02,   8.94236947e+02,
          7.07022578e+02,   1.00000000e-02]),
 10.524875310154474,
 1117,
 1790,
 0)

#############################################################################

stpHenryCaModel = (array([  5.55958868e-02,   6.89318952e-01,   1.07527446e+00,
          1.2, 4.52551145e+01,   1.39363720e+02,   1.79134285e+01,
          5.34374094e-03]), 0.010786124085281916, 1010, 1553, 0)

stpHenryCaModel1 = (array([4.29623682e-02, 9.89689797e-01, 1.30113653e+00, 1.30115146e+00,
        1.41318889e+02, 4.75062273e+02, 8.20629282e+00, 9.35691197e-03]),
 0.006932448507854171,
 1364,
 2033,
 0)


HenryCaModel = (array([  4.34915697e-02,   7.60879350e-01,   5.40444149e-01,
          1.20000193e+00,   1.02991787e+02,   3.43294617e+02,
          2.24381827e+01,   5.44177570e-03]),
 0.0025496328288422401,
 1710,
 2607,
 0)

HenryCaModel2 = (array([  4.29332026e-02,   6.89399376e-01,   5.13738624e-01,
          1.10501587e+00,   1.81391042e+02,   5.54263112e+02,
          6.49272005e+01,   9.03305205e-03]),
 0.00029582039748969732,
 3249,
 4777,
 0)

##############################################################################

stpJesperCaModel = (array([  3.96510510e-02,   1.99545551e+00,   9.82700197e-01,
          1.53428149e+00,   1.79933101e+02,   9.23424261e+02,
          4.05025178e+02,   1.00000000e-02]),
 0.056369708030680699,
 1198,
 1819,
 0)

JesperCaModel = (array([  3.30575851e-02,   7.72419466e-01,   9.99999939e-01,
          1.58918925e+00,   1.61342338e+02,   8.20904263e+02,
          3.12655733e+02,   1.00000000e-02]),
 0.021564784713236817,
 4474,
 6568,
 0)

##############################################################################

tmmFacilitationJesperCaModel = (array([  3.96510510e-02,   1.99545551e+00,   9.82700197e-01,
          1.53428149e+00,   1.79933101e+02,   9.23424261e+02,
          4.05025178e+02,   1.00000000e-02]),
 0.056369708030680699,
 1198,
 1819,
 0)

##############################################################################

tmmfacilitationDepressionJesperCaModel = (array([  3.96510510e-02,   1.99545551e+00,   9.82700197e-01,
          1.53428149e+00,   1.79933101e+02,   9.23424261e+02,
          4.05025178e+02,   1.00000000e-02]),
 0.056369708030680699,
 1198,
 1819,
 0)
