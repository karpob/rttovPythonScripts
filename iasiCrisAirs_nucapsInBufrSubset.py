#!/usr/bin/env python3
import  configparser
pathInfo = configparser.ConfigParser()
pathInfo.read('rttov.cfg')
rttovPath = pathInfo['RTTOV']['rttovPath']
import matplotlib
matplotlib.use('Agg')
import pyrttov
import example_data as ex
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
import matplotlib
import matplotlib.cm as mplcm
import matplotlib.colors as colors
from cycler import cycler
import h5py

pyRttovPath = os.path.join(rttovPath,'wrapper')
if not pyRttovPath in sys.path:
    sys.path.append(pyRttovPath)

wn_ozoneNucaps = {}
wn_ozoneNucaps['cris-fsr'] = [1019.3750, 1020.6250, 1021.8750, 1023.1250, 1024.3750, 1025.6250, 1026.8750, 1031.2500, 1033.1250, 1034.3750, 1035.6250, 1036.8750, 1038.1250, 1040.6250, 1048.1250, 1053.1250, 1054.3750, 1056.8750, 1061.2500]

wn_ozoneNucaps['cris'] = [ 996.8750, 997.5000, 998.1250, 999.3750, 1000.0000, 1000.6250, 1001.8750, 1002.5000, 1003.1250, 1005.0000, 1007.5000, 1008.1250, 1010.0000, 1011.8750, 1012.5000, 1014.3750, 1015.0000, 1016.2500, 1018.1250, 1018.7500, 1020.6250, 1022.5000, 1023.1250, 1025.0000, 1026.8750, 1028.7500, 1031.2500, 1033.1250, 1034.3750, 1035.0000, 1036.2500, 1038.1250, 1040.6250, 1041.8750, 1045.6250, 1047.5000, 1048.1250, 1049.3750, 1050.0000, 1050.6250, 1051.8750, 1053.1250, 1053.7500, 1055.6250, 1056.2500, 1056.8750, 1058.1250, 1058.7500, 1060.0000, 1061.2500, 1063.7500, 1066.2500, 1068.1250]

wn_ozoneNucaps['iasi'] = [ 997.0000, 998.2500, 999.7500, 1000.7500, 1001.5000, 1002.2500, 1003.2500, 1004.7500, 1005.2500, 1006.0000, 1007.2500, 1008.2500, 1009.7500, 1010.5000, 1012.0000, 1013.2500, 1014.5000, 1015.5000, 1016.5000, 1018.2500, 1018.7500, 1020.2500, 1021.0000, 1022.2500, 1023.0000, 1024.2500, 1025.0000, 1026.2500, 1027.0000, 1027.7500, 1029.0000, 1030.0000, 1031.0000, 1031.7500, 1033.0000, 1034.7500, 1036.7500, 1038.2500, 1039.5000, 1040.5000, 1041.5000, 1046.2500, 1051.2500, 1054.5000, 1055.5000, 1057.7500, 1059.5000, 1061.2500, 1062.5000, 1063.5000, 1065.0000, 1068.2500, 1069.0000]

wn_ozoneNucaps['airs'] = [ 1001.3840, 1005.2627, 1008.2999, 1010.4804, 1013.1092, 1016.6352, 1030.5279, 1061.3320, 1061.8124, 1063.2559, 1065.1864, 1068.5809 ]

def calculateWeightingFunctions(chan_list, rttovInstance, myProfiles):
    wf = []
    nlevels = np.asarray(rttovInstance.TauLevels).shape[2] 
    for c in range(1,rttovInstance.Nchannels+1):
        if c in chan_list:
            num = rttovInstance.TauLevels[0, c-1, 1::] - rttovInstance.TauLevels[0, c-1, 0:nlevels-1]
            den = np.log(myProfiles.P[0, 0:nlevels-1]) - np.log(myProfiles.P[0, 1:nlevels])
            wf.append(num/den)
    return wf
 
def plotWeightingFunctions(chan_list, profiles, wf, instrument):
    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = len(chan_list)
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    for item in dir(plt.gca()):
        print(item)
    plt.gca().set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
    #plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    for f in wf:
        plt.plot( f[0:99], profiles.P[0,0:99] )
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    print(len(chan_list),len(wn_ozoneNucaps[instrument]))
    for i,c in enumerate(chan_list):
        legendList.append('{} {:4.3f} {:4.3f}'.format(c, wn_ozoneNucaps[instrument][i], 10000.0/wn_ozoneNucaps[instrument][i]))
    plt.legend(legendList, fontsize=6,  ncol=3)
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Weighting Function')
    plt.savefig(instrument+'_weighting_functions.png')
    
    if(len(chan_list)>30): matplotlib.rc('xtick', labelsize=6) 
    else: matplotlib.rc('xtick', labelsize=10) 
 

    plt.figure()
    plt.title('Weighting Function by Instrument Channel')
    wf = np.asarray(wf)
    plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T, norm = LogNorm( vmin = wf[:,0:99].min(), vmax = wf[:,0:99].max() ) )
    plt.colorbar()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Instrument Channel')
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.xticks(np.arange(len(chan_list)), chan_list, rotation='vertical')
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    plt.tight_layout() 
    plt.savefig(instrument+'_weighting_functions_pcolor.png')



if __name__ == '__main__':

    # This example program simulates two profiles for each of three instruments
    # The example profile data are defined in example_data

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nlevels = len(ex.p_ex)
    nprofiles = 2
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    # Associate the profiles and other data from example_data.h with myProfiles
    # Note that the simplecloud, clwscheme, icecloud and zeeman data are not mandatory and
    # are omitted here

    def expand2nprofiles(n, nprof):
        # Transform 1D array to a [nprof, nlevels] array
        outp = np.empty((nprof, len(n)), dtype=n.dtype)
        for i in range(nprof):
            outp[i, :] = n[:]
        return outp

    myProfiles.GasUnits = ex.gas_units
    myProfiles.P = expand2nprofiles(ex.p_ex, nprofiles)
    myProfiles.T = expand2nprofiles(ex.t_ex, nprofiles)
    # Modify the temperature of the second profile
    myProfiles.T[1, :] += 2
    myProfiles.Q = expand2nprofiles(ex.q_ex, nprofiles)
    myProfiles.CO2 = expand2nprofiles(ex.co2_ex, nprofiles)
    myProfiles.O3 = expand2nprofiles(ex.o3_ex, nprofiles)
    myProfiles.Angles = ex.angles.transpose()
    myProfiles.S2m = ex.s2m.transpose()
    myProfiles.Skin = ex.skin.transpose()
    myProfiles.SurfType = ex.surftype.transpose()
    myProfiles.SurfGeom = ex.surfgeom.transpose()
    myProfiles.DateTimes = ex.datetimes.transpose()

    # ------------------------------------------------------------------------
    # Set up Rttov instances for each instrument
    # ------------------------------------------------------------------------

    # Create three Rttov objects for Four instruments
    iasiRttov = pyrttov.Rttov()
    crisRttov = pyrttov.Rttov()
    crisFsrRttov = pyrttov.Rttov()
    airsRttov = pyrttov.Rttov()

    chan_list_fsr  = [ 592, 594, 596, 598, 600, 602, 604, 611, 614, 616, 618, 620, 622, 626, 638, 646, 648, 652, 659]
    chan_list_cris =  [556, 557, 558, 560, 561, 562, 564, 565, 566, 569, 573, 574, 577, 580, 581, 584, 585, 587, 590, 591, 594,\
                                       597, 598, 601, 604, 607, 611, 614, 616, 617, 619, 622, 626, 628, 634, 637, 638, 640, 641, 642, 644, 646,\
                                       647, 650, 651, 652, 654, 655, 657, 659, 663, 667, 670]
    chan_list_iasi  = [1409, 1414, 1420, 1424, 1427, 1430, 1434, 1440, 1442, 1445, 1450, 1454, 1460, 1463, 1469, 1474, 1479, 1483, 1487, 1494,\
                                      1496, 1502, 1505, 1510, 1513, 1518, 1521, 1526, 1529, 1532, 1537, 1541, 1545, 1548, 1553, 1560, 1568, 1574, 1579, 1583,\
                                      1587, 1606, 1626, 1639, 1643, 1652, 1659, 1666, 1671, 1675, 1681, 1694, 1697]
    chan_list_airs =[ 1003, 1012, 1019, 1024, 1030, 1038, 1069, 1115, 1116, 1119, 1123 , 1130]
    #chan_list_airs = [1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082, 1088, 1090, 1092, 1104, 1111, 1115, 1116, 1119, 1120, 1123]
    #chan_list_iasi = [1479, 1509, 1513, 1521, 1536, 1574, 1578, 1579, 1585, 1587, 1626, 1639, 1643, 1652, 1658, 1671] 

    # Set the options for each Rttov instance:
    # - the path to the coefficient file must always be specified
    # - turn RTTOV interpolation on (because input pressure levels differ from
    #   coefficient file levels)
    # - set the verbose_wrapper flag to true so the wrapper provides more
    #   information
    # - enable solar simulations for iasi
    # - enable CO2 simulations for cris (the CO2 profiles are ignored for
    #   the iasi and airs simulations)
    # - enable the store_trans wrapper option for airs to provide access to
    #   RTTOV transmission structure

    iasiRttov.FileCoef = '{}/{}'.format(rttovPath,
                                           "rtcoef_rttov12/rttov9pred101L/rtcoef_metop_2_iasi_so2.H5")
#                                            "rtcoef_rttov12/rttov8pred101L/rtcoef_metop_2_iasi.H5")
    iasiRttov.Options.AddInterp = True
    iasiRttov.Options.AddSolar = True
    iasiRttov.Options.CO2Data = True
    iasiRttov.Options.OzoneData = True
    iasiRttov.Options.StoreTrans = True
    iasiRttov.Options.VerboseWrapper = True
    
    crisRttov.FileCoef = '{}/{}'.format(rttovPath,
                                        "rtcoef_rttov12/rttov9pred101L/rtcoef_jpss_0_cris_so2.H5")
    crisRttov.Options.AddInterp = True
    crisRttov.Options.AddSolar = True
    crisRttov.Options.CO2Data = True
    crisRttov.Options.OzoneData = True
    crisRttov.Options.StoreTrans = True
    crisRttov.Options.VerboseWrapper = True

    crisFsrRttov.FileCoef = '{}/{}'.format(rttovPath,
                                        "rtcoef_rttov12/rttov9pred101L/rtcoef_jpss_0_cris-fsr_so2.H5")
    crisFsrRttov.Options.AddInterp = True
    crisFsrRttov.Options.AddSolar = True
    crisFsrRttov.Options.CO2Data = True
    crisFsrRttov.Options.OzoneData = True
    crisFsrRttov.Options.StoreTrans = True
    crisFsrRttov.Options.VerboseWrapper = True

    airsRttov.FileCoef = '{}/{}'.format(rttovPath,
                                       "rtcoef_rttov12/rttov9pred101L/rtcoef_eos_2_airs_so2.H5")
    airsRttov.Options.AddInterp = True
    airsRttov.Options.StoreTrans = True
    airsRttov.Options.VerboseWrapper = True

    airsRttov.Options.AddInterp = True
    airsRttov.Options.AddSolar = True
    airsRttov.Options.CO2Data = True
    airsRttov.Options.OzoneData = False
    airsRttov.Options.StoreTrans = True
    airsRttov.Options.VerboseWrapper = True

    # Load the instruments: for cris and airs do not supply a channel list and
    # so read all channels
    try:
        iasiRttov.loadInst()
        crisRttov.loadInst()
        crisFsrRttov.loadInst()
        airsRttov.loadInst()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    # Associate the profiles with each Rttov instance
    iasiRttov.Profiles = myProfiles
    crisRttov.Profiles = myProfiles
    crisFsrRttov.Profiles = myProfiles
    airsRttov.Profiles = myProfiles
    try:
        #iasiRttov.printOptions()
        iasiRttov.runDirect()
        crisRttov.runDirect()
        crisFsrRttov.runDirect()
        airsRttov.runDirect()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)
    print("calc weighting functions?")
    iasi_wf    = calculateWeightingFunctions(chan_list_iasi, iasiRttov,    myProfiles)
    cris_wf    = calculateWeightingFunctions(chan_list_cris, crisRttov,    myProfiles)
    crisFsr_wf = calculateWeightingFunctions(chan_list_fsr,  crisFsrRttov, myProfiles)
    airs_wf    = calculateWeightingFunctions(chan_list_airs, airsRttov,    myProfiles)
    print("Plot weighting Functions?")    
    plotWeightingFunctions(chan_list_iasi, myProfiles, iasi_wf,    'iasi')
    plotWeightingFunctions(chan_list_cris, myProfiles, cris_wf,    'cris')
    plotWeightingFunctions(chan_list_fsr,  myProfiles, crisFsr_wf, 'cris-fsr')
    plotWeightingFunctions(chan_list_airs, myProfiles, airs_wf,    'airs')

