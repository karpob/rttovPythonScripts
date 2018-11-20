#!/usr/bin/env python2.7
# sorry, 2.7 for now, you can run 2to3 on rttov, but I need to get my scripts together to automate that a bit.
rttovPath = "/discover/nobackup/bkarpowi/rt/rttov12_gcc7.2/"
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

def calculateWeightingFunctions(chan_list, rttovInstance, myProfiles):
    wf = []
    nlevels = np.asarray(rttovInstance.TauLevels).shape[2] 
    for c in range(1,rttovInstance.Nchannels+1):
        if c in chan_list:
            num = rttovInstance.TauLevels[0, c-1, 1::] - rttovInstance.TauLevels[0, c-1, 0:nlevels-1]
            den = np.log(myProfiles.P[0, 0:nlevels-1]) - np.log(myProfiles.P[0, 1:nlevels])
            wf.append(num/den)
    return wf
 
def plotWeightingFunctions(chan_list, profiles, wf, instrument, wavenumbers):
    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = len(chan_list)
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    for f in wf:
        plt.plot( f[0:99], profiles.P[0,0:99] )
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    for i,c in enumerate(chan_list):
        legendList.append('{} {:4.3f} {:4.3f}'.format(c, wavenumbers[instrument][c-1], 10000.0/wavenumbers[instrument][c-1]))
    plt.legend(legendList, fontsize=6,  ncol=3)
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Weighting Function')
    plt.savefig(instrument+'_weighting_functions.png')
    
    if(len(chan_list)>30): matplotlib.rc('xtick', labelsize=6) 
    else: matplotlib.rc('xtick', labelsize=10) 
 

    plt.figure()
    plt.title('Weighting Function by Instrument Channel')
    wf = np.asarray(wf)
    #plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T, norm = LogNorm( vmin = wf[:,0:99].min(), vmax = 0.3 ))#vmax = wf[:,0:99].max() ) )
    plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T,vmin=wf[:,0:99].min(), vmax=0.3 )#vmax = wf[:,0:99].max() ) )
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
    fsrH5  = h5py.File('etc/cris-fsr_wavenumbers.h5','r')
    crisH5 = h5py.File('etc/cris_wavenumbers.h5','r')
    iasiH5 = h5py.File('etc/iasi_wavenumbers.h5','r')
    airsH5 = h5py.File('etc/airs_wavenumbers.h5','r')

    chan_list_fsr  = fsrH5['idxEcmwfOzoneInInstrument']
    chan_list_cris = crisH5['idxEcmwfOzoneInInstrument'] 
    chan_list_iasi = iasiH5['idxEcmwfOzoneInInstrument']
    chan_list_airs = airsH5['idxEcmwfOzoneInInstrument']
    wavenumbers = {}
    wavenumbers['cris-fsr'] = fsrH5['wavenumbers']
    wavenumbers['cris'] = crisH5['wavenumbers']
    wavenumbers['iasi'] = iasiH5['wavenumbers']
    wavenumbers['airs'] = airsH5['wavenumbers']
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
    plotWeightingFunctions(chan_list_iasi, myProfiles, iasi_wf,    'iasi', wavenumbers)
    plotWeightingFunctions(chan_list_cris, myProfiles, cris_wf,    'cris',wavenumbers)
    plotWeightingFunctions(chan_list_fsr,  myProfiles, crisFsr_wf, 'cris-fsr',wavenumbers)
    plotWeightingFunctions(chan_list_airs, myProfiles, airs_wf,    'airs',wavenumbers)

