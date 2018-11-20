
# Example of using the Rttov class to call RTTOV for multiple instruments
# with the emissivity and BRDF atlases

# Three Rttov instances are created representing three instruments
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

def calculateWeightingFunctions(chan_list, rttovInstance, myProfiles):
    
    nlevels = np.asarray(rttovInstance.TauLevels).shape[2]
    nprofiles = myProfiles.P.shape[0]
    wf = np.zeros([nprofiles, nlevels-1, len(chan_list)])
    for p in range(myProfiles.P.shape[0]): 
        iList = 0
        for c in range(0,rttovInstance.Nchannels):
            if c+1 in chan_list:
                num = rttovInstance.TauLevels[p, c, 1::] - rttovInstance.TauLevels[p, c, 0:nlevels-1]
                den = np.log(myProfiles.P[p, 0:nlevels-1]) - np.log(myProfiles.P[p, 1:nlevels])
                wf[p,:,iList] = num/den
                iList+=1
    return wf
 
def plotWeightingFunctions(chan_list, profiles, wf, instrument, profileNumber, wavenumbers):
    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = len(chan_list)
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    for i in range(wf.shape[1]):
        plt.plot( wf[0:99,i], profiles.P[int(profileNumber)-1,0:99] )
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.gca().set_xlim(0,0.5)
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    for i,c in enumerate(chan_list):
        legendList.append('{} {:4.3f} {:4.3f}'.format(c, wavenumbers[instrument][c-1], 10000.0/wavenumbers[instrument][c-1]))
    plt.legend(legendList, fontsize=6,  ncol=3)
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Weighting Function')
    plt.savefig(instrument+'_profile_'+profileNumber+'_weighting_functions.png')
    
    if(len(chan_list)>30): matplotlib.rc('xtick', labelsize=6) 
    else: matplotlib.rc('xtick', labelsize=10) 
    plt.close()

    plt.figure()
    plt.title('Weighting Function by Instrument Channel')
    wf = np.asarray(wf)
    #plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T, norm = LogNorm( vmin = wf[:,0:99].min(), vmax = 0.3 ))#vmax = wf[:,0:99].max() ) )
    plt.pcolor(np.arange(wf.shape[1]), myProfiles.P[int(profileNumber)-1,0:99], wf[0:99,:],vmin=wf[0:99,:].min(), vmax=0.5 )#vmax = wf[:,0:99].max() ) )
    plt.colorbar()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Instrument Channel')
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.xticks(np.arange(len(chan_list)), chan_list, rotation='vertical')
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    plt.tight_layout() 
    plt.savefig(instrument+'_profile_'+profileNumber+'_weighting_functions_pcolor.png')
    plt.close()
    print(profileNumber, myProfiles.O3[int(profileNumber)-1,:].sum())
rttov_installdir = '../'

def plotJacobians(chan_list, profiles, jacobians, instrument, profileNumber, wavenumbers, ofWhat):
    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = len(chan_list)
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    jacSub = []
    for i in range(jacobians.shape[0]):
        if(i+1 in chan_list):
            plt.plot( jacobians[i,:], profiles.P[int(profileNumber)-1,:] )
            jacSub.append(jacobians[i,:])
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    #plt.gca().set_xlim(0,0.5)
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    for i,c in enumerate(chan_list):
        legendList.append('{} {:4.3f} {:4.3f}'.format(c, wavenumbers[instrument][c-1], 10000.0/wavenumbers[instrument][c-1]))
    plt.legend(legendList, fontsize=6,  ncol=3)
    plt.ylabel('Pressure [hPa]')
    if(ofWhat.lower() =='temperature'): plt.xlabel('Jacobian [K/K]')
    else: plt.xlabel('Jacobian [K/ppmv]')
    plt.savefig(instrument+'_profile_'+profileNumber+'_jacobians_'+ofWhat+'.png')
    
    if(len(chan_list)>30): matplotlib.rc('xtick', labelsize=6) 
    else: matplotlib.rc('xtick', labelsize=10) 
    plt.close()

    plt.figure()
    if(ofWhat.lower() == 'temperature'): plt.title(ofWhat.capitalize()+' Jacobian by Instrument Channel [K/K]')
    else: plt.title(ofWhat.capitalize()+' Jacobian by Instrument Channel [K/ppmv]')
    jacSub = np.asarray(jacSub)
    print(profileNumber)
    #plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T, norm = LogNorm( vmin = wf[:,0:99].min(), vmax = 0.3 ))#vmax = wf[:,0:99].max() ) )
    plt.pcolor(np.arange(len(chan_list)), myProfiles.P[int(profileNumber)-1,:], jacSub[:,:].T ) 
    plt.colorbar()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Instrument Channel')
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.xticks(np.arange(len(chan_list)), chan_list, rotation='vertical')
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    plt.tight_layout() 
    plt.savefig(instrument+'_profile_'+profileNumber+'_jacobians_'+ofWhat+'_pcolor.png')
    plt.close()
    print(profileNumber, myProfiles.O3[int(profileNumber)-1,:].sum())


def plotProfilesO3(profiles):

    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = profiles.P.shape[0]
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    for i in range(profiles.O3.shape[0]):
        plt.plot( profiles.O3[i,0:99], profiles.P[i,0:99] )
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    plt.legend(np.arange(1,7), fontsize=6,  ncol=1)
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Ozone [ppmv]')
    plt.savefig('ozoneProfiles.png')
    
    plt.close()
def plotProfilesT(profiles):

    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = profiles.P.shape[0]
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    for i in range(profiles.O3.shape[0]):
        plt.plot( profiles.T[i,0:99], profiles.P[i,0:99] )
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    plt.legend(np.arange(1,7), fontsize=6,  ncol=1)
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Temperature [K]')
    plt.savefig('temperatureProfiles.png')
    
    plt.close()



def readProfileH5(filename,nprofiles):
    h5 = h5py.File( filename )
    nprofiles = 6 
    groups = []
    for i in range(1,nprofiles+1):
        groups.append("{:04d}".format(i))
    P = np.zeros([nprofiles,101])
    T = np.zeros([nprofiles,101])
    Q = np.zeros([nprofiles,101])
    CO2 = np.zeros([nprofiles,101])
    O3 = np.zeros([nprofiles,101])
    for i,g in enumerate(groups):
        P[i,:] = np.asarray(h5['PROFILES'][g]['P'])
        Q[i,:] = np.asarray(h5['PROFILES'][g]['Q'])
        T[i,:] = np.asarray(h5['PROFILES'][g]['T'])
        CO2[i,:] = np.asarray(h5['PROFILES'][g]['CO2'])
        O3[i,:] = np.asarray(h5['PROFILES'][g]['O3'])
        GasUnits = int(np.asarray(h5['PROFILES'][g]['GAS_UNITS']))
    
    return P, T, Q, CO2, O3, GasUnits 
if __name__ == '__main__':

    # This example program simulates two profiles for each of three instruments
    # The example profile data are defined in example_data

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nprofiles = 6
    nlevels = 101
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)
    h5ProfileFilename = '../rttov_test/profile-datasets-hdf/standard101lev_allgas.H5'
    myProfiles.P, myProfiles.T, myProfiles.Q, myProfiles.CO2, myProfiles.O3, myProfiles.GasUnits = readProfileH5(h5ProfileFilename, nprofiles)
   
    # repeat first example parameter for all profiles
    myProfiles.Angles = ex.angles.transpose()[0,:]*np.ones([nprofiles,ex.angles.transpose()[0,:].shape[0]])
    myProfiles.S2m = ex.s2m.transpose()[0,:]*np.ones([nprofiles,ex.s2m.transpose()[0,:].shape[0]])
    myProfiles.Skin = ex.skin.transpose()[0,:]*np.ones([nprofiles,ex.skin.transpose()[0,:].shape[0]])
    myProfiles.SurfType = ex.surftype.transpose()[0,:]*np.ones([nprofiles,ex.surftype.transpose()[0,:].shape[0]])
    myProfiles.SurfGeom = ex.surfgeom.transpose()[0,:]*np.ones([nprofiles,ex.surfgeom.transpose()[0,:].shape[0]])
    myProfiles.DateTimes = ex.datetimes.transpose()[0,:]*np.ones([nprofiles,ex.datetimes.transpose()[0,:].shape[0]])

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
 
    chan_list_fsr_nucaps  = fsrH5['idxNucapsOzoneInInstrument']
    chan_list_cris_nucaps = crisH5['idxNucapsOzoneInInstrument'] 
    chan_list_iasi_nucaps = iasiH5['idxNucapsOzoneInInstrument']
    chan_list_airs_nucaps = airsH5['idxNucapsOzoneInInstrument']

    chan_list_fsr_bufr  = fsrH5['idxBufrSubset']
    chan_list_cris_bufr = crisH5['idxBufrSubset'] 
    chan_list_iasi_bufr = iasiH5['idxBufrSubset']
    chan_list_airs_bufr = airsH5['idxBufrSubset']

    chan_list_fsr = []
    chan_list_cris = []
    chan_list_iasi = []
    chan_list_airs = []

    # go through and include only channels in buffer subset.
    for c in chan_list_fsr_nucaps:
        if(c in chan_list_fsr_bufr):
            chan_list_fsr.append(c)

    for c in chan_list_cris_nucaps:
        if(c in chan_list_cris_bufr):
            chan_list_cris.append(c)

    for c in chan_list_iasi_nucaps:
        if(c in chan_list_iasi_bufr):
            chan_list_iasi.append(c)

    for c in chan_list_airs_nucaps:
        if(c in chan_list_airs_bufr):
            chan_list_airs.append(c)


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

    iasiRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                           "rtcoef_rttov12/rttov9pred101L/rtcoef_metop_2_iasi_so2.H5")
#                                            "rtcoef_rttov12/rttov8pred101L/rtcoef_metop_2_iasi.H5")
    iasiRttov.Options.AddInterp = True
    iasiRttov.Options.AddSolar = True
    iasiRttov.Options.CO2Data = True
    iasiRttov.Options.OzoneData = True
    iasiRttov.Options.StoreTrans = True
    iasiRttov.Options.VerboseWrapper = True
    
    crisRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                        "rtcoef_rttov12/rttov9pred101L/rtcoef_jpss_0_cris_so2.H5")
    crisRttov.Options.AddInterp = True
    crisRttov.Options.AddSolar = True
    crisRttov.Options.CO2Data = True
    crisRttov.Options.OzoneData = True
    crisRttov.Options.StoreTrans = True
    crisRttov.Options.VerboseWrapper = True

    crisFsrRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                        "rtcoef_rttov12/rttov9pred101L/rtcoef_jpss_0_cris-fsr_so2.H5")
    crisFsrRttov.Options.AddInterp = True
    crisFsrRttov.Options.AddSolar = True
    crisFsrRttov.Options.CO2Data = True
    crisFsrRttov.Options.OzoneData = True
    crisFsrRttov.Options.StoreTrans = True
    crisFsrRttov.Options.VerboseWrapper = True

    airsRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                       "rtcoef_rttov12/rttov9pred101L/rtcoef_eos_2_airs_so2.H5")
    airsRttov.Options.AddInterp = True
    airsRttov.Options.StoreTrans = True
    airsRttov.Options.VerboseWrapper = True

    airsRttov.Options.AddInterp = True
    airsRttov.Options.AddSolar = True
    airsRttov.Options.CO2Data = True
    airsRttov.Options.OzoneData = True
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
        #iasiRttov.runDirect()
        #crisRttov.runDirect()
        #crisFsrRttov.runDirect()
        #airsRttov.runDirect()
        iasiRttov.runK()
        crisRttov.runK()
        crisFsrRttov.runK()
        airsRttov.runK()
        airsRttov.printOptions()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)
    print("calc weighting functions.")
    iasi_wf    = calculateWeightingFunctions(chan_list_iasi, iasiRttov,    myProfiles)
    cris_wf    = calculateWeightingFunctions(chan_list_cris, crisRttov,    myProfiles)
    crisFsr_wf = calculateWeightingFunctions(chan_list_fsr,  crisFsrRttov, myProfiles)
    airs_wf    = calculateWeightingFunctions(chan_list_airs, airsRttov,    myProfiles)
    print("Plot weighting Functions.")

    for i in range(nprofiles):

        print('Plotting IASI')    
        plotWeightingFunctions(chan_list_iasi, myProfiles, iasi_wf[i,:,:],    'iasi','{:04d}'.format(i+1), wavenumbers)
        plotJacobians(chan_list_iasi, myProfiles, iasiRttov.O3K[i,:,:],    'iasi', '{:04d}'.format(i+1), wavenumbers,'ozone')
        plotJacobians(chan_list_iasi, myProfiles, iasiRttov.TK[i,:,:],    'iasi', '{:04d}'.format(i+1), wavenumbers,'temperature')
        plotJacobians(chan_list_iasi, myProfiles, iasiRttov.QK[i,:,:],    'iasi', '{:04d}'.format(i+1), wavenumbers,'water')
 

        print('Plotting CrIS')
        plotWeightingFunctions(chan_list_cris, myProfiles, cris_wf[i,:,:],    'cris','{:04d}'.format(i+1),wavenumbers)
        plotJacobians(chan_list_cris, myProfiles, crisRttov.O3K[i,:,:],    'cris', '{:04d}'.format(i+1), wavenumbers,'ozone')
        plotJacobians(chan_list_cris, myProfiles, crisRttov.TK[i,:,:],    'cris', '{:04d}'.format(i+1), wavenumbers,'temperature')
        plotJacobians(chan_list_cris, myProfiles, crisRttov.QK[i,:,:],    'cris', '{:04d}'.format(i+1), wavenumbers,'water')
 

        print('Plotting CrIS-FSR')
        plotWeightingFunctions(chan_list_fsr,  myProfiles, crisFsr_wf[i,:,:], 'cris-fsr','{:04d}'.format(i+1),wavenumbers)
        plotJacobians(chan_list_fsr, myProfiles, crisFsrRttov.O3K[i,:,:],    'cris-fsr', '{:04d}'.format(i+1), wavenumbers,'ozone')
        plotJacobians(chan_list_fsr, myProfiles, crisFsrRttov.TK[i,:,:],    'cris-fsr', '{:04d}'.format(i+1), wavenumbers,'temperature')
        plotJacobians(chan_list_fsr, myProfiles, crisFsrRttov.QK[i,:,:],    'cris-fsr', '{:04d}'.format(i+1), wavenumbers,'water')
 
        print('Plotting AIRS')
        plotWeightingFunctions(chan_list_airs, myProfiles, airs_wf[i,:,:],    'airs', '{:04d}'.format(i+1),wavenumbers)
        plotJacobians(chan_list_airs, myProfiles, airsRttov.O3K[i,:,:],    'airs', '{:04d}'.format(i+1), wavenumbers,'ozone')
        plotJacobians(chan_list_airs, myProfiles, airsRttov.TK[i,:,:],    'airs', '{:04d}'.format(i+1), wavenumbers,'temperature')
        plotJacobians(chan_list_airs, myProfiles, airsRttov.QK[i,:,:],    'airs', '{:04d}'.format(i+1), wavenumbers,'water')

    plotProfilesO3(myProfiles)   
    plotProfilesT(myProfiles)   
