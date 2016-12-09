#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

# import h5py as h5
import numpy as np
import pylab as plt
from random import sample, seed
from os.path import getsize as getFileSize
import matplotlib.cm as cm

import re

# ================================================================================
# Basic variables
# ================================================================================

# Set up some basic attributes of the run

whichsimulation = 0
whichimf = 1        # 0=Slapeter; 1=Chabrier
dilute = 7500       # Number of galaxies to plot in scatter plots
sSFRcut = -11.0     # Divide quiescent from star forming galaxies (when plotmags=0)
MERGER_NUM = 10
ContinuousAccretionOn = 1

Grav_const = 6.67e-8
R = 8.31e7

matplotlib.rcdefaults()

plt.rc('axes', color_cycle=[
    'k',
    'b',
    'r',
    'g',
    'm',
    '0.5',
    ], labelsize='x-large')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

OutputDir = '' # set in main below

OutputFormat = '.png'
TRANSPARENT = False

OutputList = []


class Results:

    """ The following methods of this class generate the figures and plot them.
    """

    def __init__(self):
        """Here we set up some of the variables which will be global to this
        class."""

        if whichsimulation == 0:    # Mini-Millennium
          self.Hubble_h = 0.73
          self.BoxSize = 62.5       # Mpc/h
          self.MaxTreeFiles = 8     # FilesPerSnapshot
          self.OmegaM = 0.25
          self.OmegaL = 0.75
          self.z = 6.197
          self.aList = "../input/treefiles/millennium/millennium.a_list"

        elif whichsimulation == 1:  # Full Millennium
          self.Hubble_h = 0.73
          self.BoxSize = 500        # Mpc/h
          self.MaxTreeFiles = 512   # FilesPerSnapshot
          self.OmegaM = 0.25
          self.OmegaL = 0.75
          self.z = 6.197
          self.aList = "/lustre/projects/p014_swin/raw_data/millennium/full/millennium.a_list"

        elif whichsimulation == 2:  # Full Bolshoi
          self.Hubble_h = 0.7
          self.BoxSize = 250        # Mpc/h
          self.MaxTreeFiles = 125  # FilesPerSnapshot
          self.OmegaM = 0.27
          self.OmegaL = 0.73
          self.z = 6.197
          self.aList = "/lustre/projects/p004_swin/msinha/tao/data_products/input/bolshoi-planck/lhalotree/run1/bolshoi-planck.a_list"

        else:
          print "Please pick a valid simulation!"
          exit(1)



    def read_gals(self, model_name, first_file, last_file):

        # The input galaxy structure:
        Galdesc_full = [
            ('SnapNum'                      , np.int32),
            ('Type'                         , np.int32),
            ('GalaxyIndex'                  , np.int64),
            ('CentralGalaxyIndex'           , np.int64),
            ('SAGEHaloIndex'                , np.int32),
            ('SAGETreeIndex'                , np.int32),
            ('SimulationHaloIndex'          , np.int64),
            ('mergeType'                    , np.int32),
            ('mergeIntoID'                  , np.int32),
            ('mergeIntoSnapNum'             , np.int32),
            ('dT'                           , np.float32),
            ('Pos'                          , (np.float32, 3)),
            ('Vel'                          , (np.float32, 3)),
            ('Spin'                         , (np.float32, 3)),
            ('Len'                          , np.int32),
            ('Mvir'                         , np.float32),
            ('CentralMvir'                  , np.float32),
            ('Rvir'                         , np.float32),
            ('Vvir'                         , np.float32),
            ('Vmax'                         , np.float32),
            ('VelDisp'                      , np.float32),
            ('ColdGas'                      , np.float32),
            ('StellarMass'                  , np.float32),
            ('BulgeMass'                    , np.float32),
            ('HotGas'                       , np.float32),
            ('EjectedMass'                  , np.float32),
            ('BlackHoleMass'                , np.float32),
            ('IntraClusterStars'            , np.float32),
            ('MetalsColdGas'                , np.float32),
            ('MetalsStellarMass'            , np.float32),
            ('MetalsBulgeMass'              , np.float32),
            ('MetalsHotGas'                 , np.float32),
            ('MetalsEjectedMass'            , np.float32),
            ('MetalsIntraClusterStars'      , np.float32),
            ('SfrDisk'                      , np.float32),
            ('SfrBulge'                     , np.float32),
            ('SfrDiskZ'                     , np.float32),
            ('SfrBulgeZ'                    , np.float32),
            ('DiskRadius'                   , np.float32),
            ('Cooling'                      , np.float32),
            ('Heating'                      , np.float32),
            ('QuasarModeBHaccretionMass'    , np.float32),
            ('TimeOfLastMajorMerger'        , np.float32),
            ('TimeOfLastMinorMerger'        , np.float32),
            ('OutflowRate'                  , np.float32),
            ('infallMvir'                   , np.float32),
            ('infallVvir'                   , np.float32),
            ('infallVmax'                   , np.float32),
            ('QSOBHaccrete'                 , (np.float32, MERGER_NUM)),
            ('QSOmergeAge'                  , (np.float32, MERGER_NUM)),
            ('QSOmergeTime'                 , (np.float32, MERGER_NUM)),
            ('QSOBH'                        , (np.float32, MERGER_NUM)),
            ('MergSnap'                     , np.int32),
            ('QSOmergSnap'                  , (np.int32, MERGER_NUM)),
            ('MergNum'                      , np.int32),
            ('QSOmergeType'                 , (np.int32, MERGER_NUM)),
            ('QSOBHaccretionRate'           , np.float32),
            ('QSOBHaccretionMass'           , np.float32),
            ('QSOLuminosity'                , np.float32),
            ('ColdGasToAccrete'             , np.float32)
            ]
        names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
        formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
        Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)


        # Initialize variables.
        TotNTrees = 0
        TotNGals = 0
        FileIndexRanges = []

        print "Determining array storage requirements."

        # Read each file and determine the total number of galaxies to be read in
        goodfiles = 0
        for fnr in xrange(first_file,last_file+1):
            fname = model_name+'_'+str(fnr)  # Complete filename

            if not os.path.isfile(fname):
              # print "File\t%s  \tdoes not exist!  Skipping..." % (fname)
              continue

            if getFileSize(fname) == 0:
                print "File\t%s  \tis empty!  Skipping..." % (fname)
                continue

            fin = open(fname, 'rb')  # Open the file
            Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
            NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
            TotNTrees = TotNTrees + Ntrees  # Update total sim trees number
            TotNGals = TotNGals + NtotGals  # Update total sim gals number
            goodfiles = goodfiles + 1  # Update number of files read for volume calculation
            fin.close()

        print
        print "Input files contain:\t%d trees ;\t%d galaxies ." % (TotNTrees, TotNGals)
        print

        # Initialize the storage array
        G = np.empty(TotNGals, dtype=Galdesc)

        offset = 0  # Offset index for storage array

        # Open each file in turn and read in the preamble variables and structure.
        print "Reading in files."
        for fnr in xrange(first_file,last_file+1):
            fname = model_name+'_'+str(fnr)  # Complete filename

            if not os.path.isfile(fname):
              continue

            if getFileSize(fname) == 0:
                continue

            fin = open(fname, 'rb')  # Open the file
            Ntrees = np.fromfile(fin, np.dtype(np.int32), 1)  # Read number of trees in file
            NtotGals = np.fromfile(fin, np.dtype(np.int32), 1)[0]  # Read number of gals in file.
            GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
            print ":   Reading N=", NtotGals, "   \tgalaxies from file: ", fname
            GG = np.fromfile(fin, Galdesc, NtotGals)  # Read in the galaxy structures

            FileIndexRanges.append((offset,offset+NtotGals))

            # Slice the file array into the global array
            # N.B. the copy() part is required otherwise we simply point to
            # the GG data which changes from file to file
            # NOTE THE WAY PYTHON WORKS WITH THESE INDICES!
            G[offset:offset+NtotGals]=GG[0:NtotGals].copy()

            del(GG)
            offset = offset + NtotGals  # Update the offset position for the global array

            fin.close()  # Close the file


        print
        print "Total galaxies considered:", TotNGals

        # Convert the Galaxy array into a recarray
        G = G.view(np.recarray)

        w = np.where(G.StellarMass > 1.0)[0]
        print "Galaxies more massive than 10^10Msun/h:", len(w)
        print

        # Calculate the volume given the first_file and last_file
        self.volume = self.BoxSize**3.0 * goodfiles / self.MaxTreeFiles

        return G


# --------------------------------------------------------

    def QSOBHaccrete_distribution(self, G):

        print 'Plotting the QSO BHaccrete distribution of all galaxies'

        seed(2222)

        plt.figure()
        ax = plt.subplot(111)  # 1 plot on the figure

        w = np.where((G.Mvir  > 0.0) & (G.StellarMass > 0.))[0]
        # if(len(w) > dilute): w = sample(w, dilute)

        thisPlot = plt.scatter(np.log10(G.Mvir[w] / self.Hubble_h * 1.e10), np.log10(G.QSOBHaccretionRate[w] / self.Hubble_h), marker='o', s=20., alpha=0.8, label='Stars')
        thisPlot.set_clim(vmin=0,vmax=1)

        plt.ylabel(r'$\log \dot{M_{BH}}$    $[\mathrm{M}_{\odot}]$')  # Set the y...
        plt.xlabel(r'$\log M_{vir}$    $[\mathrm{M}_{\odot}]$')  # and the x-axis labels

        outputFile = OutputDir + '16.QSOBHaccreteDistribution_z' + str(G.z) + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)

# --------------------------------------------------------

    def QSOBHaccrete_function(self, G):

        print 'Plotting the QSO BHaccrete function of all galaxies'

        seed(2222)

        plt.figure()


        w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccretionRate > 0.))[0]

        BHaccrete = np.zeros(len(w))

        BHaccrete = BHaccrete + G.QSOBHaccretionRate[w] / self.Hubble_h

        ax = plt.subplot(111)  # 1 plot on the figure

        mi = -8#round(np.min(np.log10(BHaccrete)))
        ma = 4#round(np.max(np.log10(BHaccrete)))

        binwidth = 0.25
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(np.log10(BHaccrete), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        ax.set_yscale('log')
        plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='los-velocity')

        plt.ylabel(r'$\mathrm{number}\ \mathrm{density}$    $[\mathrm{Mpc}^{-3}]$')  # Set the y...
        plt.xlabel(r'$\log \dot{M_{BH}}$    $[\mathrm{M}_{\odot}]$')  # and the x-axis labels

        outputFile = OutputDir + '17.QSOBHaccreteFunction_z' + str(G.z) + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)

# --------------------------------------------------------

    def QSOluminosity_function(self, G):

        print 'Plotting the QSO luminosity distribution of all galaxies'

        seed(2222)

        Msun = 2.e33
        c = 3.e10
        yr = 365.*24.*3600.
        km_cm = 1.e5
        pc_cm = 3.086e18
        Mpc_cm = 3.086e24
        Lsol = 3.9e33

        efficiency = 0.1
        gamma = 5./3.
        mu = 1.22
        Temp = 6.e4

        H = 100.
        AgeRedshift = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM/(1.+G.z)**3)**0.5)
        AgeUniverse = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM)**0.5)

        plt.figure()

        w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccretionRate > 0.))[0]

        # -----------------------------------------------------------------------------
        ax = plt.subplot(111)  # 1 plot on the figure

        mi = 43
        ma = 50#round(np.max(np.log10(luminosity)))
        binwidth = 0.25
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(np.log10(G.QSOLuminosity[w])+np.log10(Msun)-np.log10(yr), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        ax.set_yscale('log')
        if(ContinuousAccretionOn == 1):
            plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='Eddington limit')
        elif(ContinuousAccretionOn == 2):
            plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='$f_{Edd} (z)$')
        elif(ContinuousAccretionOn == 3):
            plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='Eddington-limited growth & decay')

        # -----------------------------------------------------------------------------
        # Observational Fits

        xaxes = np.arange(mi,ma,0.1)
        if(G.z < 0.5):
            phi = -5.45
            Lstar = 11.94 +np.log10(3.9e33)
            gamma1 = 0.868
            gamma2 = 1.97
            dphi = 0.28
            dLstar = 0.21
            dgamma1 = 0.050
            dgamma2 = 0.17
        elif(G.z < 1.5):
            phi = -4.63
            Lstar = 12.59 +np.log10(3.9e33)
            gamma1 = 0.412
            gamma2 = 2.23
            dphi = 0.15
            dLstar = 0.11
            dgamma1 = 0.122
            dgamma2 = 0.15
        elif(G.z < 2.5):
            phi = -4.83
            Lstar = 13.10 +np.log10(3.9e33)
            gamma1 = 0.320
            gamma2 = 2.39
            dphi = 0.05
            dLstar = 0.04
            dgamma1 = 0.046
            dgamma2 = 0.07
        elif(G.z < 3.5):
            phi = -5.23
            Lstar = 13.17 +np.log10(3.9e33)
            gamma1 = 0.395
            gamma2 = 2.10
            dphi = 0.12
            dLstar = 0.10
            dgamma1 = 0.060
            dgamma2 = 0.12
        elif(G.z < 4.5):
            phi = -4.66
            Lstar = 12.39 +np.log10(3.9e33)
            gamma1 = 0.254
            gamma2 = 1.69
            dphi = 0.37
            dLstar = 0.32
            dgamma1 = 0.736
            dgamma2 = 0.18
        elif(G.z < 5.5):
            phi = -5.38
            Lstar = 12.46 +np.log10(3.9e33)
            gamma1 = 0.497
            gamma2 = 1.57
            dphi = 1.19
            dLstar = 1.10
            dgamma1 = 0.458
            dgamma2 = 0.41
        elif(G.z < 6.5):
            phi = -5.13
            Lstar = 11.00 +np.log10(3.9e33)
            gamma1 = 0.
            gamma2 = 1.11
            dphi = 0.38
            dLstar = 0.
            dgamma1 = 0.
            dgamma2 = 0.13
        QLF = 10.**phi / ((10.**(xaxes-Lstar))**gamma1 + (10.**(xaxes-Lstar))**gamma2)
        QLF_min = 10.**(phi-dphi) / ((10.**(xaxes-(Lstar-dLstar)))**(gamma1-dgamma1) + (10.**(xaxes-(Lstar-dLstar)))**(gamma2-dgamma2))
        QLF_max = 10.**(phi+dphi) / ((10.**(xaxes-(Lstar+dLstar)))**(gamma1+dgamma1) + (10.**(xaxes-(Lstar+dLstar)))**(gamma2+dgamma2))

        plt.plot(xaxes, QLF, color='grey', label='Hopkins 2007')
        ax.fill_between(xaxes, QLF_min, QLF_max, alpha=0.1, facecolor='grey')
        # -----------------------------------------------------------------------------


        plt.ylabel(r'$\mathrm{number}\ \mathrm{density}$    $[\mathrm{Mpc}^{-3}]$')  # Set the y...
        plt.xlabel(r'$\log L$    $[\mathrm{erg}\ \mathrm{s}^{-1}]$')  # and the x-axis labels

        plt.axis((mi,ma,1.e-9,1.e-3))
        leg = plt.legend(loc='upper right', prop={'size':12})
        leg.draw_frame(False)  # Don't want a box frame

        outputFile = OutputDir + '19.QSOLuminosityFunction_z' + str(G.z) + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


        print 'Plotting the QSO luminosity function in terms of M1450'

        seed(2222)

        alpha1 = -1.76
        alpha2 = -0.44
        lambdaBreak = 1050.
        nuBreak = c/lambdaBreak*1.e8
        nu1450 = c/1450.*1.e8

        tmp = (1. + alpha1) * (1. + alpha2) / (alpha1 - alpha2)
        M1450luminosity = G.QSOLuminosity[w] * tmp * (nu1450/nuBreak)**alpha2 / nuBreak * Msun / yr

        zBins = 100000
        thisSum=0.
        for i in range(zBins):
            red = self.z/np.float64(zBins)*np.float64(i)
            delta_z = self.z/np.float64(zBins)
            thisSum = thisSum + c/(self.Hubble_h*km_cm*1.e2/Mpc_cm)/(self.OmegaM*(1.+red)**3+self.OmegaL)**0.5*delta_z
        thisSum = (1.+self.z)*thisSum
        lumDistance = thisSum

        MUV1450 = -2.5*np.log10(M1450luminosity/(4.*np.pi*lumDistance**2))-48.6-5.*np.log10(lumDistance/(10.*pc_cm))

        print np.max(M1450luminosity)
        print np.min(MUV1450)

        if(G.z > 5.5):
            Jiang16_x,Jiang16_y,Jiang16_ymax,Jiang16_ymin,Jiang16_xmin,Jiang16_xmax = np.loadtxt('/lustre/projects/p004_swin/ahutter/QSO_sage/observations/Jiang2016_z6.dat', unpack='True', usecols=(0,1,2,3,4,5))
            Kash15_x,Kash15_y,Kash15_ymax,Kash15_ymin = np.loadtxt('/lustre/projects/p004_swin/ahutter/QSO_sage/observations/Kashikawa2015_z6.dat', unpack='True', usecols=(0,1,2,3))
            Giall15_x,Giall15_y,Giall15_ymax,Giall15_ymin = np.loadtxt('/lustre/projects/p004_swin/ahutter/QSO_sage/observations/Giallongo2015_z6.dat', unpack='True', usecols=(0,1,2,3))
        elif(G.z > 4.5):
            Giall15_x,Giall15_y,Giall15_ymax,Giall15_ymin = np.loadtxt('/lustre/projects/p004_swin/ahutter/QSO_sage/observations/Giallongo2015_z5.dat', unpack='True', usecols=(0,1,2,3))
            Greer15_x,Greer15_y,Greer15_ymax,Greer15_ymin = np.loadtxt('/lustre/projects/p004_swin/ahutter/QSO_sage/observations/McGreer2013_z5.dat', unpack='True', usecols=(0,1,2,3))

        ax = plt.subplot(111)  # 1 plot on the figure

        mi = -29.
        ma = -18.
        binwidth = 0.25
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(MUV1450, range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        ax.set_yscale('log')
        if(ContinuousAccretionOn == 1):
            plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='Eddington limit')
        elif(ContinuousAccretionOn == 2):
            plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='$f_{Edd} (z)$')
        elif(ContinuousAccretionOn == 3):
            plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='Eddington-limited growth & decay')

        if(G.z > 5.5):
            plt.errorbar(Jiang16_x, Jiang16_y*1.e-9, yerr=[Jiang16_y*1.e-9-Jiang16_ymin*1.e-9, Jiang16_ymax*1.e-9-Jiang16_y*1.e-9], xerr=[Jiang16_x-Jiang16_xmin, Jiang16_xmax-Jiang16_x], fmt='o', color='red', label='Jiang 2016')
            plt.errorbar(Kash15_x, Kash15_y, yerr=[Kash15_y-Kash15_ymin, Kash15_ymax-Kash15_y], fmt='o', color='blue', label='Kashikawa 2015')
            plt.errorbar(Giall15_x, 10.**Giall15_y, yerr=[10.**Giall15_y-10.**Giall15_ymin, 10.**Giall15_ymax-10.**Giall15_y], fmt='o', color='green', label='Giallongo 2015')
        elif(G.z > 4.5):
            plt.errorbar(Giall15_x, 10.**Giall15_y, yerr=[10.**Giall15_y-10.**Giall15_ymin, 10.**Giall15_ymax-10.**Giall15_y], fmt='o', color='green', label='Giallongo 2015')
            plt.errorbar(Greer15_x, Greer15_y, yerr=[Greer15_y-Greer15_ymin, Greer15_ymax-Greer15_y], fmt='o', color='blue', label='McGreer 2013')

        plt.ylabel(r'$\mathrm{number}\ \mathrm{density}$    $[\mathrm{Mpc}^{-3}]$')  # Set the y...
        plt.xlabel(r'$\mathrm{M}_{1450}$')  # and the x-axis labels

        plt.axis((mi,ma,1.e-11,1.e-4))
        leg = plt.legend(loc='upper left', prop={'size':12})
        leg.draw_frame(False)  # Don't want a box frame

        outputFile = OutputDir + '19a.QSOLuminosityFunction1450_z' + str(G.z) + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)

# # --------------------------------------------------------
#
#     def BHaccretion_time(self, G):
#
#         print 'Plotting the accretion onto BH over time for all galaxies'
#
#         seed(2222)
#
#         H = 100.
#         AgeRedshift = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM/(1.+self.z)**3)**0.5)
#         AgeUniverse = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM)**0.5)
#
#         a = np.loadtxt(self.aList, unpack='True')
#
#         z_snap = 1./a-1.
#
#         plt.figure()
#         ax = plt.subplot(111)  # 1 plot on the figure
#
#         BH = np.empty(0)
#         BHaccrete = np.empty(0)
#         redshiftmerger = np.empty(0)
#
#         for i in range(MERGER_NUM):
#             w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,i] > 0.))[0]
#             if(len(w) > dilute): w = sample(w, dilute)
#
#             time = AgeUniverse - G.QSOmergeAge[w,i]
#             redshift = (self.OmegaL/self.OmegaM/(np.sinh(1.5*self.OmegaL**0.5*H*time))**2)**(1./3.)-1.
#
#             BH = np.append(BH, G.BlackHoleMass[w] / self.Hubble_h)
#             BHaccrete = np.append(BHaccrete, G.QSOBHaccrete[w,i] / self.Hubble_h)
#             redshiftmerger = np.append(redshiftmerger, redshift)
#
#         BH_index = np.int32((np.log10(BH) - np.log10(np.min(BH))) * 2)#/ (round(np.log10(np.min(BH))) - round(np.log10(np.max(BH))))
#         redshift_index = np.int32((redshiftmerger - round(G.z)) * 2)
#
#         BH_index_max = np.max(BH_index) + 1
#         redshift_index_max = np.max(redshift_index) + 1
#         redshift_plot = np.arange(redshift_index_max)*0.5 + round(G.z)
#
#         hist_median = np.zeros((redshift_index_max, BH_index_max))
#         hist_mean = np.zeros((redshift_index_max, BH_index_max))
#         hist_std = np.zeros((redshift_index_max, BH_index_max))
#         hist_25 = np.zeros((redshift_index_max, BH_index_max))
#         hist_75 = np.zeros((redshift_index_max, BH_index_max))
#
#         i_min = 8
#         i_max = 13
#
#         cmap = plt.get_cmap('jet_r')
#         colors = [cmap(i) for i in np.linspace(0, 1, 6)]
#
#         for i in range(BH_index_max):
#             w = np.where(BH_index == i)
#             for j in range(redshift_index_max):
#                 v = np.where((BH_index == i) & (redshift_index == j))[0]
#                 if(len(v)>0):
#                     hist_median[j,i] = np.median(BHaccrete[v]/BH[v])
#                     hist_mean[j,i] = np.mean(BHaccrete[v]/BH[v])
#                     hist_std[j,i] = np.std(BHaccrete[v]/BH[v])
#                     hist_25[j,i] = np.percentile(BHaccrete[v]/BH[v],25)
#                     hist_75[j,i] = np.percentile(BHaccrete[v]/BH[v],75)
#                 else:
#                     hist_std[j,i] = 1.e-4
#                     hist_25[j,i] = 1.e-4
#                     hist_75[j,i] = 1.e-4
#
#             if(i >= i_min and i <= i_max):
#                 k = i*0.5 + np.log10(np.min(BH))+10.
#                 plt.plot(redshift_plot, hist_median[:,i], c=colors[i-9],  label='$M=10^{k}$'.format(k=k))
#                 ax.fill_between(redshift_plot, hist_25[:,i], hist_75[:,i], alpha=0.1, facecolor=colors[i-9])
#
#         plt.axis([G.z, 12, 0.01, 1.])
#         ax.set_yscale('log')
#         plt.ylabel(r'$\Delta\mathrm{M_{BH,Q}} / \mathrm{M_{BH}}$')  # Set the y...
#         plt.xlabel(r'$z_{\mathrm{merge}}$')  # and the x-axis labels
#         plt.legend(loc='best')
#
#         outputFile = OutputDir + '21.BHaccretion_time_z' + str(G.z) + OutputFormat
#         plt.savefig(outputFile)  # Save the figure
#         print 'Saved file to', outputFile
#         plt.close()
#
#         # Add this plot to our output list
#         OutputList.append(outputFile)

# =================================================================


#  'Main' section of code.  This if statement executes if the code is run from the
#   shell command line, i.e. with 'python allresults.py'

if __name__ == '__main__':

    from optparse import OptionParser
    import os

    parser = OptionParser()
    parser.add_option(
        '-d',
        '--dir_name',
        dest='DirName',
        default='./results/',
        help='input directory name (default: ./results/)',
        metavar='DIR',
        )
    parser.add_option(
        '-f',
        '--file_base',
        dest='FileName',
        default='model_z6.197',
        help='filename base (default: model_z0.000)',
        metavar='FILE',
        )
    parser.add_option(
        '-n',
        '--file_range',
        type='int',
        nargs=2,
        dest='FileRange',
        default=(0, 511),
        help='first and last filenumbers (default: 0 7)',
        metavar='FIRST LAST',
        )

    (opt, args) = parser.parse_args()

    if opt.DirName[-1] != '/':
        opt.DirName += '/'

    OutputDir = opt.DirName + 'plots/'

    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    res = Results()

    print 'Running allresults...'

    FirstFile = opt.FileRange[0]
    LastFile = opt.FileRange[1]

    fin_base = opt.DirName + opt.FileName

    G = res.read_gals(fin_base, FirstFile, LastFile)
    G.z = np.float32(re.findall('\d+\.\d+',opt.FileName)[0])

    print 'Anaylzing snapshot', G.SnapNum[0], "at redshift", G.z

    res.QSOBHaccrete_distribution(G)
    res.QSOBHaccrete_function(G)
    res.QSOluminosity_function(G)
    # res.BHaccretion_time(G)
