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
          self.BoxSize = 157.5#62.5       # Mpc/h
          self.MaxTreeFiles = 8     # FilesPerSnapshot
          self.OmegaM = 0.25
          self.OmegaL = 0.75
          self.z = 0.

        elif whichsimulation == 1:  # Full Millennium
          self.Hubble_h = 0.73
          self.BoxSize = 500        # Mpc/h
          self.MaxTreeFiles = 512   # FilesPerSnapshot
          self.OmegaM = 0.25
          self.OmegaL = 0.75
          self.z = 0.

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
            ('MergNum'                      , np.int32)
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

        for i in range(MERGER_NUM):
            w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,i] > 0.))[0]
            if(len(w) > dilute): w = sample(w, dilute)
            print i
            print G.QSOBHaccrete[w,i]
            t = i*np.ones(len(G.QSOBHaccrete[w,i]))/MERGER_NUM
            thisPlot = plt.scatter(np.log10(G.Mvir[w]), np.log10(G.QSOBHaccrete[w,i]), marker='o', s=20., c=t, cmap = cm.jet_r, alpha=0.8, label='Stars')
            thisPlot.set_clim(vmin=0,vmax=1)

        plt.ylabel(r'$\log \dot{M_{BH}}$    $[\mathrm{M}_{\odot} \mathrm{yr}^{-1}]$')  # Set the y...
        plt.xlabel(r'$\log M_{vir}$    $[\mathrm{M}_{\odot}]$')  # and the x-axis labels

        outputFile = OutputDir + '16.QSOBHaccreteDistribution' + OutputFormat
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


        w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,0] > 0.))[0]
        if(len(w) > dilute): w = sample(w, dilute)

        BHaccrete = np.zeros(len(w))

        for i in range(MERGER_NUM):
            BHaccrete = BHaccrete+G.QSOBHaccrete[w,i]

        print len(BHaccrete)

        ax = plt.subplot(111)  # 1 plot on the figure

        mi = round(np.min(np.log10(BHaccrete)))
        ma = round(np.max(np.log10(BHaccrete)))

        binwidth = 0.25
        NB = (ma - mi) / binwidth

        print BHaccrete
        (counts, binedges) = np.histogram(np.log10(BHaccrete), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        ax.set_yscale('log')
        plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='los-velocity')

        plt.ylabel(r'$\mathrm{number}\ \mathrm{density}$    $[\mathrm{Mpc}^{-3}]$')  # Set the y...
        plt.xlabel(r'$\log \dot{M_{BH}}$    $[\mathrm{M}_{\odot} \mathrm{yr}^{-1}]$')  # and the x-axis labels

        outputFile = OutputDir + '17.QSOBHaccreteFunction' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)

# --------------------------------------------------------

    def QSO_distribution(self, G):

        print 'Plotting the QSO distribution of all galaxies'

        seed(2222)

        plt.figure()
        ax = plt.subplot(111)  # 1 plot on the figure

        for i in range(MERGER_NUM):
            w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,i] > 0.))[0]
            if(len(w) > dilute): w = sample(w, dilute)

            H = 100.#(self.Hubble_h*1.e7/3.086e24)

            AgeUniverse = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM)**0.5)
            time = AgeUniverse - G.QSOmergeAge[w,i]
            redshift = (self.OmegaL/self.OmegaM/(np.sinh(1.5*self.OmegaL**0.5*H*time))**2)**(1./3.)-1.

            t = i*np.ones(len(G.QSOBHaccrete[w,i]))/MERGER_NUM
            thisPlot = plt.scatter(redshift, np.log10(G.QSOBHaccrete[w,i]/G.QSOBHaccrete[w,0]), marker='o', s=10., c=t, cmap = cm.jet_r, alpha=0.8, label='Stars')
            thisPlot.set_clim(vmin=0,vmax=1)

        # ax.set_xscale('log')
        plt.axis([G.z, 10.0, -6., 0.5])
        plt.ylabel(r'$\log \dot{M}_{BH}^{(i)} / \dot{M}_{BH}^{(0)}$')  # Set the y...
        plt.xlabel(r'$z$')  # and the x-axis labels

        outputFile = OutputDir + '18.QSODistribution' + OutputFormat
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
        Mpc_cm = 3.086e24
        Lsol = 3.9e33

        efficiency = 0.1
        gamma = 5./3.
        mu = 1.22
        Temp = 6.e4

        H = 100.
        AgeRedshift = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM/(1.+G.z)**3)**0.5)
        AgeUniverse = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM)**0.5)

        print AgeRedshift, AgeUniverse

        plt.figure()

        w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,0] > 0.))[0]
        if(len(w) > dilute): w = sample(w, dilute)

        # MODEL 1 ---------------------------------------------------------------------
        rBondi = np.empty(len(w), dtype=np.float32)
        BHaccrete = np.empty(len(w), dtype=np.float64)
        mergeAge = np.empty(len(w))
        mergeTime = np.empty(len(w))
        luminosity = np.zeros(len(w), dtype=np.float64)

        # Compute Bondi radius
        cs = (gamma*R*Temp/mu)**0.5
        rBondi = Grav_const/cs**2*G.BlackHoleMass[w]*2.e33*1.e10

        # Compute accretion time
        tacc = 1.*rBondi/cs*(self.Hubble_h*km_cm/Mpc_cm)
        # tacc = 0.01*G.Rvir[w]*3.086e24/(G.Vvir[w]*1.e3)*(self.Hubble_h*km_cm/Mpc_cm)
        # Compute Eddington luminosity

        for i in range(MERGER_NUM):
            BH = np.float64(G.QSOBH[w,i] * 1.e10 / self.Hubble_h)
            BHaccrete = np.float64(G.QSOBHaccrete[w,i] * 1.e10 / self.Hubble_h)
            mergeAge = np.float64(G.QSOmergeAge[w,i])
            mergeTime = G.QSOmergeTime[w,i]

            time = np.where(BHaccrete>0. , mergeAge-(AgeUniverse-AgeRedshift), 0.)
            time_factor = np.where(BHaccrete>0., np.exp(-time/tacc), 0.)
            luminosity = luminosity + BHaccrete*time_factor

            print time[0], AgeUniverse, AgeRedshift, mergeAge[0], mergeAge[1], mergeAge[2]
        luminosity = efficiency*Msun/yr*c**2*luminosity*1.e-7
        # -----------------------------------------------------------------------------
        # MODEL 2 ---------------------------------------------------------------------
        luminosity2 = np.zeros(len(w), dtype=np.float64)

        tEdd = 0.45e9*yr*(self.Hubble_h*km_cm/Mpc_cm)
        tEdd_yr = 0.45e9

        print "--------------------------------------------"
        for i in range(MERGER_NUM):
            BH = np.float64(G.QSOBH[w,i] * 1.e10 / self.Hubble_h)
            BHaccrete = np.float64(G.QSOBHaccrete[w,i] * 1.e10 / self.Hubble_h)
            mergeAge = np.float64(G.QSOmergeAge[w,i])
            mergeTime = G.QSOmergeTime[w,i]

            print AgeRedshift
            time = np.where(BHaccrete>0. , mergeAge-(AgeUniverse-AgeRedshift), 0.)
            # if(i==0):
            #     for j in range(len(w)):
            #         print np.log10(BH[j]), np.log10(BHaccrete[j]), np.log10(time[j]), np.log10(mergeAge[j]), np.log10(AgeRedshift), G.z
            BHadd = np.where(BHaccrete>0. ,(BH+1.e2)/tEdd_yr*np.exp(time/tEdd*(1.-efficiency)/efficiency), 0.)
            luminosity2 = luminosity2 + np.where(BHadd <= BHaccrete, BHadd, 0.)

        luminosity2 = Msun/yr*c**2*luminosity2*1.e-7
        # -----------------------------------------------------------------------------
        # MODEL 3 ---------------------------------------------------------------------
        luminosity3 = np.zeros(len(w), dtype = np.float64)
        F = 0.1
        tEdd = 0.45e9*yr*(self.Hubble_h*km_cm/Mpc_cm)

        print "--------------------------------------------"
        for i in range(MERGER_NUM):
            BH = np.float64(G.QSOBH[w,i] * 1.e10 / self.Hubble_h)
            BHaccrete = np.float64(G.QSOBHaccrete[w,i] * 1.e10 / self.Hubble_h)
            mergeAge = np.float64(G.QSOmergeAge[w,i])
            mergeTime = G.QSOmergeTime[w,i]

            BHpeak = BH + F * BHaccrete * (1. - efficiency)
            Lpeak = 3.3e4*Lsol * BHpeak
            alpha = -0.95 + 0.32 * np.log10(Lpeak / (1.e12*Lsol))

            time = np.where(BHaccrete>0. ,(mergeAge-(AgeUniverse-AgeRedshift)) / (yr*(self.Hubble_h*km_cm/Mpc_cm)), 0.)
            lum = Lpeak * ( 1. + time/1.e9 * (Lpeak / (1.e9*Lsol))**-alpha  )**(1./alpha)

            luminosity3 = luminosity3 + lum

        # -----------------------------------------------------------------------------
        ax = plt.subplot(111)  # 1 plot on the figure

        mi = 43
        ma = 48#round(np.max(np.log10(luminosity)))
        binwidth = 0.125
        NB = (ma - mi) / binwidth

        (counts, binedges) = np.histogram(np.log10(luminosity), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth

        (counts2, binedges2) = np.histogram(np.log10(luminosity2), range=(mi, ma), bins=NB)
        xaxeshisto2 = binedges2[:-1] + 0.5 * binwidth

        (counts3, binedges3) = np.histogram(np.log10(luminosity3), range=(mi, ma), bins=NB)
        xaxeshisto3 = binedges3[:-1] + 0.5 * binwidth

        ax.set_yscale('log')
        plt.step(xaxeshisto, counts / binwidth / (self.BoxSize/self.Hubble_h)**3, 'k-', label='dynamic time')
        plt.step(xaxeshisto2, counts2 / binwidth / (self.BoxSize/self.Hubble_h)**3, 'r-', label='Eddington limit')
        plt.step(xaxeshisto3, counts3 / binwidth / (self.BoxSize/self.Hubble_h)**3, 'b-', label='Hopkins')

        # -----------------------------------------------------------------------------
        # Observational Fits

        xaxes = np.arange(43.,48.,0.1)
        if(G.z < 0.5):
            phi = -5.45
            Lstar = 11.94 +np.log10(3.9e33)
            gamma1 = 0.868
            gamma2 = 1.97
        elif(G.z < 1.5):
            phi = -4.63
            Lstar = 12.59 +np.log10(3.9e33)
            gamma1 = 0.412
            gamma2 = 2.23
        elif(G.z < 2.5):
            phi = -4.83
            Lstar = 13.10 +np.log10(3.9e33)
            gamma1 = 0.320
            gamma2 = 2.39
        elif(G.z < 3.5):
            phi = -5.23
            Lstar = 13.17 +np.log10(3.9e33)
            gamma1 = 0.395
            gamma2 = 2.10
        elif(G.z < 4.5):
            phi = -4.66
            Lstar = 12.39 +np.log10(3.9e33)
            gamma1 = 0.254
            gamma2 = 1.69
        elif(G.z < 5.5):
            phi = -5.38
            Lstar = 12.46 +np.log10(3.9e33)
            gamma1 = 0.497
            gamma2 = 1.57
        elif(G.z < 6.5):
            phi = -5.13
            Lstar = 11.00 +np.log10(3.9e33)
            gamma1 = 0.
            gamma2 = 1.11
        QLF = 10.**phi / ((10.**(xaxes-Lstar))**gamma1 + (10.**(xaxes-Lstar))**gamma2)

        plt.plot(xaxes, QLF)
        # -----------------------------------------------------------------------------


        plt.ylabel(r'$\mathrm{number}\ \mathrm{density}$    $[\mathrm{Mpc}^{-3}]$')  # Set the y...
        plt.xlabel(r'$\log L$    $[\mathrm{erg}\ \mathrm{s}^{-1}]$')  # and the x-axis labels

        plt.axis((43,48,1.e-8,1.e-3))
        leg = plt.legend(loc='upper right')
        leg.draw_frame(False)  # Don't want a box frame

        outputFile = OutputDir + '19.QSOLuminosityFunction_z' + str(G.z) + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def check_merger_snapnum(self, G):

        print 'Plotting the QSO luminosity distribution of all galaxies'

        seed(2222)

        H = 100.
        AgeRedshift = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM/(1.+self.z)**3)**0.5)
        AgeUniverse = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM)**0.5)

        a = np.loadtxt("../input/treefiles/millennium/millennium.a_list", unpack='True')

        z_snap = 1./a-1.

        plt.figure()
        ax = plt.subplot(111)  # 1 plot on the figure

        for i in range(MERGER_NUM):
            w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.))[0]
            #print len(w)
            w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,i] > 0.))[0]
            #print len(w)
            if(len(w) > dilute): w = sample(w, dilute)

            time = AgeUniverse - G.QSOmergeAge[w,i]
            redshift = (self.OmegaL/self.OmegaM/(np.sinh(1.5*self.OmegaL**0.5*H*time))**2)**(1./3.)-1.

            t = i*np.ones(len(G.QSOBHaccrete[w,i]))/MERGER_NUM
            thisPlot = plt.scatter(redshift, G.QSOmergSnap[w,i], marker='o', s=10., c=t, cmap = cm.jet_r, alpha=0.8, label='Stars')
            thisPlot.set_clim(vmin=0,vmax=1)

        x = np.arange(len(a))
        y = 1./a[x]-1.
        plt.plot(y,x)
        plt.axis([G.z-0.1, 11.0, 0., G.SnapNum[0]])
        #ax.set_yscale('log')
        plt.ylabel(r'$\mathrm{snapshot\ when\ estimate\ merge\ time}$')  # Set the y...
        plt.xlabel(r'$z_{\mathrm{merge}}$')  # and the x-axis labels

        #plt.axis([0., 10.5, 15., 64.])

        #leg = plt.legend(loc='upper right')
        #leg.draw_frame(False)  # Don't want a box frame

        outputFile = OutputDir + '20.CheckMergerSnapNum' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)


# --------------------------------------------------------

    def QSOnumMergers_redshift(self, G):

        print 'Plotting the QSO luminosity distribution of all galaxies'

        seed(2222)

        H = 100.
        AgeRedshift = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM/(1.+self.z)**3)**0.5)
        AgeUniverse = 2./(3.*H*self.OmegaL**0.5)*np.arcsinh((self.OmegaL/self.OmegaM)**0.5)

        a = np.loadtxt("../input/treefiles/millennium/millennium.a_list", unpack='True')

        z_snap = 1./a-1.

        plt.figure()
        ax = plt.subplot(111)  # 1 plot on the figure
        #
        # for i in range(MERGER_NUM):
        #     w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.) & (G.QSOBHaccrete[:,i] > 0.))[0]
        #     if(len(w) > dilute): w = sample(w, dilute)
        #
        #     time = AgeUniverse - G.QSOmergeAge[w,i]
        #     redshift = (self.OmegaL/self.OmegaM/(np.sinh(1.5*self.OmegaL**0.5*H*time))**2)**(1./3.)-1.
        #
        #     t = i*np.ones(len(G.QSOBHaccrete[w,i]))/MERGER_NUM
        #     thisPlot = plt.scatter(redshift, G.MergNum[w], marker='o', s=10., c=t, cmap = cm.jet_r, alpha=0.8, label='Stars')
        #     thisPlot.set_clim(vmin=0,vmax=1)

        w = np.where((G.Mvir > 0.0) & (G.StellarMass > 0.))[0]
        if(len(w) > dilute): w = sample(w, dilute)
        plt.scatter(np.log10(G.BlackHoleMass[w]), G.MergNum[w], marker='o', s=10., alpha=0.8, label='Stars')
        plt.axis([-10, 0, 0., 50])
        # ax.set_xscale('log')
        plt.ylabel(r'$\mathrm{snapshot\ when\ estimate\ merge\ time}$')  # Set the y...
        plt.xlabel(r'$z_{\mathrm{merge}}$')  # and the x-axis labels

        outputFile = OutputDir + '21.QSOnumMergers_redshift' + OutputFormat
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

        # Add this plot to our output list
        OutputList.append(outputFile)

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
        default='./results/millennium_QSO/',
        help='input directory name (default: ./results/millennium_QSO/)',
        metavar='DIR',
        )
    parser.add_option(
        '-f',
        '--file_base',
        dest='FileName',
        default='model_z0.000',
        help='filename base (default: model_z0.000)',
        metavar='FILE',
        )
    parser.add_option(
        '-n',
        '--file_range',
        type='int',
        nargs=2,
        dest='FileRange',
        default=(0, 7),
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

    print 'Anayzing snapshot', G.SnapNum[0], "at redshift", G.z

    res.QSOBHaccrete_distribution(G)
    res.QSOBHaccrete_function(G)
    res.QSO_distribution(G)
    res.QSOluminosity_function(G)
    res.check_merger_snapnum(G)
    res.QSOnumMergers_redshift(G)
