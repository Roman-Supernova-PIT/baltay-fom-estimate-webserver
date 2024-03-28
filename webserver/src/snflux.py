# Program to calculate supernova fluxes starting with the Hsiao typical Type 1a spectrum

import numpy as np
import math
import scipy
import sys

def zint(z, omegam, omegak, omegade, w):
    eofz = math.sqrt(omegam*(1.0+z)**3+omegak*(1.0+z)**2+omegade*(1.0+z)**(3.0*(1.0+w)))
    zint = 1./eofz
    return zint

def snflux( inparams ):
        # Instrument
    params = {'telarea': 0.78*(math.pi*(120.0**2)),
              'thruput': 0.7, # thruput includes detector quantum efficiency
        # Astronomy
              'Mb' : -19.7, # SN B band absolute magnitude
        # Cosmology
              'h0': 70,
              'omegam' : 0.28,
              'omegak' : 0.,
              'omegade' : 0.72,
              'w' : -1.
             }
    params.update( inparams )

    # Set parameters
    c = scipy.constants.c
    h0 = params['h0']
    hp = scipy.constants.h
    telarea = params['telarea']
    thruput = params['thruput']
    
    Mb = params['Mb']
    anorm = (10**(-0.4*Mb))

    # Set cosmology
    omegam = params['omegam']
    omegak = params['omegak']
    omegade = params['omegade']
    w = params['w']

    try:
        # Setup output
        f1 = open('snfluxes.d', 'w')
        f2 = open('hsiao.txt', 'r')
        f3 = open('fluxout.txt', 'w')

        f3.write('omegam\tomegak\tomegade\tw\n')
        f3.write('\t'.join([str(x) for x in [omegam, omegak, omegade, w]])+'\n')

        # Read in Hsiao spectrum in rest frame
        wl = [0]*1800
        flux = [0]*1800
        for i in range(1800):
            a, wl[i], flux[i] = [float(x) for x in f2.readline().split()]
            flux[i] *= anorm

            # Converting ergs/cm2/sec/A to photons/cm2/sec/10Abin
            egamma = ((hp*3.0e8)/(wl[i]*1.0e-10))*1.0e7
            flux[i] *= 10./egamma

        # Set min and max wavelength range on Hsiao spectrum
        filtmin = [
            0.760,
            0.927,
            1.131,
            1.380,
            0.500]
        filtmax = [
            0.977,
            1.192,
            1.454,
            1.774,
            0.600]

        rfmin = np.empty((5, 20))
        rfmax = np.empty((5, 20))
        imin = np.empty((5, 20), dtype=int)
        imax = np.empty((5, 20), dtype=int)

        z = np.array([i*0.1 + 0.05 for i in range(20)])
        for j in range(4):
            rfmin[j] = filtmin[j]/(1.+z)
            rfmax[j] = filtmax[j]/(1.+z)
            
            imin[j] = (1000*rfmin[j]-100).astype(int)
            imax[j] = (1000*rfmax[j]-100).astype(int)

        rfmin[4] = filtmin[4]
        rfmax[4] = filtmax[4]

        imin[4] = (1000*rfmin[4]-100).astype(int)
        imax[4] = (1000*rfmax[4]-100).astype(int)

        # Add up flux between limits in rest frame
        signal = np.zeros((5, 20))
        for j in range(5):
            signal[j] = np.array([np.sum(flux[imin[j, i]:imax[j, i]+1]) for i in range(20)])

        # Signal is in photons/cm2/sec/1000A band in the supernova rest frame
        # convert flux from rest frame to observer frame

        for j in range(5):
            for i in range(20):
                z = i*0.1 + 0.05
        # calculate luminosity distance dlum
        # Calculate integral of 1/E(z)
                dist = scipy.integrate.quad(zint, 0, z, args=(omegam, omegak, omegade, w))[0]
                dlum = (c/h0)*(1.0+z)*dist
                distmod=5.0*math.log10(dlum)+25.0

                signal[j,i] = signal[j,i]/((dlum/1.0e-5)**2)
                signal[j,i] = signal[j,i]*telarea*thruput

        # signal[j,i] is in counts/sec/restframe band in observer frame
        for i in range(20):
            f1.write('\t'.join([f'{x:.3f}' for x in signal[:,i]])+'\n')

    except Exception as ex:
        sys.stderr.write( f"Exception running snflux; params is {params}\n" )
        raise

if __name__ == "__main__":
    snflux( {} )