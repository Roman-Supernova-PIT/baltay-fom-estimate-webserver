import sys
import math

def prism( inparams ):
        # Survey
    params = {'length1': 1.9,   # Number of years of survey
              'length2': 1.9,
              'length3': 1.9,
              'sqdegim1': 19.04,  # Square Degrees of imaging survey
              'sqdegim2': 4.20,
              'sqdegim3': 0.0,
              'tfixim1': 115.0,   # Exposure Time in seconds of imaging survey
              'tfixim2': 450.,
              'tfixim3': 0,
              'sqdegspec1': 3.36,   # Square degrees of spectroscopic survey
              'sqdegspec2': 1.12,
              'sqdegspec3': 0.0,
              'tfixspec1': 900.0,   # Exposure time of spectroscopic survey
              'tfixspec2': 3600.0,
              'tfixspec3': 0,
              'z1': 2.0,      # Redshift limit -- what does this mean?  Complete?
              'z2': 2.0,
              'z3': 0.0,
              'nlim': 20,     # Number of redshift bins from snflux ; don't change this!
        # Instrument
              'eff': 0.9,     # Efficiency of... WHAT?
        # Errors
              'constsysim': 0.015,  # Systematic error at redshift given by divim
              'divim': 1.8,         # 1+z of where sys error is constsysim ; scaled by (1+z)/divim
              'stnim': 10.0,        # S/N (of what? PSF? 1FWHM radius?) for imaging
              'constsysspec': 0.01,  # spectroscopic equivalent of constsysim and divim
              'divspec': 1.8,
              'stnspec': 20.0        # S/N of what? for spectra (one resolution element?  20 is a lot!)
              }
    params.update( inparams )

    try:
        f2 = open('prismoutSept15.txt', 'w')
        f3 = open('sninvarsp.txt', 'w')
        f5 = open('sninvarim.txt', 'w')
        f4 = open('snfluxes.d', 'r')
        f7 = open('errout312.txt', 'w')

        years=[params['length1'], params['length2'], params['length3']]

        sqdegminim= params['sqdegim1']
        sqdegmidim= params['sqdegim2']
        sqdegmaxim= params['sqdegim3']


        sqdgyrminim=sqdegminim*years[0]
        sqdgyrmidim=sqdegmidim*years[1]
        sqdgyrmaxim=sqdegmaxim*years[2]

        tfixim = [params['tfixim1'], params['tfixim2'], params['tfixim3']]

        zmin=params['z1']
        zmid= params['z2']
        zmax=params['z3']
        nlim=params['nlim']

        eff=params['eff']

        constsysim=params['constsysim']
        divim=params['divim']

        stnimag=params['stnim']

        # Spectroscopy

        sqdegminsp = params['sqdegspec1']
        sqdegmidsp = params['sqdegspec2']
        sqdegmaxsp = params['sqdegspec3']

        sqdgyrminsp = sqdegminsp*years[0]
        sqdgyrmidsp = sqdegmidsp*years[1]
        sqdgyrmaxsp = sqdegmaxsp*years[2]

        tfixspec = [params['tfixspec1'], params['tfixspec2'], params['tfixspec3']]


        constsyssp=params['constsysspec']
        divsp=params['divspec']

        stnspec=params['stnspec']

        # Instrument parameters

        pix=14.6
        dc=0.015
        pixim=8.97
        dcim=0.015
        zodiim=0.14
        rnoise=5.0
        res=100.0
        resifu=80.0
        wifu=3.0
        lmin=1.0
        lmax=2.0
        rlamb=4.0
        telarea=0.78*(math.pi*(120.0**2))
        telarea=telarea*0.70
        pixarea=0.11**2

        rifu = [0,0,0,0,
            580.0,
            370.0,
            224.0,
            161.0,
            129.0,
            109.0,
            97.0,
            89.0,
            83.0,
            80.0,
            80.0,
            80.0,
            83.0,
            89.0,
            95.0,
            100.0]

        filtmin= [0.760,
            0.927,
            1.131,
            1.380,
            0.500]

        filtmax= [0.977,
            1.192,
            1.454,
            1.774,
            0.600]

        sr = [0.1,
            1.0,
            3.0,
            5.8,
            9.7,
            14.2,
            18.8,
            23.4,
            28.8,
            33.3,
            36.7,
            38.4,
            38.6,
            37.5,
            34.8,
            31.2,
            27.4,
            23.7,
            20.0,
            16.7]

        for i in range(3):
            if tfixim[i] == 0:
                tfixim[i] = float('nan')
            if tfixspec[i] == 0:
                tfixspec[i] = float('nan')

        f2.write('\n')
        f2.write('Input supernova signals in five bands\n')
        f2.write('\n')

        signal = []
        for i in range(nlim):
            signal += [[float(x) for x in f4.readline().split()]]

        for i in range(nlim):
            f2.write('\t'.join([str(x) for x in signal[i]]))
            f2.write('\n')

        hrswide=((sqdegminsp/0.28)*(tfixspec[0]+70.0))/(3600.0)
        hrsdeep=((sqdegmidsp/0.28)*(tfixspec[1]+70.0))/(3600.0)
        hrstotal=hrswide+hrsdeep

        f2.write('\n')
        f2.write('SPECTROSCOPIC SURVEY\n')
        f2.write('\n')
        f2.write('survey parameters Area, Exp time, wide, deep\n')
        f2.write('\n')
        f2.write(f'{sqdegminsp}\t{sqdegmidsp}\t{tfixspec[0]}\t{tfixspec[1]}\n')
        f2.write('\n')
        f2.write('Hours/visit for spectroscopy, wide, deep, total\n')
        f2.write('\n')
        f2.write(f'{hrswide}\t{hrsdeep}\t{hrstotal}\n')
        f2.write('\n')


        snesp=[0.0]*nlim
        nsnesp=[0.0]*nlim
        snmin=[0.0]*nlim
        snmid=[0.0]*nlim
        snmax=[0.0]*nlim

        snsp = [[0.]*nlim for x in range(3)]
        totalsne = [0.]*nlim

        for i in range(1, nlim):
            z=i*0.1+0.05
            if z <= zmin:
                snmin[i]=sr[i]*sqdgyrminsp*eff
            if z <= zmid:
                snmid[i]=sr[i]*sqdgyrmidsp*eff
            if z <= zmax:
                snmax[i]=sr[i]*sqdgyrmaxsp*eff

            sntot=snmin[i]+snmid[i]+snmax[i]

            snsp[0][i]=snmin[i]
            snsp[1][i]=snmid[i] 
            snsp[2][i]=snmax[i]

            totalsne[i]=sntot
            snesp[i]=sntot

        f2.write('\n')
        f2.write('zodi backgrounds in the rest frame V band\n')
        f2.write('\n')

        zodi = [0]*nlim

        for i in range(1, nlim):
            z=i*0.10+0.05

            alambdasl=1.2
            dlambdasl=1.05

            fzodi=(0.0247*alambdasl)*math.exp(-1.68*alambdasl)
            zodi[i]=fzodi*telarea*pixarea*dlambdasl

        f2.write('\t'.join([f'{x:.3f}' for x in zodi[1:]]))

        stnsp = [[0.]*nlim for x in range(3)]
        tfixcomb = [0.]*3
        for k in range(3):
            for i in range(1, nlim):
                z=i*0.1+0.05
                nspec=2.0*(1.0+z)+1.0
                specno=nspec
                tfixcomb[k]=specno*tfixspec[k]

                readpairs=(tfixcomb[k]/specno)/5.6
                rn=(20.0/math.sqrt(readpairs))*1.732
                rnoise=math.sqrt(rn**2+5.0**2)

                rnoise=rnoise*math.sqrt(specno)
                centlamb=((filtmax[4]+filtmin[4])/2.0)*(1.0+z)
                dellamb=(filtmax[4]-filtmin[4])*(1.0+z)

                il=math.floor(centlamb*10.0)
                resifu=rifu[il]
                delpix=centlamb/(2.0*resifu)
                pixno=(dellamb/delpix)*wifu

                if z <= 0.5:
                    st=signal[i][0]*tfixcomb[k]
                else:
                    st=signal[i][4]*tfixcomb[k]

                anoise=pixno*(((zodi[i]+dc)*tfixcomb[k])+(rnoise**2))
                stnsp[k][i]=st/(math.sqrt(st+2.5*st+anoise))

                if stnsp[k][i] < stnspec:
                    snsp[k][i]=0.0

        for i in range(1, nlim):
            snesp[i] = snsp[0][i] + snsp[1][i] + snsp[2][i]


        f2.write('\n\n')
        f2.write('z, prism S/N in three tiers\n')
        f2.write('\n')

        for i in range(1, nlim):
            z=i*0.1+0.05

            f2.write(f'{z:.2f}\t{stnsp[0][i]:.2f}\t{stnsp[1][i]:.2f}\t{stnsp[2][i]:.2f}\n')

        f2.write('\n')    
        f2.write('The spectroscopic survey\n')
        f2.write('\n')
        f2.write('z,nlow,nmid,nhigh,ntot,staterr,sigstat,sigsys,sigtotal\n')
        f2.write('\n')

        nsummin=0
        nsummid=0
        nsummax=0
        nsum=0

        sigstat = [0.]*nlim
        sigsys = [0.]*nlim
        sigout = [0.]*nlim
        sigdist = [0.]*nlim
        for i in range(1, nlim):
            z=i*0.1+0.05

            sigmeas=0.08
            sigint=0.08
            siglens=0.07*z
            sigstat[i]=math.sqrt(sigmeas**2+sigint**2+siglens**2)
            if snesp[i] == 0:
                statsig=float('inf')
            else:
                statsig=sigstat[i]/math.sqrt(snesp[i])
            sigsys[i]=constsyssp*((1.0+z)/divsp)
            s=math.sqrt(statsig**2+sigsys[i]**2)
            sigout[i]=s
            sigdist[i]=s/2.0

            nmin=snsp[0][i]
            nmid=snsp[1][i]
            nmax=snsp[2][i]
            ntotal=snesp[i]

            nsummin=nsummin+nmin
            nsummid=nsummid+nmid
            nsummax=nsummax+nmax
            nsum=nsum+ntotal

            f2.write(f'{z:.2f}\t{int(nmin)}\t{int(nmid)}\t{int(nmax)}\t{int(ntotal)}\t{sigstat[i]:.3f}\t{statsig:.3f}\t{sigsys[i]:.3f}\t{s:.3f}\n')

            s1=statsig/2.0
            s2=sigsys[i]/2.0
            s3=s/2.0

        f2.write('\n')
        f2.write('total number of supernovae\n')
        f2.write(f'{int(nsummin)}\t{int(nsummid)}\t{int(nsummax)}\t{int(nsum)}\n')
        f2.write('\n')

        f2.write('\n')
        f2.write('z,sigdist\n')
        f2.write('\n')

        sigout[0]=0.006
        sigdist[0]=0.003

        for i in range(nlim):
            z=i*0.1+0.05
            f2.write(f'{z:.2f}\t{sigdist[i]:.4f}\n')
            f7.write(f'{z:.2f}\t{sigout[i]:.4f}\n')

        om=0.28
        nw0=-1.0
        nwa=0.0
        nsamestat=0
        nsamesys=0
        nyescmb=1

        snesp[0]=800.0
        snminsp=[0.]*nlim
        snmidsp=[0.]*nlim
        snmaxsp=[0.]*nlim

        nsne=[0.]*nlim

        for i in range(nlim):
            z=i*0.1+0.05
            nsne[i]=snesp[i]

        sigstat[0]=0.17
        sigsys[0]=0.0

        f3.write(f'{om}\t{nw0}\t{nwa}\n')
        f3.write('\t'.join([str(int(x)) for x in nsne])+'\n')
        f3.write(f'{nsamestat}\t{nsamesys}\t{nyescmb}\n')
        f3.write('\t'.join([f'{x:.4f}' for x in sigstat])+'\n')
        f3.write('\t'.join([f'{x:.4f}' for x in sigsys])+'\n')

        # Imaging Survey.  Xxxxxxxxxxxxxxxxxxxxxx

        hrswide=((sqdegminim/0.28)*(tfixim[0]+70.0)*4.0)/(3600.0)
        hrsdeep=((sqdegmidim/0.28)*(tfixim[1]+70.0)*4.0)/(3600.0)
        hrstotal=hrswide+hrsdeep

        f2.write('\n')
        f2.write('IMAGING ONLY SURVEY\n')
        f2.write('\n')
        f2.write('survey parameters Area, Exp time, wide, deep\n')
        f2.write('\n')
        f2.write(f'{sqdegminim}\t{sqdegmidim}\t{tfixim[0]}\t{tfixim[1]}\n')
        f2.write('\n')
        f2.write('Hours/visit for spectroscopy, wide, deep, total\n')
        f2.write('\n')
        f2.write(f'{hrswide}\t{hrsdeep}\t{hrstotal}\n')
        f2.write('\n')


        sneim=[0.0]*nlim
        nsneim=[0]*nlim
        snminsp=[0.0]*nlim
        snmidsp=[0.0]*nlim
        snmaxsp=[0.0]*nlim
        snim=[[0.0]*nlim for x in range(3)]
        for i in range(1, nlim):
            z=i*0.1+0.05

            if z <= zmin:
                snmin[i]=sr[i]*sqdgyrminim*eff
            if z <= zmid:
                snmid[i]=sr[i]*sqdgyrmidim*eff
            if z <= zmax:
                snmax[i]=sr[i]*sqdgyrmaxim*eff

            sntot=snminsp[i]+snmidsp[i]+snmaxsp[i]

            snim[0][i]=snmin[i]
            snim[1][i]=snmid[i]
            snim[2][i]=snmax[i]

            totalsne[i]=sntot
            sneim[i]=sntot

        stnim=[[0.0]*nlim for x in range(3)]
        for k in range(3):
            for i in range(1, nlim):
                z=i*0.1+0.05
                j=1

                readpairs=tfixim[k]/5.6
                rn=(20.0/math.sqrt(readpairs))*1.732
                rnoise=math.sqrt(rn**2+5.0**2)

                st=signal[i][j]*tfixim[k]

                anoise=pixim*(((zodiim+dc)*tfixim[k])+(rnoise**2))

                stnim[k][i]=st/(math.sqrt(st+(0.5)*st+anoise))

                if stnim[k][i] < stnimag:
                    snim[k][i]=0.0

        f2.write('\n')
        f2.write('z, imaging S/N in three tiers\n')
        f2.write('\n')

        for i in range(1, nlim):
            z=i*0.1+0.05

            f2.write(f'{z:.2f}\t{stnim[0][i]:.2f}\t{stnim[1][i]:.2f}\t{stnim[2][i]:.2f}\n')

        for k in range(3):
            for i in range(nlim):
                snim[k][i]=snim[k][i]-snsp[k][i]

        for i in range(1, nlim):
            sneim[i]=snim[0][i]+snim[1][i]+snim[2][i]

        f2.write('\n')
        f2.write('The imaging only survey\n')
        f2.write('\n')
        f2.write('z,nlow,nmid,nhigh,ntot,staterr,sigstat,sigsys,sigtotal\n')
        f2.write('\n')

        nsummin=0
        nsummid=0
        nsummax=0
        nsum=0

        for i in range(1, nlim):
            z=i*0.1+0.05

            sigmeas=0.08
            sigint=0.11+0.033*z
            siglens=0.07*z
            sigstat[i]=math.sqrt(sigmeas**2+sigint**2+siglens**2)
            if sneim[i] < 0:
                raise Exception( f"sneim is {sneim[i]:.02f} in bin {i}"
                                 f"(redshift {i*0.1+0.05:.02f}); probably an spectroscopic survey"
                                 f"is too deep compared to its imaging counterpart" )
            if sneim[i] == 0:
                statsig = float('inf')
            else:
                statsig=sigstat[i]/math.sqrt(sneim[i])
            sigsys[i]=constsysim*((1.0+z)/divim)
            s=math.sqrt(statsig**2+sigsys[i]**2)
            sigout[i]=s
            sigdist[i]=s/2.0

            nmin=snim[0][i]
            nmid=snim[1][i]
            nmax=snim[2][i]
            ntotal=sneim[i]


            nsummin=nsummin+nmin
            nsummid=nsummid+nmid
            nsummax=nsummax+nmax
            nsum=nsum+ntotal


            f2.write(f'{z:.2f}\t{int(nmin)}\t{int(nmid)}\t{int(nmax)}\t{int(ntotal)}\t{sigstat[i]:.3f}\t{statsig:.3f}\t{sigsys[i]:.3f}\t{s:.3f}\n')

            s1=statsig/2.0
            s2=sigsys[i]/2.0
            s3=s/2.0

        f2.write('\n')
        f2.write('total number of supernovae\n')
        f2.write(f'{int(nsummin)}\t{int(nsummid)}\t{int(nsummax)}\t{int(nsum)}\n')
        f2.write('\n')


        f2.write('\n')
        f2.write('z,sigdist\n')
        f2.write('\n')

        sigout[0]=0.006
        sigdist[0]=0.003

        for i in range(17):
            z=i*0.1+0.05
            f2.write(f'{z:.2f}\t{sigdist[i]:.4f}\n')
            f7.write(f'{z}\t{sigout[i]}\n')

        om=0.28
        nw0=-1.0
        nwa=0.0
        nsamestat=0
        nsamesys=0
        nyescmb=1


        sneim[0]=800.0

        for i in range(nlim):
            z=i*0.1+0.05
            nsne[i]=sneim[i]

        sigstat[0]=0.17
        sigsys[0]=0.0

        f5.write(f'{om}\t{nw0}\t{nwa}\n')
        f5.write('\t'.join([str(int(x)) for x in nsne])+'\n')
        f5.write(f'{nsamestat}\t{nsamesys}\t{nyescmb}\n')
        f5.write('\t'.join([f'{x:.4f}' for x in sigstat])+'\n')
        f5.write('\t'.join([f'{x:.4f}' for x in sigsys])+'\n')

        f2.close()
        f3.close()
        f5.close()
        f4.close()
        f7.close()

    except Exception as ex:
        sys.stderr.write( f"Exception running prism; params is {params}\n" )
        raise

# ======================================================================

if __name__ == "__main__":
    prism( {} )
    
