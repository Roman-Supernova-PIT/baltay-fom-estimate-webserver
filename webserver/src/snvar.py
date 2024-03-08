# This script adapted from snvar.f by Eric Linder 10/25/02

import scipy
import numpy as np
import math
    

f1 = open('sninvar.txt', 'r')
f2 = open('snoutput.txt', 'w')
f3 = open('snerrout.txt', 'w')

# Input line 1: matter density Omega_m, w0, wa (fiducial 0.28, -1, 0)
om, w0, wa = [float(x) for x in f1.readline().split()]
# Line 2: number of SN per 0.1 bin in redshift, put 0 for empty bins, put z<0.1 SN in 1st bin
sn_num = [int(x) for x in f1.readline().split()]
# Line 3: Toggle (1 yes, 0 no) 3 values: 1) Is intrinsic uncertainty of 1 SN the same for all redshifts? 2) Is systematic error level given by dm=sigsys*(1+z)/2.7 with sigsys the same for all redshifts? 3) Do you want to include CMB Fisher matrix (see below)?
samesigm0, samesigsys, yescmb = [int(x) for x in f1.readline().split()]
# Line 4: intrinsic uncertainty on 1 SN: either 1 number if constant with redshift, or 20 numbers giving sigm0 for each 0.1 bin in redshift. Standard is all sigm0=0.1
if samesigm0 == 1:
    sigm0 = float(f1.readline())
    sigm = [sigm0]*20
else:
    sigm = [float(x) for x in f1.readline().split()]
# Line 5: systematic error level: either 1 number for sigsys if follow dm=sigsys*(1+z)/2.7, or 20 numbers giving dm (not sigsys!) for each 0.1 bin in redshift.
# Systematic error is added in quadrature with intrinsic uncertainty.
# Standard is all sigsys=0.02 (so dm ~ 1+z)
# CMB Fisher matrix approximates Planck as 0.2% accuracy on shift parameter R (here called d). This is a good approximation only when combined with SN
if samesigsys == 1:
    sigsys = float(f1.readline())
    dm = [sigsys*(1.+0.1*i - 0.05)/2.7 for i in range(1, 21)]
else:
    dm = [float(x) for x in f1.readline().split()]

# Write input to output file for a check
f2.write('Reprint input to check\n')
f2.write('Omega_m, w0, wa\n')
f2.write(f'{om}\t{w0}\t{wa}\n')
f2.write('\nNumber of SN in each 0.1 bin of redshift\n')
f2.write('\t'.join([str(x) for x in sn_num])+'\n')
f2.write('\nIntrinsic error constant? Systematic error constant? Include CMB?\n')
f2.write(f'{samesigm0}\t{samesigsys}\t{yescmb}\n')
f2.write('\nIntrinsic error in each 0.1 bin of redshift\n')
f2.write('\t'.join([f'{x:.3f}' for x in sigm])+'\n')
f2.write('\nSystematic error dm in each 0.1 bin of redshift\n') 
f2.write('\t'.join([f'{x:.3f}' for x in dm])+'\n')

q = 1.-om # Dark energy density - assumes flatness
u = w0
v = wa
s = 4 # Dimension of Fisher matrix (scriptM, Omega_m, w0, wa)

# fnr calculates conformal distance
# fndX calculates analytic derivative of distance wrt parameter X

def fnr(x):
    powarg = 2.*(1. + u + v)
    exparg = -3. * v * (1. - 1./x)
    ewfull = x**powarg*math.exp(exparg)
    fnr = 1./math.sqrt(om*x**3+q*ewfull)
    return fnr

def fndo(x):
    powarg=3.*(1.+u+v)
    exparg=-3.*v*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    fndo=(x**3.-ewfull)/(om*x**3+q*ewfull)**1.5
    fndo=-fndo/2.
    return fndo

def fndu(x):
    powarg=3.*(1.+u+v)
    exparg=-3.*v*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    fndu=ewfull*math.log(x)/(om*x**3+q*ewfull)**1.5
    fndu=-3.*q*fndu/2.
    return fndu

def fndv(x):
    powarg=3.*(1.+u+v)
    exparg=-3.*v*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    fndv=ewfull/(om*x**3+q*ewfull)**1.5
    fndv=fndv*(math.log(x)-1.+1./x)
    fndv=-3.*q*fndv/2.
    return fndv

# dfn calculates reduced distance to last scattering (shift param)
# ddXfn calculates analytic derivative of distance wrt parameter X

def dfn(x):
    w0 = u
    wa = v
    powarg = 3.*(1.+w0+wa)
    exparg=-3.*wa*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    dfn=1./math.sqrt(x**3.+q/om*ewfull)
    return dfn

def ddofn(x):
    w0=u
    wa=v
    powarg=3.*(1.+w0+wa)
    exparg=-3.*wa*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    ddofn=ewfull/(x**3.+q/om*ewfull)**1.5
    return ddofn

def ddw0fn(x):
    w0=u
    wa=v
    powarg=3.*(1.+w0+wa)
    exparg=-3.*wa*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    ddw0fn=ewfull*math.log(x)/(x**3.+q/om*ewfull)**1.5
    return ddw0fn

def ddwafn(x):
    w0=u
    wa=v
    powarg=3.*(1.+w0+wa)
    exparg=-3.*wa*(1.-1./x)
    ewfull=x**powarg*math.exp(exparg)
    ddwafn=ewfull*(math.log(x)-1.+1./x)/(x**3.+q/om*ewfull)**1.5
    return ddwafn

# Calculate CMB Fisher matrix as reduced distance to last scatter

botc = 1.
topc = 1090

d = scipy.integrate.quad(dfn, botc, topc)[0]
sigd=0.002*d # 0.2% precision on shift parameter d(=R)

ddo = scipy.integrate.quad(ddofn, botc, topc)[0]
ddw0 = scipy.integrate.quad(ddw0fn, botc, topc)[0]
ddwa = scipy.integrate.quad(ddwafn, botc, topc)[0]
ddo = 0.5/om/om*ddo
ddw0 = -1.5*q/om*ddw0
ddwa = -1.5*q/om*ddwa

fcoo=yescmb*ddo*ddo/sigd**2
fcou=yescmb*ddo*ddw0/sigd**2
fcov=yescmb*ddo*ddwa/sigd**2
fcuu=yescmb*ddw0*ddw0/sigd**2
fcuv=yescmb*ddw0*ddwa/sigd**2
fcvv=yescmb*ddwa*ddwa/sigd**2
foo=fcoo
fou=fcou
fov=fcov
fuu=fcuu
fuv=fcuv
fvv=fcvv

# Begin SN Fisher calculation
bot = 1.
sigm2 = [0]*20
sigbin = [0]*20

fmm = 0
fmo = 0
fmu = 0
fmv = 0
for i in range(20):
    z = 0.1*i+0.05
    y = 1 + z
    sigm2[i] = sigm[i]*sigm[i]+sn_num[i]*dm[i]*dm[i]

    if sn_num[i] != 0:
        sigbin[i] = math.sqrt(sigm2[i]/sn_num[i])
    else:
        sigbin[i] = float('inf')
    f3.write(f'{z}\t{sigbin[i]}\n')

    top = y
    ansr = scipy.integrate.quad(fnr, bot, top)[0]
    ansdo = scipy.integrate.quad(fndo, bot, top)[0]
    ansdu = scipy.integrate.quad(fndu, bot, top)[0]
    ansdv = scipy.integrate.quad(fndv, bot, top)[0]
    pre = 5./math.log(10)/ansr
    dmdo=pre*ansdo
    dmdu=pre*ansdu
    dmdv=pre*ansdv
# Derivatives with respect to scriptM labeled by m, to Omega_m by o, to w0 by u, to wa by v
    fmm+=sn_num[i]/sigm2[i]
    fmo+=sn_num[i]*dmdo/sigm2[i]
    fmu+=sn_num[i]*dmdu/sigm2[i]
    fmv+=sn_num[i]*dmdv/sigm2[i]
    foo+=sn_num[i]*dmdo*dmdo/sigm2[i]
    fou+=sn_num[i]*dmdo*dmdu/sigm2[i]
    fov+=sn_num[i]*dmdo*dmdv/sigm2[i]
    fuu+=sn_num[i]*dmdu*dmdu/sigm2[i]
    fuv+=sn_num[i]*dmdu*dmdv/sigm2[i]
    fvv+=sn_num[i]*dmdv*dmdv/sigm2[i]

# Populate symmetric Fisher matrix matin(4, 4)
matin = np.empty((4, 4))
matin[0, 0] = fmm
matin[0, 1] = fmo
matin[0, 2] = fmu
matin[0, 3] = fmv
matin[1, 1] = foo
matin[1, 2] = fou
matin[1, 3] = fov
matin[2, 2] = fuu
matin[2, 3] = fuv
matin[3, 3] = fvv
matin[1, 0] = matin[0, 1]
matin[2, 0] = matin[0, 2]
matin[3, 0] = matin[0, 3]
matin[2, 1] = matin[1, 2]
matin[3, 1] = matin[1, 3]
matin[3, 2] = matin[2, 3]

# Output Fisher matrix
f2.write('\n')
f2.write('OUTPUT:\n')
f2.write('Fisher matrix in scriptM, Omega_m, w0, wa\n')
for ri in range(len(matin)):
    f2.write('\t'.join([f'{x:.4E}' for x in matin[ri]])+'\n')

# Invert to get covariance matrix
det = scipy.linalg.det(matin)
matin = scipy.linalg.inv(matin)

# Calculate pivot redshift zp and minimum w(a) uncertainty
# sigwp also uncertainty on constant w (fixing wa=0)
ap = 1.+matin[2, 3]/matin[3, 3]
zp = 1./ap-1.
sigwp = math.sqrt(matin[2, 2] - matin[2, 3]**2/matin[3, 3])

# Output parameter uncertainties: scriptM, Omega_m, w0, wa, wp
f2.write('\nUncertainties on scriptM, Omega_m, w0, wa, wp\n')
f2.write('\t'.join([f'{x:.4E}' for x in np.sqrt(np.diag(matin))])+'\n')

# Marginalize over all parameters except w0, wa, and reinvert to get 2D Fisher matrix for plotting contour ellipse
matsm = np.empty((2, 2))
matsm[0, 0] = matin[2, 2]
matsm[0, 1] = matin[2, 3]
matsm[1, 0] = matsm[0, 1]
matsm[1, 1] = matin[3, 3]

# Calculate error ellipse parameters
sigxsq=matsm[0, 0]
sigysq=matsm[1,1]
sigxy=matsm[0, 1]
ellsum=(sigxsq+sigysq)/2.
elldiff=math.sqrt(((sigxsq-sigysq)/2.)**2+(sigxy)**2)
ellmajor=math.sqrt(ellsum+elldiff)
ellminor=math.sqrt(ellsum-elldiff)
elltan=(2.*sigxy)/(sigxsq-sigysq)
elltheta=(math.atan(elltan)/2.)*57.296+90

# Write out error ellipse parameters
f2.write('\nerror ellipse major axis, minor axis, theta\n')
f2.write(f'{ellmajor:.3f}\t{ellminor:.3f}\t{elltheta:.3f}\n')

# Invert 2x2 error matrix
det = scipy.linalg.det(matsm)
matsm = scipy.linalg.inv(matsm)

# Write area Figure of Merit = 1/COV[w0,wa] = 1/sigma(wp)*sigma(wa)

f2.write(f'\nFOM = {math.sqrt(1./det):.2f}\n')

# Output quantities needed for plotting contour ellipse of w0-wa

f2.write('\nValues of w0, wa, reduced Fisher matrix F11, F12, F22\n')
f2.write('\t'.join([str(x) for x in [w0, wa, matsm[0, 0], matsm[0, 1], matsm[1, 1]]])+'\n')

f1.close()
f2.close()
f3.close()