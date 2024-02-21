c  snvar--Fisher method on supernovae for flat, w0-wa cosmologies
c	 version with errors varying with redshift
c                          -- Eric Linder  25 Oct 02
c
c     om - matter density (omega)
c     q  - quintessence (q=1-om)
c     w  - w=w0+wa*(1-a)
c
      double precision om,q,u,v,det,w0,wa,n(20),sigm0,sigsys
      double precision yescmb,bot,top,z,y,dfn,d,sigd,ddofn
      double precision fmm,fmo,fmu,fmv,foo,fou,fov,fuu,fuv,fvv
      double precision botc,topc,ddo,ddw0,ddwa,ddw0fn,ddwafn
      double precision fcoo,fcuu,fcvv,fcov,fcou,fcuv,sigm2(20)
      double precision fnr,fndo,fndu,fndv,ansr,ansdo,ansdu
      double precision ansdv,pre,dmdo,dmdu,dmdv,matin(4,4)
      double precision work1(4),work2(4),ap,zp,sigwp,matsm(2,2)
      double precision samesigm0,samesigsys,dm(20),sigm(20)
      double precision sigxsq,sigysq,sigxy,ellsum,elldiff
      double precision ellmajor,ellminor,elltan
      double precision sigbin(20)
      integer s,i,ss
      external fnr,fndo,fndu,fndv,ddofn,ddw0fn,ddwafn,dfn
      common om,q,u,v

      open(1,file='sninvar.txt')
      open(2,file='snoutput.txt')
      open(3,file='snerrout.txt')

  100 format(20f5.0)
  200 format(20f6.3)
  300 format(20f7.3)
   94 format(2f12.4)

c
c Input values:
c Line 1 - matter density Omega_m, w0, wa (fiducial 0.28, -1, 0)
c Line 2 - number of SN per 0.1 bin in redshift, put 0 for
c   empty bins, put z<0.1 SN in 1st bin
c   SNAP fiducial is 300,35,64,95,124,150,171,183,179,170,
c   155,142,130,119,107,94,80,0,0,0
c Line 3 -Toggle (1 yes, 0 no) 3 values: Is intrinsic uncertainty
c   of 1 SN the same for all redshifts? Is systematic error level
c   given by dm=sigsys*(1+z)/2.7 with sigsys the same for all redshifts?
c   Do you want to include CMB Fisher matrix (see below)?
c Line 4 - intrinsic uncertainty on 1 SN: either 1 number if constant
c   with redshift, or 20 numbers giving sigm0 for each 0.1 bin in redshift
c   Standard is all sigm0=0.1
c Line 5 - systematic error level: either 1 number for sigsys if follow
c   dm=sigsys*(1+z)/2.7, or 20 numbers giving dm (not sigsys!) for each
c   0.1 bin in redshift.
c   Systematic error is added in quadrature with intrinsic uncertainty.
c   Standard is all sigsys=0.02 (so dm ~ 1+z)
c CMB Fisher matrix approximates Planck as 0.2% accuracy on shift
c   parameter R (here called d)
c   This is a good approximation only when combined with SN
	read(1,*)om,w0,wa
	read(1,*)(n(i),i=1,20)
        read(1,*)samesigm0,samesigsys,yescmb
	if(samesigm0.eq.1.d0) then
           read(1,*)sigm0
           do i=1,20
              sigm(i)=sigm0
           enddo
        else
           read(1,*)(sigm(i),i=1,20)
        endif
	if(samesigsys.eq.1.d0) then
           read(1,*)sigsys
           do i=1,20
              dm(i)=sigsys*(1.d0+0.1d0*i-0.05d0)/2.7d0
           enddo
        else
           read(1,*)(dm(i),i=1,20)
        endif
c Write input to output file for a check
        write(2,*)'Reprint input to check'
        write(2,*)'Omega_m, w0, wa'
        write(2,*)om,w0,wa
        write(2,*)'Number of SN in each 0.1 bin of redshift'
        write(2,100)(n(i),i=1,20)
        write(2,*)'Intrinsic error constant? Systematic error constant?'
        write(2,*)'Include CMB?'
        write(2,*)samesigm0,samesigsys,yescmb
        write(2,*)'Intrinsic error in each 0.1 bin of redshift'
        write(2,200)(sigm(i),i=1,20)
        write(2,*)'Systematic error dm in each 0.1 bin of redshift'
        write(2,300)(dm(i),i=1,20)
	q=1.d0-om   ! dark energy density - assumes flatness
        u=w0
        v=wa
      s=4   ! dimension of Fisher matrix (scriptM,Omega_m,w0,wa)
c Initialize
      det=0.d0
      fmm=0.d0
      fmo=0.d0
	fmu=0.d0
	fmv=0.d0
      foo=0.d0
	fou=0.d0
	fov=0.d0
	fuu=0.d0
      fuv=0.d0
	fvv=0.d0
c  Calculate CMB Fisher matrix as reduced distance to last scatter
      botc=1.d0
      topc=1090.d0
	call dqromb(dfn,botc,topc,d)
	sigd=0.002d0*d  ! 0.2% precision on shift parameter d(=R)
         call dqromb(ddofn,botc,topc,ddo)
         call dqromb(ddw0fn,botc,topc,ddw0)
         call dqromb(ddwafn,botc,topc,ddwa)
	ddo=0.5d0/om/om*ddo
	ddw0=-1.5d0*q/om*ddw0
	ddwa=-1.5d0*q/om*ddwa
         fcoo=yescmb*ddo*ddo/sigd**2.
         fcou=yescmb*ddo*ddw0/sigd**2.
	   fcov=yescmb*ddo*ddwa/sigd**2.
	   fcuu=yescmb*ddw0*ddw0/sigd**2.
	   fcuv=yescmb*ddw0*ddwa/sigd**2.
	   fcvv=yescmb*ddwa*ddwa/sigd**2.
      foo=foo+fcoo
      fou=fou+fcou
      fov=fov+fcov
      fuu=fuu+fcuu
      fuv=fuv+fcuv
      fvv=fvv+fcvv
c  Begin SN Fisher calculation
      bot=1.d0
      do i=1,20
         z=0.1d0*i-0.05d0   ! all SN at redshift bin center
         y=1.d0+z
         sigm2(i)=sigm(i)*sigm(i)+n(i)*dm(i)*dm(i)

         sigbin(i)=dsqrt(sigm2(i)/n(i))
         write(3,94) z,sigbin(i)

         top=y
         call dqromb(fnr,bot,top,ansr)
         call dqromb(fndo,bot,top,ansdo)
         call dqromb(fndu,bot,top,ansdu)
         call dqromb(fndv,bot,top,ansdv)
         pre=5.d0/dlog(10.d0)/ansr   ! convert distance to mag
         dmdo=pre*ansdo
	   dmdu=pre*ansdu
	   dmdv=pre*ansdv
c Derivatives with respect to scriptM labeled by m,
c to Omega_m by o, to w0 by u, to wa by v
         fmm=fmm+n(i)/sigm2(i)
         fmo=fmo+n(i)*dmdo/sigm2(i)
         fmu=fmu+n(i)*dmdu/sigm2(i)
         fmv=fmv+n(i)*dmdv/sigm2(i)
         foo=foo+n(i)*dmdo*dmdo/sigm2(i)
         fou=fou+n(i)*dmdo*dmdu/sigm2(i)
         fov=fov+n(i)*dmdo*dmdv/sigm2(i)
         fuu=fuu+n(i)*dmdu*dmdu/sigm2(i)
         fuv=fuv+n(i)*dmdu*dmdv/sigm2(i)
         fvv=fvv+n(i)*dmdv*dmdv/sigm2(i)
      enddo
c Populate symmetric Fisher matrix matin(4,4)
      matin(1,1)=fmm
      matin(1,2)=fmo
      matin(1,3)=fmu
      matin(1,4)=fmv
      matin(2,2)=foo
      matin(2,3)=fou
      matin(2,4)=fov
      matin(3,3)=fuu
      matin(3,4)=fuv
      matin(4,4)=fvv
      matin(2,1)=matin(1,2)
      matin(3,1)=matin(1,3)
      matin(3,2)=matin(2,3)
      matin(4,1)=matin(1,4)
      matin(4,2)=matin(2,4)
      matin(4,3)=matin(3,4)
c Output Fisher matrix
      write(2,*)' '
      write(2,*)'OUTPUT:'
      write(2,*)'Fisher matrix in scriptM, Omega_m, w0, wa'
      write(2,90) (matin(1,j),j=1,s)
      write(2,90) (matin(2,j),j=1,s)
      write(2,90) (matin(3,j),j=1,s)
      write(2,90) (matin(4,j),j=1,s)
c Invert to get covariance matrix
      call matinv(s,matin,det,work1,work2)
c Calculate pivot redshift zp and minimum w(a) uncertainty
c sigwp also uncertainty on constant w (fixing wa=0)
      ap=1.d0+matin(3,4)/matin(4,4)
      zp=1.d0/ap-1.d0
      sigwp=dsqrt(matin(3,3)-(matin(3,4))**2./matin(4,4))
c Output parameter uncertainties: scriptM, Omega_m, w0, wa, wp
      write(2,*)'Uncertainties on scriptM, Omega_m, w0, wa, wp'
      write(2,90) (dsqrt(matin(i,i)),i=1,s),sigwp
c Marginalize over all parameters except w0, wa and reinvert
c to get 2D Fisher matrix for plotting contour ellipse
      matsm(1,1)=matin(3,3)
      matsm(1,2)=matin(3,4)
      matsm(2,1)=matsm(1,2)
      matsm(2,2)=matin(4,4)

c calculate error ellipse parameters
      sigxsq=matsm(1,1)
      sigysq=matsm(2,2)
      sigxy=matsm(1,2)
      ellsum=(sigxsq+sigysq)/2.d0
      elldiff=dsqrt(((sigxsq-sigysq)/2.d0)**2.+(sigxy)**2.)
      ellmajor=dsqrt(ellsum+elldiff)
      ellminor=dsqrt(ellsum-elldiff)
      elltan=(2.d0*sigxy)/(sigxsq-sigysq)
      elltheta=(atan(elltan)/2.d0)*57.296d0+90.d0
c write out error ellipse parameters
      write(2,*)
      write(2,*) 'error ellipse major axis, minor axis, theta'
      write(2,93) ellmajor,ellminor,elltheta

c  invert 2x2 error matrix
      ss=2
      call matinv(ss,matsm,det,work1,work2)
c Write area Figure of Merit = 1/COV[w0,wa] = 1/sigma(wp)*sigma(wa)
      write(2,91) 'FOM = ',dsqrt(1.d0/det)
c Output quantities needed for plotting contour ellipse of w0-wa
      write(2,*)'Values of w0, wa, reduced Fisher matrix F11, F12, F22'
      write(2,*) w0,wa,matsm(1,1),matsm(1,2),matsm(2,2)
 90   format(5e12.4)
 91   format(A,f8.2)
 93   format(3f8.3)



      stop
      end
c
c fnr calculates conformal distance
c fndX calculates analytic derivative of distance wrt parameter X
      function fnr(x)
        double precision fnr,x,om,q,u,v,ewfull,powarg,exparg
        common om,q,u,v
        powarg=3.d0*(1.d0+u+v)
        exparg=-3.d0*v*(1.d0-1.d0/x)
        ewfull=x**powarg*dexp(exparg)
        fnr=1.d0/dsqrt(om*x**3.+q*ewfull)
        return
        end

      function fndo(x)
        double precision om,q,u,v,w1,x,ewfull,fndo,powarg,exparg
        common om,q,u,v
        powarg=3.d0*(1.d0+u+v)
        exparg=-3.d0*v*(1.d0-1.d0/x)
        ewfull=x**powarg*dexp(exparg)
        fndo=(x**3.-ewfull)/(om*x**3.+q*ewfull)**1.5
        fndo=-fndo/2.d0
        return
        end

        function fndu(x)
        double precision om,q,u,v,x,ewfull,fndu
        double precision powarg,exparg
        common om,q,u,v
        powarg=3.d0*(1.d0+u+v)
        exparg=-3.d0*v*(1.d0-1.d0/x)
        ewfull=x**powarg*dexp(exparg)
        fndu=ewfull*dlog(x)/(om*x**3.+q*ewfull)**1.5
        fndu=-3.d0*q*fndu/2.d0
        return
        end
c
        function fndv(x)
        double precision om,q,u,v,x,ewfull,fndv
        double precision powarg,exparg
        common om,q,u,v
        powarg=3.d0*(1.d0+u+v)
        exparg=-3.d0*v*(1.d0-1.d0/x)
        ewfull=x**powarg*dexp(exparg)
        fndv=ewfull/(om*x**3.+q*ewfull)**1.5
        fndv=fndv*(dlog(x)-1.d0+1.d0/x)
        fndv=-3.d0*q*fndv/2.d0
        return
        end
c
c dfn calculates reduced distance to last scattering (shift param)
c ddXfn calculates analytic derivative of distance wrt parameter X

      function dfn(x)
      double precision om,dfn,x,q,w0,wa,powarg,exparg,ewfull,u,v
      common om,q,u,v
      w0=u
      wa=v
      powarg=3.d0*(1.d0+w0+wa)
      exparg=-3.d0*wa*(1.d0-1.d0/x)
      ewfull=x**powarg*dexp(exparg)
	dfn=1.d0/dsqrt(x**3.+q/om*ewfull)
      return
      end

      function ddofn(x)
      double precision om,q,x,w0,wa,ddofn,powarg,exparg,ewfull,u,v
      common om,q,u,v
      w0=u
      wa=v
      powarg=3.d0*(1.d0+w0+wa)
      exparg=-3.d0*wa*(1.d0-1.d0/x)
      ewfull=x**powarg*dexp(exparg)
      ddofn=ewfull/(x**3.+q/om*ewfull)**1.5
      return
      end

      function ddw0fn(x)
      double precision om,x,q,ddw0fn,w0,wa,powarg,exparg,ewfull,u,v
      common om,q,u,v
      w0=u
      wa=v
      powarg=3.d0*(1.d0+w0+wa)
      exparg=-3.d0*wa*(1.d0-1.d0/x)
      ewfull=x**powarg*dexp(exparg)
      ddw0fn=ewfull*dlog(x)/(x**3.+q/om*ewfull)**1.5
      return
      end

      function ddwafn(x)
      double precision om,x,q,ddwafn,w0,wa,powarg,exparg,ewfull,u,v
      common om,q,u,v
      w0=u
      wa=v
      powarg=3.d0*(1.d0+w0+wa)
      exparg=-3.d0*wa*(1.d0-1.d0/x)
      ewfull=x**powarg*dexp(exparg)
      ddwafn=ewfull*(dlog(x)-1.d0+1.d0/x)/(x**3.+q/om*ewfull)**1.5
      return
      end


      SUBROUTINE dqromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-10, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      pause 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                                      C
C  Invert a symmetric matrix and calculate its determinant.            C
C                                                                      C
C                                                                      C
C  To call:      CALL MATINV(M,ARRAY,DET,W1,W2)   where                C
C                                                                      C
C                                                                      C
C  M       : dimension of ...                                          C
C  ARRAY   : input matrix which is replaced by its inverse.            C
C  NORDER  : degree of matrix (order of determinant)                   C
C  DET     : determinant of input matrix.                              C
C  W1, W2  : work vectors of dimension M.                              C
C                                                                      C
C                                                                      C
C  Reference: Philip B Bevington, "Data Reduction and Error Analysis   C
C             for the Physical Sciences", McGraw-Hill, New York, 1969, C
C             pp. 300-303.                                             C
C                                                                      C
C                                                                      C
C  F. Murtagh, ST-ECF, Garching-bei-Muenchen, January 1986.            C
C                                                                      C
C----------------------------------------------------------------------C
        SUBROUTINE MATINV(M,ARRAY,DET,IK,JK)
        double precision    ARRAY(M,M), IK(M), JK(M), det
C
   10   DET = 1.d0
   11   DO 100 K = 1, M
C       Find largest element ARRAY(I,J) in rest of matrix.
        AMAX = 0.d0
   21      DO 30 I = K, M
              DO 30 J = K, M
   23            IF (ABS(AMAX)-ABS(ARRAY(I,J))) 24,24,30
   24            AMAX = ARRAY(I,J)
                 IK(K) = I
                 JK(K) = J
   30      CONTINUE
C          Interchange rows and columns to put AMAX in ARRAY(K,K).
   31      IF (AMAX) 41,32,41
   32      DET = 0.d0
           GOTO 140
   41      I = IK(K)
           IF (I-K) 21,51,43
   43      DO 50 J = 1, M
              SAVE = ARRAY(K,J)
              ARRAY(K,J) = ARRAY(I,J)
   50      ARRAY(I,J) = -SAVE
   51      J = JK(K)
           IF (J-K) 21,61,53
   53      DO 60 I = 1, M
              SAVE = ARRAY(I,K)
              ARRAY(I,K) = ARRAY(I,J)
   60      ARRAY(I,J) = -SAVE
C          Accumulate elements of inverse matrix.
   61      DO 70 I = 1, M
              IF (I-K) 63,70,63
   63         ARRAY(I,K) = -ARRAY(I,K)/AMAX
   70      CONTINUE
   71      DO 80 I = 1, M
              DO 80 J = 1, M
                 IF (I-K) 74,80,74
   74            IF (J-K) 75,80,75
   75            ARRAY(I,J) = ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
   80      CONTINUE
   81      DO 90 J = 1, M
              IF (J-K) 83,90,83
   83         ARRAY(K,J) = ARRAY(K,J)/AMAX
   90      CONTINUE
           ARRAY(K,K) = 1.d0/AMAX
  100   DET = DET * AMAX
C       Restore ordering of matrix.
  101   DO 130 L = 1, M
           K = M - L + 1
           J = IK(K)
           IF (J-K) 111,111,105
  105      DO 110 I = 1, M
              SAVE = ARRAY(I,K)
              ARRAY(I,K) = -ARRAY(I,J)
  110      ARRAY(I,J) = SAVE
  111      I = JK(K)
           IF (I-K) 130,130,113
  113      DO 120 J = 1, M
              SAVE = ARRAY(K,J)
              ARRAY(K,J) = -ARRAY(I,J)
  120      ARRAY(I,J) = SAVE
  130   CONTINUE
  140   RETURN
        END

