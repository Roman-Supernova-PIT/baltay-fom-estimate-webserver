      program snflux

c program to calculate supernova fluxes
c starting with the Hsiao typical Type 1a spectrum

      common/parms/omegam,omegak,omegade,h0,w,w0,w1
      external zint

      real*8 wl(2000),flux(2000),sig(20)
      real*8 filtmin(5),filtmax(5)
      real*8 rfmin(5,20),rfmax(5,20)
      integer imin(5,20),imax(5,20)
      real*8 signal(5,20)
     

c set parameters

      c=3.0e5
      h0=70.0
      signal=0.0
      hp=6.67e-34

      telarea=0.78*(3.14*(120.0**2))
      thruput=0.70
      

c thruput includes detector quantum efficiency


c normalize Hsiao template spectra B=0 to B=-19.7

      anorm=0.759e8


c set cosmology

      omegam=0.28
      omegak=0.0
      omegade=0.72
      w=-1.0

c setup output

      open(1,file='snfluxes.d')
      open(2,file='hsiao.txt')
      open(3,file='fluxout.txt')

  100 format(4f10.2)
  101 format(5f12.3)
  102 format(f10.2,3f10.4,i5)
  103 format(5f14.3)

      write(3,*) '   '
      write(3,*) '   '
      write(3,*) '    omegam    omegak   omegade       w  '
      write(3,100) omegam,omegak,omegade,w
      write(3,*) '   '


c read in Hsiao spectrum in rest frame

      do 1 i=1,1800

      read(2,*)a,wl(i),flux(i)
      flux(i)=flux(i)*anorm
     

c converting ergs/cm2/sec/A to photons/cm2/sec/10Abin

      egamma=((hp*3.0e8)/(wl(i)*1.0e-10))*1.0e7
      flux(i)=(flux(i)/egamma)*10.0
   1  continue
     
c set min and max wavelength range on Hsiao spectrum

      filtmin(1)=0.760
      filtmin(2)=0.927
      filtmin(3)=1.131
      filtmin(4)=1.380
      filtmin(5)=0.500

      filtmax(1)=0.977
      filtmax(2)=1.192
      filtmax(3)=1.454
      filtmax(4)=1.774
      filtmax(5)=0.600

      do 30 j=1,4
      do 31 i=1,20

      z=i*0.1-0.05

      rfmin(j,i)=filtmin(j)/(1.0+z)
      rfmax(j,i)=filtmax(j)/(1.0+z)

      imin(j,i)=1000.0*rfmin(j,i)-100.0
      imax(j,i)=1000.0*rfmax(j,i)-100.0

         
      rfmin(5,i)=filtmin(5)
      rfmax(5,i)=filtmax(5)

      imin(5,i)=1000.0*rfmin(5,i)-100.0
      imax(5,i)=1000.0*rfmax(5,i)-100.0

  31  continue
  30  continue

    
c  add up flux between limits in restframe

      do 3 j=1,5
      
      do 2 i=1,20

      signal(j,i)=0

      ia=imin(j,i)
      ib=imax(j,i)

      do 5 k=ia,ib 

      signal(j,i)=signal(j,i)+flux(k)

   5  continue

   2  continue

   3  continue

 
c signal is in photons/cm2/sec/1000A band
c in the supernova rest frame
c convert flux from rest frame to observer frame

c loop over z

      do 11 j=1,5

      do 10 i=1,20

      zsne=i*0.10-0.05

c calculate luminosity distance dlum

c calculate integral of 1/E(z)

      sum=0.0
      dz=zsne/100
      absm=-19.7

      do 20 k=1,100

      ak=float(k)
      z=(ak-0.5)*dz
      sum=sum+zint(z)*dz

   20 continue

      dist=sum

      dlum=(c/h0)*(1.0+zsne)*dist
      distmod=5.0*log10(dlum)+25.0


      signal(j,i)=signal(j,i)/((dlum/1.0e-5)**2)
      signal(j,i)=signal(j,i)*telarea*thruput

c signal(j,i) is in counts/sec/restframe band in observer frame


   10 continue

   11 continue

      do 77 i=1,20

      write(1,101)(signal(j,i),j=1,5)

 
  77  continue

  78  continue



      end



      function zint(z)

      common/parms/omegam,omegak,omegade,h0,w,w0,w1


      eofz=sqrt(omegam*(1.0+z)**3+omegak*(1.0+z)**2+omegade*(1.0+z)**
     &(3.0*(1.0+w)))

      zint=1.0/eofz

      end

