c     program  hom4_3qc(2a),   06.09.98
c     .........................................................................
c     .  Light Scattering by Spheroidal Particles
c     .  by N.V. Voshchinnikov and V.G. Farafonov
c     .  copyright (c) 1996/98 Astronomical Institute, St. Petersburg Univ.
c     .........................................................................
c     .........................................................................
c     .  Calculations of extinction, scattering and absorption efficiency
c     .  factors for prolate or oblate homogeneous spheroids.
c     .
c     .  Description of method, equations, and some results see:
c     .  N.V. Voshchinnikov and V.G. Farafonov. Optical properties of
c     .  spheroidal particles, Ap. Space Sci. 204, 19-86, 1993.
c     .........................................................................
c     .  Parameters:
c     .      nterms: number of terms in Eqs. for efficiency factors,
c     .     nalphas: number of angles for incident radiation,
c     .........................................................................in
c     .  Input data:          file:   homn.dat
c     .            k: parameter describing type of particle,
c     .                           = 0 for prolate spheroids,
c     .                           = 1 for oblate spheroids,
c     . ri = n + k*i: complex index of refraction,
c     .          eps: relative error for calculations of eigenvalues,
c     .         eps1: relative error for calculations of spheroidal functions,
c     .         eps2: relative error for calculations of efficiency factors,
c     .          wvl: wavelength of radiation (in mkm),
c     .          icx: parameter describing size of spheroids,
c     .               = 1 - Rv (radius of equivolume sphere, in mkm),
c     .               = 2 - Xv (size parameter of equivolume sphere),
c     .               = 3 - c  (size parameter proportional to focus distance),
c     .               = 4 - 2*pi*a/wvl (=bcksi),
c     .               = 5 - Xc (size parameter for cylinders),
c     .          NRv: number of sizes of spheroids,
c     .         aRv0: initial size,
c     .          dRv: step over sizes,
c     .         nalf: number of angles for incident radiation 'alpha',
c     .         ial0: initial angle 'alpha',
c     .         idal: step over angles 'alpha',
c     .          a/b: aspect ratio for spheroid (value 1)
c     .          a/b: aspect ratio for spheroid (value 2)
c     .          ........................................
c     .........................................................................
c     .  Output data:         file:   hom.res
c     .   Cext/(Pi*Rv**2) : normalized extinction cross-sections,
c     .   Csca/(Pi*Rv**2) : normalized scattering cross-sections,
c     .   Cabs/(Pi*Rv**2) : normalized absorption cross-sections.
c     .....................
c     .                       file:   h.dat
c     .   Qext(TM,TE) : extinction efficiency factors,
c     .   Qsca(TM,TE) : scattering efficiency factors,
c     .   Qabs(TM,TE) : absorption efficiency factors.
c     .....................
c     .                       file:   h1.dat
c     .        (non-polarized incident radiation, C=C/(Pi*Rv**2))
c     .   Cext=1/2*(Cext(TM)+Cext(TE)) : extinction cross-section,
c     .   Csca=1/2*(Csca(TM)+Csca(TE)) : scattering cross-section,
c     .   Csca=1/2*(Cabs(TM)+Cabs(TE)) : absorption cross-section,
c     .   Alb                          : albedo,
c     .   Cpol                         : polarization cross-section,
c     .   Pol                          : polarization (in %).
c     .........................................................................
c     .  Intermediate data:   file:  hom.out
c     .   Qext: extinction efficiency factors,
c     .   Qsca: scattering efficiency factors,
c     .   Qabs: absorption efficiency factors.
*     .........................................................................
*     .   Functions and subroutines used (alphabetically):
*     .      functions: AKAPPA, COFF,
*     .                 OMEGA,  TAU,
*     .                 dreal (Attention ---> for Lahey only)
*     .      subroutines: asig11,    BESSJJ,    cdcof3,
*     .                   CDLAMn,    CDLAMo_l,  CDLAMo_s,
*     .                   CDLAMp,    CDRB12,    CDRf12,
*     .                   CESSEL,    CMLM1,     cmlm2,
*     .                   CXSC11,    DASWF1,    delta21,
*     .                   DLEGF00,   DLEGF1,    DLEGF2,
*     .                   DRF222,    DRF2k02,
*     .                   DRKF45,    DRKFS,
*     .                   FEHL,      FUNLEG,    FUNLEG1,
*     .                   gamkap1,   gj11,      homfunq,
*     .                   INT3,      INVMAT,
*     .                   Lambda,    sphtm,     sphte
c     .........................................................................
      parameter (nterms=100, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z), COMPLEX*16 (R-S)
      REAL*8 KSI0
      COMPLEX*16 C2

      DIMENSION Ialpha(nalphas), Jalpha(nalphas), Malpha(nalphas),
     *          ASIN(nalphas), ACOS(nalphas), FACT(170),
     *          QXM(nalphas),  QSM(nalphas),
     *          QXMM(nalphas), QSMM(nalphas)
     *         ,QXE(nalphas),  QSE(nalphas),
     *          QXEM(nalphas), QSEM(nalphas)

      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /K2/ Ialpha, Jalpha
      COMMON /EPS1/ EPS1
      COMMON /EPS3/ EPS3
      COMMON /icc/ i1
      COMMON /FF/ ASIN,ACOS
      COMMON /FACT/ FACT
      COMMON /PI/ PI
  100 FORMAT(2F8.4)
  101 FORMAT(I4/(F8.4))
  102 FORMAT(I4)
  200 FORMAT(1X,'Refractive index:',2F8.4,1X,'I'/
     *       1x,'Rv=',F6.3,1X,'mkm',2X,'Lambda=',F6.3,1X,'mkm',
     *       2x,'Xv=',F7.3/1X,'a/b =',F7.3,4X,'Ksi0=',F6.3)
 1200 FORMAT(1X,'Refractive index:',2F8.4,1X,'i',5x,
     *       'Lambda=',F6.3,1X,'mkm',5x,'a/b =',F7.3)
  201 FORMAT(37X,4(D16.8,2X),D10.3,/)
  202 FORMAT(D12.4)
  203 FORMAT(1X,'SPHEROID: Prolate') 
  204 FORMAT(1X,78('*'))
  205 FORMAT(1X,78('='))
 1205 FORMAT(1X,96('='))
  206 FORMAT(11X,'C',4X,'cksi0',2X,'Nmax',1X,'alpha',5X,'Qext',
     *       10X,'Qsca',10X,'Qabs')
 1206 FORMAT(4X,'c',4X,'x_V',1X,'2*pi*a/w',1X,'alpha',2X,
     *       'Qext(TM)',4X,'Qsca(TM)',4X,'Qabs(TM)',4x,
     *       'Qext(TE)',4X,'Qsca(TE)',4X,'Qabs(TE)')
 1207 FORMAT(4X,'c',4X,'x_V',1X,'2*pi*a/w',1X,'alpha',1X,
     *       'Cext/Pi/Rv**2',2X,'Csca/Pi/Rv**2',2X,'Cabs/Pi/Rv**2',
     *        2x,'Alb',5X,'Cpol',4X,'Pol(%)')
  207 FORMAT(1X,'M=',I2,2X,2(F6.2,1X),2(I3,2X),1X,3(1pD13.6,1x),
     *       1X,'TM')
  277 FORMAT(32X,3(1pD13.6,1x))
  208 FORMAT(1X,'SPHEROID: Oblate') 
  209 FORMAT(1X,78('.'))
  211 FORMAT(1X,2('+'),'K=',I4,2X,'RI,C1,C2=',5F10.4)
  212 FORMAT(1X,'alpha:',I4)
  213 FORMAT(1X,78('-'))
  214 FORMAT(1X,'M=',I2,2X,2(F6.2,1X),2(I3,2X),1X,3(1pD13.6,1x),
     *       1X,'TE')
  219 FORMAT(2X,' ********* N>4*nterms **** NK=',I8)
  220 FORMAT(2X,' ******* NAL>nalphas *** NAL=',I8)
  221 FORMAT(1X,78('.'))
  222 FORMAT(2X,' *******   M>15  **********  M=',I8)
  223 FORMAT(1X,'C1=',F6.3,2X,'C2=',2F6.3,2X,'Cksi=',F6.3,2X,
     *       'Nmax=',I3,2X,'Nmxe=',I3)
  224 FORMAT(1X,'alpha',4x,'m',6X,'Qext',
     *       11X,'Qsca',11X,'Qabs') 
  225 FORMAT(1X,I3,4X,I2,3(1pD14.6,1X),1X,'TM')
 1224 FORMAT(1X,'alpha',3x,'m',1X,'Cext/(Pi*Rv**2)',
     *       1X,'Csca/(Pi*Rv**2)',1X,'Cabs/(Pi*Rv**2)')
 1225 FORMAT(1X,3f6.2,i6,6(1pD12.4))
12252 FORMAT(1X,3f6.2,i6,1x,1pD12.4,2(3x,d12.4),
     *       1X,0pf8.2,1pD12.4,0pf6.2)
12251 FORMAT(1X,f5.2,6(1pD13.5))
  226 FORMAT(1X,I3,4X,I2,3(1pD14.6,1X),1X,'TE')
  227 FORMAT(10X,I4,10X,5(1pD12.3,1X),3X,2(D12.3,7X))
  228 FORMAT(1X,'*** Nmax > nterms *** (',I3,')'/
     *       1x,'Continue ? (0 - Yes, 1 - No)  ')
  516 FORMAT(1X,78('*'))
c     .........................................................................
                   print *,'Hom_SPH v.4_3q2a'
c     .........................................................................
      open (unit=05,file='homn.dat',status='old',access='sequential')
      open (unit=07,file='hom.out',status='unknown',access='append')
      open (unit=09,file='hom.res',status='unknown',access='append')
      open (unit=10,file='h.dat',status='unknown',access='append')
      open (unit=11,file='h1.dat',status='unknown',access='append')
c      open (unit=11,file='hom.t1',status='unknown',
c     *      access='sequential')
c      open (unit=12,file='hom.t2',status='unknown',
c     *      access='sequential')
c     .........................................................................
c     .   set program constants
c     .........................................................................
       S=(0D0,0D0)
       S1=(0D0,1D0)
       PI=4d0 * datan(1d0)
       DPI=2D0*PI
       PI180=PI/180D0
       eps3 = 1d-80
*       eps3 = 1d-180
       FACT(1)=1D0
       FACT(2)=1D0
       DO 130 I=3,170
  130  FACT(I)=FACT(I-1)*(I-1D0)
c     ..........................................................................
c     .   input data
c     .........................................................................
      READ(5,102) K
      READ(5,100) RI
      READ(5,202) EPS,EPS1,eps2
      READ(5,100) WVL
      WRITE(7,202) EPS,EPS1,eps2
      READ (5,102) ICX
      GO TO (300,301,302,303,404),ICX
  300       READ (5,101) NC,ARV0,DRV
            WRITE (7,101) NC,ARV0,DRV,WVL
       GO TO 304
  301       READ (5,101) NC, XV0,DXV
            WRITE (7,101) NC, XV0,DXV,WVL
       GO TO 304
  302       READ (5,101) NC,  C0,DC
            WRITE (7,101) NC,  C0,DC ,WVL
       GO TO 304
  303       READ (5,101) NC,CKSI0,DCKSI
            WRITE (7,101) NC,CKSI0,DCKSI,WVL
       GO TO 304
  404       READ (5,101) NC,XC0,DXC
            WRITE (7,101) NC,  XC0,DXC  ,WVL
  304 READ(5,102) NALf,IAL0,IDAL 
      WRITE(7,102) NALF,IAL0,IDAL
      NAL=NALF
      IF(IAL0.EQ.0) NAL= NALF-1
  999 READ(5,100,END=1000) AB
      WRITE(7,101) K,AB
c     .........................................................................
c     .   calculate program constants
c     .........................................................................
      IF(K) 2,3,2
    2 CONTINUE
      KSI0=1D0/DSQRT(AB**2-1D0)
c                                    ksi0=0.2d0
c                                    ksi0=0d0
      AKSI=KSI0**2+1D0
      AK = -1D0
      GO TO 4
    3 CONTINUE
      KSI0=DSQRT(1.0D0+1.0D0/(AB**2-1.0D0))
      AKSI=KSI0**2-1D0
      AK = 1D0
    4 CONTINUE
      WRITE(7,204)
      IF(NAL.GT.nalphas) WRITE(*,220) NAL
      IF(NAL.GT.nalphas) WRITE(7,220) NAL
      IF(NAL.GT.nalphas) STOP
      IAL=IAL0
      IF(NAL.EQ.0) GO TO 666
      IF(IAL.EQ.0) GO TO 66
      IAL=IAL0-IDAL
   66 CONTINUE
      DO 6 I=1,NAL
      IAL=IAL+IDAL
      ASIN(I)=DSIN(PI180*IAL)
      ACOS(I)=DCOS(PI180*IAL)
      IF(IAL.EQ.90) ACOS(I)=0D0
    6 Ialpha(I)=IAL
  666 CONTINUE
      if(k.eq.0) WRITE(10,203)
      if(k.eq.1) WRITE(10,208)
      WRITE(10,1200) RI,wvl,AB
      WRITE(10,1205)
      WRITE(10,1206)
      WRITE(10,1205)
      if(k.eq.0) WRITE(11,203)
      if(k.eq.1) WRITE(11,208)
      WRITE(11,1200) RI,wvl,AB
      WRITE(11,1205)
      WRITE(11,1207)
      WRITE(11,1205)
c     .........................................................................
c     .   start of calculations over size parameter
c     .........................................................................
      GO TO (305,306,307,308,409),ICX
  305       ARV=ARV0-DRV
            GO TO 309
  306       XV= XV0-DXV
            GO TO 309
  307       C1=  C0-DC
            GO TO 309
  308       CKSI=CKSI0-DCKSI
            dC=dCKSI/KSI0
            IF(K.EQ.1) dC=dC/AB
            GO TO 309
  409       XC  =  XC0-DXC
  309 CONTINUE
      DO 1 I1=1,Nc

*--------------------------------------------------------
c     timer

c      character*11 itime0,itime1
c      call time(itime0)
c      write(7,*) ' System time: ',itime0
c      call timer(itty0)
*--------------------------------------------------------

c            rewind 11
c            rewind 12
      GO TO (310,311,312,313,414),ICX
  310       ARV=ARV+DRV
            XV=DPI*ARV/WVL
            GO TO 314
  311       XV=XV+DXV
            ARV=XV*WVL/DPI
  314       CONTINUE
            IF(K.EQ.0) CKSI=XV*(AB)**(2D0/3D0)
            IF(K.EQ.1) CKSI=XV/(AB)**(2D0/3D0)
            C1=CKSI/KSI0
            XC=CKSI/AB
            GO TO 315
  312       C1 =C1 + DC
            CKSI=C1*KSI0
            XC=CKSI/AB
            GO TO 316
  313       CKSI=CKSI+DCksi
            C1=CKSI/KSI0
            IF(K.EQ.1) C1=C1/AB
            XC=CKSI/AB
            GO TO 316
  414       XC=XC+DXC
            CKSI=XC*AB
            C1=CKSI/KSI0
            IF(K.EQ.1) C1=C1/AB
  316       CONTINUE
            IF(K.EQ.0)   XV=CKSI/(AB)**(2D0/3D0)
            IF(K.EQ.1)   XV=C1*KSI0*(AB)**(2D0/3D0)
            ARV=XV*WVL/DPI
  315 C2=RI*C1
      if(k.eq.0) bcksi = cksi
      if(k.eq.1) bcksi = c1 * ksi0 * ab
c      NMAX=1.2*CDABS(C2)+3
      NMAX1=1.2*CDABS(C2) + 7

       NMAX = idint(bcksi) + 5
       if(nmax1.gt.nmax) nmax = nmax1

       IF(k.eq.0.and.KSI0.GT.5D0) NMAX = nmax + 10
****      IF(k.eq.0.and.KSI0.GT.1.8D0) NMAX=idint(xv)+6
****      IF(k.eq.1.and.KSI0.GT.0.8D0) NMAX=idint(xv)+6
****      IF(k.eq.1.and.ksi0.gt.0.5D0.and.c1.gt.2d0) NMAX=nmax+6
*                                 NmAX=nmax+10
*                                 NmAX=nmax+10
*                                 NmAX=nmax+10
      IF(MOD(NMAX,2)) 33,34,33
   33 NMAX=NMAX+1
   34 CONTINUE
*       nmax=nterms
        if(nmax.gt.nterms) nmax=nterms

                    M=1

c      NK=1.5*NMAX+10
      NK=NMAX+40
      IF(NK.LT.60) NK=60
      IF(Ksi0.gt.1.5d0) NK=nk+40
c                    NK=4*nterms
      IF(NK.GT.4*nterms) NK=4*nterms
      IF(K.EQ.0) WRITE(7,203)
      IF(K.EQ.1) WRITE(7,208)
      WRITE(7,200)  RI,ARV,WVL,XV,AB,KSI0
      WRITE(7,205)
      WRITE(7,206)
      WRITE(7,213)
      IF(K.EQ.0) WRITE(*,203)
      IF(K.EQ.1) WRITE(*,208)
      WRITE(*,200)  RI,ARV,WVL,XV,AB,KSI0
      WRITE(*,205)
      WRITE(*,206)
      WRITE(*,213)
      IF(NMAX.le.nterms) go to 557
      IF(NMAX.GT.nterms) NMAX=nterms
  557 CONTINUE
      NMXE=NMAX+4
      IF(CDABS(C2).LT.4D0) NMXE=NMXE+2
c      IF(AB.GT.4D0.and.k.eq.1) NMXe=NMXe+4
      IF(AB.GT.4D0.and.k.eq.1) NMXe=NMXe+8
      IF(c1.lT.4D0.and.k.eq.1) NMXe=NMXe+10
      IF(NMXE.GT.nterms) NMXE=nterms
      if(ksi0.gt.2d0)  nmxe = nmax + 8
      if(cksi.lt.3d0.and.ksi0.gt.2d0)  nmxe=nmax
      if(ksi0.gt.8d0)  nmxe=nmax

c           nmxe=nmax

      Jalpha(1)=0
      IF(IAL0.EQ.0) GO TO 50
   54 CONTINUE
      DO 7 I=1,NAL
      QXE(I)=0D0
      QXM(I)=0D0
      QSE(I)=0D0
      QSM(I)=0D0
      Jalpha(I)=0
    7 CONTINUE
      JAM=0
c     .........................................................................
c     .   alpha =/ 0 deg.
c     .   calculate efficiency factors
c     .........................................................................
   11 CONTINUE
            WRITE(*,516)
c----------
c  TM mode
c----------
      CALL homfunq(1,m,C1,C2,KSI0,EPS,NMXE,*12)
      CALL sPHtM(1,m,RI,C1,KSI0,NMXE,NMAX)
      CALL CXSC11(1,M,C1,NMXE,NMAX,QXMM,QSMM,QEXT,QSCA)

      JM=JAM+1
         DO  I = JM, NAL

      QXM(I)=QXM(I)+QXMM(I)
      QSM(I)=QSM(I)+QSMM(I)
      QABS=QXM(I)-QSM(I)
      QABSm=QXMm(I)-QSMm(I)

       WRITE(7,207) M,C1,bCKSI,NMAX,Ialpha(I),QXM(I),QSM(I),QABS

            WRITE(7,277) QXmM(I),QSmM(I),QABSm

c       WRITE(10,207) M,C1,bCKSI,NMAX,Ialpha(I),QXM(I),QSM(I),QABS

c            WRITE(10,277) QXmM(I),QSmM(I),QABSm
            WRITE(*,277) QXmM(I),QSmM(I),QABSm

       WRITE(*,207) M,C1,bCKSI,NMAX,Ialpha(I),QXM(I),QSM(I),QABS
         end do

c----------
c  TE mode
c----------
      CALL sPHTE(M,C1,KSI0,NMXE)
      CALL CXSC11(1,M,C1,NMXE,NMXE,QXEM,QSEM,QEXT,QSCA)

         DO  I = JM, NAL

      QXE(I) = QXE(I) + QXEM(I)
      QSE(I) = QSE(I) + QSEM(I)
      QABS = QXE(I) - QSE(I)
      QABSm = QXEm(I) - QSEm(I)
      WRITE(7,214) M,C1,bCKSI,NMXE,Ialpha(I),QXE(I),QSE(I),QABS
      WRITE(7,277) QXeM(I),QSeM(I),QABSm
c            WRITE(10,277) QXeM(I),QSeM(I),QABSm
c            WRITE(10,221)
            WRITE(*,277) QXeM(I),QSeM(I),QABSm

        WRITE(*,214) M,C1,bCKSI,NMXE,Ialpha(I),QXE(I),QSE(I),QABS

      IF(DABS(QXMM(I)/QXM(I)).LE.eps2.AND.DABS(QXEM(I)/QXE(I)).LE.eps2)
c      IF(DABS(QXMM(I)/QXM(I)).LE.eps2)
     *          Jalpha(I)=1
      JAM=JAM+Jalpha(I)
                     Malpha(I)=M
      IF(JAM.EQ.NAL) GO TO 13

            end do

      WRITE(7,221)
      M=M+1
c      M=M+3

c            if(m.gt.25) go to 13
c            if(m.gt.25) print *, '     m > 25 !!!'

      GO TO 11
   13 WRITE(7,213)
      WRITE(*,213)
      GO TO 53
   50 CONTINUE
c     .........................................................................
c     .   alpha = 0 deg.
c     .........................................................................
c        nmax = nmax + 10
      QEX1=1D0
  558 CONTINUE
      IF(NMax.gt.nterms) NMax=nterms
c  612 CONTINUE
   55 CONTINUE   
      if(nal.eq.0) NMXE=NMAX
c     .........................................................................
c     .   calculate efficiency factors
c     .........................................................................
      CALL homfunq(0,1,C1,C2,KSI0,EPS,NMXE,*12)
      CALL sPHtM(0,1,RI,C1,KSI0,NMXE,NMAX)
      CALL CXSC11(0,1,C1,NMXE,NMAX,QXMM,QSMM,QEXT,QSCA)
      QAB=QEXT-QSCA

       WRITE(7,207) M,C1,bCKSI,NMAX,IAL0,QEXT,QSCA,QAB
       WRITE(*,207) M,C1,bCKSI,NMAX,IAL0,QEXT,QSCA,QAB

      IF(DABS(QEXT/QEX1-1D0).LE.eps2.and.qAB.gt.-0.1d0) GO TO 52
      IF(DABS(QAB).lt.eps2) GO TO 52
      IF(qAB.lt.-1d0) GO TO 1
      QEX1=QEXT
      NMAX=NMAX+2
      NMXE=NMAX+4
      IF(c1.lT.4D0.and.k.eq.1) NMXe=NMXe+30
      IF(NMXE.GT.nterms) NMXE=nterms
      IF(NMaX.GT.nterms) go to 52
****      GO TO 55
   52 CONTINUE
      WRITE(7,205)
      WRITE(*,205)
      IF(NAL.GT.0) GO TO 54
   53 CONTINUE
c     .........................................................................
c     .   print final results
c     .........................................................................
      IF(K.EQ.0) WRITE(9,203)                                          
      IF(K.EQ.1) WRITE(9,208)                                          
      WRITE(9,200) RI,ARV,WVL,XV,AB,KSI0
      WRITE(9,205)                                                     
      WRITE(9,223) C1,C2,bCKSI,NMAX,NMXE
      WRITE(9,213)                                                     
      WRITE(9,1224)
      WRITE(9,213)                                                     
      M0=1                                                             
      IF(k.EQ.0) gn0=1d0/ab**(2d0/3d0)
      IF(k.EQ.1) gn0=ab**(2d0/3d0)
      v1=qext*gn0
      v2=qsca*gn0
      v3=qab*gn0
c      IF(IAL0.EQ.0) WRITE(9,225) IAL0,M0,QEXT,QSCA,QAB
      IF(IAL0.EQ.0) WRITE(9,225) IAL0,M0,v1,v2,v3
      IF(IAL0.EQ.0) WRITE(10,1225) c1,xv,bcksi,ial0,QEXT,QSCA,QAB
      IF(IAL0.EQ.0) WRITE(11,12252) c1,xv,bcksi,ial0,v1,v2,v3,v2/v1
      IF(NAL.EQ.0) GO TO 114                                           
      DO 112 IA=1,NAL                                                  
      IF(k.EQ.0) gn=dsqrt(ab**2*asin(ia)**2+acos(ia)**2)/ab**(2d0/3d0)
      IF(k.EQ.1) gn=dsqrt(ab**2*acos(ia)**2+asin(ia)**2)/ab**(1d0/3d0)
      QABS=QXM(IA)-QSM(IA)                                             
      v1=qxm(ia)*gn
      v2=qsm(ia)*gn
      v3=qabs*gn
c      WRITE(9,225) Ialpha(IA),Malpha(IA),QXM(IA),QSM(IA),QABS
      WRITE(9,225) Ialpha(IA),Malpha(IA),v1,v2,v3
      WRITE(10,1225) c1,xv,bcksi,Ialpha(IA),QXM(IA),QSM(IA),QABS,
     *                   QXe(IA),QSe(IA),QXe(IA)-QSe(IA)
      qxm(ia)=v1
      qsm(ia)=v2
  112 CONTINUE                                                         
      WRITE(9,221)                                                     
      DO 412 IA=1,NAL                                                  
      IF(k.EQ.0) gn=dsqrt(ab**2*asin(ia)**2+acos(ia)**2)/ab**(2d0/3d0)
      IF(k.EQ.1) gn=dsqrt(ab**2*acos(ia)**2+asin(ia)**2)/ab**(1d0/3d0)
      QABS=QXE(IA)-QSE(IA)                                             
      v1=qxe(ia)*gn
      v2=qse(ia)*gn
      v3=qabs*gn
c      WRITE(9,226) Ialpha(IA),Malpha(IA),QXE(IA),QSE(IA),QABS
      WRITE(9,226) Ialpha(IA),Malpha(IA),v1,v2,v3
c         WRITE(10,1225) cksi,QXe(IA),QSe(IA),QABS
ccc      WRITE(10,1225) c1,xv,bcksi,Ialpha(IA),QXM(IA),QSM(IA),
ccc     *               qxm(ia)-qsm(ia), v1, v2, v3
         v4 = (qxm(ia)+v1)/2d0
         v44 = (qsm(ia)+v2)/2d0
         v5 = (qxm(ia)-v1)/2d0
      WRITE(11,12252) c1,xv,bcksi,Ialpha(IA),v4, v44,
     *               (qxm(ia)-qsm(ia)+v3)/2d0, v44/v4, v5, v5/v4*100d0
  412 CONTINUE                                                         
  114 WRITE(9,205)                                                   
      GO TO 5
   12 WRITE(7,211)  K, RI, C1, C2
      WRITE(10,211) K, RI, C1, C2
    5 CONTINUE

*--------------------------------------------------------
c     timer

c      call time(itime1)
c      write(7,*) ' System time: ',itime1
c
c      call timer(itty1)
c
c      full_time_sec=(itty1-itty0)/100.
c      i_time_hour=(full_time_sec)/3600.
c      i_time_min=(full_time_sec-i_time_hour*3600)/60.
c
c      arest_sec=full_time_sec-(i_time_hour*3600.+i_time_min*60.)
c
c      write (*,12121)  i_time_hour, i_time_min, arest_sec
c      write (7,12121)  i_time_hour, i_time_min, arest_sec
c12121 format(1x,'Running time: ',i2,' hrs',1x,i2,' min',1x,
c     * f4.1,' sec')
*--------------------------------------------------------
    1 CONTINUE
      WRITE(10,1205)
      WRITE(11,1205)
      WRITE(7,205)
      WRITE(*,205)
      GO TO 999
 1000 STOP
      END

c================================================
c                 FUNCTIONS
c================================================

c************************************************
C AKAPPA
C
      FUNCTION AKAPPA(M,L,N,NN,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      COMMON /EPS3/ EPS3
      IF(MOD(IABS(L-N),2)) 20,20,9
   20 AKAPPA =0D0
      RETURN
    9 CONTINUE
      M2=2*M
      IF(MOD(L-M,2))1,1,2
    1 I1=0
      I2=1
      K=M2
      E=A(1)*(-B(1)*M/(M2+3D0))*2D0
      GO TO 3
    2 I1=1
      I2=0
      K=M2+1
      E=A(1)*(B(1)*(M+2D0)/(M2+1D0)-B(2)*(M+1D0)*(M2+2D0)/
     *  (M2+5D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NN-1
      J=2*I-1-I2
      JM=J+M
      JM2=2*(J+M)
      F=1D0
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=A(I)*(B(I-I2)*J*(JM+1D0)/(JM2-1D0)-B(I+I1)*JM*(M2+J+1D0)/
     *  (JM2+3D0))*2D0/(JM2+1D0)*F
c      IF(DABS(G).LE.1D-160.and.i.gt.l+n.or.i.eq.nn-1) GO TO 6
      IF(DABS(G).LE.eps3.and.i.gt.l+n.or.i.eq.nn-1) GO TO 6
    5 E=E+G
    6 AKAPPA =E/C/D
      RETURN
      END

c************************************************
C COFF
C
      FUNCTION COFF(L,N,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N)
      IF(MOD(L-1,2))1,1,2
    1 I2=1
      GO TO 3
    2 I2=0
    3 D=0D0
      DO 4 I=1,N
      J=2*I-1-I2
    4 D=D+(J+1D0)*(J+2D0)*A(I)
      COFF=D
      RETURN
      END

c************************************************
C OMEGA
C
      FUNCTION OMEGA(M,L,N,NN,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      COMMON /EPS3/ EPS3
      IF(MOD(IABS(L-N),2)) 9,9,20
   20 OMEGA  =0D0
      RETURN
    9 CONTINUE
      M2=2*M
      IF(MOD(L-M,2))1,1,2
    1 I2=1
      K=M2
      E=A(1)*(B(2)*(M2+2D0)*(M2+1D0)/(M2+5D0)/(M2+3D0)-
     *  B(1)*2D0*(M+1D0)*(M2-1D0)/(M2+3D0)/(M2-1D0))*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E=A(1)*(B(2)*(M2+3D0)*(M2+2D0)/(M2+7D0)/(M2+5D0)-B(1)*2D0*
     *  (M+1)*(M2+1D0)/(M2+5D0)/(M2+1D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NN-1
      J=2*I-1-I2
      JM=J+M
      JM2=2*(J+M)
      F=1D0
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=A(I)*(B(I-1)*J*(J-1D0)/(JM2-1D0)/(JM2-3D0)+B(I+1)*(J+M2+2D0)*
     *  (J+M2+1D0)/(2D0*JM+5D0)/(2D0*JM+3D0)-B(I)*2D0*(JM*(JM+1D0)+M*M
     *  -1D0)/(2D0*JM+3D0)/(2D0*JM-1D0))*2D0/(JM2+1D0)*F
c      IF(DABS(G).LE.1D-160.and.i.gt.l+n.or.i.eq.nn-1) GO TO 6
      IF(DABS(G).LE.eps3.and.i.gt.l+n.or.i.eq.nn-1) GO TO 6
    5 E=E+G
    6 OMEGA =-E/C/D
      RETURN
      END

c************************************************
C TAU
C
      FUNCTION TAU(M,L,N,NN,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      COMMON /EPS3/ EPS3
      IF(MOD(IABS(L-N),2)) 9,9,20
   20 TAU    =0D0
      RETURN
    9 CONTINUE
      M2=2*M
      IF(MOD(L-M,2))1,1,2
    1 I2=1
      K=M2
      E=A(1)*B(1)*M*(M+1D0)*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E=A(1)*B(1)*(M+1D0)*(M+2D0)*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NN
      J=2*I-1-I2
      JM=J+M
      F=1D0
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=A(I)*B(I)*JM*(JM+1D0)*2D0/(2D0*JM+1D0)*F
c      IF(DABS(G).LE.1D-160.and.i.gt.l+n.or.i.eq.nn-1) GO TO 6
      IF(DABS(G).LE.eps3.and.i.gt.l+n.or.i.eq.nn-1) GO TO 6
    5 E=E+G
    6 TAU    =E/C/D
      RETURN
      END

c************************************************
C dreal ---- Attention ---> for Lahey (only)
C
c      FUNCTION dreal(r)
c      REAL*8 dreal
c      complex*16 r
c      dreal = dfloat(r)
c      RETURN
c      END

c================================================
c                 SUBROUTINES
c================================================

c************************************************
c aSIGma11
c
      SUBROUTINE aSIGma11(M,NN,SIG)
      parameter (nterms=100)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 s, s1
      DIMENSION A(4*nterms), B(4*nterms),
     *          aCOF(nterms,4*nterms), VV1(nterms), SIG(nterms,nterms)
      COMMON /EPS3/ EPS3
      COMMON /Kc1/ aCOF,VV1
      COMMON /K1/ S, S1, AKSI, AK, kkK, NK, nal
       M2=2*M
      DO 10 N=1,NN
        NM=N+M-1
      DO 12 I=1,NK
   12 B(I)=aCOF(N,I)
      C=VV1(N)
      DO 11 L=1,NN
        LM=L+M-1
      IF(MOD(IABS(LM-NM),2)) 9,9,20
    9 CONTINUE
      DO 13 I=1,NK
   13 A(I)=aCOF(L,I)
      D=VV1(L)
      IF(MOD(LM-M,2))1,1,2
    1 I2=1
      K=M2
      E=A(1)*(-B(2)*(M+1D0)*(M2+1D0)*(M2+2D0)/((M2+3D0)*(M2+5D0))
     *+B(1)*(2D0*M*M+3D0*M-2D0)/((M2-1D0)*(M2+3D0)))*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E=A(1)*(-B(2)*(M+2D0)*(M2+2D0)*(M2+3D0)/((M2+5D0)*(M2+7D0))+
     *B(1)*(2D0*M*M+9D0*M+4D0)/((M2+1D0)*(M2+5D0)))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NK
      J=2*I-1-I2
      JM=J+M
      JM2=2*JM
      F=1D0
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=A(I)*(B(I-1)*J*(J-1D0)*JM/(JM2-3D0)/(JM2-1D0)-B(I+1)*(JM+1D0)*
     *(J+M2+1D0)*(J+M2+2D0)/(JM2+3D0)/(JM2+5D0)+B(I)*(JM*3D0*(JM+1D0)-
     *M*M-2D0)/(JM2-1D0)/(JM2+3D0))*2D0/(JM2+1D0)*F
c      IF(CDABS(G).LE.1D-140.and.i.gt.nn) GO TO 6
      IF(DABS(G).LE.eps3.and.i.gt.nn) GO TO 6
    5 E=E+G
    6 SIG(N,L)=E/C/D
      GO TO 11
   20 SIG(N,L)=(0D0,0D0)
   11 CONTINUE
   10 CONTINUE
      RETURN
      END

c************************************************
C BESSJJ
c
      SUBROUTINE BESSJJ(A,NUM,BESJ,BESY)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1),BESY(NUM+1) 
      BESJ(NUM+1)=0D0
      BESJ(NUM)=1D-260
      N=2*NUM+1
      NUM1=NUM-1
      DO 11 I=1,NUM1
      N=N-2
      I1=NUM-I
   11 BESJ(I1)=N*A*BESJ(I1+1)-BESJ(I1+2)
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
      BESY(I)=0D0
   13 BESJ(I)=C*BESJ(I)
      RETURN
      END

c************************************************
C cdcof3
C
      SUBROUTINE cdcof3(VL,AD,IA,BD,IB,ID,VN,IN,M,N,C,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,AD,VN,C
      DIMENSION AD(IA),BD(IB),RL(550)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      common /FACT/ FACT(170)
      AC=C
      IER=0
C* èPOBEPKA BBOÑàMõX áHAóEHàâ èAPAMETPOB
      IF(M.LT.0.OR.N.LT.M) IER=1
      IF(AC.LT.0D0) IER=3
      IF(ID.EQ.1.AND.IA.LE.((N-M)/2+2)) IER = 4
      IF(ID.EQ.0.AND.IN.NE.0) IER = 6
      IF(IER.NE.0) GO TO 900
      M2=M*M
      MM=2*M
      S2=(1D0,0D0)
      AM = 4D0*M2-1D0
      RC2=C*C
      RC4=RC2*RC2
      N1=N-M
      IF(ID.EQ.0) GO TO 900
C**** BõóàCãEHàE K-TOB PAáãOÜEHàü
      IF(MOD(N1+1,2)) 21,20,21
   20 I11=1
      GO TO 22
   21 I11=0
   22 CONTINUE
      L1=(N1-I11)/2
      L11=L1+1
      RL(L11)=S2
      IF(L1.EQ.0) GO TO 24
C***  BõóàCãEHàE K-TOB èPà R<N-M
      DO 23 J=1,L1
      JJ1=N1-2*J
      RU=S2
      RB=(M+JJ1)*(M+JJ1+1D0)+AK*RC2/2D0*(1D0-AM/(MM+2D0*JJ1-1D0)/
     *   (MM+2D0*JJ1+3D0))        -VL
      RV=S2/RB
      RW=RV
      JJ2=JJ1
      IF(JJ2.LE.1) GO TO 103
  102 CONTINUE
      II=MM+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(M+JJ2-2D0)*(M+JJ2-1D0)+AK*RC2/2D0*(1D0-AM/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
*      IF(CDABS(RV).LE.1d-150) GO TO 103
c      IF(CDABS(RV).LE.1d-250) GO TO 103
      IF(CDABS(RV).LE.eps3) GO TO 103
      JJ2=JJ2-2
      IF(JJ2.LE.1) GO TO 103
      GO TO 102
  103 RL(L11-J)=-(MM+JJ1+1D0)*(MM+JJ1+2D0)*AK*RC2/(MM+2D0*JJ1+3D0)/
     *           (MM+2D0*JJ1+5D0)*RW
   23 RL(L11-J)=RL(L11-J)*RL(L11+1-J)
   24 CONTINUE
C***  BõóàCãEHàE K-TOB èPà R>N-M
      IA1=IA-1
      DO 25 J=L11,IA1
      JJ1=N1+2*(J-L11)+2
      II=MM+JJ1
      II1=II+JJ1
      RU=S2
      RB=(M+JJ1)*(M+JJ1+1D0)+AK*RC2/2D0*(1D0-AM/(II1-1D0)/(II1+3D0))-VL
      RV=-JJ1*(JJ1-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RW=RV
      JJ2=JJ1
  104 CONTINUE
      JJ2=JJ2+2
      II=MM+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(M+JJ2)*(M+JJ2+1D0)+AK*RC2/2D0*(1D0-AM/(II1-1D0)/(II1+3D0))
     *          -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
       IF(CDABS(RV).LE.eps3) GO TO 105
c       IF(CDABS(RV).LE.1d-250) GO TO 105
      GO TO 104
  105 RL(J+1)=RW*(MM+2D0*JJ1-1D0)*(MM+2D0*JJ1+1D0)/((MM+JJ1-1D0)*
     *      (MM+JJ1)*AK*RC2)
      RL(J+1)=RL(J)*RL(J+1)
       IF(J.GT.L11) RL(J+1)=RL(J+1)/RL(L11)
       IF(CDABS(RL(J+1)).LT.1D-160) GO TO  125
c       IF(CDABS(RL(J+1)).LT.eps3) GO TO  125
   25 CONTINUE
      GO TO 126
  125 IJ=J+2
         DO 127 I=IJ,IA
  127    RL(I)=S
  126 CONTINUE
C *                  BõóàCãEHàE K-TA D(N-M)
      L4=(N+M+I11)/2
      F2=1D0/Fact(l1+1)
      DO 27 I=1,L4
   27 F2=F2*(L4+I)
      D1=(-1D0)**L1*F2/2D0**N1
      SNM=S
      DO 28 J=1,IA
      J1=2*(J-1)+I11
      J2=(J1-I11)/2
      J3=(J1+MM+I11)/2
      rF2=RL(J)/2D0**J1
      DO 30 I=1,J3
   30 rF2=rF2*(J3+I)
      rF1=rF2
      DO 29 I=1,J2
   29 rF1=rF1/I
       SG=(-1D0)**J2*rF1
       IF(CDABS(SG).LT.1d-150.AND.J.GT.(N-M)) GO TO 128
c       IF(CDABS(SG).LT.eps3.AND.J.GT.ia/2) GO TO 128
   28 SNM=SNM+SG
  128 AD(L11)=D1/SNM
      DO 31 I=1,IA
   31 AD(I)=RL(I)*AD(L11)
      IF(ID.NE.2) GO TO 901
C*       BõóàCãEHàE KOùîîàñàEHTOB PAáãOÜEHàü èPà  R<0
      BD(1)=AD(1)
      DO 223 J=2,IB
      JJ1=I11-2*(J-1)
      RU=S2
      RB=(M+JJ1)*(M+JJ1+1D0)+AK*RC2/2D0*(1D0-AM/(MM+2D0*JJ1-1D0)/
     *   (MM+2D0*JJ1+3D0))        -VL
      RV=1D0/RB
      RW=RV
      JJ2=JJ1
  224 CONTINUE
      II=MM+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(M+JJ2-2D0)*(M+JJ2-1D0)+AK*RC2/2D0*(1D0-AM/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
      IF(CDABS(RV).LE.eps3) GO TO 225
c      IF(CDABS(RV).LE.1d-250) GO TO 225
      JJ2=JJ2-2
      GO TO 224
  225 CONTINUE
      IF(JJ1.EQ.-MM-2.OR.JJ1.EQ.-MM-1) GO TO 226
      BD(J)=-(MM+JJ1+1D0)*(MM+JJ1+2D0)*AK*RC2/(MM+2D0*JJ1+3D0)/
     *           (MM+2D0*JJ1+5D0)*RW
      GO TO 227
  226 BD(J)= AK*RC2/(MM+2D0*JJ1+3D0)/(MM+2D0*JJ1+5D0)*RW
      IF(I11.EQ.1) BD(J)=-BD(J)
  227 BD(J)=BD(J)*BD(J-1)
          IF(DABS(BD(J)).LT.eps3) GO TO 1223
c          IF(DABS(BD(J)).LT.1D-260) GO TO 1223
  223 CONTINUE
        GO TO 901
 1223 IJ=J+1
        DO 1224 I=IJ,IB
 1224   BD(I)=0D0
  901 IF(IN.EQ.0) GO TO 900
C*****  BõóàCãEHàE HOPMàPYûôEÉO MHOÜàTEãü *
      VN=S
      DO 32 J=1,IA
      JJ1=2*(J-1)+I11
      rF1=AD(J)**2
      IF(M.EQ.0) GO TO 333
      DO 33 I=1,MM
   33 rF1=rF1*(JJ1+I)
  333 RVN=rF1/(MM+2D0*JJ1+1D0)
c      IF(CDABS(RVN).LT.1D-250.AND.J.GT.ia/2) GO TO 34
      IF(CDABS(RVN).LT.eps3.AND.J.GT.ia/2) GO TO 34
   32 VN=VN+RVN
   34 VN=2D0*VN
  900 RETURN
      END

c************************************************
c cdrb12
c
      SUBROUTINE CDRB12(R,DR,AR,ADR,ID,M,N,C,X,AD,IA,IER)
      parameter (nterm=200)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,S1,BESJ,BESY,AC,SIG,RK,R,DR,F2,RS,RD,S,BR,ABR
      DIMENSION AD(IA),BESJ(4*nterm),BESY(4*nterm)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /pi/ pi
      IER=0
      IF(X.LT.1D-10) IER=10
      IF(IER.NE.0) GO TO 900
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 5
    2 JK=0
    5 CONTINUE
      AM=M
      AM=AM/2D0
      PI2 = pi / 2d0
      X2=X*x
      AC=C*X
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      F2=AD(I)
      DO 8 I1=1,M2
    8 F2=F2*(J+I1)
    7 SIG=SIG+F2
      RK=CDSQRT(PI2/AC)/SIG*((X2-AK)/X2)**AM
      NUM=3*CDABS(AC)+10
      if(NUM.lt.60) num=60
      if(NUM-1.ge.4*nterm) write(*,*) num, 4*nterm
      if(NUM-1.ge.4*nterm) stop
      IF(ID.EQ.1) CALL CESSEL(1D0/AC,NUM,BESJ,BESY)
      IF(ID.EQ.0) CALL BESSJJ(1D0/AC,NUM,BESJ,BESY)
      R=S
      DR=S
      AR=0D0
      ADR=0D0
      F1=AK/X*M/(X2-AK)
      NM2=MIN0(IA,NUM/2-2)
      DO 13 I=1,NM2
      J=2*(I-1)+JK
      F2=AD(I)*S1**(J+M-N)
      DO 18 I1=1,M2
   18 F2=F2*(J+I1)
      IF(I.EQ.1) GO TO 21
c      IF(CDABS(RD).LT.1D-160.AND.CDABS(RS).LT.1D-160.AND.I.GT.ia/2
      IF(CDABS(RD).LT.eps3.AND.CDABS(RS).LT.eps3.AND.I.GT.(M+5)
     *   .AND.ID.EQ.0)   GO TO 19
   21 RS=F2*BESJ(J+M+1)

*          write(*,*) i,ad(i)
*          write(*,*) f2
*          write(*,*) j+m+1,besj(j+m+1)
*          write(*,*) rs
*          write(*,*) r

c          write(7,*) i,ad(i)
c          write(7,*) f2
c          write(7,*) j+m+1,besj(j+m+1)
c          write(7,*) rs
c          write(7,*) r


      RD=F2*((F1+(M+J)/X)*BESJ(J+M+1)-C*BESJ(J+M+2))
      R=R+RS
      DR=DR+RD
   20 CONTINUE
      IF(ID.EQ.0) GO TO 13
      IF(I.EQ.1) GO TO 22
c      IF(CDABS(BR).LT.1D-160.AND.CDABS(ABR).LT.1D-160.AND.I.GT.(M+5))
      IF(CDABS(BR).LT.eps3.AND.CDABS(ABR).LT.eps3.AND.I.GT.(M+5))
     *   GO TO 14
   22 BR=F2*BESY(J+M+1)
                     AR  = AR  + BR
      ABR=F2*((F1+(M+J)/X)*BESY(J+M+1)-C*BESY(J+M+2))
                     ADR = ADR + ABR
   13 CONTINUE
   14 AR= AR* RK
      ADR= ADR* RK

c          write(7,*) r,1d0/rk
*          write(*,*) r,1d0/rk

   19 R=R*RK
      DR=DR*RK
  900 RETURN
      END

c************************************************
c CDRF12
c
      SUBROUTINE CDRF12(R,DR,AR,ADR,I12,M,N,C,AD,IA,BD,IB,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 AD,C,DR
      DIMENSION AD(IA ),BD(IB+1),P(250),PD(250),Q(250),QD(250),FACT(170)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P, PD, Q, QD
      COMMON /FACT/FACT
      ier=0
      M2=2*M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE
      IF(MOD(N,2)) 3,4,3
    3  CONTINUE
      Jn=1
      GO TO 5
    4  CONTINUE
      Jn=0
    5 CONTINUE
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      if(m2+j.gt.169) go to 291
      rf2=AD(I)*fact(m2+j+1)/fact(j+1)
      go to 292
 291  CONTINUE
      rf2=ad(i)
      do 293 jj=1,m2
      rf2=rf2*(jj+j)
 293  CONTINUE
 292  CONTINUE
c      if(cdabs(rf2).lt.1d-160) go to 8
      if(cdabs(rf2).lt.eps3) go to 8
    7 SIG=SIG+rf2
    8  CONTINUE
      RK1=(m2+1d0+jk*2d0)/2d0**nm2*fact(nm2+1+jk)/
     *    AD(1)/fact(m+1)*sig/fact((nm1-jk)/2+1)/
     *    fact((nm2+jk)/2+1)/C**(M+JK)
      if(k.eq.1) rk1=rk1/(-s1)**(M+JK)*(-s1)**jn
      R=S
      DR=S
      IF(I12.EQ.0) GO TO 70
      RK2=2d0**nm1/(m2-1d0)*fact(m2+1)/
     *    fact(m+1)*fact((nm1-jk)/2+1)*fact((nm2+jk)/2+1)/
     *    fact(nm2+jk+1)*bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(m2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(m2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(m2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      AR=0D0
      ADR=0D0
      DO 12 I=1,M
        J2=2*(M-I)+1+JK
      ADR=ADR+BD(I+1)*QD(J2)
   12 AR=AR+BD(I+1)*Q(J2)           
   70 CONTINUE
      DO 13 I=1,IA
        J1=2*(I-1)+1+JK
        J2=2*(M+I)-1+JK
      RS=AD(I)*P(J1)
      RD=AD(I)*PD(J1)
      R=R+RS
      DR=DR+RD
      IF(I12.EQ.0) GO TO 17
c      IF(DABS(BR).LT.1D-160.AND.DABS(ABR).LT.1D-160.AND.I.GT.(M+5))
      IF(DABS(BR).LT.eps3.AND.DABS(ABR).LT.eps3.AND.I.GT.(M+5))
     *   GO TO 17
      BR=AD(I)*Q(J2)
      ABR=AD(I)*QD(J2)
      AR=AR+BR
      ADR=ADR+ABR                                  
   17 CONTINUE
c      IF(CDABS(RS).LT.1D-160.AND.CDABS(RD).LT.1D-160) GO TO 14
      IF(CDABS(RS).LT.eps3.AND.CDABS(RD).LT.eps3) GO TO 14
   13 CONTINUE
   14 CONTINUE
      IF(I12.EQ.0) GO TO 71
      MR=M+1
      DO 15 I=MR,IB
        J2=2*(I-M)-JK
      if (k.eq.1.and.mod(m,2).eq.0) go to 116
      BR=BD(I+1)*P(J2)
      ABR=BD(I+1)*PD(J2)
      go to 115
  116 CONTINUE
      BR=-BD(I+1)*P(J2)
      ABR=-BD(I+1)*PD(J2)
  115 CONTINUE
      AR=AR+BR
      ADR=ADR+ABR
ccc                     write(7,*) i,ar,br,p(j2),bd(i+1)
c      IF(DABS(BR).LT.1D-160.AND.DABS(ABR).LT.1D-160.AND.I.GT.(M+5))
      IF(DABS(BR).LT.eps3.AND.DABS(ABR).LT.eps3.AND.I.GT.(M+5))
     *   GO TO 16
   15 CONTINUE
   16 AR=AR/RK2
      ADR=ADR/RK2
   71 R=R/RK1
      DR=DR/RK1
  900 RETURN
      END

C**********************************************************************
C CDLAMn  - ... eigenvalues
C
      SUBROUTINE CDLAMn(vl0,VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,C,dx,vl0(nn)
      DIMENSION VL(NN),RL(550),SL(550),RDA(550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
      if(ac.lt.32+m-3.or.ac.gt.47.4d0) SC=-sc
      IER=0
      in0 = 0
C*  Parameters check-up
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s
      M2=M*M
      MM=2*M
      RC2=C*C
             DO 111 IN=1,NN
             dx=dcmplx(0.1d0,0.1d0*af)
             n1001=0
c             if(dreal(sc).lt.0d0) sc = -sc
 1113 continue
             JN=0
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
      IF(NMAX.GT.550) WRITE(7,2000) NMAX
 2000 FORMAT(14X,'******* CDLAMn ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE

          if(cdabs(vl0(in)).lt.1d-10) then
                 in0 = in
                 vl0(in)=VL0(IN-1)+dabs(dreal(VL0(IN-1)-VL0(IN-2)))

      IF(CC.le.(N+3D0)) then
      N2=N+M
      N3=2*N-1
      N4=N3+4
      AM=4D0*M2-1D0
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      vl0(in) = N*(N+1D0)-AK*RC2/2D0*(AM/N3/N4-1D0)+RC2**2/2D0*
     *          ((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *          (N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))
      end if

                    print *,in,vl0(in)

                  end if

                 r0=vl0(in)
                 rl(1)=r0
                 go to 2

  400 CONTINUE
         JN=JN+1
         RL(1)=R0-JN*dx*SC*AK

             IF(JN.EQ.1001.and.in.gt.2) GO TO 1112
             IF(JN.EQ.1001.and.in.le.2) GO TO 1131

    2 CONTINUE
               RDEL=dx
      NUM=1
  100 CONTINUE
c -------------------------------------
C*            Iteration: 'from below to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1 = 4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2D0*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*            Iteration: 'from top to below'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-  2
      DO 12 I=  2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB
      NUM=NUM+1

      if(num.gt.540) eps1 = eps1 * 10d0
*      if(num.gt.540) eps1 = eps1 * 10d0 ** (num - 540)
*      if(num.gt.545) eps1 = 1d-08
      if(num.gt.540)  print *, num, eps1

            if(num.eq.550) then
            ier = 8
            go to 900
            end if

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE
      IF(NUM.gE.1009) GO TO 15
      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)

cc                    print *,in,vl0(in)-vl(in)

      IF(aC.LT.2.5D0) GO TO 110
        Cc2=dreal(r0)
        Cc1=dreal(VL(IN))
        CFF=dabs(cc1-cc2)
        cn=ac
        an=dfloat(n)
        if(ac.gt.13d0) an = an + 1
        if(an.lt.cn) cn = an
        if(an.gt.cn.and.ac.lt.10d0) cn = an

      IF(in.le.2) GO TO 1121
      IF(in0.ne.0) GO TO 110

         cbb=dreal(VL(IN))-dreal(VL(IN-1))
                   IF(k.eq.0.or.k.eq.1.and.
     *                MOD(in,2).eq.0.and.dimag(c).gt.0.5d0) then
        IF(cbb.lt.0d0)  go to 110
                    end if
        IF(CFF.gt.cn.or.cbb.le.0d0) go to 400
 1112  CONTINUE
 3001  FORMAT(I2,3D15.8,1x,2d10.3)
 3002  FORMAT(1x,'*',2f6.2,2I5,2x,3(f10.3,f10.3),'*')
        ch=r0
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        ck=vl(in)
       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,crn,cc0,ck
       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,crn,cc0,ck
            if(jn.lt.1001) go to 110
            sc=-sc
            n1001=n1001+1
            eps1 = eps1 * 10d0

cccc            if(n1001.ge.8) then
            if(n1001.ge.4) then
            ier = 8
            go to 900
            end if

c            if(n1001.ge.28) go to 110
            if(n1001.le.3) go to 119
             jn=0
             dx=dx/2d0
c            if(mod(n1001,2).eq.0) dx=dx/5d0
            go to 400
  119  CONTINUE
            if(mod(n1001,2).eq.0) dx=dx/1.5d0
            go to 1113
 1121  CONTINUE
        IF(CFF.gt.ac) go to 400
 1131  CONTINUE
        ch=r0
        cr=rl(1)
        ck=vl(in)
       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,ck
       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,ck
            if(jn.lt.1001) go to 110
            sc=-sc
            n1001=n1001+1
            eps1 = eps1 * 10d0

ccc*            if(n1001.ge.8) then
            if(n1001.ge.4) then
            ier = 8
            go to 900
            end if

c            if(n1001.ge.28) go to 110
            if(n1001.le.3) go to 119
             jn=0
             dx=dx/2d0
c            if(mod(n1001,2).eq.0) dx=dx/5d0
            go to 400
  110  CONTINUE

  111  CONTINUE
  900 CONTINUE
      return
      END

C**********************************************************************
C CDLAMo_l - Oblate eigenvalues (large values)
C**********************************************************************
C
      SUBROUTINE CDLAMo_l(VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
c      COMPLEX*16 VL,C,dx,alam0
      COMPLEX*16 VL,C,dx
      DIMENSION VL(NN),RL(550),SL(550),RDA(550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT(170)

*****************************************
c      COMMON /r0/ alam0(500)
*****************************************
      inull = 0
 991  continue
      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
      IER=0
c -------------------------------------
C**  Parameters check-up
c -------------------------------------
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s
      DO 113 IN = 1,550
          rl(in) = s
          sl(in) = s
 113      rda(in) = s

      sdelta=s

      M2=M*M
      MM=2*M
      AM=4D0*M2-1D0
      RC2=C*C
      RC4=RC2*RC2
cc      NC=idint(aC)
             jcn=0
             DO 111 IN=1,NN

              inout = 4
c              if(ac.gt.75d0)  inout = 8
              nt = 0
              iv = 0

*mmm             if(ac.gt.27d0.and.m.gt.10) sc=-sc
c             dx=dcmplx(0.1d0,0.1d0*af)
             dx=dcmplx(0.25d0,0.25d0*af)
c             if(ac.gt.70d0)  dx = dcmplx(0.3d0,0.3d0*af)
c             if(ac.gt.25d0) dx=dcmplx(0.25d0,0.25d0*af)
             if(ac*m.gt.400d0)
     *          dx=dcmplx(0.35d0+(ac*m-400d0)/1000d0,0.35d0*af)

             if(in.gt.6.and.dreal(vl(in-1)).lt.-4*ac) then
             dx=dcmplx(0.1d0,0.1d0*af)
             if(ac.gt.60d0) dx=dcmplx(0.25d0,0.25d0*af)
             if(ac.gt.70d0) dx=dcmplx(0.3d0,0.3d0*af)
             end if

             if(in.gt.6.and.ac.gt.40d0.and.dreal(vl(in-2)).lt.0d0.and.
     *          dreal(vl(in-1)).gt.-4*ac)
     *          dx=dcmplx(0.35d0+(ac-40d0)/200d0,0.35d0*af)

             if(in.gt.6.and.ac.gt.70d0.and.
     *          dabs(dreal(vl(in-2)) + dreal(vl(in-1))).lt.cn)
     *          dx=dcmplx(0.35d0,0.35d0*af)

             n1001=0
             if(dreal(sc).lt.0d0.and.n1001.eq.0) sc = -sc
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
 1113 continue
             JN=0
      IF(NMAX.GT.500) WRITE(7,2000) NMAX
 2000 FORMAT(14X,'******* CDLAMo_l ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE
      IF(aC.LT.5.50D0) GO TO 401

      IF(n.lT.7.or.in.le.3) GO TO 401
      IF(in.lT.7) GO TO 401
      IF(dreal(vl(in-1)).lt.0d0.and.m.ge.10) GO TO 401

      IF(dreal(vl(in-1)).lt.3d0*ac.and.m.ge.15) then
      if(dreal(sc).lt.0d0.and.n1001.eq.0) sc = - sc
      GO TO 401
      end if

        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an

                   IF(dabs(dreal(vl(in-1))).lT.cn.
     *                and.dreal(vl(in-2)).lT.cn.
     *                or.dabs(dreal(vl(in-2)+vl(in-1))).lT.ac+m) then
           IF(dreal(vl(in-1)).lT.0d0)  RL(1) = S
           IF(dreal(vl(in-1)).gT.0d0)  RL(1) = 2d0*vl(in-1)
           IF(dreal(vl(in-1)).gT.cn)   RL(1) = vl(in-1) + ac/10d0
             r0=rl(1)
c             if(mod(n-m+1,2).ne.0) nt = 1
             if(n1001.gt.1.and.m.lt.5) dx = dx / 1.25d0
             if(n1001.ne.0.and.m.ge.5)  dx = 2d0 * dx
ccccc             if(n1001.eq.0.and.dreal(dx).lt.ac/100d0)  dx = ac / 100d0
ccccc             if(dreal(sc).gt.0d0.and.ac.lt.60d0.and.m.lt.5) sc = - sc

             if(n1001.eq.0) then
             if(dreal(dx).lt.ac/100d0)  dx = ac / 100d0
             if(dreal(sc).gt.0d0.and.ac.lt.60d0.and.m.lt.5) sc = - sc
             end if

c             write(*,*) in,r0

             print *, n, '   0'
c             write(7,*) n, '   0'
             iv = 100

             GO TO 2
                    end if

         IF(dreal(vl(in-1)).lt.cn.
     *      and.dabs(dreal(vl(in-2)+vl(in-1))).lT.cn.
     *      or.dreal(vl(in-1)).lt.cn.
     *      and.dabs(dreal(vl(in-2)+vl(in-3))).lT.cn) then
            nt = 2
            GO TO 401
         end if


c~~~~~~~~~~~~~~~~~
         IF(dreal(vl(in-1)).gt.cn.and.dreal(vl(in-1)+vl(in-2)).gt.ac+m.
     *      or.dreal(vl(in-2)).gt.0d0) then

c         IF(dreal(dx).gt.0.1d0.and.dabs(VL(IN-1)-VL(IN-2)).gt.cn+2*m+5)

                 IF(dreal(dx).gt.0.1d0.and.
c     *           cdabs(dreal(VL(IN-1)-VL(IN-2))).gt.ac+2*m+5)
     *           cdabs(VL(IN-1)-VL(IN-2)).gt.ac+2*m+5)
     *           dx = dcmplx(0.1d0,0.1d0*af)


c===============
       if(cdabs(VL(IN-2)-VL(IN-3)).lt.cn.or.
     *    cdabs(VL(IN-1)-VL(IN-2)).lt.cn) then

       if(mod(n-m+1,2).eq.0) r0 = VL(IN-1) + cdabs(VL(IN-2)-VL(IN-3))
       if(mod(n-m+1,2).ne.0) r0 = VL(IN-2) + 2d0*ac
             print *, n, '   ++++'
c            write (7,*) n, '   ++++'
             iv = 1111
c----------------
             if(n1001.eq.0.or.n1001.eq.1) then
             dx = dcmplx(0.25d0+ac*m/10000d0,0.25d0*af)
             if(ac.gt.60d0) dx = dcmplx(0.3d0+ac*m/10000d0,0.3d0*af)
             if(ac.gt.75d0) dx = dcmplx(0.4d0+ac*m/10000d0,0.4d0*af)
             end if
c----------------

       go to 223
       end if
c===============

       r0=VL(IN-1) + cdabs(VL(IN-1)-VL(IN-2))

c        write(*,*) in,r0

c^^^^^^^^^^^^^^^
      if(cdabs(VL(IN-1)-VL(IN-2)).lt.10d0-m/15d0) then

c----------------
         if(ac.ge.40d0) then
         if(ac*m.ge.600d0.or.dreal(vl(in-1)).gt.cn+2*m)
     *      r0 = r0 + VL(IN-2)/2d0
         if(ac*m.lt.600d0.or.dreal(vl(in-1)).lt.cn+2*m)
     *      r0 = r0 + VL(IN-2)
         end if
c----------------

c----------------
          if(ac.lt.40d0.and.m.le.15) then
          if(dreal(vl(in-1)).gt.ac+m) r0 = r0 + VL(IN-2)/2d0
          if(dreal(vl(in-1)).le.ac+m) r0 = r0 + VL(IN-2)
          end if
c----------------

      end if
c^^^^^^^^^^^^^^^

 223    continue
      if(ac.ge.30d0.and.m.gt.3) r0=r0+m/5d0
      if(ac.ge.16d0.and.ac.lt.29d0) r0=r0+2.086d0
      if(ac.ge.30d0.and.ac.lt.40d0) r0=r0+(ac-30d0)/2d0
      if(ac.ge.30d0.and.ac.lt.40d0.and.dreal(vl(in-2)).lt.cn) r0=r0-m
      if(ac.ge.40d0.and.ac.lt.50d0)  r0=r0+n/10d0-(ac-40d0)/3d0-m/2.5d0
c      if(ac.ge.40d0.and.ac.lt.50d0)  r0=r0+n/10d0-(ac-40d0)/3d0
      if(ac.ge.50d0.and.ac.le.60d0) r0=r0+(ac-50d0)/3d0
      if(ac.gt.60d0.and.ac.le.70d0) r0=r0+(ac-60d0)/3d0
      if(ac.gt.70d0.and.ac.le.80d0) r0=r0+(ac-70d0)/4d0
      if(ac.gt.80d0.and.ac.le.90d0) r0=r0-(ac-80d0)/8d0+3.5d0
      if(ac.gt.90d0.and.ac.le.100d0) r0=r0-(ac-90d0)/8d0+3.5d0
      if(ac.gt.100d0) r0=r0-(ac-100d0)/8d0+3.5d0
      if(ac.gt.110d0) r0=r0+(ac-110d0)/8d0+0.5d0

      if(ac.gt.50d0.and.m.gt.15)  r0=r0-n/40d0
             jcn=jcn+1

c        write(*,*) in,r0

c             print *, n, '   +'
c            write(7,*) n, '   +'
             iv = 11
      GO TO 2
          end if
c~~~~~~~~~~~~~~~~~

                    nt = 0
cccc                    nt = 3
                    GO TO 401

*             write(*,*) in,r0

c                   GO TO 2
  400 CONTINUE

cccccc         sdelta=s

             JN=JN+1
 3000        FORMAT(1X,'CDLAMo_l N,RL(1)=',I5,5D20.10)
             RL(1)=R0-JN*dx*SC
             IF(JN.EQ.1001) GO TO 1112
c            GO TO 2
             GO TO 22
  401 CONTINUE

c             print *, n, '   401'
c             write(7,*) n, '   401'

             iv = 401

      RL(1)=S
      r0=rl(1)
      IF(CC.LE.0.5D0) GO TO 2
      N3=2*N-1
      N4=N3+4
      IF(CC.GT.(N+3D0)) GO TO 1
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      N2=N+M
      RL(1)=N*(N+1D0)-AK*RC2/2D0*(AM/N3/N4-1D0)+RC4/2D0*
     *((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *(N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))
      r0=rl(1)
c                                 WRITE(7,3000) N,RL(1)
c                                 WRITE(*,3000) N,RL(1)
      GO TO 2
    1 CONTINUE
c -------------------------------------
C** Initial approximation (large 'c')
c -------------------------------------
      IF(MOD(N1,2)) 3,4,3
    3 q=N
      as = 1d0
      GO TO 16
    4 q=N+1d0
      as = - 1d0
   16 continue
      q2=q*q
      q4=q2*q2
      q6=q4*q2
      d=q2+1d0-M2
      RL(1)=-RC2+2D0*C*q-d/2D0-q*d/8D0/C
     *      -(4D0*q2*(d+1D0)+d**2)/64D0/RC2
     *      -q*(33d0*q4+114d0*q2+37d0-2d0*m2*(23d0*q2+25d0)
     *      +13*m2**2)/512d0/rc2/c
     *      -(63d0*q6+340d0*q4+239d0*q2+14d0
     *      -10d0*m2*(10d0*q4+23d0*q2+3d0)+m2*m2*(39d0*q2-18d0)
     *      -2d0*m2**3)/1024d0/rc4
     *      -q*(527d0*q6+4139d0*q4+5221d0*q2+1009d0
     *      -m2*(939d0*q4+3750d0*q2+1591d0)+m2*m2*(465d0*q2+635d0)
     *      -53d0*m2**3)/8192d0/rc4/c
c     *      +as*(4d0*c)**(q+1d0)/cdexp(2d0*c)/fact(n+1)/fact(n+m+1)
c     *      *(1d0-(3d0*q2+4d0*q+1d0-m2)/8d0/c
c     *      +((9d0*q4+4d0*q2*q-18d0*q2-12d0*q-7d0)
c     *      -m2*(6d0*q2-4d0*q-6d0)+m2*m2)/128d0/rc2
c     *      -((23d0*q6-44d0*q4*q-171d0*q4+400d0*q2*q+257d0*q2+444d0*q
c     *      +51d0)-m2*(19d0*q4-56d0*q2*q-50d0*q2+240d0*q+as*55d0)
c     *      +m2*m2*(5d0*q2-12d0*q+5d0)-m2*m2*m2)/3d0/1024d0/rc2/c)

               IF(dreal(vl(in-1)).lT.0d0.and.dreal(vl(in-2)).lT.0d0.
     *             and.dreal(rl(1)).gT.0d0) then

               if(inull.eq.0.or.ac.lt.83d0.or.dreal(rl(1)).ge.ac/2d0)
     *            then
                 rl(1) = s
                 go to 121
                 end if

               if(inull.ne.0) then
               if(ac.ge.83d0.and.dreal(rl(1)).lT.ac/2d0)
     *            rl(1) = vl(in-1) / 2d0
               sdelta = s
               end if

  121          CONTINUE
               end if

      if(ac.gt.24.92d0.and.ac.lt.24.93d0.and.m.eq.1.and.n.eq.15) then
      rl(1)=-2.78d0
      inout = 32
      end if

      if(ac.gt.55.32d0.and.ac.lt.55.33d0.and.m.eq.1.and.n.eq.34) then
      rl(1)=-52.247d0
      inout = 32
      sdelta = s
      end if

      if(ac.gt.63.69d0.and.ac.lt.63.70d0.and.m.eq.1.and.n.eq.37) then
      rl(1)=-144.2d0
      inout = 32
      sdelta = s
      end if

      r0=rl(1)
    2 CONTINUE


*****************************************

cc      IF(JN.EQ.0.and.n1001.eq.0) then
      IF(JN.EQ.0) then


      delta=sdelta
c      if(dreal(r0).le.0d0.or.nt.ne.0.or.dreal(vl(in-1)).gt.5d0*ac)

c         write(*,*) nt, r0

c       if(dreal(r0).le.0d0.or.nt.ne.0)
       if(dreal(r0).lt.0d0.or.nt.ne.0)
     *                   r0 = r0 + delta

c         write(*,*) in, r0

         if(m.gt.6.and.dabs(dreal(r0)).lt.cn.
     *      and.dreal(dx).lt.0.7d0) dx = 0.7d0

         if(m.gt.6.and.n1001.eq.0)  then
         if(dabs(dreal(r0)).lt.cn.or.dabs(dreal(vl(in-1))).lt.cn) then
c     *                   dx = 0.5d0
                        dx = 0.7d0
             if(m.gt.11.and.ac*m.gt.400d0.and.mod(n-m+1,2).ne.0) then
c             if(m.gt.11.and.ac*m.gt.600d0) then
             if(dreal(sc).gt.0d0) sc = - sc
             dx = 2d0*(ac+m-11)/10d0
             end if
             end if
         end if

      if(in.gt.1.and.dreal(r0).lt.dreal(vl(in-1)))  r0 = vl(in-1)
      rl(1) = r0

c         write(*,*) r0

      end if

c      alam0(in)=r0

c         write(7,*) in, r0


*****************************************
   22 CONTINUE
c      RDEL=dx
c      if(dreal(dx).gt.0.2d0) RDEL = dcmplx(0.2d0,0.2d0*af)
        RDEL0 = dcmplx(0.1d0,0.1d0*af)
        RDEL = rdel0
  222 CONTINUE
      NUM=1

  100 CONTINUE
c -------------------------------------
C*       Iteration: `from below to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1=4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*       Iteration: `from top to below'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-2
      DO 12 I=2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB

c      if(in.eq.1)   write(7,*) 'sl  ', num, sl(num)

      NUM=NUM+1

      if(num.gt.540) eps1 = 1d-10
*      if(num.gt.500) eps1 = 1d-08
      if(num.gt.540)  print *, num

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE

*         print *,num,num-2

      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))


c      if(in.eq.1)   print *, num, rl(num)
c      if(in.eq.1)   write(7,*) num, rl(num)
*      if(in.eq.1)   print *, num, rdel


      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL

c      if(in.eq.1)   print *, num, rl(num)
c      if(in.eq.1)   write(7,*) num, rl(num)
c      if(in.eq.1)   write(7,*) rdel
c      if(in.eq.1)   print *,  rdel

      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)
      IF(aC.LT.2.5D0) GO TO 110
      IF(in.le.2) GO TO 110

        Cc2=dreal(rl(1))

***        if(m.ge.10.and.dreal(vl(in-1)).gt.0d0)  Cc2=dreal(r0)

        Cc1=dreal(VL(IN))
        CFF=dabs(cc1-cc2)
c                if(cff.gt.10d0*ac.and.cc1.lt.10d0*ac) then
                if(cff.gt.20d0*ac.and.cc1.lt.-ac) then
                rdel = rdel/2d0
                rl(1) = r0
                write(*,*) in, rdel
                write(*,*) cff, cc1
                go to 222
                end if
        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an
                   IF(n1001.eq.2.
     *                and.dreal(vl(in)).gT.0d0.
     *                and.dreal(vl(in-1)).lT.0d0.
     *                and.dreal(vl(in-2)).lT.0d0) then
                cn=ac+an
                end if
         cbb=dreal(VL(IN))-dreal(VL(IN-1))

        IF(CFF.gt.cn) go to 400
        IF(cbb.le.0d0.and.dreal(vl(in)).gT.-5*ac) go to 400
        IF(cbb.lt.0d0.and.dreal(vl(in)).le.-5*ac) go to 400

 1112  CONTINUE

 3001  FORMAT(I2,3D15.8,1x,2d10.3)
 3002  FORMAT(1x,'L',2f7.3,2I5,2x,3(f10.3,f10.3))

        ch=r0
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        ck=vl(in)
       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,crn,cc0,ck
       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,crn,cc0,ck
            if(jn.lt.1001) go to 110
            sc = -sc
            n1001=n1001+1

c            if(n1001.ge.4) then
            if(n1001.ge.inout) then
       write(7,*) ac
       write(*,*) ac
       write(7,*) '**   Problems with computations of eigenvalues!'
       write(*,*) '**   Problems with computations of eigenvalues!'
            ier = 8
            go to 900
            end if

  109  CONTINUE
ccc            if(mod(n1001,2).eq.0) dx=dx/1.5d0
            if(mod(n1001,2).eq.0.and.nt.eq.0) dx=dx/3d0
            go to 1113
  110  CONTINUE
c        sdelta=vl(in)-r0
       if(dreal(vl(in)).gt.0d0) sdelta=vl(in)-r0
       if(dreal(vl(in)).lt.0d0.and.mod(in,2).ne.0) sdelta=vl(in)-r0

                if(cdabs(sdelta).gt.5d0*ac.and.in.lt.3) then
c                rdel = rdel/200d0
                rdel = rdel/2d0
                eps1 =  eps1 * 10d0
                write(*,*) in, rdel
                write(*,*) r0, vl(in)
                go to 222
                end if


         IF(dreal(vl(in)-vl(in-1)).gt.in.or.
c         IF(dreal(vl(in)-vl(in-1)).gt.in.and.
c     *      dreal(vl(in-1)-vl(in-2)).gt.in) go to 111
     *      dreal(vl(in-1)-vl(in-2)).gt.in-1-ac/10d0) go to 111

          IF(CC.GT.N+3D0.and.iv.ne.401) then

      IF(MOD(N1,2)) 83,84,83
   83 q=N
      as = 1d0
      GO TO 816
   84 q=N+1d0
      as = - 1d0
  816 continue
      q2=q*q
      q4=q2*q2
      q6=q4*q2
      d=q2+1d0-M2
      RL(1)=-RC2+2D0*C*q-d/2D0-q*d/8D0/C
     *      -(4D0*q2*(d+1D0)+d**2)/64D0/RC2
     *      -q*(33d0*q4+114d0*q2+37d0-2d0*m2*(23d0*q2+25d0)
     *      +13*m2**2)/512d0/rc2/c
     *      -(63d0*q6+340d0*q4+239d0*q2+14d0
     *      -10d0*m2*(10d0*q4+23d0*q2+3d0)+m2*m2*(39d0*q2-18d0)
     *      -2d0*m2**3)/1024d0/rc4
     *      -q*(527d0*q6+4139d0*q4+5221d0*q2+1009d0
     *      -m2*(939d0*q4+3750d0*q2+1591d0)+m2*m2*(465d0*q2+635d0)
     *      -53d0*m2**3)/8192d0/rc4/c
     *      +as*(4d0*c)**(q+1d0)/cdexp(2d0*c)/fact(n+1)/fact(n+m+1)
     *      *(1d0-(3d0*q2+4d0*q+1d0-m2)/8d0/c
     *      +((9d0*q4+4d0*q2*q-18d0*q2-12d0*q-7d0)
     *      -m2*(6d0*q2-4d0*q-6d0)+m2*m2)/128d0/rc2
     *      -((23d0*q6-44d0*q4*q-171d0*q4+400d0*q2*q+257d0*q2+444d0*q
     *      +51d0)-m2*(19d0*q4-56d0*q2*q-50d0*q2+240d0*q+as*55d0)
     *      +m2*m2*(5d0*q2-12d0*q+5d0)-m2*m2*m2)/3d0/1024d0/rc2/c)

               IF(dreal(vl(in-1)).lT.0d0.and.dreal(vl(in-2)).lT.0d0.
     *             and.dreal(rl(1)).gT.0d0) then
                           rl(1) = s
               end if

          r0=rl(1)
          IF(dabs(dreal(vl(in)-r0)).gT.cn) then
          iv = 401
          if(dreal(sc).lt.0d0) sc = - sc
          go to 2
          end if
          end if
c             write(*,*) in,r0,vl(in)

  111  CONTINUE
      N = nn
      N = M + nN - 1
      N1 = N - M
      N3=2*N-1
      N4=N3+4
      N2=N+M
      R_cont = N*(N+1D0)-AK*RC2/2D0*(AM/N3/N4-1D0)+RC4/2D0*
     *((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *(N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))

      if(cdabs(r_cont - vl(nn)).gt.3d0*ac) then
      inull = inull + 1
      write(7,*) '**   ', r_cont
      write(7,*) '**   ', vl(nn)
      write(*,*) '**   ', r_cont
      write(*,*) '**   ', vl(nn)
      if(inull.eq.2) then
      ier = 8
      go to 900
      end if
      go to 991
      end if

  900 CONTINUE
      return
      END

C**********************************************************************
C CDLAMo_s  - Oblate eigenvalues (small values)
C**********************************************************************
C
      SUBROUTINE CDLAMo_s(VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
c      COMPLEX*16 VL,C,dx,alam0
      COMPLEX*16 VL,C,dx
      DIMENSION VL(NN),RL(550),SL(550),RDA(550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal

*****************************************
c      COMMON /r0/ alam0(500)
*****************************************

      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
      IER=0
c -------------------------------------
C**  Parameters check-up
c -------------------------------------
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s

      DO 113 IN = 1,550
          rl(in) = s
          sl(in) = s
 113      rda(in) = s
c      sdelta=s

      M2=M*M
      MM=2*M
      AM=4D0*M2-1D0
      RC2=C*C
      RC4=RC2*RC2
      NC=idint(aC)
             jcn=0
             DO 111 IN=1,NN

              inout = 8

             dx=dcmplx(0.1d0,0.1d0*af)
             IF(aC.gt.20D0.and.m.gt.10)
     *             dx = dcmplx(0.1d0+(m-10d0)/100d0,0.1d0*af)
             n1001=0
 1113 continue
             JN=0
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
      IF(NMAX.GT.500) WRITE(7,2000) NMAX
 2000 FORMAT(14X,'******* CDLAMo_s ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE
      IF(aC.LT.5.50D0) GO TO 401
c      IF(in.lT.3) GO TO 401
      IF(in.lT.7) GO TO 401
ccccc      IF(in.lT.7.and.m.lt.12) GO TO 401
      IF(n.lT.Nc-nc/3) GO TO 401
      IF(dreal(vl(in-1)).lT.0d0) GO TO 401
c      IF(dreal(vl(in-1)-vl(in-2)).lT.n.and.mod(n-m,2).ne.0) GO TO 401
      r0=VL(IN-1)+2d0*N+1d0
c        IF(aC.gt.22D0) r0 = r0 - 2.5d0
        IF(aC.gt.18D0) r0 = r0 - 2.532d0
             jcn=jcn+1

ccc          write(*,*) n,r0

         IF(m.eq.1.and.aC.gt.15.998D0.and.aC.lt.16.003D0.and.in.eq.13)
     *      then
            r0 = 66.1923d0 - (ac - 16d0) * 10d0
            dx = dx / 100d0
         end if

         IF(m.eq.7.and.aC.gt.18.690D0.and.aC.lt.18.695D0.and.in.eq.10)
     *      then
            r0 = 134.161d0 - (ac - 18.691d0) * 13.69744d0
            dx = dx / 5.5d0
         end if
         IF(m.eq.8.and.aC.gt.18.535D0.and.aC.lt.18.540D0.and.in.eq.8)
     *      then
            r0 = 116.313d0 - (ac - 18.535d0) * 13.69744d0
            dx = dx / 5.5d0
         end if

      GO TO 2
  400 CONTINUE
         JN=JN+1
           RL(1)=R0-JN*dx*SC
         IF(JN.EQ.1001) GO TO 1112
         GO TO 2
  401 CONTINUE
      RL(1)=S
      IF(CC.LE.0.5D0) GO TO 2
      N3=2*N-1
      N4=N3+4
      IF(CC.GT.(N+3D0)) GO TO 1
c      IF(CC.GT.(N+3D0).and.cc.le.21d0) GO TO 1
c      IF(dabs(CC-(N+3D0)).gt.0.5d0.and.cc.gt.21d0) GO TO 1
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      N2=N+M
      RL(1)=N*(N+1D0)-AK*RC2/2D0*(AM/N3/N4-1D0)+RC4/2D0*
     *((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *(N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))
      r0=rl(1)
      GO TO 2
    1 CONTINUE
c -------------------------------------
C** Initial approximation (large 'c')
c -------------------------------------
      IF(MOD(N1,2)) 3,4,3
    3 q=N
      GO TO 16
    4 q=N+1d0
   16 continue
      q2=q*q
      q4=q2*q2
      q6=q4*q2
      d=q2+1d0-M2
      RL(1)=-RC2+2D0*C*q-d/2D0-q*d/8D0/C
     *      -(4D0*q2*(d+1D0)+d**2)/64D0/RC2
     *      -q*(33d0*q4+114d0*q2+37d0-2d0*m2*(23d0*q2+25d0)
     *      +13*m2**2)/512d0/rc2/c
     *      -(63d0*q6+340d0*q4+239d0*q2+14d0
     *      -10d0*m2*(10d0*q4+23d0*q2+3d0)+m2*m2*(39d0*q2-18d0)
     *      -2d0*m2**3)/1024d0/rc4
     *      -q*(527d0*q6+4139d0*q4+5221d0*q2+1009d0
     *      -m2*(939d0*q4+3750d0*q2+1591d0)+m2*m2*(465d0*q2+635d0)
     *      -53d0*m2**3)/8192d0/rc4/c
      r0=rl(1)
    2 CONTINUE
      IF(dreal(vl(in-1)).lT.0d0.and.jn.eq.1) r0=r0-2d0
      IF(dreal(vl(in-1)).lT.0d0.and.jn.eq.1.and.dreal(r0).gt.0d0)
     *                                      r0=r0-2d0
               RDEL=dx
      NUM=1

ccc          write(*,*) n,r0

c -------------------------------------
c      alam0(in)=r0
c -------------------------------------

  100 CONTINUE
c -------------------------------------
C*       Iteration: 'from below to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1=4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2D0*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*       Iteration: 'from top to below'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-2
      DO 12 I=2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB
      NUM=NUM+1

      if(num.gt.540) eps1 = 1d-10
c      if(num.gt.540) eps1 = eps1 * 10d0
      if(num.gt.540)  print *, num, eps1

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE
      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)
      IF(aC.LT.2.5D0) GO TO 111
      IF(in.le.2) GO TO 111
        Cc2=dreal(rl(1))
        Cc1=dreal(VL(IN))
        CFF=dabs(cc1-cc2)
        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an
        if(an.gt.cn.and.ac.lt.10d0) cn=an
         cbb=dreal(VL(IN))-dreal(VL(IN-1))
c        IF(CFF.gt.cn.or.cbb.lt.0d0) go to 400
        IF(CFF.gt.cn.or.cbb.le.0d0) go to 400
 1112  CONTINUE
c 3001  FORMAT(I2,3D15.8,1x,2d10.3)
 3002  FORMAT(1x,'s',2f7.3,2I5,2x,3(f10.3,f10.3))
        ch=r0
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        ck=vl(in)
       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,crn,cc0,ck
       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,crn,cc0,ck

            if(jn.lt.1001) then
            go to 111
            end if

            sc=-sc
            n1001=n1001+1

            if(n1001.ge.inout) then
       write(7,*) '**   Problems with computations of eigenvalues!'
       write(*,*) '**   Problems with computations of eigenvalues!'
            ier = 8
            go to 900
            end if

  878 continue

  119  CONTINUE
            if(mod(n1001,2).eq.0) dx=dx/3d0
            go to 1113

  111  CONTINUE
  900 CONTINUE
      return
      END

C**********************************************************************
C CDLAMP  - Prolate eigenvalues
C
      SUBROUTINE CDLAMp(VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,C,dx
      DIMENSION VL(NN),RL(550),SL(550),RDA(550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal

*****************************************
*      COMMON /r0/ alam0(500)
*****************************************

      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
*****      if(ac.lt.32.5+m-3.or.ac.gt.47.4d0) SC=-sc
      if(ac.lt.29d0) SC=-sc
      IER=0
c -------------------------------------
C*  Parameters check-up
c -------------------------------------
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s
      nt=0
      sdelta=s
      M2=M*M
      m4=m2*m2
      MM=2*M
      AM = 4D0*M2-1D0
      RC2=C*C
      RC4=RC2*RC2
      NC=idint(aC)
             jcn=0
             DO 111 IN=1,NN

              inout = 4
****             write(*,*) in,nt

             dx=dcmplx(0.1d0,0.1d0*af)
             if(ac.gt.60d0) dx=dcmplx(0.2d0,0.2d0*af)
             n1001=0
 1113 continue
             JN=0
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
      IF(NMAX.GT.1090) WRITE(7,2000) NMAX
 2000 FORMAT(14X,'******* CDLAMp ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE
      IF(aC.LT.11.90D0) GO TO 401
      IF(in.lT.7) GO TO 401
**      IF(n.lT.Nc-nc/3-1.and.nt.eq.0) GO TO 401
      IF(n.lT.Nc-nc/3-m.and.nt.eq.0) GO TO 401
       r0=VL(IN-1)+dabs(dreal(VL(IN-1))-dreal(VL(IN-2)))
*       if(ac.ge.30d0.and.m.gt.3) r0=r0+m/2d0
       if(ac.ge.30d0.and.m.gt.3) r0=r0+m/5d0
      if(ac.ge.16d0.and.ac.lt.29d0) r0=r0+1.5d0
      if(ac.ge.30d0.and.ac.lt.40d0) r0=r0+(ac-30d0)/2d0
      if(ac.ge.40d0.and.ac.lt.50d0) r0=r0+n/10d0-(ac-40d0)/3d0
      if(ac.ge.50d0.and.ac.le.60d0) r0=r0+(ac-50d0)/3d0
      if(ac.gt.60d0.and.ac.le.70d0) r0=r0+(ac-60d0)/3d0
      if(ac.gt.70d0.and.ac.le.80d0) r0=r0+(ac-70d0)/4d0
      if(ac.gt.80d0.and.ac.le.90d0) r0=r0-(ac-80d0)/8d0+3.5d0
      if(ac.gt.90d0.and.ac.le.100d0) r0=r0-(ac-90d0)/8d0+4d0
      if(ac.gt.100d0) r0=r0-(ac-100d0)/8d0+4.5d0
      if(ac.gt.110d0) r0=r0+(ac-110d0)/8d0+1.5d0

**********      sdelta=s

             jcn=jcn+1
      GO TO 2
  400 CONTINUE

         sdelta=s

         JN=JN+1
         RL(1)=R0-JN*dx*SC*AK
         IF(JN.EQ.1001) GO TO 1112
         GO TO 2
  401 CONTINUE
      RL(1)=S
      R0=S
      IF(CC.LE.4.3D0) GO TO 2
      IF(CC.GT.(N+3D0)) GO TO 1
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      an3=2d0*N-1d0
      an4=an3+4d0
      N2=N+M
      RL(1)=N*(N+1D0)-AK*RC2/2D0*(AM/an3/an4-1D0)+RC4/2D0
     *      *((N1-1D0)*N1*(N2-1D0)*N2/((an3+2D0)*an3**3*an4)-(N1+1D0)
     *      *(N1+2D0)*(N2+1D0)*(N2+2D0)/((an4+2D0)*an4**3*(an3+2D0)))
     *      +rc4*rc2*(4d0*m2-1d0)*((N1+1D0)*(N1+2D0)*(N2+1D0)*(N2+2D0)
     *      /(an3*(2d0*N+1D0)*an4**5*(an4+2D0)*(an4+4D0))
     *      -(N1-1D0)*N1*(N2-1D0)*N2
     *      /((an3-4D0)*(an3-2D0)*an3**5*(an3+2D0)*an4))
      r0=rl(1)
      GO TO 2
    1 CONTINUE
c -------------------------------------
C** Initial approximation (large 'c')
c -------------------------------------
       q=2d0*N1+1d0
       q2=q*q
       q3=q2*q
       q4=q3*q
      RL(1)=C*q+M2-(q2+5D0)/8D0-q*(q2+11D0-32D0*M2)/C/64D0
     *      -(5D0*(q2**2+26D0*q2+21D0)-384D0*M2*(q2+1D0))/1024D0/RC2
     *      -1d0/rc2/c*(1d0/128d0**2*(33d0*q4*q+1594d0*q3
     *      +5621d0*q)-m2/128d0*(37d0*q3+167d0*q)+m4/8d0*q)
     *      -1d0/rc2/rc2*(1d0/256d0**2*(63d0*q4*q2+4940d0*q4
     *      +43327d0*q2+22470d0)-m2/512d0*(115d0*q4
     *      +1310d0*q2+735d0)+3d0*m4/8d0*(q2+1d0))
     *      -1d0/rc2/rc2/c*(1d0/1024d0**2*(527d0*q4*q3+61529d0*q4*q
     *      +1043961d0*q3+2241599d0*q)-m2/32d0/1024d0*(5739d0*q4*q
     *      +127550d0*q3+298951d0*q)+m4/512d0*(355d0*q3+1505d0*q)
     *      -m4*m2*q/16d0)
      r0=rl(1)
    2 CONTINUE
               RDEL=0.1D0*SC
      NUM=1

*****************************************
      delta=sdelta
      IF(delta.gT.1d0) delta=1d0
      IF(delta.lT.-1d0) delta=-1d0
**      if(n1001.eq.0)  r0=r0+sdelta
      if(n1001.eq.0)  r0=r0+delta
ccc         r0=r0+delta
*      alam0(in)=r0
c      write(*,*) in, r0,sdelta

        IF(cdabs(sdelta).gT.0.8d0.and.nt.eq.0) then
        sdelta=s
        nt=1
        end if
*****************************************

  100 CONTINUE
c -------------------------------------
C*            Iteration: 'from bottom to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1 = 4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2D0*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*            Iteration: 'from top to bottom'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-  2
      DO 12 I=  2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB
      NUM=NUM+1

      if(num.gt.540) eps1 = 1d-10
*      if(num.gt.500) eps1 = 1d-08
      if(num.gt.540)  print *, num

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE
      IF(NUM.gE.1009) GO TO 15
      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)
       sdelta=vl(in)-r0
      IF(aC.LT.11.90D0) GO TO 111
      IF(in.lT.3) GO TO 111
        Cc2=r0
        Cc1=VL(IN)
        CFF=dabs(cc1-cc2)
        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an
        IF(CFF.ge.cn) go to 400
 1112  CONTINUE
 3002  FORMAT(1x,'p',2f6.2,2I5,2x,3(f10.3,f10.3))
        if(jn.gt.1) then
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        WRITE(7,3002) c,N,jn,cc2,cr,crn,cc0,cc1
        WRITE(*,3002) c,N,jn,cc2,cr,crn,cc0,cc1
        end if
            if(jn.lt.1001) go to 111
            sc=-sc
            n1001=n1001+1
            if(mod(n1001,2).eq.0) dx=dx/1.5d0

            if(n1001.ge.inout) then
       write(7,*) ac
       write(*,*) ac
       write(7,*) '**   Problems with computations of eigenvalues!'
       write(*,*) '**   Problems with computations of eigenvalues!'
            ier = 8
            go to 900
            end if

c            if(mod(n1001,48).eq.0) go to 900
cyy            if(mod(n1001,2).eq.0) go to 900
            go to 1113
  111  CONTINUE
  900 CONTINUE
      return
      END

c************************************************
C CESSEL
c
      SUBROUTINE CESSEL(A,NUM,BESJ,BESY)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1),BESY(NUM+1)
      BESJ(NUM+1)=(0D0,0D0)
      BESJ(NUM)=(1D-260,1D-260)
      N=2*NUM+1
      NUM1=NUM-1
      DO 11 I=1,NUM1
      N=N-2
      I1=NUM-I
   11 BESJ(I1)=N*A*BESJ(I1+1)-BESJ(I1+2)
      B1=A*BESJ(1)-BESJ(2)
      B2=-A*B1-BESJ(1)
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
   13 BESJ(I)=C*BESJ(I)
      BESY(1)=-C*B1
      BESY(2)=C*B2
      DO 14 I=3,NUM
      I2=I-2
   14 BESY(I)=(2.0D0*I2+1.0D0)*A*BESY(I-1)-BESY(I2)
      RETURN
      END

c************************************************
C CMLM1
C
      SUBROUTINE CMLM1(N,A1,A2,A3)
      parameter (nterms=100)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      real*8 A1(nterms,nterms)
      DIMENSION A2(nterms,nterms),A3(nterms,nterms),A4(nterms,nterms)
      DO 1 I=1,N
      DO 1 J=1,N
      B=(0D0,0D0)
      DO 2 K=1,N
    2 B=B+A1(I,K)*A2(K,J)
      A4(I,J)=B
    1 CONTINUE
      DO 3 I=1,N
      DO 3 J=1,N
    3 A3(I,J)=A4(I,J)
      RETURN
      END

c************************************************
C CMLM2
C
      SUBROUTINE CMLM2(N,A1,A2,A3)
      parameter (nterms=100)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      real*8 A2(nterms,nterms)
      DIMENSION A1(nterms,nterms),A3(nterms,nterms),A4(nterms,nterms)
      DO 1 I=1,N
      DO 1 J=1,N
      B=(0D0,0D0)
      DO 2 K=1,N
    2 B=B+A1(I,K)*A2(K,J)
      A4(I,J)=B
    1 CONTINUE
      DO 3 I=1,N
      DO 3 J=1,N
    3 A3(I,J)=A4(I,J)
      RETURN
      END

c************************************************
C CXSC11
C
      SUBROUTINE CXSC11(II,M,C,NE,NN,QX,QS,QEXT,QSCA)
      parameter (nterms=100, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      REAL*8      KAP
      DIMENSION B0(nterms),ACOF(nterms,4*nterms),VV(nterms),
     *          OMG(nterms,nterms),TU(nterms,nterms),
     *          BD(nalphas,2*nterms),KAP(nterms,nterms),IALFA(nalphas),
     *          JALFA(nalphas),RB0(nalphas,nterms),
     *          RA(nalphas,nterms),RB (nalphas,nterms),ASIN(nalphas),
     *          QX(nalphas),QS(nalphas)
      COMMON /EPS3/ EPS3
       COMMON /kc1/ ACOF,VV
       COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
       COMMON /K2/ IALFA, JALFa
       COMMON /K10/ B0
       COMMON /K11/ BD
       COMMON /FF/ ASIN,acos(nalphas)
       COMMON /RES/ RB0, RA, RB

       CALL INT3(M,NE,KAP,OMG,TU)

       IF(II.EQ.0) GO TO 101
       CK=C**2*DSQRT(AKSI)
       DO 20 IA=1,NAL
       DK=CK*DSQRT(AKSI+AK*ASIN(IA)**2)
       R=S
       D=0D0
       E=0D0
       IF(JALFA(IA).EQ.1) GO TO 20
       IF(M.GT.1) GO TO 23
       DO 21 L=1,NN

c         write(7,*) l,rb0(1,l)

       V=VV(L)**2
       R1=S1**(-L)*RB0(IA,L)*BD(IA,L)*V
       D1=RB0(IA,L)*DCONJG(RB0(IA,L))*V
       R=R+R1
       D=D+D1
c      IF(CDABS(R1).LT.1D-130.AND.DABS(D1).LT.1D-130.
      IF(CDABS(R1).LT.eps3.AND.DABS(D1).LT.eps3.
     *   AND.IALFA(IA).NE.90)    GO TO 22
   21 CONTINUE
   22 D=2D0*D

ccc         go to 26

   23 CONTINUE
       DO 24 L=1,NN
       LM=L+M-1
       R1=-S1**(-LM)*(S1*RA(IA,L)*BD(IA,L)*VV(L)**2-
     *    RB(IA,L)*BD(IA,L+NE))*ASIN(IA)
       R=R+R1
       DO 25 N=1,NN
       NM=N+M-1
        RAN=DCONJG(RA(IA,N))
        RBN=DCONJG(RB(IA,N))
      D1=S1**(NM-LM)*(RA(IA,L)*RAN*OMG(L,N)+
     *  S1*(RB(IA,L)*RAN* KAP(L,N)-
     *  RA(IA,L)*RBN* KAP(N ,L))+
     *  RB(IA,L)*RBN*TU(L,N))*VV(N)*VV(L)
       D=D+D1
c       IF(DABS(D1).LT.1D-130.AND.IALFA(IA).NE.90) GO TO 27
       IF(DABS(D1).LT.eps3.AND.IALFA(IA).NE.90) GO TO 27
   25 CONTINUE
c   27  IF(CDABS(R1).LT.1D-130.AND.IALFA(IA).NE.90) GO TO 26
   27  IF(CDABS(R1).LT.eps3.AND.IALFA(IA).NE.90) GO TO 26
   24 CONTINUE
   26 CONTINUE
      R = R * 4D0/DK
      QX(IA)=R
      QS(IA)=D/DK
   20 CONTINUE
        IF(II.nE.0) GO TO 100
c.................................................................
c     alpha = 0 deg.
c.................................................................
  101 CONTINUE
      D=0D0
      E=0D0
      CK=C**2*AKSI
      DO 6 L=1,NN

c         write(7,*) l,ra(1,l)
c         write(7,*) rb(1,l)
c         write(7,*) b0(l)

      DO 5 N =1,NN
        RAN=DCONJG(RA(1,N))
        RBN=DCONJG(RB(1,N))
      G=S1**(N-L)*(RA(1,L)*RAN*OMG(L,N)+S1*(RB(1,L)*RAN*KAP(L,N)-
     *  RA(1,L)*RBN*KAP(N,L))+RB(1,L)*RBN*TU(L,N))*VV(L)*VV(N)
c      IF(DABS(G).LE.1D-120.AND.L.GT.3) GO TO 7
      IF(DABS(G).LE.eps3.AND.L.GT.3) GO TO 7
    5 E=E+G
    7 G1=S1**(-L)*RB(1,L)*B0(L)
c      IF(DABS(G1).LE.1D-120.AND.L.GT.3) GO TO 11
      IF(DABS(G1).LE.eps3.AND.L.GT.3) GO TO 11
      D=D+G1
    6 CONTINUE
   11 QEXT=-2D0*D/CK
      QSCA=E/CK    
  100 RETURN
      END

c************************************************
c  daswf1
C                                                                
      subroutine DASWF1(ar,adr,is,M,n,x,icos,AD,ia,IER)
      parameter (nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      DIMENSION AD(IA),P(nalphas,250),PD(nalphas,250)
      COMMON /EPS3/ EPS3
      COMMON /LEG1/ P,PD
      ier=0
      ar=0d0                                                           
      aDR=0d0
      if(dabs(x-1d0).lt.1d-6) return
      if(dabs(x).lt.1d-6) go to 14
      IF(MOD((N-M),2)) 1,2,1                                          
    1 JK=1                                                
      GO TO 5                                                         
    2 JK=0                                               
    5 CONTINUE
      DO 13 I=1,IA                                                 
      J1=2*(I-1)+1+JK
      aRS=AD(I)*P(icos,J1)
      aRD=AD(I)*PD(icos,J1)
      aR=aR+aRS
      aDR=aDR+aRD                                                  
c      IF(DABS(aRS).LT.1D-160.AND.DABS(aRD).LT.1D-160) GO TO 900
      IF(DABS(aRS).LT.eps3.AND.DABS(aRD).LT.eps3) GO TO 900
   13 CONTINUE                   
      go to 900                
c                          x = 0
   14 CONTINUE                        
      CALL DLEGF00(ar,adr,n,M,is,IER1)
  900 continue
      RETURN
      END                    

c************************************************
c DELta21
c
      SUBROUTINE DELta21(M,NN,DEL)
      parameter (nalphas=19,nterms=100)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 BCOF,B,A,VV2,DEL,C,E1,G1,D,s,s1
      DIMENSION A(4*nterms),B(4*nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms),
     *          ACOF(nterms,4*nterms), VV1(nterms),
     *          DEL(nterms,nterms)
      COMMON /EPS3/ EPS3
      COMMON /Kc1/  ACOF, VV1
      COMMON /Kc2/  BCOF, VV2
      COMMON /K1/ S, S1, AKSI, AK, kkK, NK, nal
       M2=2*M
      DO 10 N=1,NN
        NM=N+M-1

      DO I = 1, NK
      B(I) = BCOF(N,I)
      end do
      C = VV2(N)

      DO 11 L=1,NN
        LM=L+M-1
        IF(MOD(IABS(LM-NM),2)) 9,9,20
    9 CONTINUE

      DO I = 1, NK
      a(I) = aCOF(l,I)
      end do
      d = VV1(l)

      IF(MOD(LM-M,2))1,1,2
    1 I2=1
      K=M2
      E1=A(1)*B(1)*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E1=A(1)*B(1)*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E1=E1*F
      DO 5 I=2,NK
      J=2*I-1-I2
      JM=J+M
      JM2=2*JM
      F=1D0
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G1=A(I)*B(I)*2D0/(JM2+1D0)*F
c      IF(CDABS(G1).LE.1D-140.and.i.gt.nn) GO TO 6
      IF(CDABS(G1).LE.eps3.and.i.gt.nn) GO TO 6
      E1=E1+G1
    5 continue
    6 DEL(N,L)=E1/C/D
      GO TO 11
   20 DEL(N,L)=(0D0,0D0)
   11 CONTINUE
   10 CONTINUE
      RETURN
      END

c************************************************
C DLEGf00
C
      subroutine DLEGF00(res,red,n,M,id,iER)
      IMPLICIT REAL*8(A-H,O-Z)                                        
      COMMON  /FACT/ FACT(170)
      ier=0
      res=0d0
      red=0d0
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(ier.ne.0) go to 100
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(N-M),2)) 50,1,50
    1 continue
      k1=(m+n)/2
      e=1d0
      do 2 kk=1,k1
    2 e=e*(k1+kk)
      k2=(n-m)/2
      res=(-1)**k2*e/fact(k2+1)/2d0**n
      if(id) 50,100,50
   50 continue
c=========================             <<     Ø‡Æ®ß¢Æ§≠†Ô    >>
      if(id.eq.0) go to 100
c************************************  << x1 = 0 >>
      IF(MOD(iabs(N-M),2)) 71,100,71
   71 continue
      k1=(m+n+1)/2
      e=1d0
      do 72 kk=1,k1
   72 e=e*(k1+kk)
      k2=(n-m-1)/2
      red=(-1)**k2*e/fact(k2+1)/2d0**n   
  100 CONTINUE                                                      
      RETURN                                                        
      END      

c************************************************
C DLEGf1
C   
      subroutine DLEGF1(rs,rd,X1,x2,x3,n,M,ic,id,iER)
      IMPLICIT REAL*8(A-H,O-Z) 
      complex*16 s, s1
      COMMON /FACT/ FACT(170)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      ier=0
      res=0d0
      red=0d0
      xx2=dabs(x2)
      nm=n-m
      j=mod(n,4)
      jm=mod(m,4)
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(id.ne.0.and.dabs(x1-1d0).lt.1d-6.and.ic.eq.0) ier=2
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6)  go to 20
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(nm),2)) 100,1,100
    1 continue
      k1=(m+n)/2
      e=1d0
      do 2 k=1,k1
    2 e=e*(k1+k)
      res=e/fact(nm/2+1)/2d0**n
      if(id) 50,102,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.gt.n) go to 100
      if(nm.ge.2.and.x1.gt.1d0.and.ic.eq.0) go to 48
      z=1d0/(2d0**n)
      k1=(nm+2)/2
      d=0d0
      if(ic.ne.0) go to 21
c.....................................   << ic = 0 >>
      if(nm.gt.169) go to 291
      do 22 k=1,k1
      k3=nm-2*(k-1)
      h=1d0/fact(k3+1)
      k4=n-k+1
      do 27 i=1,k4
   27 h=h*(k4+i)
   22 d=d+(-1d0)**(k-1)*h/fact(k)*x1**k3
      go to 292
  291 continue
      do 293 k=1,k1
      k3=nm-2*(k-1)
      k5=k3-169
      h=1d0/fact(169)
      k4=n-k+1
      k6=k4/2
      do 294 i=1,k6
  294 h=h*(k4+i)
      do 295 i=170,k5
  295 h=h/i
      h=h/fact(k)
      do 296 i=k6+1,k4
  296 h=h*(k4+i)
  293 d=d+(-1d0)**(k-1)*h*x1**k3
  292 continue
       res=d*z*dsqrt(dabs(1d0-x1**2)**m)
      if(id) 50,102,50     
c<<<<<<<<<++++++++++++++++++++++++++++   << ic /= 0 >>
   21 continue
      if(nm.gt.169) go to 29
      do 28 k=1,k1
      k3=nm-2*(k-1)
      h=1d0/fact(k3+1)
      k4=n-k+1
      do 33 i=1,k4
   33 h=h*(k4+i)
   28 d=d+h/fact(k)*x1**k3
      go to 31
   29 continue
      do 30 k=1,k1
      k3=nm-2*(k-1)
      k5=k3-169
      h=1d0/fact(169)
      k4=n-k+1
      k6=k4/2
      do 34 i=1,k6
   34 h=h*(k4+i)
      do 35 i=170,k5
   35 h=h/i
      h=h/fact(k)
      do 36 i=k6+1,k4
   36 h=h*(k4+i)
   30 d=d+h*x1**k3
   31 continue
      res=d*z*dsqrt((1+x1**2)**m)
      if(id) 50,102,50     
   48 continue
      res=((2d0*n-1d0)*x1*x2-(n+m-1d0)*x3)/nm
      if(id) 50,102,50     
   50 continue
c=========================             <<  derivature  >>
c*****************************************   << x1 /= + - 1 >>
      if(ic.ne.0) go to 101
      red=((nm)*res-x1*n*x2)/(x1**2-1d0)
       rs=res
       rd=red
       return
  101 continue
      red=((nm)*res-x1*n*xx2)/(x1**2-ak)
  100 CONTINUE        
      if(jm.eq.2.or.jm.eq.3) go to 103
      rs=res
      if(j.eq.2.or.j.eq.3) rs=-dabs(res)
      go to 104
  103 continue
      rs=-dabs(res)
      if(j.eq.0.or.j.eq.1) rs=-rs
  104 continue
      j=mod(n-1,4)                          
      rd=red
      if(j.eq.2.or.j.eq.3) rd=-dabs(red)
      RETURN      
  102 continue
       rs=res
      if(ic.eq.0.and.jm.eq.3.and.x1.gt.1d0) rs=dabs(res)
      if(ic.ne.0.and.jm.eq.2) rs=-dabs(res)
      if(ic.ne.0.and.jm.eq.3) rs=-dabs(res)
      RETURN                                               
      END

c************************************************
C DLEGf2
C   
      subroutine DLEGF2(res,red,X1,x2,n,M,ic,id,iER) 
      IMPLICIT REAL*8(A-H,O-Z)  
      complex*16 s1, s
      dimension r(500),q(500)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      COMMON /FACT/ FACT(170)
      COMMON /PI/ PI
      ier=0
      res=0d0
      red=0d0
      if(n+m.lt.0) ier=1
      if(ic.eq.0.and.dabs(x1).le.1d0) ier=2
      if(ic.ne.0.and.dabs(x1).lt.0d0) ier=2
      if(n+m-1.lt.0.and.id.ne.0) ier=3
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6.or.ic.eq.0)  go to 20
c**************************************   << x1 = 0 & ic /= 0 >>
      if(m.gt.n) go to 3
c                             m  < = n
      k1=n-m
      k2=n+m-1
      d=1d0
      e=1d0
      IF(MOD(k1,2)) 2,1,2
    1 continue
      if(k1.eq.0) go to 6
      do 4 i=1,k2,2
    4 d=d*i
      do 5 i=2,k1,2
    5 e=e*i
    6 res=-pi/2d0*d/e
      if(id) 50,100,50
    2 continue
      do 7 i=2,k2,2
    7 d=d*i
      do 8 i=1,k1,2
    8 e=e*i
      res=d/e
      if(id) 50,100,50
c                                m  >  n 
    3 continue
      IF(MOD(iabs(N-M),2)) 10,9,10
    9 go to 100
   10 continue
      k1=(n-m-1)/2
      k2=(m+n-1)/2
      d=1d0
      do 11 i=1,k1
   11 d=d*(k1+i)
      res=2d0**n*d*fact(k1+1)
      if(id) 50,100,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.eq.0)  go to 61
      if(m.eq.1)  go to 71
c *********************************     m  > 1
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      r(1)=-ak/(dsqrt(x1**2-ak))
      r(2)=-2d0*x1*ak/(dsqrt(x1**2-ak))*r(1)
      if(m.eq.2) go to 221
      do 21 i=3,m
   21 r(i)=-2d0*(i-1d0)*x1*ak/(dsqrt(x1**2-ak))*r(i-1)+
     *     (2d0-i)*(i-1d0)*r(i-2)
  221 continue
      q0=r(m)
      if(n.eq.0) res=q0
      if(n.eq.0) go to 26
      j=n
      if(n.lt.0) j=-n
      nmax=30*x1+j+50
      if(nmax.gt.500) nmax=500
      r(nmax)=x1-dsqrt(x1**2-ak)
      xer=dabs(x1-1d0/dsqrt(3d0))
      do 22 i=nmax-1,1,-1
   22 r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
      q(1)=r(1)*q0*ak
      if(xer.lt.1d-10.and.m.eq.3) q(1)=3d0
      if(xer.lt.1d-10.and.m.eq.6) q(1)=-144d0
      if(xer.lt.1d-10.and.m.eq.9) q(1) = 4.536d4
      if(xer.lt.1d-10.and.m.eq.12) q(1) = -4.35456d+07
      if(xer.lt.1d-10.and.m.eq.15) q(1) = 9.340529643126767d+10
      if(xer.lt.1d-10.and.m.eq.18) q(1) = -3.76610155207262d+14
      if(xer.lt.1d-10.and.m.eq.21) q(1) = 2.554546682742534d+18
      if(xer.lt.1d-10.and.mod(m,3).eq.0) q0=0d0
      if(n.eq.1) go to 24
      if(n.lt.0) go to 225
      do 23 i=1,j-1
   23 q(i+1)=r(i+1)*q(i)*ak**(i+1)
       go to 24
c******************************************      n < 0
  225 continue 
      if(n.eq.-1) res=(ak*x1*q0+(m-1d0)*q(1))/m
      if(n.eq.-1) go to 26
      r(1)=q0
      r(2)=(ak*x1*q0+(m-1d0)*q(1))/m
      do 25 i=-2,n,-1
      r(1-i)=(ak**i*(2d0*i+3d0)*x1*r(-i)-
     *       (i-m+2d0)*r(-i-1))/(i+m+1d0)
   25 continue 
      res=r(-n+1)
      go to 26
   24 continue
      res=q(n)
   26 continue
      if(id) 50,100,50
c***********************************    m = 0
   61 continue
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      if(n.eq.0) go to 65
      q0=res
      nmax=1.5*x1+n+10
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 62 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   62 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 64
      do 63 i=1,n-1
      q(i+1)=r(i+1)*q(i)*ak**(i+1)
   63 continue
   64 continue
      res=q(n)
   65 if(id) 50,100,50
c***********************************    m = 1
   71 continue
      res=-ak/(dsqrt(x1**2-ak))
      if(n.eq.0) go to 75
      q0=res
      if(n.eq.-1) res=ak*x1*q0
      if(n.eq.-1) go to 75
      nmax=30*x1+n+50
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 72 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   72 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 74
      do 73 i=1,n-1
   73 q(i+1)=r(i+1)*q(i)*ak**(i+1)
   74 continue
      res=q(n)
   75 if(id) 50,100,50
   50 continue
c=========================             <<     Ø‡Æ®ß¢Æ§≠†Ô    >>
      if(id.eq.0) go to 100
      if(ic.eq.0.or.mod(n,2).eq.0)  go to 174
      red=-((n-m)*res+x1*n*x2)/(x1**2-ak)
      go to 100
  174 continue
      red=((n-m)*res-x1*n*x2)/(x1**2-ak)
  100 CONTINUE                                                      
      RETURN                                                        
      END

c************************************************
C DRF222
c
      SUBROUTINE DRF222(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,RL
      DIMENSION Y(4),YP(4)
      COMMON /F2/ RL,C,M,i12
       AKSI=X**2+1D0

      if(i12.eq.0.or.i12.eq.1) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2)+(RL-(C*X)**2-M**2/AKSI)*Y(1))/AKSI
      end if

      if(i12.eq.2) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2)+(RL-(C*X)**2-M**2/AKSI)*Y(1))/AKSI
       YP(3)=Y(4)
       YP(4)=(-2D0*X*Y(4)+(RL-(C*X)**2-M**2/AKSI)*Y(3))/AKSI
      end if

      RETURN
      END

c************************************************
C DRF2K02
C
      SUBROUTINE DRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,R1,S1,r,s,s2,t,r2
      DIMENSION AD(IA),BD(IB),FACT(170),Y(4)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT
      COMMON /pi/ pi
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 3
    2 JK=0
    3 CONTINUE
      S2 = s
      DO 4 I=1,IA
      J=2*(I-1)+JK
      F=1D0
      DO 5 JJ=1,M2
    5 F=F*(J+JJ)
      T=AD(I)*F
      S2=S2+T
c      IF(cDABS(T).LT.1D-130.and.i.gt.ia/2) GO TO 6
      IF(cDABS(T).LT.eps3.and.i.gt.ia/2) GO TO 6
    4 CONTINUE
    6 CONTINUE
      MF=N+M+1+JK
      IF(MF.LE.170) F=FACT(MF)
      IF(MF.LE.170) GO TO 7
      F=FACT(170)
      DO 8 I=171,MF
    8 F=F*I
    7 CONTINUE

      if(i12.eq.1.or.i12.eq.2) then
      R=S1**(N-M-JK)*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
     *  FACT(M2+1)/BD(M+1)/S2*
     *  (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1))**2
      R1=S1**(M-N-JK)*S2/2D0**M/C**(M+1+JK)/AD(1)/FACT(M+1)
      end if


      if(i12.eq.0.or.i12.eq.2)
     *            R2=S1**(n-m-JK)*2D0**M*C**(M+JK)*AD(1)*FACT(M+1)/s2


      IF(JK.EQ.1) GO TO 10

      if(i12.eq.0) then
      Y(1)=R2/(m2+1d0)
      Y(2)=0D0
      end if

      if(i12.eq.1) then
      Y(1)=R
      Y(2)=R1*(M2+1D0)
      end if

      if(i12.eq.2) then
      Y(1)=R2/(m2+1d0)
      Y(2)=0D0
      Y(3)=R
      Y(4)=R1*(M2+1D0)
      end if

      GO TO 900

   10 CONTINUE

      if(i12.eq.0) then
      Y(1)=0D0
      Y(2)=R2/(M2+3D0)
      end if

      if(i12.eq.1) then
      Y(1)=R1*(M2+3D0)
      Y(2)=R*(M2-3D0)
      end if

      if(i12.eq.2) then
      Y(1)=0D0
      Y(2)=R2/(M2+3D0)
      Y(3)=R1*(M2+3D0)
      Y(4)=R*(M2-3D0)
      end if

  900 CONTINUE

c      write(*,*) n
c      write(*,22) y(1)
c      write(*,22) y(2)
c      write(7,22) y(1)
c      write(7,22) y(2)
c      write(*,22) y(3)
c      write(*,22) y(4)
c      write(7,22) y(3)
c      write(7,22) y(4)
c  22  format(1x,2d24.15)

      RETURN
      END

c************************************************
C drkf45
c
      SUBROUTINE DRKF45 (F,NEQN,Y,T,TOUT,
     1   RELERR,ABSERR,IFLAG,WORK,IWORK)
      INTEGER IWORK(5),NEQN
      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,
     1   WORK(27)
      EXTERNAL F
      K1M = NEQN + 1
      K1 = K1M +1
      K2 = K1 + NEQN
      K3 = K2 + NEQN
      K4 = K3 + NEQN
      K5 = K4 + NEQN
      K6=K5+NEQN
      CALL DRKFS (F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),
     1   WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
     2   WORK(K6),WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),
     3   IWORK(4),IWORK(5))
      RETURN
      END

c************************************************
C drkfs
c
      SUBROUTINE DRKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,
     1   IFLAG,YP,H,F1,F2,F3,F4,F5,SAVRE,SAVAE,
     2   NFE,KOP,INIT,JFLAG,KFLAG)
      LOGICAL  HFAILD, OUTPUT
      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,H,YP(NEQN),F1(NEQN),
     1   F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN),SAVRE,SAVAE
      EXTERNAL F
      REAL*8  A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,REMIN,
     1   RER,S,SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK
      DATA REMIN /1.0D-16/
      DATA MAXNFE /500000/
      IF(NEQN.LT.1) GO TO 10
      IF((RELERR.LT.0D0).OR.(ABSERR.LT.0D0)) GO TO 10
      MFLAG = IABS(IFLAG)
      IF((MFLAG.EQ.0).OR.(MFLAG.GT.8))  GO TO 10
      IF(MFLAG.NE.1) GO TO 20
      EPS = 1.0D0
    5 EPS = EPS/2.0D0
      EPSP1 = EPS + 1.0D0
      IF(EPSP1.GT.1.0D0) GO TO 5
      U26 = 26.0D0 * EPS
      GO TO 50
   10 IFLAG =8
      RETURN
   20 IF((T.EQ.TOUT).AND.(KFLAG.NE.3)) GO TO 10
      IF(MFLAG.NE.2) GO TO 25
      IF((KFLAG.EQ.3).OR.(INIT.EQ.0)) GO TO 45
      IF(KFLAG.EQ.4) GO TO 40
      IF((KFLAG.EQ.5).AND.(ABSERR.EQ.0.0D0)) GO TO 30
      IF((KFLAG.EQ.6).AND.(RELERR.LE.SAVRE).AND.
     1   (ABSERR.LE.SAVAE)) GO TO 30
      GO TO 50
   25 IF(IFLAG.EQ.3) GO TO 45
      IF(IFLAG.EQ.4) GO TO 40
      IF((IFLAG.EQ.5).AND.(ABSERR.GT.0.0D0))  GO TO 45
   30 STOP
   40 NFE = 0
      IF(MFLAG.EQ.2) GO TO 50
   45 IFLAG = JFLAG
      IF(KFLAG.EQ.3)  MFLAG = IABS(IFLAG)
   50 JFLAG = IFLAG
      KFLAG = 0
      SAVRE = RELERR
      SAVAE = ABSERR
      RER = 2.0D0*EPS + REMIN
      IF(RELERR.GE.RER)  GO TO 55
      RELERR = RER
      IFLAG = 3
      KFLAG = 3
      RETURN
   55 DT = TOUT -T
      IF(MFLAG.EQ.1)  GO TO 60
      IF(INIT.EQ.0)  GO TO 65
      GO TO 80
   60 INIT = 0
      KOP = 0
      A = T
      CALL F(A,Y,YP)
      NFE = 1
      IF(T.NE.TOUT) GO TO 65
      IFLAG = 2
      RETURN
   65 INIT = 1
      H = DABS(DT)
      TOLN = 0.0D0
      DO 70 K = 1,NEQN
         TOL=RELERR*DABS(Y(K)) + ABSERR
          IF(TOL.LE.0.0D0)  GO TO 70
         TOLN = TOL
         YPK = DABS(YP(K))
         IF(YPK *H**5.GT.TOL)  H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF(TOLN.LE.0.0D0)  H = 0.0D0
      H=1D-4
      H= DMAX1(H, U26*DMAX1(DABS(T),DABS(DT)))
      JFLAG = ISIGN(2,IFLAG)
   80 H =DSIGN(H,DT)
      IF(DABS(H).GE.2.0D0*DABS(DT))  KOP = KOP+1
      IF(KOP.NE.100) GO TO 85
      KOP = 0
      IFLAG = 7
      RETURN
   85 IF(DABS(DT).GT.U26*DABS(T))  GO TO 95
      DO 90  K = 1,NEQN
   90 Y(K) = Y(K) + DT*YP(K)
      A = TOUT
      CALL F(A,Y,YP)
      NFE = NFE+ 1
      GO TO  300
   95 OUTPUT=.FALSE.
      SCALE= 2D0/RELERR
      AE = SCALE * ABSERR
  100 HFAILD = .FALSE.
      HMIN = U26 * DABS(T)
      DT = TOUT -T
      IF(DABS(DT).GE.2.0D0*DABS(H))  GO TO 200
      IF(DABS(DT).GT.DABS(H))  GO TO 150
      OUTPUT= .TRUE.
      H = DT
      GO TO  200
  150 H = 0.5D0 * DT
  200 IF(NFE.LE.MAXNFE)  GO TO 220
      IFLAG = 4
      KFLAG = 4
      RETURN
  220 CALL FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE = NFE + 5
      EEOET=0D0
      DO 250  K = 1,NEQN
         ET = DABS(Y(K))+ DABS(F1(K)) + AE
         IF(ET.GT.0.0D0)  GO TO 240
               IFLAG =5
               RETURN
  240    EE = DABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  250    EEOET = DMAX1(EEOET,EE/ET)
      ESTTOL = DABS(H)*EEOET*SCALE/7.524D5
      IF (ESTTOL.LE.1.0D0)  GO TO 260
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1D0
      IF(ESTTOL.LT.5.9049D4)  S = 0.9D0/ESTTOL**0.2D0
      H = S * H
      IF(DABS(H).GT.HMIN)  GO TO  200
      IFLAG = 6
      KFLAG = 6
      RETURN
  260 T=T+H
      DO 270  K=1,NEQN
  270    Y(K) =F1(K)
      A = T
      CALL F(A,Y,YP)
      NFE = NFE + 1
      S = 5.0D0
      IF(ESTTOL.GT.1.889568D-4) S = 0.9D0/ESTTOL**0.2D0
      IF (HFAILD)  S = DMIN1(S,1.0D0)
      H= DSIGN( DMAX1(S*DABS(H),HMIN),H)
      IF(OUTPUT)  GO TO 300
      IF(IFLAG.GT.0) GO TO 100
      IFLAG = -2
      RETURN
  300 T = TOUT
      IFLAG = 2
      RETURN
      END

c************************************************
C DRSF212
C
      SUBROUTINE DRSF212(R,DR,ar,adr,i12,m,N,c,X,rl,AD,IA,BD,IB,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,AD(IA),RL,cc,rll
      DIMENSION WORK(27),IWORK(5),Y(4),BD(IB)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /F2/ RLl,Cc,Mm,i12a
      COMMON /FACT/ FACT(170)
      EXTERNAL  DRF222
   11 FORMAT(1X,I10,2D10.3,4D20.8)
   31 FORMAT(17H TOLERANCES RESET,2D12.3)
   41 FORMAT(11H MANY STEPS)
   71 FORMAT(12H MUCH OUTPUT)
   81 FORMAT(14H IMPROPER CALL)

          mm = m
          cc = c
          rll = rl
          i12a = i12

      NEQN=2
      if(i12.eq.2) NEQN=4
       T=0D0
cccc       ELERR=1D-5
       ELERR=1D-15
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL DRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
   10 CALL DRKF45(DRF222,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)

      if(i12.eq.0) then
      R=Y(1)
      DR=Y(2)
      end if

      if(i12.eq.1) then
      aR=Y(1)
      aDR=Y(2)
      end if

      if(i12.eq.2) then
      R=Y(1)
      DR=Y(2)
      aR=Y(3)
      aDR=Y(4)
      end if

c          write(7,*) c
c          write(7,*) ad(1)
c          write(7,*) bd(1)
c          write(7,*) r
c          write(7,*) dr
c          write(7,*) ar
c          write(7,*) adr
      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T.LT.TFINAL) GO TO 10
      IER=0
      GO TO  111
   30 WRITE(7,31)  ELERR,ABSERR
      IER=3
      GO TO 10
   40 WRITE(7,41)
      IER = 4
      GO TO 10
   50 ABSERR=1D-9
      WRITE(7,31)  ELERR,ABSERR
      IER=5
      GO TO 10
   60  ELERR=10D0 * ELERR
      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
      GO TO 10
   70 WRITE(7,71)
      IFLAG=2
      IER=7
      GO TO 10
   80 WRITE(7,81)
      IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
      END

c************************************************
C fehl
c
      SUBROUTINE FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
      REAL*8  Y(NEQN),T,H,YP(NEQN),F1(NEQN),F2(NEQN),
     1   F3(NEQN),F4(NEQN),F5(NEQN),S(NEQN),CH
      CH = H/4.0D0
      DO 221  K = 1,NEQN
  221    F5(K) = Y(K) + CH * YP(K)
      CALL F(T+CH,F5,F1)
      CH = 3.0D0 * H/3.2D1
      DO 222  K = 1,NEQN
  222    F5(K) = Y(K) + CH*(YP(K)+3.0D0*F1(K))
      CALL F(T+3.0D0*H/8.0D0,F5,F2)
      CH = H/2.197D3
      DO 223  K =1,NEQN
  223    F5(K) = Y(K)+CH*(1.932D3*YP(K)+(7.296D3*F2(K)-7.2D3*F1(K)))
      CALL F(T+1.2D1*H/1.3D1,F5,F3)
      CH = H/4.104D3
      DO 224  K= 1,NEQN
  224    F5(K)=Y(K)+CH*((8.341D3*YP(K)-8.45D2*F3(K))
     1         +(2.944D4*F2(K) - 3.2832D4*F1(K)))
      CALL F(T+H,F5,F4)
      CH = H/2.052D4
      DO 225 K= 1,NEQN
  225    F1(K)=Y(K)+CH*((-6.08D3*YP(K)+(9.295D3 *F3(K)-
     1         5.643D3*F4(K)))+(4.104D4*F1(K)-2.8352D4 *F2(K)))
      CALL F(T+H/2.0D0,F1,F5)
      CH = H/7.61805D6
      DO 230  K= 1,NEQN
  230    S(K)=Y(K)+CH*((9.0288D5*YP(K)+(3.855735D6 *F3(K) -
     1        1.371249D6*F4(K)))+(3.953664D6*F2(K)+2.7702D5*F5(K)))
      RETURN
      end

c************************************************
C FUNLEG
C
      SUBROUTINE FUNLEG(M,X,NF)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(250),PD(250),Q(250),QD(250)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P, PD, Q, QD
      IF(NF.GT.250) WRITE(7,1) NF
    1 FORMAT(1X,'FUNLEG  NF > 250, NF=',2I10)
      IF(NF.GT.250) GO TO 900
         DO I  = 1, 250
         P (I) = 0D0
         PD(I) = 0D0
         Q (I) = 0D0
         QD(I) = 0D0
         end do
      CALL DLEGF1(P(1),D,X,D,d,M,M,K,0,IER1)
      CALL DLEGF2(Q(1),D,X,D,-M,M,K,0,IER2)
       IF((IER1+IER2).NE.0) WRITE(7,2) M,IER1,IER2
    2 FORMAT(1X,'FUNLEG  M,I,IER1,2,3,4=',6I10)
       DO 3 I=2,NF
       if(i.gt.2) d=P(I-2)
      CALL DLEGF1(P(I),PD(I-1),X,P(I-1),d,M+I-1,M,K,1,IER3)
      CALL DLEGF2(Q(I),QD(I-1),X,Q(I-1),-M+I-1,M,K,1,IER4)
       IF((IER3+IER4).NE.0) WRITE(7,2) M,I,IER1,IER2,IER3,IER4
    3 CONTINUE
  900 RETURN
      END

c************************************************
C FUNLEG1
C
      SUBROUTINE FUNLEG1(M,NF)
      parameter (nalphas=19)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(nalphas,250),PD(nalphas,250)
      COMMON /LEG1/ P,PD
      COMMON /K1/ S, S1, AKSI, AK, kkK, NK, nal
c      COMMON /K2/ Ialpha(nalphas), Jalpha(nalphas)
      COMMON /FF/ ASIN(nalphas),ACOS(nalphas)
      IF(NF.GT.250) WRITE(7,1) NF
    1 FORMAT(1X,'FUNLEG1 NF > 250, NF=',2I10)
      IF(NF.GT.250) GO TO 900
         DO j  = 1, nalphas
         DO I  = 1, 250
         P (j,I) = 0D0
         PD(j,I) = 0D0
         end do
         end do
      do 6 ia=1,nal
      x=acos(ia)
      CALL DLEGF1(P(ia,1),D,X,D,d,M,M,0,0,IER1)
       IF((IER1).NE.0) WRITE(7,2) M,IER1
    2 FORMAT(1X,'FUNLEG1 M,I,IER1,2,3,4=',6I10)
       DO 3 I=2,NF
       if(i.gt.2) d=p(ia,i-2)
      CALL DLEGF1(P(ia,I),PD(ia,I-1),X,P(ia,I-1),d,M+I-1,M,0,1,IER3)
       IF((IER3).NE.0) WRITE(7,2) M,I,IER3
    3 CONTINUE
    6 CONTINUE
  900 RETURN
      END

c************************************************
c GAMKAP1
c
      SUBROUTINE GAMKAP1(M,NN,GAM,AKAP)
      parameter (nalphas=19,nterms=100)
      IMPLICIT REAL*8 (A-H,O-Z)
      complex*16 s, s1
      DIMENSION A(4*nterms), B(4*nterms),
     *          aCOF(nterms,4*nterms), VV1(nterms),
     *          GAM(nterms,nterms), AKAP(nterms,nterms)
      COMMON /EPS3/ EPS3
      COMMON /Kc1/  aCOF, VV1
      COMMON /K1/ S, S1, AKSI, AK, kkK, NK, nal

        M2=2*M
       DO 10 N=1,NN
          NM=N+M-1
       DO 12 I=1,NK
   12 B(I)=aCOF(N,I)
      C = VV1(N)
       DO 11 L=1,NN
          LM=L+M-1
      IF(MOD(IABS(LM-NM),2)) 20,20,9
    9 CONTINUE
       DO 13 I=1,NK
   13 A(I) = aCOF(L,I)
      D = VV1(L)
      IF(MOD(LM-M,2))1,1,2
    1 I1=0
      I2=1
      K=M2
      E=A(1)*(-B(1)*M/(M2+3D0))*2D0
      E1=A(1)*(B(1)/(M2+3D0))*2D0
      GO TO 3
    2 I1=1
      I2=0
      K=M2+1
      E=A(1)*(B(1)*(M+2D0)/(M2+1D0)-B(2)*(M+1D0)*(M2+2D0)/
     *(M2+5D0))*2D0/(M2+3D0)
      E1=A(1)*(B(1)/(M2+1D0)+B(2)*(M2+2D0)/(M2+5D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      E1=E1*F
      DO 5 I=2,NK
      J=2*I-1-I2
      JM=J+M
      JM2=2*JM
      F=1D0
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=A(I)*(B(I-I2)*J*(JM+1D0)/(JM2-1D0)-B(I+I1)*JM*(M2+J+1D0)/
     *(JM2+3D0))*2D0/(JM2+1D0)*F
      G1=A(I)*(B(I+I1)*(J+M2+1D0)/(JM2+3D0)+B(I-I2)*J/(JM2-1D0))
     *   *2D0/(JM2+1D0)*F
c      IF(CDABS(G).LE.1D-140.AND.CDABS(G1).LE.1D-140.and.i.gt.nn)
      IF(DABS(G).LE.eps3.AND.DABS(G1).LE.eps3.and.i.gt.nn)
     *   GO TO 6
      E1=E1+G1
    5 E=E+G
    6 GAM(N,L)=E1/C/D
      AKAP(N,L)=E/C/D
      GO TO 11
   20 GAM(N,L)=(0D0,0D0)
      AKAP(N,L)=(0D0,0D0)
   11  CONTINUE
   10  CONTINUE
      RETURN
      END

c************************************************
C gj11
c
      SUBROUTINE gj11(A,N,NP,B,M,MP)
      implicit complex*16(a-h,o-z)
      PARAMETER (nterms=100, nalphas=19)
      DIMENSION A(NP,NP),B(mP,nP),IPIV(nterms),
     *          INDXR(nterms),INDXC(nterms)
      COMMON /K2/ Ialpha(nalphas),Jalpha(nalphas)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=(0.0d0,0.0d0)
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (cdABS(A(J,K)).GE.cdabs(BIG))THEN
                  BIG=cdABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
           IF(Jalpha(l).EQ.1) GO TO 15
            DUM=B(L,IROW)
            B(L,IROW)=B(L,ICOL)
            B(L,ICOL)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (cdabs(A(ICOL,ICOL)).EQ.(0.0d0,0.0d0))
     *  PAUSE 'Singular matrix.'
        PIVINV=(1.0d0,0.0d0)/A(ICOL,ICOL)
        A(ICOL,ICOL)=(1.0d0,0.0d0)
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
           IF(Jalpha(l).EQ.1) GO TO 17
          B(L,ICOL)=B(L,ICOL)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=(0.0d0,0.0d0)
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
             IF(Jalpha(l).EQ.1) GO TO 19
              B(L,LL)=B(L,LL)-B(L,ICOL)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

c************************************************
C homfunq
C
      SUBROUTINE homfunq(II,M,C1,C2,ksi0,EPS,NE,*)
      parameter (nterms=100, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      REAL*8 ksi0
      COMPLEX*16 C2, BCOF, VV2, DEL21
      DIMENSION
     *          RLC1(nterms), RLC2(nterms),
     *          ADC1(4*nterms), bDC1(4*nterms),
     *          rDC2(4*nterms),

     *          ACOF(nterms,4*nterms), VV1(nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms),

     *          BD(nalphas,2*nterms), B0(nterms),
     *          Ialpha(nalphas),Jalpha(nalphas),
     *          ASIN(nalphas), ACOS(nalphas)

  202 FORMAT(1X,'homFUNq NN>nterms',5X,'NN=',I5)
  210 FORMAT(1X,'L=',I4,5X,'IER1,2,3,4,5=',5I5)
  212 FORMAT(1X,'I,L=',2I5,5X,'IER=',I5)
      COMMON /Kc1/  ACOF, VV1
      COMMON /Kc2/  BCOF, VV2
      COMMON /Kfun/ ar1_1(2*nterms), r3_1(2*nterms),
     *              aar1_1(nterms), rr3_1(nterms),
     *              rr1_2(nterms)

      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /K2/  Ialpha, Jalpha
      COMMON /K10/ B0
      COMMON /K11/ BD
      COMMON /FF/  ASIN, ACOS
      COMMON /INT/ DEL21(nterms,nterms),  GAM11(nterms,nterms),
     *             akap11(nterms,nterms), asig11(nterms,nterms)
      COMMON /EPS1/ EPS1
      W1=1D0/(C1*(ksi0**2-1D0+2*K))
      IF(NE-nterms)  40,40,41
   41 WRITE(7,202) NE
      WRITE(*,202) NE
      RETURN 1
   40 CONTINUE
      RC1=DCMPLX(C1,0D0)

c++++++++
c     Lambda ---> c1
c++++++++
      if(k.eq.0) CALL CDLAMp(RLC1,M,NE,RC1,EPS,IE)
c-------
         if(k.eq.1) then
      if(dreal(rc1).lt.22d0) CALL CDLAMo_s(RLC1,M,NE,RC1,EPS,IE)
      if(dreal(rc1).ge.22d0) CALL CDLAMo_l(RLC1,M,NE,RC1,EPS,IE)
         end if
         if(ie.ne.0) print *, c1, ie
         if(ie.ne.0) return 1

c            do jj = 1,ne
c            write(7,*) jj, rlc1(jj)
c            end do

c-------

c++++++++
c     Lambda ---> c2
c++++++++
      call lambda(K,M,ne,C2,EPS,rlc2,ie)
         if(ie.ne.0) return 1

c-------

       if(ii.eq.1) CALL FUNLEG1(M,2*nk+2)

      if(k.eq.0) then
      IFUN1 = 11
      IF(ksi0.GE.2D0) IFUN1 = 22
      end if

      if(k.eq.1) then
       ifun1 = 33
       if(ksi0.gt.1.4d0) ifun1 = 22
      end if

      inum1=nk
      if(nk.lt.inum1) inum1=nk+10
      if(inum1.gt.230) inum1=230
      IF(ifun1.ne.22)  CALL FUNLEG(m,ksi0,inum1)

c
c-------------- functions: c1,c2/ksi0
c
      DO 1 L=1,NE
      LM=L+M-1
      L1=L+NE

c----------------------------------------------------
c  Functions for `c1'
c----------------------------------------------------
      RL=RLC1(L)

      CALL cdcof3(RL,RDC2,NK,BDC1,NK,2,rvn,1,M,LM,RC1,IER1)

      if(ifun1.eq.11) go to 11
      if(ifun1.eq.22) go to 22
      if(ifun1.eq.33) go to 33

c----------------------------------------------------
c    LEG/LEG    ifun1 = 11
c----------------------------------------------------
   11  CONTINUE
      CALL CDRF12(R1,R2,ar3,ar4,1,M,LM,RC1,RDC2,NK,
     *            BDC1,NK,IER2)
* -----
      ar1_1(L)  = R1
      ar1_1(L1) = R2
      r3_1(l) = r1 + s1 * ar3
      r3_1(l1) = r2 + s1 * ar4
* -----
      W = r1 * ar4 - r2 * ar3
      W2a = dabs(W - W1)
      if(W2a.gt.eps1)  write(*,99999) ifun1,l,W,W1,W2a

99999 format(1x,'ifun=',i2,2x,'L=',i3,2x,'W=',d13.6,2x,'W1=',d13.6,2x,
     * 'W2=',d13.6)

      go to 5

c----------------------------------------------------
c    BES/BES    ifun1 = 22
c----------------------------------------------------
   22  CONTINUE
      CALL CDRB12(R1,R2,ar3,ar4,1,M,LM,RC1,ksi0,RDC2,NK,IER3)
* -----
      ar1_1(L)  = R1
      ar1_1(L1) = R2
      r3_1(l) = dcmplx(ar1_1(L),ar3)
      r3_1(l1) = dcmplx(ar1_1(L1),ar4)
* -----
      W = r1 * ar4 - r2 * ar3
      W2a = dabs(W - W1)
      if(W2a.gt.eps1)  write(*,99999) ifun1,l,W,W1,W2a

      go to 5

c----------------------------------------------------
c    DIFF/DIFF   ifun1 = 33
c----------------------------------------------------
   33  CONTINUE

      CALL DRSF212(aR1,aR2,ar3,ar4,2,m,LM,rc1,ksi0,rl,
     *             RDC2,NK,BDC1,m+1,IER4)
* -----
      ar1_1(L)  = aR1
      ar1_1(L1) = aR2
      r3_1(l) = dcmplx(ar1_1(L),ar3)
      r3_1(l1) = dcmplx(ar1_1(L1),ar4)
* -----
      W = ar1 * ar4 - ar2 * ar3
      W2a = dabs(W - W1)
      if(W2a.gt.eps1)  write(*,99999) ifun1,l,W,W1,W2a
c        write(*,99999) ifun1,l,W,W1,W2a

    5   continue

      IF((IER1+IER2+IER3+IER4).NE.0) WRITE(7,210) L,
     *    IER1,IER2,IER3,IER4

       DO I = 1,NK
       ADC1(I) = RDC2(I)
       ACOF(L,i) = ADC1(i)
       end do

      VC1L2 = RVN
      VV1(L) = DSQRT(VC1L2)

      IF(II.EQ.0) then
         B0(L) = COFF(L,NK,ADC1)
         GO TO 88
      end if

      DO 8 I = 1, NAL
      IF(Jalpha(I).EQ.1) GO TO 8
      CALL DASWF1(BD(I,L),BD(I,L1),1,M,LM,ACOS(I),i,ADC1,NK,IER)
      IF(IER.NE.0) WRITE(7,212) I,L,IER
      BD(I,L) = BD(I,L) / VC1L2
    8 CONTINUE

   88 CONTINUE

      ar1_1(L) = ar1_1(L) * VV1(L)
      ar1_1(L1) = ar1_1(L1) * VV1(L)
      r3_1(L) = r3_1(L) * VV1(L)
      r3_1(L1) = r3_1(L1) * VV1(L)
      aar1_1(L) = ar1_1(L1) / ar1_1(L)
      rr3_1(L) = r3_1(L1) / r3_1(L)

c----------------------------------------------------
c  Functions for `c2'/ `ksi0'
c----------------------------------------------------
      RP=RLC2(L)

      CALL cdcof3(RP,RDC2,NK,BDC1,NK,2,RVN2,1,M,LM,C2,IER5)

      if(ifun1.eq.11) ifun2 = 1
      if(ifun1.eq.22.or.ksi0.gt.1.6d0) ifun2 = 2
       go to (731, 732), ifun2
  731 CALL CDRF12(r1_2l,r1_2l1,v,v,0,M,LM,C2,RDC2,NK,
     *             BDC1,NK,IER4)

      GO TO 7

  732 CALL CDRB12(r1_2l,r1_2l1,v,v,0,M,LM,C2,ksi0,
     *             RDC2,NK,IER3)
    7 CONTINUE

c=====
c      rW1 = 1D0/(C2*(ksi0**2-1D0+2*K))
c      rW = r1_2(l) * r4 - r1_2(l1) * r3
c      rW2 = cdabs(rW - rW1)
c      if(W2a.gt.eps1) write(*,98) ifun2,l,rW,rW1,rW2

c   98 format(1x,'c2/ksi0 L=',2i2,1x,3(2d10.3,2x))
c=====

      IF((IER3+IER4+IER5).NE.0) WRITE(7,210) L,
     *    IER1,IER2,IER3,IER4,IER5
c.......................
           if(cdabs(rvn2).lt.1d20) vv2(l) = cdsqrt(rvn2)
           if(cdabs(rvn2).gt.1d20.and.cdabs(rvn2).lt.1d80)
     *        vv2(l) = cdsqrt(rvn2/1d40)*1d20
           if(cdabs(rvn2).gt.1d80) vv2(l) = cdsqrt(rvn2/1d80)*1d40
c.......................

c      r1_2(L ) = r1_2(L ) * VV2(L)
c      r1_2(L1) = r1_2(L1) * VV2(L)
      rr1_2(L) = r1_2L1 / r1_2L

      DO N = 1, NK
      BCOF(L,N)=rDC2(N)
      end do

    1 CONTINUE

        CALL DELta21(M,NE,DEL21)
        CALL GAMKAP1(M,NE,GAM11,akap11)
        CALL asigma11(M,NE,asig11)

* ----------------------------------
      RETURN
      END

c************************************************
C INT3
C
      SUBROUTINE INT3(M,NN,KAP,OMG,TU)
      parameter (nterms=100)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KAP
      complex*16 s, s1
      DIMENSION A(4*nterms),B(4*nterms),ACOF(nterms,4*nterms),
     *          VV(nterms), OMG(nterms,nterms),TU(nterms,nterms),
     *          KAP(nterms,nterms)
       COMMON /kc1/ ACOF,VV
       COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      DO 6 L=1,NN
      LM=L+M-1
       V1=VV(L)
      DO 12 I=1,NK
   12 A(I)=ACOF(L,I)
      DO 5 N =1,NN
      NM=N+M-1
       V2=VV(N)
      DO 13 I=1,NK
   13 B(I)=ACOF(N ,I)
          OMG(L,N)=OMEGA (M,LM,NM,NK,A,B,V1,V2)
          TU (L,N)=TAU   (M,LM,NM,NK,A,B,V1,V2)
          KAP(L,N)=AKAPPA(M,LM,NM,NK,A,B,V1,V2)
    5 CONTINUE
    6 CONTINUE
       RETURN
        END

c************************************************
C INVMAt
C
      SUBROUTINE INVMAT(M,NN,KSI0,AINV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KSI0
      COMPLEX*16 S, S1
      DIMENSION AINV(200,NN),FACT(170)
      COMMON /FACT/ FACT
      COMMON /K1/  S, S1, AKSI, AK, k, nk, nal
      M2=2*M
      AM=(FACT(M2+1)/(FACT(M+1)*2D0**M))**2
      B=0D0
      IM=(-1)**(M-1)
      DO 1 I=1,M
      I1=I-1
      JM=M-I1
      IF(K.EQ.1) IM=(-1)**JM
      CMK=FACT(M+1)/FACT(I1+1)/FACT(JM+1)
      D=0D0
      DO 2 J=1,JM
      J1=J-1
      CMJ=FACT(JM)/FACT(J1+1)/FACT(JM-J1)
      IF(K.EQ.0) CK=(-1)**J1*2D0
      IF(K.EQ.1) CK=1D0
      D=D+CK*KSI0**(2D0*J1)*CMJ/(2D0*(JM-J1)-1D0)
    2 CONTINUE
      B=B+IM*CMK*AKSI**I1*D
    1 CONTINUE
      IF(K.EQ.0) AINV(1,1)=AM*(B-IM*AKSI**M/KSI0
     *                     *DLOG((KSI0+1D0)/(KSI0-1D0)))
      IF(K.EQ.1) AINV(1,1)=AM*(2D0*B+AKSI**M*2D0/KSI0*DATAN(1D0/KSI0))
      AINV(2,2)=(M+0.5D0)*(KSI0**2*AINV(1,1)*(M+0.5D0)-FACT(M2+1))*AK*
     *          4D0
      DO 3 N=1,NN
      N1=N+1
      IF(N.GT.2) N1=N
      DO 4 L=N1,NN
      LM=L+M-1
      IF(MOD(IABS(L-N),2)) 5,5,6
    6 AINV(N,L)=S
      GO TO 7
    5 CONTINUE
      IL=LM-2
      I2=L-4
      L3=2*IL+3
      IF(I2.LT.1) B2=0D0
      IF(I2.LT.1) GO TO 8
      B2=-L3*(IL+M)*(IL+M-1D0)/((L3-4D0)*(IL-M+1D0)*
     *   (IL-M+2D0))*AINV(N,I2)
         ggksi=1d+303/10d0**(ksi0/3d0)
      IF(DABS(B2).GT.ggksi) AINV(N,L)=AINV(N,I2)
      IF(DABS(B2).GT.ggksi) GO TO 7
    8 CONTINUE
      AINV(N,L)=AK*(L3-2D0)*L3/((IL-M+1D0)*(IL-M+2D0))*
     *          (KSI0**2-AK*(2D0*(IL-M)*(IL+M+1D0)+M2-1D0)/
     *          (L3*(2D0*IL-1D0)))*AINV(N,L-2)+B2
      IF((IL+M+1).GT.170) GO TO 10
      IF(N.EQ.L-2) AINV(N,L)=AINV(N,L)-2D0*L3/FACT(IL-M+3)*
     *                     FACT(IL+M+1)*AK
      GO TO 11
   10 CONTINUE
       ML=M2-1
       F=1D0
       DO 12 I=1,ML
   12  F=F*(L+I-1)
      IF(N.EQ.L-2) AINV(N,L)=AINV(N,L)-2D0*L3*F
   11 CONTINUE
      IF(N.EQ.L) GO TO 4
    7 AINV(L,N)=AINV(N,L)
    4 CONTINUE
    3 CONTINUE
      RETURN
      end

c************************************************
C Lambda
C
      SUBROUTINE Lambda(K,M,ne,C2,EPS,rlc2,ie)
      parameter (nterms=100)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      REAL*8 stepr,stepi,imc2,stepi0
      COMPLEX*16 C2,cc2a
      DIMENSION RLC2(nterms),RLcC2(nterms)
      cC2 = dreal(c2)
      imC2 = dimag(c2)

c++++
          if(k.eq.0.and.imc2.lt.0.5d0) then
          CALL CDLAMp(RLC2,M,NE,C2,EPS,IE)
          if(ie.ne.0)      print *,c2, ie
          if(ie.ne.0) go to 484
          go to 443
          end if
c&-------
                if(k.eq.1) then
      if(imc2.lt.0.11d0.and.cc2.le.22d0.or.
     *   imc2.lt.0.3d0.and.cc2.gt.22d0) then
      if(cc2.lt.22d0) CALL CDLAMo_s(rlc2,M,NE,C2,EPS,IE)
      if(cc2.ge.22d0) CALL CDLAMo_l(rlc2,M,NE,C2,EPS,IE)
         if(ie.ne.0)  go to 484

c            do jj = 1,ne
c            write(7,*) jj, rlc2(jj)
c            end do

         go to 443
      end if
                end if
c&-------

 484        continue
c++++++++
      if(k.eq.0)   stepi0 = 0.032d0
      if(k.eq.1)   stepi0 = 0.091d0
      if(k.eq.1.and.cc2.gt.22d0)   stepi0 = 0.018d0
      if(k.eq.1.and.cc2.gt.32d0)   stepi0 = 0.0082d0
         if(stepi0.gt.imc2) then
         stepi0 = imc2 / 3d0
         end if

c-------------
             stepr = 0d0
             stepi = stepi0
                  ac2 = 0.091d0
       if(k.eq.1) ac2 = 0d0
       cc2a = DCMPLX(cc2,ac2)

c###
             if(k.eq.0) then
1443    continue
             CALL CDLAMp(RLcC2,M,NE,cc2a,EPS,IE)

         if(ie.ne.0) then
         cc2a = DCMPLX(cc2,0d0)
         go to 1443
         end if
              end if

c###
             if(k.eq.1) then

         if(dreal(cc2a).le.22d0) CALL CDLAMo_s(rlcc2,M,NE,cc2a,EPS,IE)
         if(dreal(cc2a).gt.22d0) CALL CDLAMo_l(rlcc2,M,NE,cc2a,EPS,IE)
             end if
c-------------

1444    continue
              cc2a = cc2a + dcmplx(stepr,stepi)
              aimc2 = imc2 - dimag(cc2a)
              if(dabs(dreal(cc2a)-cc2).lt.stepr.or.aimc2.lt.0d0)
     *        cc2a = c2

              epss = eps
              if(k.eq.0) epss = 1d-10
              CALL CDLAMn(rlcc2,rlc2,M,NE,cc2a,EPSs,IE)
              if(ie.ne.0) print *,cc2a, ie
ccc                          print *,cc2a, ie

         if(ie.ne.0) then
         cc2a = cc2a - 2d0 * dcmplx(stepr,stepi)

         if(dimag(cc2a).lt.0d0) cc2a = dcmplx(dreal(cc2a),0d0)
         stepi = stepi/2.1d0
ccc         if(stepi.lt.1d-6) stop 000
         if(stepi.lt.1d-6) return

         go to 1444
         end if

              if(cdabs(c2)-cdabs(cc2a).lt.1d-5) go to 1445

           do jj=1,ne
             rlcc2(jj)=rlc2(jj)
           end do

              go to 1444
 1445    continue
c++++++++
c++++++++
  443 CONTINUE
       return
       end

c************************************************
C sphTe
C   Te mode, homogeneous spheroid
C................................................
C   Matrix, Eq.(2.33)
C ⁄ƒƒƒ                                                                            ƒƒƒƒƒø
C ≥ ksi1*rw1             ksi1*rw1             ksi1*rw1             ksi1*rw1            ≥
C ≥ G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3 ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ¿ƒƒƒƒ                                                                           ƒƒƒƒƒŸ
C................................................

      SUBROUTINE sphte(M,C1,ksi0,nn)
      parameter (nterms=100, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16(R-S)
      REAL*8 ksi0
      COMPLEX*16 DEL21,
     *           AAA(nterms,nterms),BBB(nalphas,2*nterms),
     *           CCC(nalphas,nterms),
     *           aa3(nterms,nterms), bb3(nterms,nterms),
     *           ag3(nterms,nterms), bg3(nterms,nterms)
      DIMENSION
     *          Ialpha(nalphas),Jalpha(nalphas),
     *          q1(nterms,nterms), Ainv(200,200)

      COMMON /Kc1/  ACOF(nterms,4*nterms), VV1(nterms)
      COMMON /Kfun/ ar1_1(2*nterms), r3_1(2*nterms),
     *              aar1_1(nterms), rr3_1(nterms),
     *              rr1_2(nterms)

      common /bbb/ Ra00(nalphas,nterms), Ra0(nalphas,nterms)

      common /ires/ RA3(nterms,nterms), Rg3(nterms,nterms), Rf1

      COMMON /K1/  S, S1, AKSI, AK, k, nk, na
      COMMON /K2/  Ialpha, Jalpha
      COMMON /INT/ DEL21(nterms,nterms),  GAM11(nterms,nterms),
     *             akap11(nterms,nterms), asig11(nterms,nterms)
      COMMON /RES/ Rz0(nalphas,nterms),
     *             Rz1(nalphas,nterms),Rz2(nalphas,nterms)

      NDIM = nterms
      NP2=NN/2

         IF(M.GT.1) GO TO 20
*
*  Axisymmeric part
*
c**
c** matrix Eq.(2.19)
c**
c--------
      DO  N = 1, NN
        DO  L = 1, NN
         AAA(N ,L ) = ra3(n,l)
       if(n.eq.l) AAA(N,L) = AAA(N,L) - rr3_1(l)
        end do
      end do
c--------

       DO 23 IA = 1, NA
        DO N = 1, NN
          S2 = S
           DO i = 1, NN
            S2 = S2 - ra3(n,I) * RA00(ia,i)
            if(n.eq.i) S2 = S2 + aar1_1(i) * RA00(ia,i)
           end do
          ccc(ia,N) = S2
        end do
   23  CONTINUE

          CALL gj11(AAA,nn,NDIM,ccc,na,nalphas)

* calculation of Z0 (Eq.(2.16))

      do ia = 1, na
      do i = 1, nn
         rz0(ia,i) = ccc(ia,i) / R3_1(i)
      end do
      end do

c**
*  Non-axisymmeric part
c**
   20     CONTINUE

c--------
c  Inverse matrix
c--------
***      inmm=160+idint(c1)+idint(aksi)+m+10
      inmm=nk
      if (inmm.gt.200) inmm = 200

* Q1 -- q(c_1,ksi0)
      CALL INVMAT(M,inmm,KSI0,Ainv)
      NK1=MIN0(NK,inmm/2)

      DO  N = 1, nn
        DO  L = 1, nn
          q1(N,L) = 0d0
          IF(MOD(IABS(L-N),2)) 4,4,3
    4     CONTINUE
          IF(MOD(N,2)) 6,7,6
    6     JN=1
          GO TO 88
    7     JN=0
   88     CONTINUE
            DO I = 1, NK1
              aF = 0d0
               DO 90 I1 = 1, NK1
   90          aF = aF + aCOF(L,I1) * ainv(2*I-JN,2*I1-JN)
            q1(N,L) = q1(N,L) + aCOF(N,I) * aF
            end do
          q1(N,L) = q1(N,L) / VV1(N) / VV1(L)

    3   continue
        end do
      end do

c--------
      CALL Cmlm2(NN,rg3,gam11,ag3)

          do i=1,nn
          do j=1,nn
          bi = 0d0
          if(i.eq.j) bi = 1d0
          bg3(i,j) = bi + ksi0 * ra3(i,j)
          end do
          end do

      CALL Cmlm2(NN,bg3,gam11,bg3)

c--------
      do i = 1, nn
       do j = 1, nn
       bi = 0d0
       if(i.eq.j) bi = 1d0
       aa3(i,j) = ksi0 * (bi + ksi0 * ra3(i,j)) - ak * ag3(i,j)
       bb3(i,j) = ksi0 * rg3(i,j) - bg3(i,j)
       end do
      end do

      CALL Cmlm2(NN,aa3,q1,aa3)
      CALL Cmlm2(NN,bb3,q1,bb3)

      do i = 1, nn
       do j = 1, nn
        bb3(i,j) = bb3(i,j) +  akap11(i,j) / aksi
       end do
      end do

      CALL Cmlm2(NN,aa3,gam11,ag3)
      CALL Cmlm2(NN,bb3,gam11,bg3)

      DO 11 IA = 1, NA
      IF(Jalpha(IA).EQ.1) GO TO 11
      DO N = 1, nn
       S2=S
       S3=S
         DO I = 1, nn
          S2 = S2 - ksi0 * ((rf1 + 1d0) * Ra3(n,I)
     *            - rf1 * aa3(n,i)) * RA0(ia,i)
            if(n.eq.i) S2 = S2 +
     *                 ksi0 * (rf1 + 1d0) * aar1_1(i) * RA0(ia,i)
          S3 = S3 - ((rf1 + 1d0) * (rg3(n,i) - gam11(n,i) * aar1_1(i))
     *            - ksi0 * rf1 * bb3(n,i)) * RA0(ia,i)
         end do
       BBB(ia,N) = S2
       BBB(ia,N+nn) = S3
      end do
   11 CONTINUE

c**
c** matrix Eq.(2.33)
c**

      DO 5 IV=1,2
      IG=IV-1
      DO 8 N=1,NP2
      N1=N  + NP2
      ND=2*(N-1)+IV
      NG=2*N-IG
      DO 9 L=1,NP2
      L1=L + NP2
      LD=2*(L-1)+IV
      LG=2*L-IG
c--------
      AAA(N ,L ) = ksi0 * ((rf1 + 1d0) * ra3(nd,ld) - rf1 * aa3(nd,ld))
       if(n.eq.l) AAA(N,L) = AAA(N,L) - ksi0 * (rf1 + 1d0) * rr3_1(ld)
      AAA(N ,L1) = ak * ((rf1 + 1d0) * (rg3(nd,lg)
     *             - gam11(nd,lg) * rr3_1(lg)) - rf1 * ag3(nd,lg))
c--------
      AAA(N1,L ) = (rf1 + 1d0) * (rg3(ng,ld) - gam11(ng,ld) *
     *             rr3_1(ld)) - ksi0 * rf1 * bb3(ng,ld)
      AAA(N1,L1) = ksi0 * (rf1 + 1d0) * ra3(ng,lg)
     *             - ak * rf1 * bg3(ng,lg)
       if(n.eq.l) AAA(N1,L1) = AAA(N1,L1) - ksi0
     *                         * (rf1 + 1d0) * rr3_1(lg)
c--------
    9 CONTINUE
    8 CONTINUE


      DO 111 IA = 1, NA
      IF(Jalpha(IA).EQ.1) GO TO 111
      DO 13 N = 1, NP2
        ND = 2*(N-1)+IV
        NG = 2*N-IG
        CCC(ia,N    )   = BBB(ia,ND)
        CCC(ia,N+NP2)   = BBB(ia,NG+nn)
   13 continue
  111 CONTINUE
c---
          CALL gj11(AAA,nn,NDIM,ccC,na,nalphas)
c---

* calculations of Z1 / Z2

      do ia = 1, na
       IF(Jalpha(IA).EQ.1) GO TO 112
       do i = 1, np2
        J1=2*(I-1)+IV
        J2=2*I-IG
        rz1(ia,j1) = ccc(ia,i) / R3_1(j1)
        rz2(ia,j2) = ccc(ia,i+np2) / R3_1(j2) / c1
       end do
 112   continue
      end do
    5 CONTINUE

      RETURN
      END

c************************************************
C sphTM
C   TM mode, homogeneous spheroid
C................................................
C   Matrix, Eq.(2.34)
C ⁄ƒƒƒ                                                                            ƒƒƒƒƒø
C ≥ ksi1*rw1             ksi1*rw1             ksi1*rw1             ksi1*rw1            ≥
C ≥ G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3 ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ¿ƒƒƒƒ                                                                           ƒƒƒƒƒŸ
C................................................

      SUBROUTINE sphtm(II,M,RI,C1,ksi0,ne,nn)
      parameter (nterms=100, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16(R-S)
      REAL*8 ksi0
      COMPLEX*16 DEL21,
     *           AAA(nterms,nterms),BBB(nalphas,2*nterms),
     *           CCC(nalphas,nterms)
      DIMENSION
     *          BD(nalphas,2*nterms), B0(nterms),
     *          Ialpha(nalphas),Jalpha(nalphas),
     *          ASIN(nalphas), ACOS(nalphas)

      COMMON /Kc1/  ACOF(nterms,4*nterms), VV1(nterms)
      COMMON /Kfun/ ar1_1(2*nterms), r3_1(2*nterms),
     *              aar1_1(nterms), rr3_1(nterms),
     *              rr1_2(nterms)

      common /bbb/  Ra00(nalphas,nterms), Ra0(nalphas,nterms)

      common /ires/ RA3(nterms,nterms), Rg3(nterms,nterms), Rf1

      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /K2/  Ialpha, Jalpha
      COMMON /K10/ B0
      COMMON /K11/ BD
      COMMON /FF/  ASIN, ACOS
      COMMON /INT/ DEL21(nterms,nterms),  GAM11(nterms,nterms),
     *             akap11(nterms,nterms), asig11(nterms,nterms)
      COMMON /RES/ Rz0(nalphas,nterms),
     *             Rz1(nalphas,nterms),Rz2(nalphas,nterms)

      NDIM = nterms
      NP2=NN/2
      rf1 = ri**2 - 1d0

c------------
c ra3 --> ~R_3 (after Eq.(2.20), p.26)
c      CALL Crr(Ne,del21,rr1_2,RA3)

* multiplication of complex matrix (a1), vector (a2) and a1^*

      DO I = 1, Ne
        DO J = 1, Ne
c         a4(i,j) = a1(j,i) * a2(j)
         rg3(i,j) = del21(j,i) * rr1_2(j)
        end do
      end do

      DO I = 1, Ne
        DO J = 1, Ne
         s2 = s
          DO  Kk = 1, Ne
c   2     B = B + A4(I,K) * A1(K,j)
          s2 = s2 + rg3(I,Kk) * del21(Kk,j)
          end do
c        A3(I,J) = B
         rA3(I,J) = s2
        end do
      end do
c------------

         IF(M.GT.1.OR.II.EQ.0) GO TO 20
*
*  Axisymmeric part
*
c**
c** matrix Eq.(2.20)
c**
      DO  N = 1, NN
      n1 = n + nn
        DO  L = 1, NN
        l1 = l + nn
c--------
         AAA(N ,L ) = - aksi * ra3(n,l)
       if(n.eq.l) AAA(N,L) = AAA(N,L) +
     *            aksi * (rf1 + 1d0) * rr3_1(l) + rf1 * ksi0
c--------
        end do
      end do

       DO 23 IA = 1, NAL
         DO N = 1, Ne
           RA00(ia,n) = - 2D0 * S1**N * BD(IA,N) * ar1_1(n)
         end do

        DO N = 1, NN
          S2 = S
           DO i = 1, NN
            S2 = S2 + aksi * ra3(n,I) * RA00(ia,i)
            if(n.eq.i) S2 = S2 -
     *      (ksi0 * rf1 + aksi * (rf1 + 1d0) * aar1_1(i)) * RA00(ia,i)
           end do
          ccc(ia,N) = S2
        end do
   23  CONTINUE

          CALL gj11(AAA,nn,NDIM,ccc,nal,nalphas)

* calculation of Z0 (Eq.(2.16))

      do ia = 1, nal
      do i = 1, nn
         rz0(ia,i) = ccc(ia,i) / R3_1(i)
      end do
      end do

c**
*  Non-axisymmeric part
c**
   20     CONTINUE

      CALL Cmlm1(Ne,gam11,ra3,rg3)

      sf1 = rf1 * ksi0 / aksi

      IF(II.EQ.0) then
         DO  L = 1, Nn
            RA0(1,L) = 2D0*S1**(L-1)*B0(L)*ar1_1(L)/VV1(L)**2
         end do
      end if

      IF(II.nE.0) then
         DO IA = 1, NAl
           IF(Jalpha(IA).EQ.1) GO TO 15
            DO N = 1, Ne
              RA0(ia,n) = 4D0*S1**(N+M-2)*BD(IA,N)*ar1_1(N)/ASIN(IA)
            end do
   15      CONTINUE
         end do
      end if

      NA=NAL
      IF(II.EQ.0) na = 1
      DO 11 IA = 1, NA
      IF(Jalpha(IA).EQ.1) GO TO 11
      DO N = 1, Ne
       S2=S
       S3=S
         DO I = 1, Ne
          S2 = S2 - ksi0 * Ra3(n,I) * RA0(ia,i)
            if(n.eq.i) S2 = S2 + ksi0 * aar1_1(i) * RA0(ia,i)
          S3 = S3 - (rg3(n,i) - gam11(n,i) * (rf1 + 1d0) * aar1_1(i)
     *            - sf1 * akap11(n,i)) * RA0(ia,i)
         end do
       BBB(ia,N) = S2
       BBB(ia,N+Ne) = S3
      end do
   11 CONTINUE

c**
c** matrix Eq.(2.34)
c**

      DO 5 IV=1,2
      IG=IV-1
      DO 8 N=1,NP2
      N1=N  + NP2
      ND=2*(N-1)+IV
      NG=2*N-IG
      DO 9 L=1,NP2
      L1=L + NP2
      LD=2*(L-1)+IV
      LG=2*L-IG
c--------
      AAA(N ,L ) = ksi0 * ra3(nd,ld)
       if(n.eq.l) AAA(N,L) = AAA(N,L) - ksi0 * rr3_1(ld)
      AAA(N ,L1) = ak * (rg3(nd,lg) - gam11(nd,lg) * rr3_1(lg))
c--------
      AAA(N1,L ) = (rg3(ng,ld) - gam11(ng,ld) *
     *             (rf1 + 1d0) * rr3_1(ld)) - sf1 * akap11(ng,ld)
      AAA(N1,L1) = ksi0 * ra3(ng,lg) - ak * rf1 / aksi * asig11(ng,lg)
c       if(n.eq.l) AAA(N1,L1) = AAA(N1,L1) - ksi0 *
c     *                         (rf1 + 1d0) * rr3_1(l) - rf1
       if(n.eq.l) AAA(N1,L1) = AAA(N1,L1) + 1d0 -
     *                         (rf1 + 1d0) * (1d0 + ksi0 * rr3_1(lg))
c--------
    9 CONTINUE
    8 CONTINUE


      DO 111 IA = 1, NA
      IF(Jalpha(IA).EQ.1) GO TO 111
      DO 13 N = 1, NP2
        ND = 2*(N-1)+IV
        NG = 2*N-IG
        CCC(ia,N    )   = BBB(ia,ND)
        CCC(ia,N+NP2)   = BBB(ia,NG+Ne)
   13 continue
  111 CONTINUE
c---
          CALL gj11(AAA,nn,NDIM,ccC,na,nalphas)
c---

* calculations of Z1 / Z2

      do ia = 1, na
       IF(Jalpha(IA).EQ.1) GO TO 112
       do i = 1, np2
        J1=2*(I-1)+IV
        J2=2*I-IG
        rz1(ia,j1) = ccc(ia,i) / R3_1(j1)
        rz2(ia,j2) = ccc(ia,i+np2) / R3_1(j2) / c1
       end do
 112   continue
      end do
    5 CONTINUE

      RETURN
      END
