      subroutine zlev(z,t,q,zb,sp,ie,je,ke,ak,bk)
c     ===========================================

c     argument declaration
      integer  ie,je,ke,its,iqs
      real     z(ie,je,ke),t(ie,je,ke),getps
      real     q(ie,je,ke),zb(ie,je),sp(ie,je)
      real     ak(ke),bk(ke)

c     variable declaration
      integer  i,j,k
      real     pu, po, zo, zu, tv
      real     rdg,tzero,eps
      data     rdg,tzero,eps /29.271, 273.15, 0.6078/

c     t and q are staggered (stag(3)=-0.5; its=int(2*stag(3))
      its=-1
      iqs=-1

c     computation of height of main levels
      do i=1,ie
        do j=1,je
          zu = zb(i,j)
          pu = sp(i,j)
          tv = (t(i,j,1)+tzero)*(1.+eps*q(i,j,1))
          po = ak(1)+bk(1)*sp(i,j)
          if (po.le.0) po=0.5*pu
          zo = zu + rdg*tv*alog(pu/po)
          z(i,j,1) = zo
          pu = po
          zu = zo
          do k=2,ke
            tv = (0.5*(t(i,j,k)+t(i,j,k-1-its))+tzero)*
     >           (1.+eps*0.5*(q(i,j,k)+q(i,j,k-1-iqs)))
            po = ak(k)+bk(k)*sp(i,j)
            if (po.le.0) po=0.5*pu
            zo = zu + rdg*tv*alog(pu/po)
            z(i,j,k) = zo
            pu = po
            zu = zo
          enddo
        enddo
      enddo
      end

      subroutine flightlev(fl,p,ie,je,ke)
c     ===================================

c     Convert an array of pressures (hPa) into an array of altitudes (m)
c     or vice versa using the ICAO standard atmosphere.
c
c     See http://www.pdas.com/coesa.htm
c
c     The 7 layers of the US standard atmosphere are:
c
c     h1   h2     dT/dh    h1,h2 geopotential alt in km
c      0   11     -6.5     dT/dh in K/km
c     11   20      0.0
c     20   32      1.0
c     32   47      2.8
c     47   51      0.0
c     51   71     -2.8
c     71   84.852 -2.0

c     based upon IDL routine from Dominik Brunner
c     converted to f90: Jan 2002  Heini Wernli

      implicit none

c     argument declaration
      integer   ie,je,ke
      real      fl(ie,je,ke),p(ie,je,ke)

c     variable declaration
      integer   i,j,k,n,il
      real      limits(0:7),lrs(7)
      integer   iszero(7)
      real      pb(0:7),tb(0:7)
      real      g,r,gmr,feet
      data      g,r,feet /9.80665, 287.053, 0.3048/

c     define layer boundaries in m, lapsre rates (9 means 0) and isothermal flag
      limits(0)=0.
      limits(1)=11000.
      limits(2)=20000.
      limits(3)=32000.
      limits(4)=47000.
      limits(5)=51000.
      limits(6)=71000.
      limits(7)=84852.

      lrs(0)=-6.5/1000.
      lrs(1)=9.0/1000.
      lrs(2)=1.0/1000.
      lrs(3)=2.8/1000.
      lrs(4)=9.0/1000.
      lrs(5)=-2.8/1000.
      lrs(6)=-2.0/1000.

      iszero(0)=0
      iszero(1)=1
      iszero(2)=0
      iszero(3)=0
      iszero(4)=1
      iszero(5)=0
      iszero(6)=0

      gmr=g/r         ! Hydrostatic constant

c     calculate pressures at layer boundaries
      pB(0)=1013.25 ! pressure at surface
      TB(0)=288.15  ! Temperature at surface

c     loop over layers and get pressures and temperatures at level tops
      do i=0,6
        TB(i+1)=TB(i)+(1-iszero(i))*lrs(i)*(limits(i+1)-limits(i))
        pB(i+1)=(1-iszero(i))*pB(i)*exp(alog(TB(i)/TB(i+1))*gmr/lrs(i))+
     >          iszero(i)*PB(i)*exp(-gmr*(limits(i+1)-limits(i))/TB(i))
      enddo

c     now calculate which layer each value belongs to
c     and calculate the corresponding altitudes
      do i=1,ie
      do j=1,je
      do k=1,ke
        do n=0,6
          if ((pb(n).ge.p(i,j,k)).and.(pb(n+1).le.p(i,j,k))) then
            il=n
            goto 100
          endif
        enddo
  100   continue
        fl(i,j,k)=iszero(il)*(-ALOG(p(i,j,k)/pB(il))*TB(il)/gmr+
     >            limits(il))+(1-iszero(il))*((TB(il)/((p(i,j,k)/
     >            pB(il))**(lrs(il)/gmr))-TB(il))/lrs(il)+limits(il))
        fl(i,j,k)=fl(i,j,k)/feet/100.
      enddo
      enddo
      enddo

      end

      subroutine pottemp(pt,t,sp,ie,je,ke,ak,bk)
c     ==========================================
 
c     argument declaration
      integer   ie,je,ke
      real      pt(ie,je,ke),t(ie,je,ke),sp(ie,je),
     >     	ak(ke),bk(ke)
 
c     variable declaration
      integer   i,j,k
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
 
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf

c     computation of potential temperature
      do i=1,ie
      do j=1,je
        psrf=sp(i,j)
        do k=1,ke
c     distinction of temperature in K and deg C
          if (t(i,j,k).lt.100.) then
            pt(i,j,k)=(t(i,j,k)+tzero)*( (1000./pp(k))**rdcp )
          else
            pt(i,j,k)=t(i,j,k)*( (1000./pp(k))**rdcp )
          endif
        enddo
      enddo
      enddo
      end

      subroutine tnat(tn,t,p,ie,je,ke)
c     ================================

c     argument declaration
      integer   ie,je,ke
      real      tn(ie,je,ke),t(ie,je,ke),p(ie,je,ke)

c     variable declaration
      integer   i,j,k
      real      tzero
      data      tzero /273.15/

      real      mb_to_mm
      real      a1,a2,a3,a4,a5,c0,c1,c2,d
      real      chiH2O,chiHNO3,pH2O,pHNO3,log_H2O,log_NAT

      mb_to_mm=760./1000.
      a1 = -2.7836
      a2 = -0.00088
      a3 = 89.7674
      a4 = -26242.0
      a5 = 0.021135

      chiH2O  = 5.e-6
      chiHNO3 = 5.e-9

      do i=1,ie
      do j=1,je
      do k=1,ke
        pH2O  = p(i,j,k)*chiH2O
        pHNO3 = p(i,j,k)*chiHNO3
        log_H2O = log(pH2O*mb_to_mm)
        log_NAT = log(pHNO3*mb_to_mm)
        c0 = log_NAT-log_H2O*a1-a3
        c1 = a2*log_H2O+a5
        c2 = a4
        d = c0*c0-4.0*c1*c2
        tn(i,j,k)=t(i,j,k)+tzero-(c0+sqrt(d))/(2.*c1)
      enddo
      enddo
      enddo
      end

      subroutine read_tice(tice_arr)
c     ==============================

      integer	n
      real	tice_arr(100)

      open(10,
     >file='/home/henry/prog/tnat_tice/t_ice_from_10_every_2_hPa')
      do n=1,71
        read(10,*)tice_arr(n)
      enddo
      close(10)
      end

      subroutine tice(ti,t,p,tice_arr,ie,je,ke)
c     =========================================

c     argument declaration
      integer   ie,je,ke
      real      ti(ie,je,ke),t(ie,je,ke),p(ie,je,ke)

c     variable declaration
      integer   i,j,k
      real      tzero
      data      tzero /273.15/

      real	tice_arr(100)
      real      p0,p1,pinc,rind
      integer   ind,np

c     define the lowest p-value and the p-increment for the table
      p0=10.
      pinc=2.
      np=71
      p1=p0+real(np-1)*pinc

      do i=1,ie
      do j=1,je
      do k=1,ke

        if (p(i,j,k).lt.p0) then
          print*,'*** error: table is not prepared for p < p0 ***'
          print*,ie,je,ke,i,j,k,p(i,j,k),p0
          call exit(1)
        endif

        if (p(i,j,k).gt.p1) then
          print*,'*** error: table is not prepared for p > p1 ***'
          call exit(1)
        endif

c       calculate the real index of the given pressure level
        rind=(p(i,j,k)-p0)/pinc+1.
        ind=int(rind)
        rind=rind-real(ind)

c       interpolate tice
        ti(i,j,k)=t(i,j,k)+tzero-
     >           (tice_arr(ind)+rind*(tice_arr(ind+1)-tice_arr(ind)))

      enddo
      enddo
      enddo
      end

      subroutine temp(t,pt,p,ie,je,ke,ak,bk)
c     ======================================
 
c     argument declaration
      integer   ie,je,ke
      real      pt(ie,je,ke),t(ie,je,ke),p(ie,je,ke),
     >     	ak(ke),bk(ke)

c     variable declaration
      integer   i,j,k
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
 
c     computation of potential temperature
      do i=1,ie
         do j=1,je
            do k=1,ke
               t(i,j,k)=pt(i,j,k)*( (p(i,j,k)/1000.)**rdcp )
            enddo
         enddo
      enddo
      end

      subroutine laitpv(lpv,pv,th,ie,je,ke)
c     =====================================
 
c     argument declaration
      integer   ie,je,ke
      real      lpv(ie,je,ke),pv(ie,je,ke),th(ie,je,ke)
 
c     variable declaration
      integer   i,j,k
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
 
c     computation of Lait PV
      do i=1,ie
      do j=1,je
      do k=1,ke
        lpv(i,j,k)=pv(i,j,k)*((th(i,j,k)/420.)**(-9./2.))
      enddo
      enddo
      enddo

      end

      subroutine relhum(rh,q,t,sp,ie,je,ke,ak,bk)
c     ===========================================
 
c     argument declaration
      integer   ie,je,ke
      real      rh(ie,je,ke),t(ie,je,ke),q(ie,je,ke),
     >     sp(ie,je),ak(ke),bk(ke)
 
c     variable declaration
      integer   i,j,k
      real      rdcp,tzero
      real      b1,b2w,b3,b4w,r,rd,gqd,ge
      data      rdcp,tzero /0.286,273.15/
      data      b1,b2w,b3,b4w,r,rd /6.1078, 17.2693882, 273.16, 35.86,
     &     287.05, 461.51/
 
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
 
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke
               ge = b1*exp(b2w*t(i,j,k)/(t(i,j,k)+b3-b4w))
               gqd= r/rd*ge/(pp(k)-(1.-r/rd)*ge)
               rh(i,j,k)=100.*q(i,j,k)/gqd
            enddo
         enddo
      enddo
      end

      subroutine relhum_i(rh,q,p,ie,je,ke,ak,mdv)
c     ===========================================
 
c     argument declaration
      integer   ie,je,ke
      real      rh(ie,je,ke),p(ie,je,ke),q(ie,je,ke),
     >     	ak(ke),mdv
 
c     variable declaration
      integer   i,j,k
      real      p0,kappa
      data      p0,kappa /1000.,0.286/
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
      real      b1,b2w,b3,b4w,r,rd,gqd,ge
      data      b1,b2w,b3,b4w,r,rd /6.1078, 17.2693882, 273.16, 35.86,
     &     287.05, 461.51/
 
      do i=1,ie
      do j=1,je
      do k=1,ke
        if (p(i,j,k).eq.mdv) then
          rh(i,j,k)=mdv
        else
          tt=ak(k)*((p(i,j,k)/p0)**kappa)-tzero
          ge = b1*exp(b2w*tt/(tt+b3-b4w))
          gqd= r/rd*ge/(p(i,j,k)-(1.-r/rd)*ge)
          rh(i,j,k)=100.*q(i,j,k)/gqd
        endif
      enddo
      enddo
      enddo
      end

      subroutine gradth(gth,th,sp,levtyp,cl,ie,je,ke,ak,bk,vmin,vmax)
C     ===============================================================
 
c     argument declaration
      integer   ie,je,ke,levtyp
      real      gth(ie,je,ke),th(ie,je,ke),sp(ie,je),cl(ie,je)
      real      ak(ke),bk(ke),vmin(4),vmax(4)
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dthdx,dthdy
      integer   i,j,k,ind,ind2,stat
 
      allocate(dthdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'gradth: error allocating dthdx'
      allocate(dthdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'gradth: error allocating dthdy'
      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'gradth: error allocating dspdx'
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'gradth: error allocating dspdy'
 
      if (levtyp.eq.1) then
         call ddh2(th,dthdx,cl,'X',ie,je,1,vmin,vmax)
         call ddh2(th,dthdy,cl,'Y',ie,je,1,vmin,vmax)
      else if (levtyp.eq.3) then
         call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
         call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
         call ddh3(th,dthdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
         call ddh3(th,dthdy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      endif
 
      gth=sqrt(dthdx**2.+dthdy**2.)
 
      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dthdx)) DEALLOCATE(dthdx)
      IF (ALLOCATED(dthdy)) DEALLOCATE(dthdy)
 
      end

       subroutine calc_gradpv(gpv,pv,cl,ie,je,ke,ak,bk,vmin,vmax)
C      ===============================================================
  
c      argument declaration
       integer   ie,je,ke,levtyp
       real      gpv(ie,je,ke),pv(ie,je,ke),sp(ie,je),cl(ie,je)
       real      ak(ke),bk(ke),vmin(4),vmax(4)
  
c      variable declaration
       REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dpvdx,dpvdy
       integer   i,j,k,ind,ind2,stat
       
       allocate(dpvdx(ie,je,ke),STAT=stat)
       if (stat.ne.0) print*,'gradpv: error allocating dpvdx'
       allocate(dpvdy(ie,je,ke),STAT=stat)
       if (stat.ne.0) print*,'gradpv: error allocating dpvdy'
  
       call ddh2(pv,dpvdx,cl,'X',ie,je,1,vmin,vmax)
       call ddh2(pv,dpvdy,cl,'Y',ie,je,1,vmin,vmax)
  
       gpv=1e6*sqrt(dpvdx**2.+dpvdy**2.)
  
       IF (ALLOCATED(dpvdx)) DEALLOCATE(dpvdx)
       IF (ALLOCATED(dpvdy)) DEALLOCATE(dpvdy)
  
       end

      subroutine wetbpt(thw,the,sp,ie,je,ke,ak,bk)
c     ============================================
 
c     argument declaration
      integer   ie,je,ke
      real      thw(ie,je,ke),the(ie,je,ke),
     >     sp(ie,je),ak(ke),bk(ke)
 
c     variable declaration
      real      tsa
      integer   i,j,k
 
      do k=1,ke
         do j=1,je
            do i=1,ie
               thw(i,j,k)=tsa(the(i,j,k)-273.15,1000.)+273.15
            enddo
         enddo
      enddo
      end

      real function tsa(os,p)
c     =======================
 
C     This function returns the temperature tsa (celsius) on a
c     saturation adiabat at pressure p (millibars). os is the equivalent
c     potential temperature of the parcel (celsius). sign(a,b) replaces
c     the algebraic sign of a with that of b.
C     b is an empirical constant approximately equal to 0.001 of the
c     latent heat of vaporization for water divided by the specific
c     heat at constant pressure for dry air.
 
      real      a,b,os,p,tq,d,tqk,x,w
      integer   i
 
      data b/2.6518986/
      a=os+273.16
 
C     tq is the first guess for tsa
 
      tq=253.16
 
C     d is an initial value used in the iteration below
 
      d=120.
 
C     Iterate to obtain sufficient accuracy....see table 1, p.8
C     of Stipanuk (1973) for equation used in iteration
      do 1 i=1,12
         tqk=tq-273.16
         d=d/2.
         x=a*exp(-b*w(tqk,p)/tq)-tq*((1000./p)**.286)
         if (abs(x).lt.1e-7) go to 2
         tq=tq+sign(d,x)
 1    continue
 2    tsa=tq-273.16
      end

      real function w(t,p)
c     ====================
 
C     This function returns the mixing ratio (grams of water vapor per
C     kilogram of dry air) given the dew point (celsius) and pressure
C     (millibars). If the temperature is input instead of the
C     dew point, then saturation mixing ratio (same units) is returned.
C     The formula is found in most meteorological texts.
 
      real      t,p,tkel,x,esat
 
      tkel=t+273.16
      x=esat(tkel)              ! our function esat requires t in Kelvin
      w=622.*x/(p-x)
      end
 
 
      real function esat(t)
c     =====================
 
C     This function returns the saturation vapor pressure over water
c     (mb) given the temperature (Kelvin).
C     The algorithm is due to Nordquist, W. S. ,1973: "Numerical
C     Approximations of Selected Meteorological Parameters for Cloud
C     Physics Problems" ECOM-5475, Atmospheric Sciences Laboratory,
c     U. S. Army Electronics Command, White Sands Missile Range,
c     New Mexico 88002.
 
      real p1,p2,c1,t
 
      p1=11.344-0.0303998*t
      p2=3.49149-1302.8844/t
      c1=23.832241-5.02808*log10(t)
      esat=10.**(c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/t)
      end

      subroutine equpot(ap,t,q,sp,ie,je,ke,ak,bk)
c     ===========================================

c     calculate equivalent potential temperature
 
c     argument declaration
      integer  ie,je,ke
      real     ap(ie,je,ke),t(ie,je,ke),sp(ie,je)
      real     q(ie,je,ke),ak(ke),bk(ke)
 
c     variable declaration
      integer  i,j,k
      real     rdcp,tzero
      data     rdcp,tzero /0.286,273.15/
 
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
 
c     computation of potential temperature
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke
               ap(i,j,k) = (t(i,j,k)+tzero)*(1000./pp(k))
     +              **(0.2854*(1.0-0.28*q(i,j,k)))*exp(
     +              (3.376/(2840.0/(3.5*alog(t(i,j,k)+tzero)-alog(
     +              100.*pp(k)*max(1.0E-10,q(i,j,k))/(0.622+0.378*
     +              q(i,j,k)))-0.1998)+55.0)-0.00254)*1.0E3*
     +              max(1.0E-10,q(i,j,k))*(1.0+0.81*q(i,j,k)))
            enddo
         enddo
      enddo
      end

      subroutine diabheat(dhr,t,w,rh,sp,ie,je,ke,ak,bk)
c     =================================================
 
c     argument declaration
      integer   ie,je,ke
      real      dhr(ie,je,ke),t(ie,je,ke),w(ie,je,ke),
     &     rh(ie,je,ke),sp(ie,je),ak(ke),bk(ke)
 
c     variable declaration
      integer   i,j,k
      real      p0,kappa,tzero
      data      p0,kappa,tzero /1000.,0.286,273.15/
      real      blog10,cp,r,lw,eps
      data      blog10,cp,r,lw,eps /.08006,1004.,287.,2.5e+6,0.622/
      real      esat,c,tt
 
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
 
c     computation of diabatic heating rate
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke
               if (rh(i,j,k).lt.80.) then ! only moist air of interest
                  dhr(i,j,k)=0. ! cond. heating rate set to zero
               else if (w(i,j,k).gt.0.) then ! cond. heating only
                                             ! for ascent
                  dhr(i,j,k)=0.
               else
                  tt=t(i,j,k)*((pp(k)/p0)**kappa) ! temp. from pot.temp.
                  c=lw/cp*eps*blog10*esat(tt)/pp(k)
                  dhr(i,j,k)=21600.* ! in units K per 6 hours
     >                 (1.-exp(.2*(80.-rh(i,j,k)))) ! weighting fun.
                                                    ! for 80<RH<100
     >                 *(-c*kappa*t(i,j,k)*w(i,j,k)/(100.*pp(k)))/(1.+c)
               endif
            enddo
         enddo
      enddo
      end

      subroutine diabheat2(dhr,t,w,rh,sp,ie,je,ke,ak,bk)
c     ==================================================
 
c     argument declaration
      integer   ie,je,ke
      real      dhr(ie,je,ke),t(ie,je,ke),w(ie,je,ke),
     &     rh(ie,je,ke),sp(ie,je),ak(ke),bk(ke)
 
c     variable declaration
      integer   i,j,k
      real      p0,kappa,tzero
      data      p0,kappa,tzero /1000.,0.286,273.15/
      real      blog10,cp,r,lw,eps
      data      blog10,cp,r,lw,eps /.08006,1004.,287.,2.5e+6,0.622/
      real      esat,c,tt
 
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
 
c     computation of diabatic heating rate
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke
               if (rh(i,j,k).lt.80.) then ! only moist air of interest
                  dhr(i,j,k)=0. ! cond. heating rate set to zero
               else if (w(i,j,k).gt.0.) then ! cond. heating only
                                             ! for ascent
                  dhr(i,j,k)=0.
               else
                  c=lw/cp*eps*blog10*esat((t(i,j,k)+273.15))/pp(k)
                  tt=(t(i,j,k)+273.15)*((p0/pp(k))**kappa) ! pot.temp from temp
                  dhr(i,j,k)=21600.* ! in units K per 6 hours
     >                 (1.-exp(.2*(80.-rh(i,j,k)))) ! weighting fun.
                                                    ! for 80<RH<100
     >                 *(-c*kappa*tt*w(i,j,k)/(100.*pp(k)))/(1.+c)
               endif
            enddo
         enddo
      enddo
      end

      subroutine diabpvr(dpvr,uu,vv,dhr,sp,cl,f,ie,je,ke,ak,bk,
     >     vmin,vmax)
C     =========================================================
 
c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: dpvr(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),
     >     dhr(ie,je,ke),sp(ie,je),cl(ie,je),f(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke),vmin(4),vmax(4)
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dhrdp,dudp,dvdp,dvdx
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dhrdx,dudy,dhrdy
      integer   k,stat
 
      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dspdx'
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dspdy'
 
      allocate(dhrdp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dhrdp'
      allocate(dudp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dudp'
      allocate(dvdp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dvdp'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dvdx'
      allocate(dhrdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dhrdx'
      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dudy'
      allocate(dhrdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'diabpvr: error allocating dhrdy'
 
      call ddp(dhr,dhrdp,sp,ie,je,ke,ak,bk)
      call ddp(uu,dudp,sp,ie,je,ke,ak,bk)
      call ddp(vv,dvdp,sp,ie,je,ke,ak,bk)
      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(dhr,dhrdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(dhr,dhrdy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(uu,dudy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
 
      do k=1,ke
         dpvr(1:ie,1:je,k)=-1.E6*9.80665*(
     >        -dvdp(1:ie,1:je,k)*dhrdx(1:ie,1:je,k)
     >        +dudp(1:ie,1:je,k)*dhrdy(1:ie,1:je,k)
     >        +(-dudy(1:ie,1:je,k)+dvdx(1:ie,1:je,k)
     >        +f(1:ie,1:je))*dhrdp(1:ie,1:je,k))
      enddo
 
      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dhrdp)) DEALLOCATE(dhrdp)
      IF (ALLOCATED(dudp)) DEALLOCATE(dudp)
      IF (ALLOCATED(dvdp)) DEALLOCATE(dvdp)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dhrdx)) DEALLOCATE(dhrdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dhrdy)) DEALLOCATE(dhrdy)
 
      end

      subroutine vort(var,uu,vv,sp,cl,ie,je,ke,ak,bk,vmin,vmax)
C     =========================================================
C     calculate vorticity VORT=D[V:X]-D[U*COS:Y]/COS
 
c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),
     >     sp(ie,je),cl(ie,je)
      real      ak(ke),bk(ke),vmin(4),vmax(4)
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dvdx,dudy,uucl
      integer   stat
      real	pole
 
      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dspdx(ie,je)'
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dspdy(ie,je)'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dvdx'
      allocate(uucl(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating uucl'
      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dudy'
 
      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      do k=1,ke
         uucl(1:ie,1:je,k)=uu(1:ie,1:je,k)*cl(1:ie,1:je)
      enddo
      call ddh3(uucl,dudy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
         var(1:ie,1:je,k)=1.E4*(dvdx(1:ie,1:je,k)-
     >                          dudy(1:ie,1:je,k)/cl(1:ie,1:je))
*        var(1:ie,1:je,k)=1.E4*(dvdx(1:ie,1:je,k)-
*    >                          dudy(1:ie,1:je,k))
      enddo

c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,je-1,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif

      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,2,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(uucl)) DEALLOCATE(uucl)
 
      end

      subroutine avort(var,uu,vv,sp,cl,f,ie,je,ke,ak,bk,vmin,vmax)
C     ============================================================
C     calculate absolute vorticity AVO=F+D[V:X]-D[U*COS:Y]/COS

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),
     >                    sp(ie,je),cl(ie,je),f(ie,je)
      real      ak(ke),bk(ke),vmin(4),vmax(4)

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dvdx,dudy,uucl
      integer   stat
      real      pole

      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dspdx(ie,je)'
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dspdy(ie,je)'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dvdx'
      allocate(uucl(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating uucl'
      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dudy'

      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      do k=1,ke
         uucl(1:ie,1:je,k)=uu(1:ie,1:je,k)*cl(1:ie,1:je)
      enddo
      call ddh3(uucl,dudy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
         var(1:ie,1:je,k)=1.E4*(f(1:ie,1:je)+dvdx(1:ie,1:je,k)-
     >                          dudy(1:ie,1:je,k)/cl(1:ie,1:je))
      enddo

c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,je-1,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif

      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,2,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(uucl)) DEALLOCATE(uucl)

      end

      subroutine vort_i(var,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
C     ======================================================
C     calculate vorticity VORT=1e4*D[V:X]-D[U*COS:Y]/COS on isentropic
C     levels 
C     with misdat (mdv) treatment !
C     mark liniger 990501
 
c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),
     >     cl(ie,je)
      real      vmin(4),vmax(4),mdv
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dvdx,dudy,uucl
      integer   stat
 
      allocate(dvdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dvdx'
      allocate(uucl(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating uucl'
      allocate(dudy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dudy'

      do k=1,ke

      call ddh2m(vv(1,1,k),dvdx,cl,'X',ie,je,1,vmin,vmax,mdv)
      
      do i=1,ie
         do j=1,je
            if (uu(i,j,k).eq.mdv) then
               uucl(i,j)=mdv
            else
               uucl(i,j)=uu(i,j,k)*cl(i,j)
            endif
         enddo
      enddo


      call ddh2m(uucl,dudy,cl,'Y',ie,je,1,vmin,vmax,mdv)
 
      do i=1,ie
         do j=1,je
            if ((dvdx(i,j).eq.mdv).or.
     >           (dudy(i,j).eq.mdv).or.
     >           (cl(i,j).lt.1e-3)) then
               var(i,j,k)=mdv
            else
               var(i,j,k)=1.E4*(dvdx(i,j)-dudy(i,j)
     >              /cl(i,j))
            endif
         enddo
      enddo

      enddo
 
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(uucl)) DEALLOCATE(uucl)
 
      end

      subroutine curvo(var,uu,vv,sp,cl,ie,je,ke,ak,bk,vmin,vmax)
C     ==========================================================
C     calculate curvature vorticity
C     =u^2*d[v:x]-v^2*d[u*cos:y]/cos-u*v*(d[u:x]-d[v*cos:y]/cos))/(u^2+v^2)

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),
     >     		  sp(ie,je),cl(ie,je)
      real      ak(ke),bk(ke),vmin(4),vmax(4)

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dudx,dvdx,dudy,dvdy,
     >					     uucl,vvcl
      integer   stat

      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dspdx'
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dspdy'
      allocate(dudx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dudx'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dvdx'
      allocate(uucl(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating uucl'
      allocate(vvcl(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating vvcl'
      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dudy'
      allocate(dvdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'vort: error allocating dvdy'

      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(uu,dudx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      do k=1,ke
        uucl(1:ie,1:je,k)=uu(1:ie,1:je,k)*cl(1:ie,1:je)
        vvcl(1:ie,1:je,k)=vv(1:ie,1:je,k)*cl(1:ie,1:je)
      enddo
      call ddh3(uucl,dudy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(vvcl,dvdy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
      do j=1,je
      do i=1,ie
        var(i,j,k)=1.e4*(uu(i,j,k)**2.*dvdx(i,j,k)-
     >                   vv(i,j,k)**2.*dudy(i,j,k)/cl(i,j)-
     >         uu(i,j,k)*vv(i,j,k)*(dudx(i,j,k)-dvdy(i,j,k)/cl(i,j)))/
     >         max(uu(i,j,k)**2.+vv(i,j,k)**2.,1.e-3)
      enddo
      enddo
      enddo

c     set poles values to mdv
      if (vmax(2).eq.90.) then
        do k=1,ke
        do i=1,ie
          var(i,je,k)=-999.98999
        enddo
        enddo
      endif

      if (vmin(2).eq.-90.) then
        do k=1,ke
        do i=1,ie
          var(i,1,k)=-999.98999
        enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dudx)) DEALLOCATE(dudx)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dvdy)) DEALLOCATE(dvdy)
      IF (ALLOCATED(uucl)) DEALLOCATE(uucl)
      IF (ALLOCATED(vvcl)) DEALLOCATE(vvcl)

      end

      subroutine potvort(pv,uu,vv,th,sp,cl,f,ie,je,ke,ak,bk,
     >     vmin,vmax)
C     ======================================================
C     calculate PV
 
c     argument declaration
      integer   ie,je,ke
      real      pv(ie,je,ke),uu(ie,je,ke),vv(ie,je,ke),
     >     th(ie,je,ke),sp(ie,je),
     >     cl(ie,je),f(ie,je)
      real      ak(ke),bk(ke),vmin(4),vmax(4)
      real      pvpole
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dthdp,dudp,dvdp
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: uu2,dvdx,dudy,dthdx,dthdy
      integer   i,j,k,stat
 
      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dspdx(ie,je)'
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dspdy(ie,je)'
      allocate(dthdp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dthdp'
      allocate(dudp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dudp'
      allocate(dvdp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dvdp'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dvdx'
      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dudy'
      allocate(dthdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dthdx'
      allocate(dthdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dthdy'
      allocate(uu2(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating uu2'
 
      call ddp(th,dthdp,sp,ie,je,ke,ak,bk)
      call ddp(uu,dudp,sp,ie,je,ke,ak,bk)
      call ddp(vv,dvdp,sp,ie,je,ke,ak,bk)
      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(th,dthdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
      call ddh3(th,dthdy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
c     conversion of uu for spheric coordinates (uu*cos(phi))
      do k=1,ke
         uu2(1:ie,1:je,k)=uu(1:ie,1:je,k)*cl(1:ie,1:je)
      enddo
      call ddh3(uu2,dudy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
c     conversion of dudy for spheric coordinates (dudy/cos(phi))
      do k=1,ke
         dudy(1:ie,1:je,k)=dudy(1:ie,1:je,k)/cl(1:ie,1:je)
      enddo
 
      do k=1,ke
         pv(1:ie,1:je,k)=1.E6*9.80665*(
     >        -(-dudy(1:ie,1:je,k)+dvdx(1:ie,1:je,k)
     >        +f(1:ie,1:je))*dthdp(1:ie,1:je,k)
     >        -(dudp(1:ie,1:je,k)*dthdy(1:ie,1:je,k)
     >        -dvdp(1:ie,1:je,k)*dthdx(1:ie,1:je,k)))
      enddo
 
c     interpolate poles from neighbouring points on every level
      if (vmax(2).gt.89.5) then
        do k=1,ke
          pvpole=0.
          do i=1,ie
            pvpole=pvpole+pv(i,je-1,k)
          enddo
          pvpole=pvpole/real(ie)
          do i=1,ie
            pv(i,je,k)=pvpole
          enddo
        enddo
      endif
 
      if (vmin(2).lt.-89.5) then
        do k=1,ke
          pvpole=0.
          do i=1,ie
            pvpole=pvpole+pv(i,2,k)
          enddo
          pvpole=pvpole/real(ie)
          do i=1,ie
            pv(i,1,k)=pvpole
          enddo
        enddo
      endif
 
      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dthdp)) DEALLOCATE(dthdp)
      IF (ALLOCATED(dudp)) DEALLOCATE(dudp)
      IF (ALLOCATED(dvdp)) DEALLOCATE(dvdp)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dthdx)) DEALLOCATE(dthdx)
      IF (ALLOCATED(dthdy)) DEALLOCATE(dthdy)
      IF (ALLOCATED(uu2)) DEALLOCATE(uu2)
 
      end

      subroutine potvort_i(pv,uu,vv,p,cl,f,ie,je,ke,ak,
     >     vmin,vmax,mdv)
C     =================================================
C     calculate isentropic PV
C     !! calculating PV makes sense only
C     when theta levels are very close to each others !!
C     with misdat (mdv) treatment 

c     argument declaration
      integer   ie,je,ke
      real      pv(ie,je,ke),uu(ie,je,ke),vv(ie,je,ke),
     >     	p(ie,je,ke),cl(ie,je),f(ie,je)
      real      ak(ke),vmin(4),vmax(4),mdv
      real      pvpole
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dpdth
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: uu2,dvdx,dudy
      integer   i,j,k,stat
 
      allocate(dpdth(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort_i: error allocating dpdth'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort_i: error allocating dvdx'
      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort_i: error allocating dudy'
      allocate(uu2(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort_i: error allocating uu2'
 
c     calculate dp/dth
c     k=1
      do i=1,ie
         do j=1,je
            if ((p(i,j,2).eq.mdv).or.(p(i,j,1).eq.mdv)) then
               dpdth(i,j,1)=mdv
            else
               dpdth(i,j,1)=100.*(p(i,j,2)-p(i,j,1))/(ak(2)-ak(1))
            endif
         enddo
      enddo
      
c     k=2,ke-1
      do i=1,ie
      do j=1,je
      do k=2,ke-1
        if ((p(i,j,k+1).eq.mdv).or.(p(i,j,k-1).eq.mdv)) then
          dpdth(i,j,k)=mdv
        else
          dpdth(i,j,k)=100.*(p(i,j,k+1)-p(i,j,k-1))/(ak(k+1)-ak(k-1))
        endif
      enddo
      enddo
      enddo

c     k=ke
      do i=1,ie
      do j=1,je
        if ((p(i,j,ke).eq.mdv).or.(p(i,j,ke-1).eq.mdv)) then
          dpdth(i,j,ke)=mdv
        else
          dpdth(i,j,ke)=100.*(p(i,j,ke)-p(i,j,ke-1))/(ak(ke)-ak(ke-1))
        endif
      enddo
      enddo

      call ddh2m(vv,dvdx,cl,'X',ie,je,ke,vmin,vmax,mdv)
c     conversion of uu for spheric coordinates (uu*cos(phi))
      

      do i=1,ie
         do j=1,je
            do k=1,ke
               if (uu(i,j,k).ne.mdv) then
                  uu2(i,j,k)=uu(i,j,k)*cl(i,j)
               else
                  uu2(i,j,k)=mdv
               endif
            enddo
         enddo
      enddo
      call ddh2m(uu2,dudy,cl,'Y',ie,je,ke,vmin,vmax,mdv)
c     conversion of dudy for spheric coordinates (dudy/cos(phi))
      do i=1,ie
         do j=1,je
            do k=1,ke
               if (dudy(i,j,k).ne.mdv) then
                  dudy(i,j,k)=dudy(i,j,k)/cl(i,j)
               else
                  dudy(i,j,k)=mdv
               endif
            enddo
         enddo
      enddo
  
      do i=1,ie
      do j=1,je
      do k=1,ke
        if ((dudy(i,j,k).eq.mdv).or.(dvdx(i,j,k).eq.mdv).or.
     >      (dpdth(i,j,k).eq.mdv).or.(dpdth(i,j,k).eq.0.)) then
          pv(i,j,k)=mdv
        else
          pv(i,j,k)=-1.E6*9.80665*
     >      (-dudy(i,j,k)+dvdx(i,j,k)+f(i,j))/dpdth(i,j,k)
        endif
      enddo
      enddo
      enddo
 
c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pvpole=0.
          do i=1,ie
            pvpole=pvpole+pv(i,je-1,k)
          enddo
          pvpole=pvpole/real(ie)
          do i=1,ie
            pv(i,je,k)=pvpole
          enddo
        enddo
      endif
 
      if (vmin(2).eq.-90.) then
        do k=1,ke
          pvpole=0.
          do i=1,ie
            pvpole=pvpole+pv(i,2,k)
          enddo
          pvpole=pvpole/real(ie)
          do i=1,ie
            pv(i,1,k)=pvpole
          enddo
        enddo
      endif
 
      IF (ALLOCATED(dpdth)) DEALLOCATE(dpdth)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(uu2)) DEALLOCATE(uu2)
 
      end


      subroutine gpotvort_i2(dpv,uu,vv,p,cl,f,ie,je,ke,ak,
     >     vmin,vmax,mdv)
C     =================================================
C     calculate isentropic PV gradient length
C     GPV=SQRT(D[PV:X]^2+D[PV:Y]^2) 
C     !! it calculates PV with potvort_i, what makes sense only
C     when theta levels are very close to each others !!
C     with misdat (mdv) treatment

c     argument declaration
      integer   ie,je,ke
      real      dpv(ie,je,ke),uu(ie,je,ke),vv(ie,je,ke),
     >     p(ie,je,ke),cl(ie,je),f(ie,je)
      real      ak(ke),vmin(4),vmax(4),mdv

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: pv
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dpvdx,dpvdy
      integer   stat
      real      dpvpole

      allocate(pv(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'gpotvort_i2: error allocating pv'
      allocate(dpvdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'gpotvort_i2: error allocating dpvdx'
      allocate(dpvdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'gpotvort_i2: error allocating dpvdy'

      call potvort_i(pv,uu,vv,p,cl,f,ie,je,ke,ak,
     >        vmin,vmax,mdv)

      do k=1,ke
         call ddh2m(pv(1,1,k),dpvdx,cl,'X',ie,je,1,vmin,vmax,mdv)
         call ddh2m(pv(1,1,k),dpvdy,cl,'Y',ie,je,1,vmin,vmax,mdv)
         do j=1,je
            do i=1,ie
               if ((dpvdx(i,j).ne.mdv).and.(dpvdy(i,j).ne.mdv)) then
                  dpv(i,j,k)=1.E6*sqrt(dpvdx(i,j)**2.+dpvdy(i,j)**2.)
               else
                  dpv(i,j,k)=mdv
               endif
            enddo
         enddo
      enddo

      IF (ALLOCATED(pv)) DEALLOCATE(pv)
      IF (ALLOCATED(dpvdx)) DEALLOCATE(dpvdx)
      IF (ALLOCATED(dpvdy)) DEALLOCATE(dpvdy)

      end

      subroutine gpotvort_i(dpv,pv,cl,ie,je,ke,
     >     vmin,vmax,mdv)
C     =================================================
C     calculate isentropic PV gradient length
C     GPV=SQRT(D[PV:X]^2+D[PV:Y]^2) 
C     !! it uses PV as input, best calculated with
C     potvort (with data from P file)  !!
C     with misdat (mdv) treatment
C     990201 mark liniger

c     argument declaration
      integer   ie,je,ke
      real      dpv(ie,je,ke),pv(ie,je,ke),cl(ie,je)
      real      vmin(4),vmax(4),mdv

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dpvdx,dpvdy
      integer   stat

      allocate(dpvdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'gpotvort_i: error allocating dpvdx'
      allocate(dpvdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'gpotvort_i: error allocating dpvdy'

      do k=1,ke
         call ddh2m(pv(1,1,k),dpvdx,cl,'X',ie,je,1,vmin,vmax,mdv)
         call ddh2m(pv(1,1,k),dpvdy,cl,'Y',ie,je,1,vmin,vmax,mdv)
         do j=1,je
            do i=1,ie
               if ((dpvdx(i,j).ne.mdv).and.(dpvdy(i,j).ne.mdv)) then
                  dpv(i,j,k)=1.E6*sqrt(dpvdx(i,j)**2.+dpvdy(i,j)**2.)
               else
                  dpv(i,j,k)=mdv
               endif
            enddo
         enddo
      enddo

      IF (ALLOCATED(dpvdx)) DEALLOCATE(dpvdx)
      IF (ALLOCATED(dpvdy)) DEALLOCATE(dpvdy)

      end


      subroutine divqu(var,uu,vv,qq,cl,ie,je,ke,vmin,vmax,mdv)
C     ========================================================
C     calculate divergence DIV=D[U:X]+D[V*cos(phi):Y]/cos(phi)
C     with misdat (mdv) and pole treatment

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),qq(ie,je,ke),
     >			  cl(ie,je)
      real      vmin(4),vmax(4),mdv


c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: qv2,dqudx,dqvdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: qu,qv
      integer   stat
      real      pole

      allocate(qv2(ie,je),STAT=stat)
      if (stat.ne.0) print*,'divqu: error allocating qv2'
      allocate(dqudx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'divqu: error allocating dqudx'
      allocate(dqvdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'divqu: error allocating dqvdy'
      allocate(qu(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'divqu: error allocating qu'
      allocate(qv(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'divqu: error allocating qv'

      do k=1,ke
      do j=1,je
      do i=1,ie
        qu(i,j,k)=qq(i,j,k)*uu(i,j,k)
        qv(i,j,k)=qq(i,j,k)*vv(i,j,k)
      enddo
      enddo
      enddo

      do k=1,ke
        call ddh2m(qu(1,1,k),dqudx,cl,'X',ie,je,1,vmin,vmax,mdv)

c       conversion of qv for spheric coordinates (qv*cos(phi))
        do j=1,je
           do i=1,ie
              if (qv(i,j,k).ne.mdv) then
                 qv2(i,j)=qv(i,j,k)*cl(i,j)
              else
                 qv2(i,k)=mdv
              endif
           enddo
        enddo


        call ddh2m(qv2,dqvdy,cl,'Y',ie,je,1,vmin,vmax,mdv)

c       conversion of dqvdy for spheric coordinates (dqvdy/cos(phi))
        do j=1,je
           do i=1,ie
              if ((dqudx(i,j).ne.mdv).and.(dqvdy(i,j).ne.mdv)) then
                 var(i,j,k)=1.E6*(dqudx(i,j)+dqvdy(i,j)/cl(i,j))
              else
                 var(i,j,k)=mdv
              endif
           enddo
        enddo
      enddo

c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          counter=0.
          do i=1,ie
             if (var(i,je-1,k).ne.mdv) then
                counter=counter+1
                pole=pole+var(i,je-1,k)
             endif
          enddo
          if (counter.ne.0) then
             pole=pole/counter
          else
             pole=mdv
          endif
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif

      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          counter=0.
          do i=1,ie
             if (var(i,2,k).ne.mdv) then
                counter=counter+1
                pole=pole+var(i,2,k)
             endif
          enddo
          if (counter.ne.0) then
             pole=pole/counter
          else
             pole=mdv
          endif
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif

      IF (ALLOCATED(dqudx)) DEALLOCATE(dqudx)
      IF (ALLOCATED(dqvdy)) DEALLOCATE(dqvdy)
      IF (ALLOCATED(qu)) DEALLOCATE(qu)
      IF (ALLOCATED(qv)) DEALLOCATE(qv)
      IF (ALLOCATED(qv2)) DEALLOCATE(qv2)

      end
      subroutine div_i(var,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
C     =====================================================
C     calculate divergence DIV=D[U:X]+D[V*cos(phi):Y]/cos(phi)
C     with misdat (mdv) and pole treatment
C     990201 mark liniger
 
c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),cl(ie,je)
      real      vmin(4),vmax(4),mdv
 
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: vv2,dudx,dvdy
      integer   stat
      real      pole
 
      allocate(vv2(ie,je),STAT=stat)
      if (stat.ne.0) print*,'div_i: error allocating vv2'
      allocate(dudx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'div_i: error allocating dudx'
      allocate(dvdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'div_i: error allocating dvdy'
 
 
      do k=1,ke
        call ddh2m(uu(1,1,k),dudx,cl,'X',ie,je,1,vmin,vmax,mdv)
 
c       conversion of vv for spheric coordinates (vv*cos(phi))
        do j=1,je
           do i=1,ie
              if (vv(i,j,k).ne.mdv) then
                 vv2(i,j)=vv(i,j,k)*cl(i,j)
              else
                 vv2(i,k)=mdv
              endif
           enddo
        enddo
 
 
        call ddh2m(vv2,dvdy,cl,'Y',ie,je,1,vmin,vmax,mdv)
 
c       conversion of dvdy for spheric coordinates (dvdy/cos(phi))
        do j=1,je
           do i=1,ie
              if ((dudx(i,j).ne.mdv).and.(dvdy(i,j).ne.mdv)) then
                 var(i,j,k)=1.E6*(dudx(i,j)+dvdy(i,j)/cl(i,j))
              else
                 var(i,j,k)=mdv
              endif
           enddo
        enddo
      enddo
 
c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          counter=0.
          do i=1,ie
             if (var(i,je-1,k).ne.mdv) then
                counter=counter+1
                pole=pole+var(i,je-1,k)
             endif
          enddo
          if (counter.ne.0) then
             pole=pole/counter
          else
             pole=mdv
          endif
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif
 
      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          counter=0.
          do i=1,ie
             if (var(i,2,k).ne.mdv) then
                counter=counter+1
                pole=pole+var(i,2,k)
             endif
          enddo
          if (counter.ne.0) then
             pole=pole/counter
          else
             pole=mdv
          endif
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif
 
      IF (ALLOCATED(dudx)) DEALLOCATE(dudx)
      IF (ALLOCATED(vv2)) DEALLOCATE(vv2)
      IF (ALLOCATED(dvdy)) DEALLOCATE(dvdy)
 
      end

      subroutine def_i(var,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
C     =====================================================
C     calculate deformation (strain) DEF= 1e6*SQRT(
C       ( D[V:X]+D[U*cos(phi):Y]/cos(phi) )^2
C       ( D[U:X]-D[V*cos(phi):Y]/cos(phi) )^2 )
C     with misdat (mdv) treatment
C     990501 mark liniger

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),cl(ie,je)
      real      vmin(4),vmax(4),mdv
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: uu2,vv2,dudx,dvdx,dudy,dvdy
      integer   stat
      real      pole
 
      allocate(uu2(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating uu2'
      allocate(vv2(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating vv2'
      allocate(dudx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dudx'
      allocate(dvdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dvdx'
      allocate(dudy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dudy'
      allocate(dvdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dvdy'
 
 
      do k=1,ke
        call ddh2m(uu(1,1,k),dudx,cl,'X',ie,je,1,vmin,vmax,mdv)
        call ddh2m(vv(1,1,k),dvdx,cl,'X',ie,je,1,vmin,vmax,mdv)
 
c       conversion of uu/vv for spheric coordinates (??*cos(phi))
 
        do j=1,je
           do i=1,ie
              if (uu(i,j,k).ne.mdv) then
                 uu2(i,j)=uu(i,j,k)*cl(i,j)
              else
                 uu2(i,j)=mdv
              endif
              if (vv(i,j,k).ne.mdv) then
                 vv2(i,j)=vv(i,j,k)*cl(i,j)
              else
                 vv2(i,j)=mdv
              endif
           enddo
        enddo
 
        call ddh2m(uu2,dudy,cl,'Y',ie,je,1,vmin,vmax,mdv)
        call ddh2m(vv2,dvdy,cl,'Y',ie,je,1,vmin,vmax,mdv)
 
c       conversion of d?dy for spheric coordinates (d?dy/cos(phi))
        do j=1,je
           do i=1,ie
              if ((dudx(i,j).ne.mdv).and.
     >             (dudy(i,j).ne.mdv).and.
     >             (dvdx(i,j).ne.mdv).and.
     >             (dvdy(i,j).ne.mdv)) then
 
                 var(i,j,k)=1.e6*sqrt(
     >                (dvdx(i,j)+dudy(i,j)/cl(i,j))
     >                *(dvdx(i,j)+dudy(i,j)/cl(i,j))
     >                +(dudx(i,j)-dvdy(i,j)/cl(i,j))
     >                *(dudx(i,j)-dvdy(i,j)/cl(i,j)))
 
              else
                 var(i,j,k)=mdv
              endif
           enddo
        enddo
      enddo
 
c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,je-1,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif
 
      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,2,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif
 
      IF (ALLOCATED(uu2)) DEALLOCATE(uu2)
      IF (ALLOCATED(vv2)) DEALLOCATE(vv2)
      IF (ALLOCATED(dudx)) DEALLOCATE(dudx)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dvdy)) DEALLOCATE(dvdy)
 
      end


      subroutine defn_i(var,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
C     =====================================================
C     calculate normal deformation (normal strain) 
C     DEFN= 1e6*( D[U:X]-D[V*cos(phi):Y]/cos(phi) )
C     with misdat (mdv) treatment
C     990727 mark liniger

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),cl(ie,je)
      real      vmin(4),vmax(4),mdv
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: vv2,dudx,dvdy
      integer   stat
      real      pole
 
      allocate(vv2(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating vv2'
      allocate(dudx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dudx'
      allocate(dvdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dvdy'
 
 
      do k=1,ke
        call ddh2m(uu(1,1,k),dudx,cl,'X',ie,je,1,vmin,vmax,mdv)
 
c       conversion of uu/vv for spheric coordinates (??*cos(phi)) 
        do j=1,je
           do i=1,ie
              if (vv(i,j,k).ne.mdv) then
                 vv2(i,j)=vv(i,j,k)*cl(i,j)
              else
                 vv2(i,j)=mdv
              endif
           enddo
        enddo
 
        call ddh2m(vv2,dvdy,cl,'Y',ie,je,1,vmin,vmax,mdv)
 
c       final calculation
        do j=1,je
           do i=1,ie
              if ((dudx(i,j).ne.mdv).and.
     >             (dvdy(i,j).ne.mdv)) then 
                 var(i,j,k)=1.e6*(dudx(i,j)-dvdy(i,j)/cl(i,j))
              else
                 var(i,j,k)=mdv
              endif
           enddo
        enddo
      enddo
 
c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,je-1,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif
 
      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,2,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif
 
      IF (ALLOCATED(vv2)) DEALLOCATE(vv2)
      IF (ALLOCATED(dudx)) DEALLOCATE(dudx)
      IF (ALLOCATED(dvdy)) DEALLOCATE(dvdy)
 
      end




      subroutine defs_i(var,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
C     =====================================================
C     calculate shear deformation (shear strain) 
C     DEFS= 1e6*( D[V:X]+D[U*cos(phi):Y]/cos(phi)) 
C     with misdat (mdv) treatment
C     990727 mark liniger

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),cl(ie,je)
      real      vmin(4),vmax(4),mdv
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: uu2,dvdx,dudy
      integer   stat
      real      pole
 
      allocate(uu2(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating uu2'
      allocate(dvdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dvdx'
      allocate(dudy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'dev_i: error allocating dudy'
 
 
      do k=1,ke
        call ddh2m(vv(1,1,k),dvdx,cl,'X',ie,je,1,vmin,vmax,mdv)
 
c       conversion of uu/vv for spheric coordinates (??*cos(phi)) 
        do j=1,je
           do i=1,ie
              if (uu(i,j,k).ne.mdv) then
                 uu2(i,j)=uu(i,j,k)*cl(i,j)
              else
                 uu2(i,j)=mdv
              endif
           enddo
        enddo
 
        call ddh2m(uu2,dudy,cl,'Y',ie,je,1,vmin,vmax,mdv)
 
c       conversion of d?dy for spheric coordinates (d?dy/cos(phi))
        do j=1,je
           do i=1,ie
              if ((dudy(i,j).ne.mdv).and.
     >             (dvdx(i,j).ne.mdv)) then 
                 var(i,j,k)=1.e6*(dvdx(i,j)+dudy(i,j)/cl(i,j))
              else
                 var(i,j,k)=mdv
              endif
           enddo
        enddo
      enddo
 
c     interpolate poles from neighbouring points on every level
      if (vmax(2).eq.90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,je-1,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,je,k)=pole
          enddo
        enddo
      endif
 
      if (vmin(2).eq.-90.) then
        do k=1,ke
          pole=0.
          do i=1,ie
            pole=pole+var(i,2,k)
          enddo
          pole=pole/real(ie)
          do i=1,ie
            var(i,1,k)=pole
          enddo
        enddo
      endif
 
      IF (ALLOCATED(uu2)) DEALLOCATE(uu2)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
 
      end





      subroutine okuboweiss_i(var,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
C     =====================================================
C     calculate okubo-weiss parameter= DEF*DEF-VORT*VORT
C     with misdat (mdv) treatment
C     990720 mark liniger

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),cl(ie,je)
      real      vmin(4),vmax(4),mdv
 
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: def,vort
      integer   stat
      real      pole
 
      allocate(def(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'okuboweiss_i: error allocating def'
      allocate(vort(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'okuboweiss_i: error allocating vort'
 
      call def_i(def,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)
      call vort_i(vort,uu,vv,cl,ie,je,ke,vmin,vmax,mdv)

      do k=1,ke
         do j=1,je
            do i=1,ie
               if ((def(i,j,k).ne.mdv).and.
     >              (vort(i,j,k).ne.mdv)) then
                  var(i,j,k)=def(i,j,k)*def(i,j,k)
     >                 -vort(i,j,k)*vort(i,j,k)*1e4
               else
                  var(i,j,k)=mdv 
               endif
            enddo
         enddo
      enddo
      
      IF (ALLOCATED(def)) DEALLOCATE(def)
      IF (ALLOCATED(vort)) DEALLOCATE(vort)
 
      end

      subroutine rich(ri,sp,uu,vv,th,ie,je,ke,ak,bk)
C     ==============================================
C     calculate Richardson number

      real      R,p0,kappa
      data      R,p0,kappa /287.,100000.,0.286/

c     argument declaration
      integer   ie,je,ke
      real      ri(ie,je,ke),uu(ie,je,ke),vv(ie,je,ke),
     >          th(ie,je,ke),sp(ie,je)
      real      ak(ke),bk(ke)

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dthdp,dudp,dvdp
      integer   i,j,k,stat

      allocate(dthdp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dthdp'
*      allocate(vel(ie,je,ke),STAT=stat)
*      if (stat.ne.0) print*,'potvort: error allocating vel'
      allocate(dudp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dudp'
      allocate(dvdp(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'potvort: error allocating dvdp'


*      vel=sqrt(uu**2+vv**2)
      call ddp(th,dthdp,sp,ie,je,ke,ak,bk)
*      call ddp(vel,dveldp,sp,ie,je,ke,ak,bk)
       call ddp(uu,dudp,sp,ie,je,ke,ak,bk)
       call ddp(vv,dvdp,sp,ie,je,ke,ak,bk)

      do k=1,ke
      do j=1,je
      do i=1,ie
         pval=(ak(k)+bk(k)*sp(i,j))*100.
        ri(i,j,k)=-R*(p0**(-kappa))*(pval**(kappa-1.))*
     >        dthdp(i,j,k)/(dudp(i,j,k)**2.+dvdp(i,j,k)**2.)
*         ri(i,j,k)=-R*(p0**(-kappa))*(pval**(kappa-1.))*
*    >             dthdp(i,j,k)/(dveldp(i,j,k)**2.)
      enddo
      enddo
      enddo

      IF (ALLOCATED(dthdp)) DEALLOCATE(dthdp)
      IF (ALLOCATED(dudp)) DEALLOCATE(dudp)
      IF (ALLOCATED(dvdp)) DEALLOCATE(dvdp)

      end

      subroutine ddp(a,d,sp,ie,je,ke,ak,bk)
c-----------------------------------------------------------------------
c     Purpose: VERTICAL DERIVATIVE
c     Compute the vertical derivative without missing data checking.
c     The derivative is taken from array 'a' in the direction of 'P'.
c     The result is stored in array 'd'.
c     3 point weighted centered differencing is used.
c     The vertical level-structure of the data is of the form
c     p=ak(k)+bk(k)*ps.
c-----------------------------------------------------------------------
 
c     declaration of arguments
      integer       ie,je,ke
      real          a(ie,je,ke),d(ie,je,ke),sp(ie,je)
 
c     variable declaration
      integer       i,j,k
      real          dpu,dpl,quot,fac,psrf
      real          ak(ke),bk(ke)
 
c     specify scaling factor associated with derivation with respect
c     to pressure
      fac=0.01
 
c     compute vertical 3 point derivative
c     ---------------------------
c     3-point vertical derivative
c     ---------------------------
      do j=1,je
         do i=1,ie
c     get surface pressure at current grid-point
            psrf=sp(i,j)
c     points at k=1
            dpu=(ak(1)+bk(1)*psrf)-(ak(2)+bk(2)*psrf)
            d(i,j,1)=(a(i,j,1)-a(i,j,2))*fac/dpu
c     points at 1<k<ke
            do k=2,ke-1
               dpu=(ak(k)+bk(k)*psrf)-(ak(k+1)+bk(k+1)*psrf)
               dpl=(ak(k-1)+bk(k-1)*psrf)-(ak(k)+bk(k)*psrf)
               quot=dpu/dpl
               d(i,j,k)=(quot*(a(i,j,k-1)-a(i,j,k))
     &              +1./quot*(a(i,j,k)-a(i,j,k+1)))*fac/(dpu+dpl)
            enddo
c     points at k=ke
            dpl=(ak(ke-1)+bk(ke-1)*psrf)-(ak(ke)+bk(ke)*psrf)
            d(i,j,ke)=(a(i,j,ke-1)-a(i,j,ke))*fac/dpl
         enddo
      enddo
      end
 
 
      subroutine ddh3(a,d,ps,dps,cl,dir,ie,je,ke,datmin,datmax,ak,bk)
c-----------------------------------------------------------------------
c     Purpose: HORIZONTAL DERIVATIVE ON PRESSURE-SURFACES WITHOUT
c     MISSING DATA
c     The derivative is taken from array 'a' in the direction of 'dir',
c     where 'dir' is either 'X','Y'. The result is stored in array 'd'.
c     The routine accounts for derivatives at the pole and periodic
c     boundaries in the longitudinal direction (depending on
c     the value of datmin, datmax). If the data-set does not reach to
c     the pole, a one-sided derivative is taken. Pole-treatment is only
c     carried out if the data-set covers 360 deg in longitude, and it
c     requires that ie=4*ii+1, where ii is an integer.
c     History:
c     Daniel Luethi
c-----------------------------------------------------------------------
 
c     declaration of arguments
      integer       ie,je,ke
      real          a(ie,je,ke),d(ie,je,ke),cl(ie,je)
      real          ps(ie,je),dps(ie,je)
      real          datmin(4),datmax(4)
      character*(*) dir
 
c     variable declaration
      integer       i,j,k
      real          ak(ke),bk(ke),as(500),bs(500)
 
c     compute vertical derivatives of ak's and bk's
      do k=2,ke-1
         as(k)=(ak(k-1)-ak(k+1))/2.
         bs(k)=(bk(k-1)-bk(k+1))/2.
      enddo
      as(1 )=ak(1)-ak(2)
      bs(1 )=bk(1)-bk(2)
      as(ke)=ak(ke-1)-ak(ke)
      bs(ke)=bk(ke-1)-bk(ke)
 
c     compute horizontal derivatives on sigma surfaces
      call ddh2(a,d,cl,dir,ie,je,ke,datmin,datmax)
 
c     apply correction for horizontal derivative on p-surfaces
      do j=1,je
         do i=1,ie
            do k=2,ke-1
               d(i,j,k)=d(i,j,k)+bk(k)*dps(i,j)/2./(as(k)+
     &              bs(k)*ps(i,j))*(a(i,j,k+1)-a(i,j,k-1))
            enddo
            k=1
            d(i,j,k)=d(i,j,k)+bk(k)*dps(i,j)/(as(k)+
     &           bs(k)*ps(i,j))*(a(i,j,k+1)-a(i,j,k))
            k=ke
            d(i,j,k)=d(i,j,k)+bk(k)*dps(i,j)/(as(k)+
     &           bs(k)*ps(i,j))*(a(i,j,k)-a(i,j,k-1))
         enddo
      enddo
      end
 
 
      subroutine ddh2(a,d,cl,dir,ie,je,ke,datmin,datmax)
c-----------------------------------------------------------------------
c     Purpose: HORIZONTAL DERIVATIVE ON DATA-SURFACES WITHOUT MISSING
C     DATA
c     Compute the horizontal derivative without missing data checking.
c     The derivative is taken from array 'a' in the direction of 'dir',
c     where 'dir' is either 'X','Y'. The result is stored in array 'd'.
c     The routine accounts for derivatives at the pole and periodic
c     boundaries in the longitudinal direction (depending on
c     the value of datmin, datmax). If the data-set does not reach to
c     the pole, a one-sided derivative is taken. Pole-treatment is only
c     carried out if the data-set covers 360 deg in longitude, and it
c     requires that ie=4*ii+1, where ii is an integer.
c-----------------------------------------------------------------------
 
c     declaration of arguments
      integer       ie,je,ke
      real          a(ie,je,ke),d(ie,je,ke),cl(ie,je)
      real          datmin(4),datmax(4)
      character*(*) dir
 
c     local variable declaration
      integer       i,j,k,ip1,im1,jp1,jm1,ip,im,j1,j2
      real          dlat,dlon,coslat,dx,dy,dxr,dyr
      integer       northpl,southpl,lonper
 
c     rerd and circ are the mean radius and diameter of the earth in
c     meter
      real          rerd,circ,pi
      data          rerd,circ,pi /6.37e6,4.e7,3.141592654/
 
c     compute flags for pole and periodic treatment
      southpl=0
      northpl=0
      lonper =0
      j1=1
      j2=je
      if (abs(datmax(1)-datmin(1)-360.).lt.1.e-3) then
         lonper=1
         if (abs(datmin(2)+90.).lt.1.e-3) then
            southpl=1
            j1=2
         endif
         if (abs(datmax(2)-90.).lt.1.e-3) then
            northpl=1
            j2=je-1
         endif
      endif
 
      dlon=((datmax(1)-datmin(1))/float(ie-1)) *pi/180.
      dlat=((datmax(2)-datmin(2))/float(je-1)) *pi/180.
 
c     print *,'Computing derivative ',dir(1:1),
c     &        ' of an array dimensioned ',ie,je,ke
 
      if (dir(1:1).eq.'X') then
c     -----------------------------
c     derivation in the x-direction
c     -----------------------------
         do k=1,ke
 
c     do gridpoints at j1<=j<=j2
            do j=j1,j2
               coslat=cl(1,j)
 
c     do regular gridpoints at 1<i<ie, 1<j<je
               dx =rerd*coslat*dlon
               dxr=1./(2.*dx)
               do i=2,ie-1
                  ip1=i+1
                  im1=i-1
                  d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
               enddo            ! i-loop
c     completed regular gridpoints at 1<i<ie, 1<j<je
 
c     do gridpoints at i=1, i=ie, 1<j<je
               if (lonper.eq.1) then
c     use periodic boundaries
                  i=1
                  ip1=2
                  im1=ie-1
                  d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
                  d(ie,j,k)=d(1,j,k)
               else
c     use one-sided derivatives
                  dxr=1./dx
                  i=1
                  ip1=2
                  im1=1
                  d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
                  i=ie
                  ip1=ie
                  im1=ie-1
                  d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
               endif
c     completed gridpoints at i=1, i=ie, j1<=j<=j2
 
            enddo               ! j-loop
c     completed gridpoints at 1<j<je
 
c     do gridpoints at j=je
            if (northpl.eq.1) then
c     for these gridpoints, the derivative in the x-direction is a
c     derivative in the y-direction at another pole-gridpoint
               dy =rerd*dlat
               dyr=1./(2.*dy)
               j=je
               jp1=je-1
               jm1=je-1
               do i=1,ie
                  ip=mod(i-1+  (ie-1)/4,ie)+1
                  im=mod(i-1+3*(ie-1)/4,ie)+1
                  d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
               enddo            ! i-loop
c     completed gridpoints at j=je
            endif
c     do gridpoints at j=1
            if (southpl.eq.1) then
               dy =rerd*dlat
               dyr=1./(2.*dy)
               j=1
               jp1=2
               jm1=2
               do i=1,ie
                  ip=mod(i-1+  (ie-1)/4,ie)+1
                  im=mod(i-1+3*(ie-1)/4,ie)+1
                  d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
               enddo            ! i-loop
            endif
c     completed gridpoints at j=1
 
         enddo                  ! k-loop
 
      else if (dir(1:1).eq.'Y') then
c     -----------------------------
c     derivation in the y-direction
c     -----------------------------
         dy =dlat*rerd
         dyr=1./(2.*dy)
         do k=1,ke
            do i=1,ie
 
c     do regular gridpoints
               do j=2,je-1
                  jp1=j+1
                  jm1=j-1
                  d(i,j,k)=dyr*(a(i,jp1,k)-a(i,jm1,k))
               enddo
 
c     do gridpoints at j=je
               if (northpl.eq.1) then
c     pole-treatment
                  j=je
                  jm1=j-1
                  jp1=j-1
                  ip=mod(i-1+(ie-1)/2,ie)+1
                  im=i
                  d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
               else
c     one-sided derivative
                  j=je
                  jm1=j-1
                  jp1=j
                  d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
               endif
c     completed gridpoints at j=je
 
c     do gridpoints at j=1
               if (southpl.eq.1) then
c     pole-treatment
                  j=1
                  jm1=2
                  jp1=2
                  ip=i
                  im=mod(i-1+(ie-1)/2,ie)+1
                  d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
               else
c     one-sided derivative
                  j=1
                  jm1=1
                  jp1=2
                  d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
               endif
c     completed gridpoints at j=1
 
            enddo
         enddo
 
      endif
      end
 
      subroutine ddh2m(a,d,cl,dir,ie,je,ke,datmin,datmax,mdv)
c-----------------------------------------------------------------------
c     Purpose: HORIZONTAL DERIVATIVE ON DATA-SURFACES WITH MISSING
C     DATA
c     Compute the horizontal derivative with missing data checking.
c     The derivative is taken from array 'a' in the direction of 'dir',
c     where 'dir' is either 'X','Y'. The result is stored in array 'd'.
c     The routine accounts for derivatives at the pole and periodic
c     boundaries in the longitudinal direction (depending on
c     the value of datmin, datmax). If the data-set does not reach to
c     the pole, a one-sided derivative is taken. Pole-treatment is only
c     carried out if the data-set covers 360 deg in longitude, and it
c     requires that ie=4*ii+1, where ii is an integer.
c     mdv is the misdat parameter (real)
c-----------------------------------------------------------------------
 
c     declaration of arguments
      integer       ie,je,ke
      real          a(ie,je,ke),d(ie,je,ke),cl(ie,je)
      real          datmin(4),datmax(4),mdv
      character*(*) dir
 
c     local variable declaration
      integer       i,j,k,ip1,im1,jp1,jm1,ip,im,j1,j2
      real          dlat,dlon,coslat,dx,dy,dxr,dyr
      integer       northpl,southpl,lonper
 
c     rerd and circ are the mean radius and diameter of the earth in
c     meter
      real          rerd,circ,pi
      data          rerd,circ,pi /6.37e6,4.e7,3.141592654/
 
c     compute flags for pole and periodic treatment
      southpl=0
      northpl=0
      lonper =0
      j1=1
      j2=je
      if (abs(datmax(1)-datmin(1)-360.).lt.1.e-3) then
         lonper=1
         if (abs(datmin(2)+90.).lt.1.e-3) then
            southpl=1
            j1=2
         endif
         if (abs(datmax(2)-90.).lt.1.e-3) then
            northpl=1
            j2=je-1
         endif
      endif
 
      dlon=((datmax(1)-datmin(1))/float(ie-1)) *pi/180.
      dlat=((datmax(2)-datmin(2))/float(je-1)) *pi/180.
 
c     print *,'Computing derivative ',dir(1:1),
c     &        ' of an array dimensioned ',ie,je,ke
 
      if (dir(1:1).eq.'X') then
c     -----------------------------
c     derivation in the x-direction
c     -----------------------------
         do k=1,ke
 
c     do gridpoints at j1<=j<=j2
            do j=j1,j2
               coslat=cl(1,j)
 
c     do regular gridpoints at 1<i<ie, 1<j<je
               dx =rerd*coslat*dlon
               dxr=1./(2.*dx)
               do i=2,ie-1
                  ip1=i+1
                  im1=i-1
                  if ((a(ip1,j,k).ne.mdv).and.(a(im1,j,k).ne.mdv)) then
                     d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
                  else
                     d(i,j,k)=mdv
                  endif
               enddo            ! i-loop
c     completed regular gridpoints at 1<i<ie, 1<j<je
 
c     do gridpoints at i=1, i=ie, 1<j<je
               if (lonper.eq.1) then
c     use periodic boundaries
                  i=1
                  ip1=2
                  im1=ie-1
                  if ((a(ip1,j,k).ne.mdv).and.(a(im1,j,k).ne.mdv)) then
                     d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
                  else
                     d(i,j,k)=mdv
                  endif
                  d(ie,j,k)=d(1,j,k)
               else
c     use one-sided derivatives
                  dxr=1./dx
                  i=1
                  ip1=2
                  im1=1
                  if ((a(ip1,j,k).ne.mdv).and.(a(im1,j,k).ne.mdv)) then
                     d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
                  else
                     d(i,j,k)=mdv
                  endif
                  i=ie
                  ip1=ie
                  im1=ie-1
                  if ((a(ip1,j,k).ne.mdv).and.(a(im1,j,k).ne.mdv)) then
                     d(i,j,k)=dxr*(a(ip1,j,k)-a(im1,j,k))
                  else
                     d(i,j,k)=mdv
                  endif
               endif
c     completed gridpoints at i=1, i=ie, j1<=j<=j2
 
            enddo               ! j-loop
c     completed gridpoints at 1<j<je
 
c     do gridpoints at j=je
            if (northpl.eq.1) then
c     for these gridpoints, the derivative in the x-direction is a
c     derivative in the y-direction at another pole-gridpoint
               dy =rerd*dlat
               dyr=1./(2.*dy)
               j=je
               jp1=je-1
               jm1=je-1
               do i=1,ie
                  ip=mod(i-1+  (ie-1)/4,ie)+1
                  im=mod(i-1+3*(ie-1)/4,ie)+1
                  if ((a(ip,jp1,k).ne.mdv).and.(a(im,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               enddo            ! i-loop
c     completed gridpoints at j=je
            endif
c     do gridpoints at j=1
            if (southpl.eq.1) then
               dy =rerd*dlat
               dyr=1./(2.*dy)
               j=1
               jp1=2
               jm1=2
               do i=1,ie
                  ip=mod(i-1+  (ie-1)/4,ie)+1
                  im=mod(i-1+3*(ie-1)/4,ie)+1
                  if ((a(ip,jp1,k).ne.mdv).and.(a(im,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               enddo            ! i-loop
            endif
c     completed gridpoints at j=1
 
         enddo                  ! k-loop
 
      else if (dir(1:1).eq.'Y') then
c     -----------------------------
c     derivation in the y-direction
c     -----------------------------
         dy =dlat*rerd
         dyr=1./(2.*dy)
         do k=1,ke
            do i=1,ie
 
c     do regular gridpoints
               do j=2,je-1
                  jp1=j+1
                  jm1=j-1
                  if ((a(i,jp1,k).ne.mdv).and.(a(i,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=dyr*(a(i,jp1,k)-a(i,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               enddo
 
c     do gridpoints at j=je
               if (northpl.eq.1) then
c     pole-treatment
                  j=je
                  jm1=j-1
                  jp1=j-1
                  ip=mod(i-1+(ie-1)/2,ie)+1
                  im=i
                  if ((a(ip,jp1,k).ne.mdv).and.(a(im,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               else
c     one-sided derivative
                  j=je
                  jm1=j-1
                  jp1=j
                  if ((a(i,jp1,k).ne.mdv).and.(a(i,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               endif
c     completed gridpoints at j=je
 
c     do gridpoints at j=1
               if (southpl.eq.1) then
c     pole-treatment
                  j=1
                  jm1=2
                  jp1=2
                  ip=i
                  im=mod(i-1+(ie-1)/2,ie)+1
                  if ((a(ip,jp1,k).ne.mdv).and.(a(im,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=dyr*(a(ip,jp1,k)-a(im,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               else
c     one-sided derivative
                  j=1
                  jm1=1
                  jp1=2
                  if ((a(i,jp1,k).ne.mdv).and.(a(i,jm1,k).ne.mdv))
     >                 then
                     d(i,j,k)=2.*dyr*(a(i,jp1,k)-a(i,jm1,k))
                  else
                     d(i,j,k)=mdv
                  endif
               endif
c     completed gridpoints at j=1
 
            enddo
         enddo
 
      endif
      end
 
      real function cosd(arg)
c-----------------------------------------------------------------------
c     compute Cos of an argument in Degree instead of Radian
c-----------------------------------------------------------------------
      real,intent(IN) :: arg
      real cos
      real,parameter :: grad2rad=3.1415926/360.
      cosd=cos(arg*grad2rad)
      return
      end

      subroutine pressure(pr,sp,stag3,ie,je,ke,aklev,bklev,aklay,bklay)
c     =================================================================
c     argument declaration
      integer  ie,je,ke
      real,intent(OUT) :: pr(ie,je,ke)
      real,intent(IN)  :: sp(ie,je),stag3
      real,intent(IN)  :: aklev(ke),bklev(ke),aklay(ke),bklay(ke)
 
c     variable declaration
      integer  i,j,k
      real     psrf
 
c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf
 
c     computation pressure
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke
               if (stag3.eq.0.) then
                  pr(i,j,k)=prlev(k)
               else
                  pr(i,j,k)=prlay(k)
               endif
            enddo
         enddo
      enddo
      end

      subroutine pres(pr,sp,ie,je,ke,ak,bk)
c     =====================================
c     argument declaration
      integer  ie,je,ke
      real,intent(OUT) :: pr(ie,je,ke)
      real,intent(IN)  :: sp(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke)

c     variable declaration
      integer  i,j,k

c     computation pressure
      do i=1,ie
      do j=1,je
      do k=1,ke
        pr(i,j,k)=ak(k)+bk(k)*sp(i,j)
      enddo
      enddo
      enddo
      end

      subroutine coslat(cl,pollon,pollat,vmin,dx,dy,ie,je)
c     ====================================================
c     argument declaration
      integer  ie,je
      real,intent(OUT) :: cl(ie,je)
 
c     variable declaration
      real	pollon,pollat,vmin(3),vmax(3)
      real	lon,lat,rlon,rlat
      integer	i,j

      real	phstoph

      real	pi
      data	pi /3.141592654/

C     Calculate cos(latitude) array and the coriolis parameter
 
      if ((pollon.ne.0.).or.(pollat.ne.90.)) then
         do j=1,je
            rlat=vmin(2)+(j-1)*dy
            do i=1,ie
               rlon=vmin(1)+(i-1)*dx
               yphys=phstoph(rlat,rlon,pollat,pollon)
C              if I use sind(lat in deg): troubles at the N-pole
               cl(i,j)=cos(rlat*pi/180.)
            enddo
         enddo
      else
         do j=1,je
            lat=vmin(2)+(j-1)*dy
            lat=pi*lat/180.
            do i=1,ie
               cl(i,j)=cos(lat)
            enddo
         enddo
      endif

      return
      end

      subroutine corpar(f,pollon,pollat,vmin,dx,dy,ie,je)
c     ====================================================
c     argument declaration
      integer  ie,je
      real,intent(OUT) :: f(ie,je)
 
c     variable declaration
      real      pollon,pollat,vmin(3),vmax(3)
      real      lon,lat,rlon,rlat
      integer   i,j
 
      real      phstoph
 
      real      pi
      data      pi /3.141592654/
 
C     Calculate cos(latitude) array and the coriolis parameter
 
      if ((pollon.ne.0.).or.(pollat.ne.90.)) then
         do j=1,je
            rlat=vmin(2)+(j-1)*dy
            do i=1,ie
               rlon=vmin(1)+(i-1)*dx
               yphys=phstoph(rlat,rlon,pollat,pollon)
C              if I use sind(lat in deg): troubles at the N-pole
               lat=2.*pi*yphys/360.
               f(i,j)=0.000145444*sin(lat)
            enddo
         enddo
      else
         do j=1,je
            lat=vmin(2)+(j-1)*dy
            lat=pi*lat/180.
            do i=1,ie
               f(i,j)=0.000145444*sin(lat)
            enddo
         enddo
      endif
 
      return
      end

      subroutine vint(intfield,field,sp,pmin,pmax,ie,je,ke,ak,bk)
C     ===========================================================

c     argument declaration
      integer   ie,je,ke
      real      intfield(ie,je),field(ie,je,ke),sp(ie,je)
      real      ak(ke),bk(ke)
      real      pmin,pmax

c     variable declaration
      integer   i,j,k,kmin,kmax
      real      coeff,pval

c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf

C     ====================================================================
C     ===  begin of main part of this subroutine  ========================
C     ====================================================================

      do i=1,ie
      do j=1,je
        intfield(i,j)=0.
        psrf=sp(i,j)
        if (psrf.lt.pmin) then                     
          intfield(i,j)=-999.98999
          goto 556
        endif
        kmin=1
        kmax=1
        do k=2,ke
          if ((pp(k).lt.pmax).and.(pp(k-1).gt.pmax)) kmin=k
          if ((pp(k).lt.pmin).and.(pp(k-1).gt.pmin)) kmax=k-1
        enddo
c       check if field is not equal mdv
        do k=kmin,kmax+1
          if (field(i,j,k).eq.-999.98999) then
            intfield(i,j)=-999.98999
            goto 556
          endif
        enddo
c       pmin and pmax are inbetween same pair of layers
c       interpolate on intermediate value (pmin+pmax)/2.
        if (kmin.eq.kmax+1) then
          pval=(pmin+pmax)/2.
          coeff=(pval-pp(kmin-1))/(pp(kmin)-pp(kmin-1))
          intfield(i,j)=(pmax-pmin)*
     &      (field(i,j,kmin)*coeff+field(i,j,kmax)*(1.-coeff))
          goto 555
        endif
c       one layer is inbetween pmin and pmax
        if (kmin.eq.kmax) then
          if (psrf.lt.pmax) then
            intfield(i,j)=field(i,j,kmin)*(psrf-pmin)
          else
            intfield(i,j)=field(i,j,kmin)*(pmax-pmin)
          endif
          goto 555
        endif
c       loop for the interior levels
        do k=kmin+1,kmax-1
          intfield(i,j)=intfield(i,j)+field(i,j,k)*
     &              0.5*(pp(k-1)-pp(k+1))
        enddo
c       special treatment of the bounding levels
        if (kmin.eq.1) then
          if (psrf.lt.pmax) then
            intfield(i,j)=intfield(i,j)+
     &          field(i,j,1)*(0.5*(pp(1)-pp(2))+psrf-pp(1))
          else
            intfield(i,j)=intfield(i,j)+
     &          field(i,j,1)*(0.5*(pp(1)-pp(2))+pmax-pp(1))
          endif
        else
          coeff=(pmax-pp(kmin-1))/(pp(kmin)-pp(kmin-1))
            intfield(i,j)=intfield(i,j)+
     &      field(i,j,kmin)*0.5*(pmax-pp(kmin+1))+
     &      (field(i,j,kmin)*coeff+field(i,j,kmin-1)*(1.-coeff))*
     &      0.5*(pmax-pp(kmin))
        endif
        if (kmax.eq.ke) then          
          intfield(i,j)=intfield(i,j)+
     &           field(i,j,ke)*0.5*(pp(ke-1)-pp(ke))
        else
          coeff=(pmin-pp(kmax+1))/(pp(kmax)-pp(kmax+1))
          intfield(i,j)=intfield(i,j)+
     &      field(i,j,kmax)*0.5*(pp(kmax-1)-pmin)+
     &      (field(i,j,kmax)*coeff+field(i,j,kmax+1)*(1.-coeff))*
     &      0.5*(pp(kmax)-pmin)
        endif
  555   continue
C       Calculate mean value
        if (psrf.lt.pmax) then
          intfield(i,j)=intfield(i,j)/(psrf-pmin)
        else
          intfield(i,j)=intfield(i,j)/(pmax-pmin)
        endif
  556   continue
      enddo
      enddo

      return
      end

      subroutine nvint(intfield,field,sp,pmi,pma,
     >                 ie,je,ke,ak,bk,iflag)
C     ===========================================
c
c     iflag     input   =0 average (like SR vint); =1 integration

c     argument declaration
      integer   ie,je,ke,iflag
      real      intfield(ie,je),field(ie,je,ke),sp(ie,je)
      real      ak(ke),bk(ke)
      real      pmin,pmax,pmi,pma

c     variable declaration
      integer   ilev
      integer   i,j,k,kmin,kmax,kshift
      real      mdv,coeff,pval

c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf

C     ====================================================================
C     ===  begin of main part of this subroutine  ========================
C     ====================================================================

C     set misdat value
      mdv=-999.98999

C     determine whether data is on model levels (ilev=0)
C     or pressure levels (ilev=1)
      ilev=1
      do k=1,ke
        if (bk(k).ne.0.) ilev=0
      enddo

      do i=1,ie
      do j=1,je
        intfield(i,j)=0.
        psrf=sp(i,j)
        kmin=1
        kmax=1
        pmin=pmi
        pmax=pma
        if (ilev.eq.0) then
          if (psrf.lt.pmin) then
            intfield(i,j)=mdv
            goto 456
          endif
        endif
        do k=2,ke
          if ((pp(k).le.pmax).and.(pp(k-1).gt.pmax)) kmin=k
          if ((pp(k).le.pmin).and.(pp(k-1).gt.pmin)) kmax=k-1
        enddo
c       if model levels: set vint=mdv if one level is equal mdv
        if (ilev.eq.0) then
          do k=kmin,kmax+1
            if (field(i,j,k).eq.mdv) then
              intfield(i,j)=mdv
              goto 456
            endif
          enddo
        endif
c       if pressure levels: check whether lowest levels are not mdv
c       if they are: adjust pmax
        if (ilev.eq.1) then
          kshift=0
          do k=kmin,kmax
            if (field(i,j,k).eq.mdv) then
              kshift=kshift+1
              pmax=0.5*(pp(k)+pp(k+1))
            endif
          enddo
c         if mdv occur go even one level higher
          if (kshift.gt.0) kshift=kshift+1
          kmin=kmin+kshift
          if (kmin.gt.kmax+1) then
            intfield(i,j)=mdv
            goto 456
          endif
        endif
c       pmin and pmax are inbetween same pair of layers
c       interpolate on intermediate value (pmin+pmax)/2.
        if (kmin.eq.kmax+1) then
          pval=(pmin+pmax)/2.
          coeff=(pval-pp(kmin-1))/(pp(kmin)-pp(kmin-1))
          intfield(i,j)=(pmax-pmin)*
     &      (field(i,j,kmin)*coeff+field(i,j,kmax)*(1.-coeff))
          goto 455
        endif
c       one layer is inbetween pmin and pmax
        if (kmin.eq.kmax) then
          if (psrf.lt.pmax) then
            intfield(i,j)=field(i,j,kmin)*(psrf-pmin)
          else
            intfield(i,j)=field(i,j,kmin)*(pmax-pmin)
          endif
          goto 455
        endif
c       loop for the interior levels
        do k=kmin+1,kmax-1
          intfield(i,j)=intfield(i,j)+field(i,j,k)*
     &              0.5*(pp(k-1)-pp(k+1))
        enddo
c       special treatment of the bounding levels
        if (kmin.eq.1) then
          if (psrf.lt.pmax) then
            intfield(i,j)=intfield(i,j)+
     &          field(i,j,1)*(0.5*(pp(1)-pp(2))+psrf-pp(1))
          else
            intfield(i,j)=intfield(i,j)+
     &          field(i,j,1)*(0.5*(pp(1)-pp(2))+pmax-pp(1))
          endif
        else
          coeff=(pmax-pp(kmin-1))/(pp(kmin)-pp(kmin-1))
            intfield(i,j)=intfield(i,j)+
     &      field(i,j,kmin)*0.5*(pmax-pp(kmin+1))+
     &      (field(i,j,kmin)*coeff+field(i,j,kmin-1)*(1.-coeff))*
     &      0.5*(pmax-pp(kmin))
        endif
        if (kmax.eq.ke) then
          intfield(i,j)=intfield(i,j)+
     &           field(i,j,ke)*0.5*(pp(ke-1)-pp(ke))
        else
          coeff=(pmin-pp(kmax+1))/(pp(kmax)-pp(kmax+1))
          intfield(i,j)=intfield(i,j)+
     &      field(i,j,kmax)*0.5*(pp(kmax-1)-pmin)+
     &      (field(i,j,kmax)*coeff+field(i,j,kmax+1)*(1.-coeff))*
     &      0.5*(pp(kmax)-pmin)
        endif
  455   continue
C       Calculate mean or integrated value
        if (iflag.eq.0) then
          if ((ilev.eq.0).and.(psrf.lt.pmax)) then
            intfield(i,j)=intfield(i,j)/(psrf-pmin)
          else
            intfield(i,j)=intfield(i,j)/(pmax-pmin)
          endif
        else
          intfield(i,j)=100./9.8*intfield(i,j)
        endif
  456   continue
      enddo
      enddo

      return
      end

C     =================================================================
      subroutine z2hno3(z,nmbl,zval,prof,hno3)
C     =================================================================
C     Input is a hno3-profile ("prof") in function of "nmbl" Z-levels
C     stored in zval.
C     Output is the "hno3"-value on this "z"

      implicit none
      integer       nmbl

C     Argument declaration
      real      zval(nmbl),prof(nmbl),hno3,z

C     Further variables declaration
      real      re
      integer   i

C     Check if z is inside profile values
      if (z.gt.zval(nmbl)) then
        print*,'*** error: z is over the range of HNO3-profile:'
        print*,z,zval(nmbl)
        stop
      endif

      if (z.lt.zval(1)) then
        print*,'*** error: z is below the range of HNO3-profile:'
        print*,z,zval(1)
        stop
      endif

C     Interpolate z of table to find the HNO3 value
      do i=1,nmbl-1
        if (z.ge.zval(i) .and. z.lt.zval(i+1)) then
          re=(z-zval(i))/(zval(i+1)-zval(i))
          goto 22
        endif
        if (i.eq.nmbl-1) then
          print*,'***Error: z:',z,' outside HNO3-Profile...'
          stop
        endif
      enddo

 22   hno3=prof(i)+re*(prof(i+1)-prof(i))
*     print*,zi,zval(i),zval(i+1),re
      return

      end

C     =================================================================
      subroutine p2h2o(p,xih2o)
C     =================================================================

      implicit none

      integer      nn
      parameter    (nn=12)
      real         wprof(nn),pprof(nn)

c     variables declaration
      real     p,xih2o,k
      integer  n
      data     pprof/150.,130.,100.,90.,80.,70.,60.,50.,40.,30.,20.,10./
      data     wprof/4.57,4.4,4.25,4.28,4.38,4.53,4.7,4.97,5.4,
     >     5.9,6.6,6.85/

      do n=1,nn-1
         if (p.le.pprof(n).and.p.ge.pprof(n+1)) then
            k=(p-pprof(n))/(pprof(n+1)-pprof(n))
            goto 100
         endif
         if (n.eq.nn-1) then
           print*,p,pprof(1),pprof(nn)
           stop'Table too small'
         endif
      enddo

 100  continue

      xih2o=wprof(n)+k*(wprof(n+1)-wprof(n))
      return
      end

C     =================================================================
      subroutine sg(temp,pres,rad,as,xHNO3,xH2O,sedi,growth)
C     =================================================================
C     Input of sg are a Temperature[K], a pressure[hPa], a Radius[um],
C     a particle shape specification stored in "as" with s for
C     spherical particles, a for aspherical a mixing ratio of HNO3.
C     Output is the sedimentation factor and a growth rate for the
C     specified particle.
      implicit none

C     Constants declaration
      real       deltarho, g, cunn_a, cunn_b, cunn_c
      real       ff1, ff2
      real       gas_diffusivity_corr, mb2mm
      real       aspect_ratio, d2, k, pi, R

      parameter  (g=9.81,  deltarho=1.62e3)   ! deltarho=(NAT Density-Air Density)
      parameter  (cunn_a=1.257, cunn_b=0.4, cunn_c=1.1)!cunn stays for Cunningham...
      parameter  (ff1=1.12, ff2=0.58248) ! Form factors for asphericity 1:3
      parameter  (mb2mm=760./1000.)
C     For aspher. part., changing the next value implies changing others as formfactors
      parameter  (aspect_ratio=1./3.)
      parameter  (d2=3.e-10, k=1.380661e-23, pi=3.1415927, R=8.31441)

C     Arguments declaration
      double precision sedi,growth
      real             temp,pres,rad,xHNO3,xH2O
      character*(1)    as

C     Functions declaration
      double precision cunn_corr,grown

C     Further variables
      double precision cc,wStokes,lamda,knudsen,capacity,deltapress
      real       eta,p_HNO3,p_H2O

      eta = 6.45e-8 * temp

      wStokes = 2.* rad**2 * g * deltarho / (9. * eta)

      lamda= k * temp / ( 2.**0.5 * 100. * pres * pi * d2**2 )

      knudsen= lamda / rad

      if (as.eq.'s') then
         cc=cunn_corr(knudsen,1.,1.,cunn_a,cunn_b,cunn_c)
         capacity=1.
      else
         cc=cunn_corr(knudsen,ff1,ff2,cunn_a,cunn_b,cunn_c)
         capacity=sqrt(1.-aspect_ratio**2) / (aspect_ratio**(2./3.) *
     >            log( (1. + sqrt(1.-aspect_ratio**2)) / aspect_ratio) )
      endif

      sedi= wStokes * cc / rad**2

      p_H2O= pres * xH2O * 1.e-6
      p_HNO3= pres * xHNO3 * 1.e-9
      deltapress= p_HNO3 - ( exp( (-2.7836 - 0.00088 * temp) *
     >            log (p_H2O * mb2mm) + 89.7674 - 26242. / temp +
     >            0.021135 * temp ) / mb2mm )

      growth= capacity * grown(temp,pres,rad,sedi,deltapress,eta)
C      print*,'T,p,r,shape,xHNO3 -> sedi,growth,',temp,pres,rad,as,
C     >        xHNO3,sedi,growth
      end

C     =================================================================
      double precision function grown(temp,pres,rad,sedi,deltapress,eta)
C     =================================================================
      implicit none

      real      molmass_g,molmass_s,rho_s
      real      p0,T0,R,pi,gas_diffusivity_corr

      parameter (molmass_g=63.e-3 ,molmass_s=117.e-3 ,rho_s=1.62e3 )
      parameter (p0=1013.25,T0=273.15,pi=3.1415927,R=8.3144)
      parameter  (gas_diffusivity_corr=1.)

      double precision sedi,deltapress,gasdiffu,gasdiffcorr,prod
      real      temp,pres,rad,eta
      real      meanvel,Reynold,Schmidt,ventil

C      if ( temp.lt.233 .or. temp.gt.313 ) print*,'*** Warning:
C     >Temperature',temp,'is out of validity for diffusivity'

      gasdiffu= sqrt( 0.018 / molmass_g ) * 0.211e-4 *
     >                   ( temp / T0 )**1.94 * (P0 / pres)
      meanvel= sqrt( 8.* R * temp / (molmass_g * pi) )
      gasdiffcorr= gasdiffu / ( 1. + 4. * gasdiffu /
     >            ( gas_diffusivity_corr * rad * meanvel ) )
      Reynold=  sedi * rad**2 * 2. * rad / eta
      Schmidt= eta / gasdiffu
      prod= Schmidt**(1./3.) * sqrt(Reynold)

      if ( prod.lt.1.4) then
         ventil= 1. + 0.108 * prod**2
      else
         ventil= 0.78 + 0.308 * prod
      endif
      grown= gasdiffcorr * ventil * molmass_s * deltapress * 100. /
     >       ( rho_s * R * temp )
      return
      end

C     =================================================================
      double precision function cunn_corr(knudsen,ff1,ff2,a,b,c)
C     =================================================================
      double precision knudsen
      real      ff1,ff2,a,b,c

      cunn_corr=(1. + (ff2/ff1) * knudsen * (a+b*exp(-c/knudsen)))/ff1
      return
      end
