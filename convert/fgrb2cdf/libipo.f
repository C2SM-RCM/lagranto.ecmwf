      real function int2d(ar,n1,n2,rid,rjd)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 2d-array to an arbitrary
c        location within the grid.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2)
c        n1,n2	   int   input   dimensions of ar
c        ri,rj     real  input   grid location to be interpolated to
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2
      real      ar(n1,n2), rid,rjd

c     local declarations
      integer   i,j,ip1,jp1,ih,jh
      real      frac0i,frac0j,frac1i,frac1j,ri,rj

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      ih=nint(ri)
      jh=nint(rj)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif

      if ((i.eq.ip1).and.(j.eq.jp1)) then
c        no interpolation at all
         int2d=ar(i,j)
      else
         frac0i=ri-float(i)
         frac0j=rj-float(j)
         frac1i=1.-frac0i
         frac1j=1.-frac0j
         int2d = ar(i  ,j  ) * frac1i * frac1j
     &         + ar(i  ,jp1) * frac1i * frac0j
     &         + ar(ip1,j  ) * frac0i * frac1j
     &         + ar(ip1,jp1) * frac0i * frac0j
      endif
      end
      real function int2dm(ar,n1,n2,rid,rjd,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 2d-array to an arbitrary
c        location within the grid. The interpolation includes the
c        testing of the missing data flag 'misdat'.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2)
c        n1,n2	   int   input   dimensions of ar
c        ri,rj     real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested
c     History:
c-----------------------------------------------------------------------
 
c     argument declarations
      integer   n1,n2
      real      ar(n1,n2), rid,rjd, misdat
 
c     local declarations
      integer   i,j,ip1,jp1,ih,jh
      real      frac0i,frac0j,frac1i,frac1j,ri,rj,int2d
 
c     check if routine without missing data checking can be called instead
      if (misdat.eq.0.) then
        int2dm=int2d(ar,n1,n2,rid,rjd)
        return
      endif
 
c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      ih=nint(ri)
      jh=nint(rj)
 
c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif
 
c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif
 
      if ((i.eq.ip1).and.(j.eq.jp1)) then
c        no interpolation at all
         int2dm=ar(i,j)
      else
         if ((misdat.eq.ar(i  ,j  )).or.
     &       (misdat.eq.ar(i  ,jp1)).or.
     &       (misdat.eq.ar(ip1,j  )).or.
     &       (misdat.eq.ar(ip1,jp1))) then
           int2dm=misdat
         else
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int2dm = ar(i  ,j  ) * frac1i * frac1j
     &            + ar(i  ,jp1) * frac1i * frac0j
     &            + ar(ip1,j  ) * frac0i * frac1j
     &            + ar(ip1,jp1) * frac0i * frac0j
         endif
      endif
      end
      real function int2dp(ar,n1,n2,rid,rjd)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 2d-array to an arbitrary
c        location within the grid. The 2d-array is periodic in the
c        i-direction: ar(0,.)=ar(n1,.) and ar(1,.)=ar(n1+1,.).
c        Therefore rid can take values in the range 0,...,n1+1
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2)
c        n1,n2     int   input   dimensions of ar
c        ri,rj     real  input   grid location to be interpolated to
c     History:
c-----------------------------------------------------------------------
 
c     argument declarations
      integer   n1,n2
      real      ar(0:n1+1,n2), rid,rjd
 
c     local declarations
      integer   i,j,ip1,jp1,ih,jh
      real      frac0i,frac0j,frac1i,frac1j,ri,rj
 
c     do linear interpolation
      ri=amax1(0.,amin1(float(n1+1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      ih=nint(ri)
      jh=nint(rj)
 
c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-5) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif
 
c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-5) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif
 
      if ((i.eq.ip1).and.(j.eq.jp1)) then
c        no interpolation at all
         int2dp=ar(i,j)
      else
         frac0i=ri-float(i)
         frac0j=rj-float(j)
         frac1i=1.-frac0i
         frac1j=1.-frac0j
         int2dp = ar(i  ,j  ) * frac1i * frac1j
     &         + ar(i  ,jp1) * frac1i * frac0j
     &         + ar(ip1,j  ) * frac0i * frac1j
     &         + ar(ip1,jp1) * frac0i * frac0j
      endif
      end
      real function cint2d(ar,n1,n2,n3,rid,rjd,ikd)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location in the horizontal. ikd specifies the level (must
c	 be an integer). A bicubic method is applied (following
c	 the Numerical Recipes).
c     Arguments:
c        ar           real  input   field, define as ar(n1,n2,n3)
c        n1,n2        int   input   dimensions of ar
c        rid,rjd      real  input   grid location to be interpolated to
c	 ikd	      int   input   level for interpolation
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3,ikd
      real      ar(n1,n2,n3), rid,rjd

c     local declarations
      integer	i,j,k
      real	y(4),y1(4),y2(4),y12(4)

c     indices of lower left corner of interpolation grid

      i=int(rid)
      j=int(rjd)
      k=ikd

c     do not interpolate if i or j are at the boundaries

      if ((i.eq.1).or.(j.eq.1).or.(i.eq.n1).or.(j.eq.n2)) then
        cint2d=ar(i,j,k)
        return
      endif

c     define the arrays y,y1,y2,y12 necessary for the bicubic
c     interpolation (following the Numerical Recipes).

      y(1)=ar(i,j,k)
      y(2)=ar(i+1,j,k)
      y(3)=ar(i+1,j+1,k)
      y(4)=ar(i,j+1,k)

      y1(1)=-(ar(i-1,j,k)-ar(i+1,j,k))/2.
      y1(2)=-(ar(i,j,k)-ar(i+2,j,k))/2.
      y1(3)=-(ar(i,j+1,k)-ar(i+2,j+1,k))/2.
      y1(4)=-(ar(i-1,j+1,k)-ar(i+1,j+1,k))/2.

      y2(1)=-(ar(i,j-1,k)-ar(i,j+1,k))/2.
      y2(2)=-(ar(i+1,j-1,k)-ar(i+1,j+1,k))/2.
      y2(3)=-(ar(i+1,j,k)-ar(i+1,j+2,k))/2.
      y2(4)=-(ar(i,j,k)-ar(i,j+2,k))/2.

      y12(1)=(ar(i+1,j+1,k)-ar(i+1,j-1,k)-ar(i-1,j+1,k)+ar(i-1,j-1,k))/4.
      y12(2)=(ar(i+2,j+1,k)-ar(i+2,j-1,k)-ar(i,j+1,k)+ar(i,j-1,k))/4.
      y12(3)=(ar(i+2,j+2,k)-ar(i+2,j,k)-ar(i,j+2,k)+ar(i,j,k))/4.
      y12(4)=(ar(i+1,j+2,k)-ar(i+1,j,k)-ar(i-1,j+2,k)+ar(i-1,j,k))/4.

      call bcuint(y,y1,y2,y12,i,j,rid,rjd,cint2d)
      return
      end
      real function int3d(ar,n1,n2,n3,rid,rjd,rkd)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
      if (abs(float(jh)-rj).lt.1.e-3) then
        j  =jh
        jp1=jh
      else
        j =min0(int(rj),n2-1)
        jp1=j+1
      endif

c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           int3d=ar(i,j,k)
c          print *,'int3d 00: ',rid,rjd,rkd,int3d
        else
c          horizontal interpolation only
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3d = ar(i  ,j  ,k  ) * frac1i * frac1j
     &           + ar(i  ,jp1,k  ) * frac1i * frac0j
     &           + ar(ip1,j  ,k  ) * frac0i * frac1j
     &           + ar(ip1,jp1,k  ) * frac0i * frac0j
c          print *,'int3d 10: ',rid,rjd,rkd,int3d
        endif
      else 
        frac0k=rk-float(k)
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           int3d = ar(i  ,j  ,k  ) * frac1k
     &           + ar(i  ,j  ,kp1) * frac0k
c          print *,'int3d 01: ',rid,rjd,rkd,int3d
        else
c          full 3d interpolation
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3d = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &           + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k 
     &           + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &           + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &           + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &           + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k 
     &           + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &           + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c          print *,'int3d 11: ',rid,rjd,rkd,int3d
        endif
      endif
      end
      real function int3dm(ar,n1,n2,n3,rid,rjd,rkd,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The interpolation includes the 
c        testing of the missing data flag 'misdat'. 
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd, misdat

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk,int3d

c     check if routine without missing data checking can be called instead
      if (misdat.eq.0.) then
        int3dm=int3d(ar,n1,n2,n3,rid,rjd,rkd)
        return
      endif

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif

c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           if (misdat.eq.ar(i,j,k)) then
             int3dm=misdat
           else
             int3dm=ar(i,j,k)
           endif
c          print *,'int3dm 00: ',rid,rjd,rkd,int3dm
        else
c          horizontal interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  ))) then
             int3dm=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dm = ar(i  ,j  ,k  ) * frac1i * frac1j
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j
c            print *,'int3dm 10: ',rid,rjd,rkd,int3dm
           endif
        endif
      else 
        frac0k=rk-float(k)
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1))) then
             int3dm=misdat
           else
             int3dm = ar(i  ,j  ,k  ) * frac1k
     &              + ar(i  ,j  ,kp1) * frac0k
c            print *,'int3dm 01: ',rid,rjd,rkd,int3dm
           endif
        else
c          full 3d interpolation
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1)).or.
     &         (misdat.eq.ar(i  ,jp1,kp1)).or.
     &         (misdat.eq.ar(ip1,j  ,kp1)).or.
     &         (misdat.eq.ar(ip1,jp1,kp1))) then
             int3dm=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dm = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &              + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &              + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k
     &              + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &              + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c            print *,'int3dm 11: ',rid,rjd,rkd,int3dm
           endif
        endif
      endif
      end
      real function int3dmlog(ar,n1,n2,n3,rid,rjd,rkd,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The interpolation includes the
c        testing of the missing data flag 'misdat'.
c        Prior to vertical interpolations the log is taken from the array.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd, misdat

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk,int3d

c     print*,'hallo in SR int3dmlog'
c     check if routine without missing data checking can be called instead
      if (misdat.eq.0.) then
        int3dmlog=int3d(ar,n1,n2,n3,rid,rjd,rkd)
        return
      endif

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif

c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           if (misdat.eq.ar(i,j,k)) then
             int3dmlog=misdat
           else
             int3dmlog=ar(i,j,k)
           endif
c          print *,'int3dmlog 00: ',rid,rjd,rkd,int3dmlog
        else
c          horizontal interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  ))) then
             int3dmlog=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dmlog = ar(i  ,j  ,k  ) * frac1i * frac1j
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j
c            print *,'int3dmlog 10: ',rid,rjd,rkd,int3dmlog
           endif
        endif
      else
        frac0k=rk-float(k)
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1))) then
             int3dmlog=misdat
           else
             int3dmlog = log(ar(i  ,j  ,k  )) * frac1k
     &                 + log(ar(i  ,j  ,kp1)) * frac0k
             int3dmlog = exp(int3dmlog)
c            print *,'int3dmlog 01: ',rid,rjd,rkd,int3dmlog
           endif
        else
c          full 3d interpolation
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1)).or.
     &         (misdat.eq.ar(i  ,jp1,kp1)).or.
     &         (misdat.eq.ar(ip1,j  ,kp1)).or.
     &         (misdat.eq.ar(ip1,jp1,kp1))) then
             int3dmlog=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dmlog = log(ar(i  ,j  ,k  )) * frac1i * frac1j * frac1k
     &                 + log(ar(i  ,jp1,k  )) * frac1i * frac0j * frac1k
     &                 + log(ar(ip1,j  ,k  )) * frac0i * frac1j * frac1k
     &                 + log(ar(ip1,jp1,k  )) * frac0i * frac0j * frac1k
     &                 + log(ar(i  ,j  ,kp1)) * frac1i * frac1j * frac0k
     &                 + log(ar(i  ,jp1,kp1)) * frac1i * frac0j * frac0k
     &                 + log(ar(ip1,j  ,kp1)) * frac0i * frac1j * frac0k
     &                 + log(ar(ip1,jp1,kp1)) * frac0i * frac0j * frac0k
             int3dmlog = exp(int3dmlog)
c            print *,'int3dmlog 11: ',rid,rjd,rkd,int3dmlog
           endif
        endif
      endif
      end

      real function int3dmvc(ar,n1,n2,n3,rid,rjd,rkd,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The interpolation includes the
c        testing of the missing data flag 'misdat'.
c	 In the vertical a Lagrangian cubic interpolation is
c	 performed.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd, misdat
 
c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh,klow,n
      real      frac0i,frac0j,frac1i,frac1j,ri,rj,rk
      real	int2d(4)

      real	int3dm

c     if n3 < 4 then do call linear interpolation in the vertical
      if (n3.lt.4) then
        int3dmvc=int3dm(ar,n1,n2,n3,rid,rjd,rkd,misdat)
        return
      endif
 
c     do linear interpolation in the horizontal
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)
 
c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif
 
c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif
 
c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif
 
      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           int3dmvc=ar(i,j,k)
c          print *,'int3dmvc 00: ',rid,rjd,rkd,int3dmvc
        else
c          horizontal interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  ))) then
             int3dmvc=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dmvc = ar(i  ,j  ,k  ) * frac1i * frac1j
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j
c            print *,'int3dmvc 10: ',rid,rjd,rkd,int3dmvc
           endif
        endif
      else
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1))) then
             int3dmvc=misdat
           else
             if (k-1.lt.1) then
               klow=1
             else if (k+2.gt.n3) then
               klow=n3-3
             else
               klow=k-1
             endif
             call cubint(ar(i,j,klow),ar(i,j,klow+1),ar(i,j,klow+2),
     &                   ar(i,j,klow+3),klow,rk,int3dmvc)
           endif
        else
c          full 3d interpolation
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1)).or.
     &         (misdat.eq.ar(i  ,jp1,kp1)).or.
     &         (misdat.eq.ar(ip1,j  ,kp1)).or.
     &         (misdat.eq.ar(ip1,jp1,kp1))) then
             int3dmvc=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             if (k-1.lt.1) then
               klow=1
             else if (k+2.gt.n3) then
               klow=n3-3
             else
               klow=k-1
             endif
             do n=1,4
               int2d(n) = ar(i  ,j  ,klow-1+n  ) * frac1i * frac1j
     &                  + ar(i  ,jp1,klow-1+n  ) * frac1i * frac0j
     &                  + ar(ip1,j  ,klow-1+n  ) * frac0i * frac1j
     &                  + ar(ip1,jp1,klow-1+n  ) * frac0i * frac0j
             enddo
             call cubint(int2d(1),int2d(2),int2d(3),int2d(4),
     &                   klow,rk,int3dmvc)
           endif
        endif
      endif
      end
      real function int4d(ar,n1,n2,n3,n4,rid,rjd,rkd,rld)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 4d-array to an arbitrary
c        location within the grid.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,..,n4  int   input   dimensions of ar
c        ri,..,rl  real  input   grid location to be interpolated to
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3,n4
      real      ar(n1,n2,n3,n4), rid,rjd,rkd,rld

c     local declarations
      integer   l,lp1,lh
      real      frac0l,frac1l,rl,int3d

c     do linear interpolation in l-direction
      rl=amax1(1.,amin1(float(n4),rld))
      lh=nint(rl)

c     Check for interpolation in l
*     if (abs(float(lh)-rl).lt.1.e-3) then
*       l  =lh
*       lp1=lh
*     else
        l =min0(int(rl),n4-1)
        lp1=l+1
*     endif

      if (l.eq.lp1) then
c       no interpolation in l
        int4d=int3d(ar(1,1,1,l),n1,n2,n3,rid,rjd,rkd)
      else
c       interpolation in l
        frac0l=rl-float(l)
        frac1l=1.-frac0l
        int4d = int3d(ar(1,1,1,l  ),n1,n2,n3,rid,rjd,rkd) * frac1l
     &        + int3d(ar(1,1,1,lp1),n1,n2,n3,rid,rjd,rkd) * frac0l
      endif
      end
      real function int4dm(ar,n1,n2,n3,n4,rid,rjd,rkd,rld,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 4d-array to an arbitrary
c        location within the grid. The interpolation includes the 
c        testing of the missing data flag 'misdat'. 
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,..,n4  int   input   dimensions of ar
c        ri,..,rl  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested.
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3,n4
      real      ar(n1,n2,n3,n4), rid,rjd,rkd,rld, misdat

c     local declarations
      integer   l,lp1,lh
      real      frac0l,frac1l,rl,rint0,rint1,int4d,int3dm

c     check whether missing data checking is required
      if (misdat.eq.0.) then
        int4dm=int4d(ar,n1,n2,n3,n4,rid,rjd,rkd,rld)
        return
      endif

c     do linear interpolation in l-direction
      rl=amax1(1.,amin1(float(n4),rld))
      lh=nint(rl)

c     Check for interpolation in l
*     if (abs(float(lh)-rl).lt.1.e-3) then
*       l  =lh
*       lp1=lh
*     else
        l =min0(int(rl),n4-1)
        lp1=l+1
*     endif

      if (l.eq.lp1) then
c       no interpolation in l
        int4dm = int3dm(ar(1,1,1,l),n1,n2,n3,rid,rjd,rkd,misdat)
      else
c       interpolation in l
        frac0l=rl-float(l)
        frac1l=1.-frac0l
        rint0 = int3dm(ar(1,1,1,l  ),n1,n2,n3,rid,rjd,rkd,misdat)
        rint1 = int3dm(ar(1,1,1,lp1),n1,n2,n3,rid,rjd,rkd,misdat)
        if ((rint0.eq.misdat).or.(rint1.eq.misdat)) then
          int4dm = misdat
        else
          int4dm = rint0*frac1l + rint1*frac0l
        endif
      endif
      end
      real function int3dl(ar,n1,n2,n3,levels,rid,rjd,rkd)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The vertical interpolation is linear
c	 in log(pressure).
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c	 levels	   real  input	 array contains pressure levels for ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c     History:
c	 Based on int3d 		July 93
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd
      real	levels(n3)

c     local declarations
      real	pval
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif

c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           int3dl=ar(i,j,k)
c          print *,'int3dl 00: ',rid,rjd,rkd,int3dl
        else
c          horizontal interpolation only
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3dl = ar(i  ,j  ,k  ) * frac1i * frac1j
     &            + ar(i  ,jp1,k  ) * frac1i * frac0j
     &            + ar(ip1,j  ,k  ) * frac0i * frac1j
     &            + ar(ip1,jp1,k  ) * frac0i * frac0j
c          print *,'int3dl 10: ',rid,rjd,rkd,int3dl
        endif
      else
*       frac0k=rk-float(k)
c       calculate the pressure value to be interpolated to
        pval=levels(int(rk))
     >            -(rk-aint(rk))*(levels(int(rk))-levels(int(rk)+1))
        frac0k=log(levels(int(rk))/pval)
     &        /log(levels(int(rk))/levels(int(rk)+1))
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           int3dl = ar(i  ,j  ,k  ) * frac1k
     &            + ar(i  ,j  ,kp1) * frac0k
c          print *,'int3dl 01: ',rid,rjd,rkd,int3dl
        else
c          full 3d interpolation
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3dl = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &            + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k
     &            + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &            + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &            + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &            + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k
     &            + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &            + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c          print *,'int3dl 11: ',rid,rjd,rkd,int3dl
        endif
      endif
      end
      subroutine bcucof(y,y1,y2,y12,c)
c-----------------------------------------------------------------------
c	Given arrays y,y1,y2 and y12, each of length 4, containing the
c	function, gradients, and cross derivative at the four grid points
c	of a rectangular grid cell (numbered counterclockwise from the 
c	lower left), and given d1 and d2, the length of the grid cell in
c	the 1- and 2-directions, this routine returns the table c that is
c	used by routine bcuint for biqubic interpolation.
c     Source: Numerical Recipes, Fortran Version, p.99
c-----------------------------------------------------------------------
      real 	c(4,4),y(4),y1(4),y2(4),y12(4),cl(16),x(16)
      integer	wt(16,16)

      data      wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,
     >		8*0,3,0,-9,6,-2,0,6,-4,
     >		10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,
     >		4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,
     >		10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,
     >		0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,
     >		10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,
     >		5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,
     >		10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

      real	xx
      integer	i,j,k,l

      do i=1,4			! pack a temporary vector x
        x(i)=y(i)
        x(i+4)=y1(i)
        x(i+8)=y2(i)
        x(i+12)=y12(i)
      enddo	

      do i=1,16			! matrix multiply by the stored table
        xx=0.
        do k=1,16
          xx=xx+wt(i,k)*x(k)
        enddo
        cl(i)=xx
      enddo

      l=0
      do i=1,4  		! unpack the result into the output table
      do j=1,4
        l=l+1
        c(i,j)=cl(l)
      enddo
      enddo

      return
      end
      subroutine bcuint(y,y1,y2,y12,x1l,x2l,x1,x2,ansy)
c-----------------------------------------------------------------------
c       Bicubic interpolation within a grid square. Input quantities are
c       y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower and
c       upper coordinates of the grid square in the 1-direction; x2l and
c       x2u likewise for the 2-direction; and x1,x2, the coordinates of
c       the desired point for the interpolation. The interplated function
c       value is returned as ansy. This routine calls bcucof.
c     Source: Numerical Recipes, Fortran Version, p.99/100
c     !!! changed the proposed code !!!
c-----------------------------------------------------------------------

      real      y(4),y1(4),y2(4),y12(4),c(4,4)
      real      ansy,x1,x2,t,u
      integer   i,x1l,x2l

      call bcucof(y,y1,y2,y12,c)

      t=x1-real(x1l)
      u=x2-real(x2l)

      ansy=0.

      do i=4,1,-1
        ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
      enddo

      return
      end
      subroutine cubint(ya1,ya2,ya3,ya4,k,x,y)
c-----------------------------------------------------------------------
c       Interface routine for SR polint for special case of cubic
c	interpolation in index space, with xa=k,k+1,k+2,k+3
c-----------------------------------------------------------------------
 
      integer	k
      real	ya1,ya2,ya3,ya4,x,y

      integer	n
      real	xa(4),ya(4),dy

      do n=1,4
        xa(1)=real(k)
        xa(2)=real(k+1)
        xa(3)=real(k+2)
        xa(4)=real(k+3)
  
        ya(1)=ya1
        ya(2)=ya2
        ya(3)=ya3
        ya(4)=ya4
      enddo

      call polint(xa,ya,4,x,y,dy)

      return
      end
      subroutine polint(xa,ya,n,x,y,dy)
c-----------------------------------------------------------------------
c       Given arrays xa and ya, each of length n, and given a value x, this
c	routine returns a value y, and an error estimate dy. If P(x) is the
c	polynomial of degree n-1 such that p(xa(i))=ya(i),i=1,...,n, then
c	the returned value y=p(x)
c     Source: Numerical Recipes, Fortran Version, p.82
c-----------------------------------------------------------------------

      integer	nmax,n
      parameter(nmax=10)
      real	xa(n),ya(n),x,y,dy
      integer	i,m,ns
      real	c(nmax),d(nmax),den,dif,dift,ho,hp,w

      ns=1
      dif=abs(x-xa(1))

      do i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo

      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if (2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo

      return
      end
      subroutine filt2d (a,af,f1,f2,nx,ny,fil,misdat,
     &                                     iperx,ipery,ispol,inpol)
c     =============================================================
c     Apply a conservative diffusion operator onto the 2d field a,
c     with full missing data checking.
c
c     a     real   inp  array to be filtered, dimensioned (nx,ny)
c     af    real   out  filtered array, dimensioned (nx,ny), can be
c                       equivalenced with array a in the calling routine
c     f1    real        workarray, dimensioned (nx+1,ny)
c     f2    real        workarray, dimensioned (nx,ny+1)
c     fil   real   inp  filter-coeff., 0<afil<=1. Maximum filtering with afil=1
c                       corresponds to one application of the linear filter.
c     misdat real  inp  missing-data value, a(i,j)=misdat indicates that
c                       the corresponding value is not available. The
c                       misdat-checking can be switched off with with misdat=0.
c     iperx int    inp  periodic boundaries in the x-direction (1=yes,0=no)
c     ipery int    inp  periodic boundaries in the y-direction (1=yes,0=no)
c     inpol int    inp  northpole at j=ny  (1=yes,0=no)
c     ispol int    inp  southpole at j=1   (1=yes,0=no)
c
c     Christoph Schaer, 1993
 
c     argument declaration
      integer     nx,ny
      real        a(nx,ny),af(nx,ny),f1(nx+1,ny),f2(nx,ny+1),fil,misdat
      integer     iperx,ipery,inpol,ispol
 
c     local variable declaration
      integer     i,j,is
      real        fh
 
c     compute constant fh
      fh=0.125*fil
 
c     compute fluxes in x-direction
      if (misdat.eq.0.) then
        do j=1,ny
        do i=2,nx
          f1(i,j)=a(i-1,j)-a(i,j)
        enddo
        enddo
      else
        do j=1,ny
        do i=2,nx
          if ((a(i,j).eq.misdat).or.(a(i-1,j).eq.misdat)) then
            f1(i,j)=0.
          else
            f1(i,j)=a(i-1,j)-a(i,j)
          endif
        enddo
        enddo
      endif
      if (iperx.eq.1) then
c       do periodic boundaries in the x-direction
        do j=1,ny
          f1(1,j)=f1(nx,j)
          f1(nx+1,j)=f1(2,j)
        enddo
      else
c       set boundary-fluxes to zero
        do j=1,ny
          f1(1,j)=0.
          f1(nx+1,j)=0.
        enddo
      endif
 
c     compute fluxes in y-direction
      if (misdat.eq.0.) then
        do j=2,ny
        do i=1,nx
          f2(i,j)=a(i,j-1)-a(i,j)
        enddo
        enddo
      else
        do j=2,ny
        do i=1,nx
          if ((a(i,j).eq.misdat).or.(a(i,j-1).eq.misdat)) then
            f2(i,j)=0.
          else
            f2(i,j)=a(i,j-1)-a(i,j)
          endif
        enddo
        enddo
      endif
c     set boundary-fluxes to zero
      do i=1,nx
        f2(i,1)=0.
        f2(i,ny+1)=0.
      enddo
      if (ipery.eq.1) then
c       do periodic boundaries in the x-direction
        do i=1,nx
          f2(i,1)=f2(i,ny)
          f2(i,ny+1)=f2(i,2)
        enddo
      endif
      if (iperx.eq.1) then
        if (ispol.eq.1) then
c         do south-pole
          is=(nx-1)/2
          do i=1,nx
            f2(i,1)=-f2(mod(i-1+is,nx)+1,2)
          enddo
        endif
        if (inpol.eq.1) then
c         do north-pole
          is=(nx-1)/2
          do i=1,nx
            f2(i,ny+1)=-f2(mod(i-1+is,nx)+1,ny)
          enddo
        endif
      endif
 
c     compute flux-convergence -> filter
      if (misdat.eq.0.) then
        do j=1,ny
        do i=1,nx
            af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
        enddo
        enddo
      else
        do j=1,ny
        do i=1,nx
          if (a(i,j).eq.misdat) then
            af(i,j)=misdat
          else
            af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
          endif
        enddo
        enddo
      endif
      end

      subroutine pipo(var3d,p3d,lev,var,nx,ny,nz,mdv,mode)
C     ====================================================

C     Interpolates the 3d variable var3d on the pressure surface
C     defined by lev. The interpolated field is returned as var.
C     p3d denotes the 3d pressure array.
C     mode determines the way of vertical interpolation:
C       mode=0 is for linear interpolation
C       mode=1 is for logarithmic interpolation

      integer   nx,ny,nz,mode
      real      lev,mdv
      real      var3d(nx,ny,nz),p3d(nx,ny,nz),var(nx,ny)

      integer   i,j,k
      real      kind
      real      int3dm

      do i=1,nx
      do j=1,ny

        kind=0.
        do k=1,nz-1
          if ((p3d(i,j,k).ge.lev).and.(p3d(i,j,k+1).le.lev)) then
            kind=float(k)+(p3d(i,j,k)-lev)/
     >                   (p3d(i,j,k)-p3d(i,j,k+1))
            goto 100
          endif
        enddo
 100    continue

        if (kind.eq.0.) then
          var(i,j)=mdv
        else
          var(i,j)=int3dm(var3d,nx,ny,nz,float(i),float(j),kind,mdv)
        endif

      enddo
      enddo

      return
      end

      subroutine thipo(var3d,th3d,lev,var,nx,ny,nz,mdv,mode)
C     ======================================================

C     Interpolates the 3d variable var3d on the isentropic surface
C     defined by lev. The interpolated field is returned as var.
C     th3d denotes the 3d theta array.
C     mode determines the way of vertical interpolation:
C       mode=0 is for linear interpolation
C       mode=1 is for logarithmic interpolation

      integer   nx,ny,nz,mode
      real      lev,mdv
      real      var3d(nx,ny,nz),th3d(nx,ny,nz),var(nx,ny)

      integer   i,j,k
      real      kind
      real      int3dm

      do i=1,nx
      do j=1,ny

        kind=0
        do k=1,nz-1
          if ((th3d(i,j,k).le.lev).and.(th3d(i,j,k+1).ge.lev)) then
            kind=float(k)+(th3d(i,j,k)-lev)/
     >                   (th3d(i,j,k)-th3d(i,j,k+1))
            goto 100
          endif
        enddo
 100    continue

        if (kind.eq.0) then
          var(i,j)=mdv
        else
          var(i,j)=int3dm(var3d,nx,ny,nz,float(i),float(j),kind,mdv)
        endif

      enddo
      enddo

      return
      end
