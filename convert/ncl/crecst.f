      program crecst

c     ************************************************************************************
c     * Create constants file for ECMWF data files                                       *
c     * Michael Sprenger / Autumn 2013                                                   *
c     ************************************************************************************

      implicit none

c     ------------------------------------------------------------------------------------
c     Declaration of parameters
c     ------------------------------------------------------------------------------------

      integer        iargc
      character*(80) arg
      character*80   cstname
      integer        stdate(5),datar(14)
      integer        nx,ny,nz
      real           xmin,ymin,xmax,ymax,dx,dy,pollon,pollat
      character*80   akbkfile
      real           aklay(1000),bklay(1000),aklay_t(1000),bklay_t(1000)
      real           aklev(1000),bklev(1000),aklev_t(1000),bklev_t(1000)
      integer        i
      character*11   datestr
      integer        ilev

c     ------------------------------------------------------------------------------------
c     Get arguments and list of akbk
c     ------------------------------------------------------------------------------------

c     check for sufficient requested arguments
      if (iargc().lt.14) then
        print*,'USAGE: crecst cstname xmin xmax ymin ymax ...'
        print*,'       ...nx ny nz dx dy pollon pollat stdate akbkfile '
        call exit(1)
      endif
 
c     read and transform input
      call getarg(1,arg)
      cstname=trim(arg)
 
      call getarg(2,arg)
      read(arg,*) xmin

      call getarg(3,arg)
      read(arg,*) xmax

      call getarg(4,arg)
      read(arg,*) ymin

      call getarg(5,arg)
      read(arg,*) ymax

      call getarg(6,arg)
      read(arg,*) nx

      call getarg(7,arg)
      read(arg,*) ny

      call getarg(8,arg)
      read(arg,*) nz

      call getarg(9,arg)
      read(arg,*) dx

      call getarg(10,arg)
      read(arg,*) dy

      call getarg(11,arg)
      read(arg,*) pollon

      call getarg(12,arg)
      read(arg,*) pollat

      call getarg(13,arg)
      read(arg,*) datestr

      call getarg(14,arg)
      read(arg,*) akbkfile

c     Read table from akbk file
      open(10,file=akbkfile)
       do i=1,nz
          read(10,*) ilev,aklev_t(i),bklev_t(i)
       enddo
      close(10)

c     ------------------------------------------------------------------------------------
c     Prepare fields and write constants file
c     ------------------------------------------------------------------------------------

c     Set grid parameters
      datar(1)=nx
      datar(2)=ny
      datar(3)=int(1000.*ymax)
      datar(4)=int(1000.*xmin)
      datar(5)=int(1000.*ymin)
      datar(6)=int(1000.*xmax)
      datar(7)=int(1000.*dx)
      datar(8)=int(1000.*dy)
      datar(9)=nz
      datar(10)=1
      datar(11)=1
      datar(12)=0
      datar(13)=int(1000.*pollon) 
      datar(14)=int(1000.*pollat) 
      
c     Define aklev, bklev
      do i=1,nz
         aklev(i)=0.01 * aklev_t(nz-i+1)
         bklev(i)=bklev_t(nz-i+1)
      enddo

c     Define aklay, bklay
      do i=1,nz-1
         aklay(i+1) = 0.5*(aklev(i)+aklev(i+1))
         bklay(i+1) = 0.5*(bklev(i)+bklev(i+1))
      enddo
      aklay(1)=0.5*(0. + aklev(1))
      bklay(1)=0.5*(1. + bklev(1))

c     Set starting date
      read(datestr(1:4),  *) stdate(1)
      read(datestr(5:6),  *) stdate(2)
      read(datestr(7:8),  *) stdate(3) 
      read(datestr(10:11),*) stdate(4)
      stdate(5) = 0
      
C     Write the constants files
      call wricst(cstname,datar,aklev,bklev,aklay,bklay,stdate)
      write(*,*)
      write(*,*)'*** cst-file ',trim(cstname),' created'

      end



