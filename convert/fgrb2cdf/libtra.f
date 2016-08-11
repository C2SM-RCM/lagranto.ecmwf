      subroutine copyar(array1,array2,nsize)
C     ======================================
C
C     Copy the array1 on array2.
C
C     array1	input	array to copy
C     array2	output	array to write
C     nsize     input   size of arrays

      integer   i,nsize
      real      array1(nsize),array2(nsize)

      do i=1,nsize
        array2(i)=array1(i)
      enddo

      return
      end
      subroutine clipreg2(plon,plat,varmin,varmax,nx,ny,ak,bk,nz,sp,
     >           prelevs,kcl0,kcl1,xcorner,ycorner,pmin,pmax,ind,
     >		 nxy,single,onelev,ierr)
C     ==============================================================
C
C     Clip a region to calculate trajectories.
C
C     plon      input   longitude of computational pole
C     plat      input   latitude of computational pole
C     varmin    input   array with minimum phys. coord.  of data region
C     varmax    input   array with maximum phys. coord.  of data region
C     nx        input   number of data points along latitude
C     ny        input   number of data points along longitude
C     ak        input   array contains ak values to calculate data-levels
C     bk        input   array contains bk values to calculate data-levels
C     nz        input   number of levels
C     sp        input   array contains surface pressure
C     prelevs   input   logical var (=.true. if data is on pressure levels)
C     kcl0,1    output  lowest,highest index for level of clipped region
C     xcorner	output	x-coordinates of corners of clipped region
C     ycorner	output	y-coordinates of corners of clipped region
C     pmin	output  minimum pressure value of clipped region
C     pmax	output	maximum pressure value of clipped region
C     ind	output	array with index of grid points within clipped region
C     nxy	output	number of points within clipped region
C     single	output 	logical flag for single trajectories
C     onelev    output  logical flag if clipmin(3)=clipmax(3)
C     ierr      output  error flag

      real      pmin,pmax,varmin(4),varmax(4)
      real	xcorner(4),ycorner(4),xmax,ymax,xmin,ymin,xx,yy
      real      ak(nz),bk(nz)
      real      sp(nx,ny)
      real      dx,dy
      real      plon,plat,xphys,yphys
      integer   kcl0,kcl1,i,j,k,l,nx,ny,nz,nxy,ierr
      integer	ind(400*100)
      logical   prelevs,onelev,single

      real      lmtolms,phtophs

C     Reset error flag

      ierr=0

C     Calculate grid increments

      dx=(varmax(1)-varmin(1))/(nx-1)
      dy=(varmax(2)-varmin(2))/(ny-1)

C     Read corners of region to clip from tape 9

      read(9,10)xcorner(1)
      read(9,10)ycorner(1)
      read(9,10)xcorner(2)
      read(9,10)ycorner(2)
      read(9,10)xcorner(3)
      read(9,10)ycorner(3)
      read(9,10)xcorner(4)
      read(9,10)ycorner(4)
      read(9,10)pmax
      read(9,10)pmin

   10 format(f7.3)

C     If necessary transform corners to the computational coordinates

      if ((plon.ne.0.).or.(plat.ne.90.)) then
        do i=1,4
          xphys=lmtolms(ycorner(i),xcorner(i),plat,plon)
          yphys=phtophs(ycorner(i),xcorner(i),plat,plon)
          xcorner(i)=xphys
          ycorner(i)=yphys
        enddo
      endif

C     Test if corners are within data domain

      do i=1,4
        if (xcorner(i).lt.varmin(1)) goto 980
        if (ycorner(i).lt.varmin(2)) goto 980
        if (xcorner(i).gt.varmax(1)) goto 980
        if (ycorner(i).gt.varmax(2)) goto 980
      enddo
      if (pmax.gt.varmin(3)) goto 980
      if (pmin.lt.varmax(3)) goto 980

C     Set onelev flag

      if (pmin.eq.pmax) then
        onelev=.true.
      else
        onelev=.false.
      endif

C     Check for single trajectory

      single=.false.
      if ((xcorner(1).eq.xcorner(2)).and.
     >    (xcorner(2).eq.xcorner(3)).and.
     >    (xcorner(3).eq.xcorner(4)).and.
     >    (ycorner(1).eq.ycorner(2)).and.
     >    (ycorner(2).eq.ycorner(3)).and.
     >    (ycorner(3).eq.ycorner(4))) then
        if (pmin.eq.pmax) then
          nxy=1
          single=.true.
          return
        else
          goto 980
        endif
      endif

C     Determine which grid points are within clipped region

      xmax=xcorner(1)
      xmin=xcorner(1)
      ymax=ycorner(1)
      ymin=ycorner(1)
      do i=2,4
        if (xcorner(i).lt.xmin) xmin=xcorner(i)
        if (xcorner(i).gt.xmax) xmax=xcorner(i)
        if (ycorner(i).lt.ymin) ymin=ycorner(i)
        if (ycorner(i).gt.ymax) ymax=ycorner(i)
      enddo

      nxy=0

      do i=1,nx
        xx=varmin(1)+(i-1)*dx
        do j=1,ny
          yy=varmin(2)+(j-1)*dy
          if (xx.lt.xmin) goto 970
          if (xx.gt.xmax) goto 970
          if (yy.lt.ymin) goto 970
          if (yy.gt.ymax) goto 970
          if ((xx-xcorner(1))*(ycorner(2)-ycorner(1))-
     >        (yy-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.) goto 970
          if ((xx-xcorner(2))*(ycorner(3)-ycorner(2))-
     >        (yy-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.) goto 970
          if ((xx-xcorner(3))*(ycorner(4)-ycorner(3))-
     >        (yy-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.) goto 970
          if ((xx-xcorner(4))*(ycorner(1)-ycorner(4))-
     >        (yy-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.) goto 970

          nxy=nxy+1
          ind(nxy)=i+(j-1)*nx
 970      continue
        enddo
      enddo
 
C     Calculate vertical indices for clipped region

      if (.not.onelev) then
        if (prelevs) then               ! data on pressure levels
          do k=1,nz-1
            if ((pmax.le.ak(k)).and.(pmax.ge.ak(k+1))) then
              if (ak(k)-pmax.gt.pmax-ak(k+1)) then
                kcl0=k+1
              else
                kcl0=k
              endif
            endif
            if ((pmin.le.ak(k)).and.(pmin.ge.ak(k+1))) then
              if (ak(k)-pmin.gt.pmin-ak(k+1)) then
                kcl1=k+1
              else
                kcl1=k
              endif
            endif
          enddo
        else                    ! data on model levels
          kcl0=1                        ! initialize kcl0 to lowest level
          kcl1=nz                       ! initialize kcl1 to highest level
          do k=nz-1,1,-1
          do l=1,nxy
            i=mod(mod(ind(l)-1,nx*ny),nx)+1
            j=int(mod(ind(l)-1,nx*ny)/nx)+1
            if ((pmax.le.ak(k)+bk(k)*sp(i,j)).and.
     >          (pmax.gt.ak(k+1)+bk(k+1)*sp(i,j))) kcl0=k+1
          enddo
          enddo
          do k=1,nz-1
          do l=1,nxy
            i=mod(mod(ind(l)-1,nx*ny),nx)+1
            j=int(mod(ind(l)-1,nx*ny)/nx)+1
            if ((pmin.lt.ak(k)+bk(k)*sp(i,j)).and.
     >          (pmin.ge.ak(k+1)+bk(k+1)*sp(i,j))) kcl1=k
          enddo
          enddo
        endif
      endif

C     Inform about eventually corrected boundary values

      if (.not.onelev) then
        if (prelevs) then
          write(*,*)pmin,' corrected to --> ',ak(kcl0)
          write(*,*)pmax,' corrected to --> ',ak(kcl1)
        else
          write(*,*)kcl0,' is level index for lower boundary'
          write(*,*)kcl1,' is level index for upper boundary'
        endif
      endif

C     Save eventually corrected boundary values

      if (prelevs.and.(.not.onelev)) then
        pmin=ak(kcl0)
        pmax=ak(kcl1)
      endif

      return

  980 ierr=1
      return
      end
      subroutine getnpoints(plon,plat,varmin,varmax,nx,ny,ak,bk,nz,sp,
     >           levtyp,grid,delh,delp,vcflag,xcorner,ycorner,
     >           pmin,pmax,nxy,nlev,ierr)
C     ================================================================
C
C     Determine approximate number of points within clipped region.
C
C     plon      input   longitude of computational pole
C     plat      input   latitude of computational pole
C     varmin    input   array with minimum phys. coord.  of data region
C     varmax    input   array with maximum phys. coord.  of data region
C     nx        input   number of data points along latitude
C     ny        input   number of data points along longitude
C     ak        input   array contains ak values to calculate data-levels
C     bk        input   array contains bk values to calculate data-levels
C     nz        input   number of levels
C     sp        input   array contains surface pressure
C     levtyp    input   int (=1 pressure =2 theta =3 model levels)
C     delh      input   desired horizontal resol of starting points (deg or km)
C     delp      input   desired vertical resol of starting points (hPa)
C     xcorner	output	corners of clipped region
C     pmin	output	minimum pressure value of clipped region
C     pmax	output	maximum pressure value of clipped region
C     nxy       output  number of points (horizontally) within clipped region
C     nv	output	number of vertical points within clipped region
C     ierr      output  error flag

      real      deltay,pi
      parameter (deltay=111.2)        ! distance in km between 2 lat circles
      parameter (pi=3.1415927)
 
      real      pmin,pmax,varmin(4),varmax(4),delx,dely,delh,delp
      real      xcorner(4),ycorner(4),xmax,ymax,xmin,ymin,x,y
      real      ak(nz),bk(nz)
      real      sp(nx,ny)
      real      dx,dy
      real      plon,plat,xphys,yphys
      integer   kcl0,kcl1,i,j,k,nx,ny,nz,nxy,nlev,vcflag,ierr
      integer	nnx,nny
      integer   levtyp
      logical	grid
 
      real      lmtolms,phtophs
 
C     Reset error flag
 
      ierr=0

C     Initialize nlev

      nlev=0
 
C     Calculate grid increments
 
      dx=(varmax(1)-varmin(1))/(nx-1)
      dy=(varmax(2)-varmin(2))/(ny-1)
 
C     Read vcflag (=1 if pmin/max denote pressure values;
C                  =2 if pmin/max denote layers)
      read(9,*)vcflag
 
C     Read corners of region to clip from tape 9
 
      read(9,10)xcorner(1)
      read(9,10)ycorner(1)
      read(9,10)xcorner(2)
      read(9,10)ycorner(2)
      read(9,10)xcorner(3)
      read(9,10)ycorner(3)
      read(9,10)xcorner(4)
      read(9,10)ycorner(4)
      read(9,10)pmax
      read(9,10)pmin
   10 format(f7.3)

      if ((plon.eq.0.).and.(plat.eq.90.)) then
        if (xcorner(1).eq.-999.) xcorner(1)=varmin(1)
        if (xcorner(2).eq.-999.) xcorner(2)=varmax(1)-dx
        if (xcorner(3).eq.-999.) xcorner(3)=varmax(1)-dx
        if (xcorner(4).eq.-999.) xcorner(4)=varmin(1)
        if (ycorner(1).eq.-999.) ycorner(1)=varmin(2)
        if (ycorner(2).eq.-999.) ycorner(2)=varmin(2)
        if (ycorner(3).eq.-999.) ycorner(3)=varmax(2)-dy
        if (ycorner(4).eq.-999.) ycorner(4)=varmax(2)-dy
      endif
 
C     Special treatment for rotated coordinate systems (transform corners to
C     rotated grid)
 
      if ((plon.ne.0.).or.(abs(plat-90.).gt.1.e-3)) then
        if (.not.grid) goto 975
        do i=1,4
          if ((xcorner(i).eq.-999.).and.(ycorner(i).eq.-999.)) then
            if (i.eq.1) then
              xcorner(i)=varmin(1)
              ycorner(i)=varmin(2)
            else if (i.eq.2) then
              xcorner(i)=varmax(1)-dx
              ycorner(i)=varmin(2)
            else if (i.eq.3) then
              xcorner(i)=varmax(1)-dx
              ycorner(i)=varmax(2)-dy
            else if (i.eq.4) then
              xcorner(i)=varmin(1)
              ycorner(i)=varmax(2)-dy
            endif
          else if ((xcorner(i).eq.-999.).or.(ycorner(i).eq.-999.)) then
            goto 985
          else
            xphys=lmtolms(ycorner(i),xcorner(i),plat,plon)
            yphys=phtophs(ycorner(i),xcorner(i),plat,plon)
            xcorner(i)=xphys
            ycorner(i)=yphys
          endif
        enddo
      endif
 
C     Test if corners are within data domain
 
      do i=1,4
*       if (xcorner(i).lt.varmin(1)) goto 980
*       if (ycorner(i).lt.varmin(2)) goto 980
*       if (xcorner(i).gt.varmax(1)) goto 980
*       if (ycorner(i).gt.varmax(2)) goto 980
        if (xcorner(i)+1.e-4.lt.varmin(1)) then
          print*,'1 ',xcorner(i),varmin(1)
          goto 980
        endif
        if (ycorner(i)+1.e-4.lt.varmin(2)) then
          print*,'2 ',xcorner(i),varmin(2)
          goto 980
        endif
        if (xcorner(i)-1.e-4.gt.varmax(1)) then
          print*,'3 ',xcorner(i),varmax(1)
          goto 980
        endif
        if (ycorner(i)-1.e-4.gt.varmax(2)) then
          print*,'4 ',xcorner(i),varmax(2)
          goto 980
        endif
      enddo

      if (vcflag.eq.1) then
        if (levtyp.eq.2) then
          if (pmax.lt.varmin(3)) goto 980
          if (pmin.gt.varmax(3)) goto 980
        else
          if (pmax.gt.varmin(3)) goto 980
          if (pmin.lt.varmax(3)) goto 980
        endif
      endif
 
C     Control output of clipped region
 
      write(*,*)'corners of the clipped region:'
      write(*,*)'        ',xcorner(1),ycorner(1)
      write(*,*)'        ',xcorner(2),ycorner(2)
      write(*,*)'        ',xcorner(3),ycorner(3)
      write(*,*)'        ',xcorner(4),ycorner(4)
 
C     Set nlev=1 if regions contains only one level
 
      if (pmin.eq.pmax) nlev=1
 
C     Check for single trajectory
 
      if ((xcorner(1).eq.xcorner(2)).and.
     >    (xcorner(2).eq.xcorner(3)).and.
     >    (xcorner(3).eq.xcorner(4)).and.
     >    (ycorner(1).eq.ycorner(2)).and.
     >    (ycorner(2).eq.ycorner(3)).and.
     >    (ycorner(3).eq.ycorner(4))) then
        if (nlev.eq.1) then
          nxy=1
          return
        else
          goto 980
        endif
      endif
 
C     Determine number of grid points within clipped region (horizontally)
 
      xmax=xcorner(1)
      xmin=xcorner(1)
      ymax=ycorner(1)
      ymin=ycorner(1)
      do i=2,4
        if (xcorner(i).lt.xmin) xmin=xcorner(i)
        if (xcorner(i).gt.xmax) xmax=xcorner(i)
        if (ycorner(i).lt.ymin) ymin=ycorner(i)
        if (ycorner(i).gt.ymax) ymax=ycorner(i)
      enddo
 
      nxy=0
 
      if (grid) then
        do i=1,nx
          x=varmin(1)+(i-1)*dx
          do j=1,ny
            y=varmin(2)+(j-1)*dy
            if (x.lt.xmin) goto 970
            if (x.gt.xmax) goto 970
            if (y.lt.ymin) goto 970
            if (y.gt.ymax) goto 970
            if ((x-xcorner(1))*(ycorner(2)-ycorner(1))-
     >          (y-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.) goto 970
            if ((x-xcorner(2))*(ycorner(3)-ycorner(2))-
     >          (y-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.) goto 970
            if ((x-xcorner(3))*(ycorner(4)-ycorner(3))-
     >          (y-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.) goto 970
            if ((x-xcorner(4))*(ycorner(1)-ycorner(4))-
     >          (y-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.) goto 970
   
            nxy=nxy+1
 970        continue
          enddo
        enddo
      else
        nny=nint((ymax-ymin)*deltay/delh)+1
        dely=(ymax-ymin)/real(nny-1)
        print*,ymax,ymin,delh,nny
        write(*,*)'value of dely set to ',dely,' (in deg)'
        do j=1,nny
          y=ymin+(j-1)*dely
          nnx=nint((xmax-xmin)*cos(y*pi/180.)*deltay/delh)
          if (nnx.gt.0) then
            delx=(xmax-xmin)/real(nnx)
          else
            delx=360.
            nnx=1
          endif
          write(*,*)'val of delx at lat ',y,' set to ',delx,' (in deg)'
          do i=1,nnx
            x=xmin+(i-1)*delx
            if ((x-xcorner(1))*(ycorner(2)-ycorner(1))-
     >          (y-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.01) goto 971
            if ((x-xcorner(2))*(ycorner(3)-ycorner(2))-
     >          (y-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.01) goto 971
            if ((x-xcorner(3))*(ycorner(4)-ycorner(3))-
     >          (y-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.01) goto 971
            if ((x-xcorner(4))*(ycorner(1)-ycorner(4))-
     >          (y-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.01) goto 971
  
            nxy=nxy+1
 971        continue
          enddo
        enddo
      endif

      write(*,*)'num of pts in the clipped region (horizontally): ',nxy
 
      if (.not.grid) then
        nlev=nint((pmax-pmin)/delp)+1
        write(*,*)'number of levels in the data domain: ',nlev
        return
      endif

C     Calculate vertical indices for whole data domain
 
      if (vcflag.eq.1) then
        if (nlev.ne.1) then
          if (levtyp.eq.1) then               ! data on pressure levels
            do k=1,nz-1
              if ((pmax.le.ak(k)).and.(pmax.ge.ak(k+1))) then
                if (ak(k)-pmax.gt.pmax-ak(k+1)) then
                  kcl0=k+1
                else
                  kcl0=k
                endif
              endif
              if ((pmin.le.ak(k)).and.(pmin.ge.ak(k+1))) then
                if (ak(k)-pmin.gt.pmin-ak(k+1)) then
                  kcl1=k+1
                else
                  kcl1=k
                endif
              endif
            enddo
          else if (levtyp.eq.2) then    ! data on theta levels
            do k=1,nz-1
              if ((pmax.ge.ak(k)).and.(pmax.le.ak(k+1))) then
                if (pmax-ak(k).gt.ak(k+1)-pmax) then
                  kcl0=k+1
                else
                  kcl0=k
                endif
              endif
              if ((pmin.ge.ak(k)).and.(pmin.le.ak(k+1))) then
                if (pmin-ak(k).gt.ak(k+1)-pmin) then
                  kcl1=k+1
                else
                  kcl1=k
                endif
              endif
            enddo
          else                    ! data on model levels
            kcl0=1                        ! initialize kcl0 to lowest level
            kcl1=nz                       ! initialize kcl1 to highest level
            do k=nz-1,1,-1
              do i=1,nx
              do j=1,ny
                if ((pmax.le.ak(k)+bk(k)*sp(i,j)).and.
     >              (pmax.gt.ak(k+1)+bk(k+1)*sp(i,j))) then
                  kcl0=k+1
                  goto 950
                endif
              enddo
              enddo
 950          continue
            enddo
            if (kcl0.eq.2) then
              do i=1,nx
              do j=1,ny
                if (pmax.ge.ak(1)+bk(1)*sp(i,j)) then
                  kcl0=1
                  goto 952
                endif
              enddo
              enddo
            endif
 952        continue
            do k=1,nz-1
              do i=1,nx
              do j=1,ny
                if ((pmin.lt.ak(k)+bk(k)*sp(i,j)).and.
     >              (pmin.ge.ak(k+1)+bk(k+1)*sp(i,j))) then
                  kcl1=k
                  goto 951
                endif
              enddo
              enddo
 951          continue
            enddo
          endif
          nlev=kcl1-kcl0+1
          write(*,*)'number of levels in the data domain: ',nlev
          write(*,*)'from indices',kcl0,' to ',kcl1
        endif
      else if (vcflag.eq.2) then
        kcl0=nint(pmin)
        kcl1=nint(pmax)
        nlev=kcl1-kcl0+1
        write(*,*)'number of levels in the data domain: ',nlev
        write(*,*)'from indices',kcl0,' to ',kcl1
      endif
 
      return
 
  975 ierr=3
      return
  980 ierr=1
      return
  985 ierr=2
      return
      end
      subroutine clipreg(plon,plat,varmin,varmax,nx,ny,ak,bk,nz,sp,
     >           levtyp,grid,delh,delp,vcflag,xcorner,ycorner,
     >           pmin,pmax,xx,yy,pp,ntra,ierr)
C     =============================================================
C
C     Clip a region to calculate trajectories.
C
C     plon      input   longitude of computational pole
C     plat      input   latitude of computational pole
C     varmin    input   array with minimum phys. coord.  of data region
C     varmax    input   array with maximum phys. coord.  of data region
C     nx        input   number of data points along latitude
C     ny        input   number of data points along longitude
C     ak        input   array contains ak values to calculate data-levels
C     bk        input   array contains bk values to calculate data-levels
C     nz        input   number of levels
C     sp        input   array contains surface pressure
C     levtyp    input   int (=1 pressure =2 theta =3 model levels)
C     delh	input	desired horizontal resol of starting points (km)
C     delp	input	desired vertical resol of starting points (hPa)
C     xx	output	starting coordinate of trajectories (longitude)
C     yy	output	starting coordinate of trajectories (latitude)
C     pp	output	starting coordinate of trajectories (pressure)
C     ntra	in/output number of grid points within clipped region
C     ierr      output  error flag
 
      real      deltay,pi
      parameter (deltay=111.2)        ! distance in km between 2 lat circles
      parameter (pi=3.1415927)

      real,allocatable, dimension (:) :: xval,yval,spval

      real	xx(*),yy(*),pp(*)
      real      pmin,pmax,pval,varmin(4),varmax(4),delh,delx,dely,delp
      real      xcorner(4),ycorner(4),xmax,ymax,xmin,ymin,x,y
      real      ak(nz),bk(nz)
      real      sp(nx,ny)
      real      dx,dy
      real      plon,plat,xphys,yphys
      integer   kcl0,kcl1,i,j,k,l,nx,ny,nz,nxy,nlev,ntra,vcflag,ierr
      integer   ii,jj
      integer	levtyp,stat
      logical	grid
 
      real      lmtolms,phtophs,int2d

C     Allocate memory for dynamical arrays

      allocate(xval(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array xval ***'
      allocate(yval(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array yval ***'
      allocate(spval(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spval ***'
 
C     Reset error flag
 
      ierr=0
 
C     Calculate grid increments
 
      dx=(varmax(1)-varmin(1))/(nx-1)
      dy=(varmax(2)-varmin(2))/(ny-1)

      write(*,*)'corners of the clipped region in SR clipreg:'
      write(*,*)'        ',xcorner(1),ycorner(1)
      write(*,*)'        ',xcorner(2),ycorner(2)
      write(*,*)'        ',xcorner(3),ycorner(3)
      write(*,*)'        ',xcorner(4),ycorner(4)
 
C     Set nlev=1 if regions contains only one level 
 
      if (pmin.eq.pmax) nlev=1
 
C     Check for single trajectory
 
      if ((xcorner(1).eq.xcorner(2)).and.
     >    (xcorner(2).eq.xcorner(3)).and.
     >    (xcorner(3).eq.xcorner(4)).and.
     >    (ycorner(1).eq.ycorner(2)).and.
     >    (ycorner(2).eq.ycorner(3)).and.
     >    (ycorner(3).eq.ycorner(4))) then
        if (nlev.eq.1) then
          ntra=1
          xx(1)=xcorner(1)
          yy(1)=ycorner(1)
          if (vcflag.eq.1) then
            pp(1)=pmin
          else if (vcflag.eq.2) then
            ii=nint((xx(1)-varmin(1))/dx)+1
            jj=nint((yy(1)-varmin(2))/dy)+1
            pp(1)=ak(nint(pmin))+bk(nint(pmin))*sp(ii,jj)
          endif
          return
        else
          goto 980
        endif
      endif
 
C     Determine which grid points are within clipped region (horizontally)
 
      xmax=xcorner(1)
      xmin=xcorner(1)
      ymax=ycorner(1)
      ymin=ycorner(1)
      do i=2,4
        if (xcorner(i).lt.xmin) xmin=xcorner(i)
        if (xcorner(i).gt.xmax) xmax=xcorner(i)
        if (ycorner(i).lt.ymin) ymin=ycorner(i)
        if (ycorner(i).gt.ymax) ymax=ycorner(i)
      enddo

C     Cut the last grid points of the data domain
C     (this prevents problems in subroutine trace: calls to getsdat)

*     if (xmax.gt.varmax(1)-dx) xmax=varmax(1)-dx
*     if (ymax.gt.varmax(2)-dy) ymax=varmax(2)-dy
 
      nxy=0
 
      if (grid) then
        do i=1,nx
          x=varmin(1)+(i-1)*dx
          do j=1,ny
            y=varmin(2)+(j-1)*dy
*           print*,x,y,xmin,xmax,ymin,ymax,
*    >             xcorner(1),ycorner(1),xcorner(2),ycorner(2),
*    >             xcorner(3),ycorner(3),xcorner(4),ycorner(4)
            if (x.lt.xmin) goto 970
            if (x.gt.xmax) goto 970
            if (y.lt.ymin) goto 970
            if (y.gt.ymax) goto 970
            if ((x-xcorner(1))*(ycorner(2)-ycorner(1))-
     >          (y-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.) goto 970
            if ((x-xcorner(2))*(ycorner(3)-ycorner(2))-
     >          (y-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.) goto 970
            if ((x-xcorner(3))*(ycorner(4)-ycorner(3))-
     >          (y-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.) goto 970
            if ((x-xcorner(4))*(ycorner(1)-ycorner(4))-
     >          (y-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.) goto 970
   
            nxy=nxy+1
*           print*,i,j,nxy
            xval(nxy)=x
            yval(nxy)=y
            spval(nxy)=sp(i,j)
 970        continue
          enddo
        enddo
      else
        nny=nint((ymax-ymin)*deltay/delh)+1
        dely=(ymax-ymin)/real(nny-1)
        print*,ymax,ymin,delh,nny
        write(*,*)'value of dely set to ',dely,' (in deg)'
        do j=1,nny
          y=ymin+(j-1)*dely
          nnx=nint((xmax-xmin)*cos(y*pi/180.)*deltay/delh)
          if (nnx.gt.0) then
            delx=(xmax-xmin)/real(nnx)
          else
            delx=360.
            nnx=1
          endif
          write(*,*)'val of delx at lat ',y,' set to ',delx,' (in deg)'
          do i=1,nnx
            x=xmin+(i-1)*delx
            if ((x-xcorner(1))*(ycorner(2)-ycorner(1))-
     >          (y-ycorner(1))*(xcorner(2)-xcorner(1)).gt.0.01) goto 971
            if ((x-xcorner(2))*(ycorner(3)-ycorner(2))-
     >          (y-ycorner(2))*(xcorner(3)-xcorner(2)).gt.0.01) goto 971
            if ((x-xcorner(3))*(ycorner(4)-ycorner(3))-
     >          (y-ycorner(3))*(xcorner(4)-xcorner(3)).gt.0.01) goto 971
            if ((x-xcorner(4))*(ycorner(1)-ycorner(4))-
     >          (y-ycorner(4))*(xcorner(1)-xcorner(4)).gt.0.01) goto 971
  
            nxy=nxy+1
*           print*,i,j,nxy
            xval(nxy)=x
            yval(nxy)=y
*           write(*,'(i8,2f9.2)')nxy,xval(nxy),yval(nxy)
            ri=(x-varmin(1))/dx+1.
            rj=(y-varmin(2))/dy+1.
            spval(nxy)=int2d(sp,nx,ny,ri,rj)
 971        continue
          enddo
        enddo
      endif
      write(*,*)'num of pts in the clipped region (horizontally): ',nxy

      if (.not.grid) then
        nlev=nint((pmax-pmin)/delp)+1
        delp=(pmax-pmin)/real(nlev-1)
        write(*,*)'number of levels in the data domain: ',nlev
        goto 960
      endif
 
C     Calculate vertical indices for clipped region
 
      if (vcflag.eq.1) then
        if (nlev.ne.1) then
          if (levtyp.eq.1) then               ! data on pressure levels
            do k=1,nz-1
              if ((pmax.le.ak(k)).and.(pmax.ge.ak(k+1))) then
                if (ak(k)-pmax.gt.pmax-ak(k+1)) then
                  kcl0=k+1
                else
                  kcl0=k
                endif
              endif
              if ((pmin.le.ak(k)).and.(pmin.ge.ak(k+1))) then
                if (ak(k)-pmin.gt.pmin-ak(k+1)) then
                  kcl1=k+1
                else
                  kcl1=k
                endif
              endif
            enddo
          else if (levtyp.eq.2) then	! data on theta levels
            do k=1,nz-1
              if ((pmax.ge.ak(k)).and.(pmax.le.ak(k+1))) then
                if (pmax-ak(k).gt.ak(k+1)-pmax) then
                  kcl0=k+1
                else
                  kcl0=k
                endif
              endif
              if ((pmin.ge.ak(k)).and.(pmin.le.ak(k+1))) then
                if (pmin-ak(k).gt.ak(k+1)-pmin) then
                  kcl1=k+1
                else
                  kcl1=k
                endif
              endif
            enddo
          else                    ! data on model levels
            kcl0=1                        ! initialize kcl0 to lowest level
            kcl1=nz                       ! initialize kcl1 to highest level
            do k=nz-1,1,-1
              do l=1,nxy
                if ((pmax.le.ak(k)+bk(k)*spval(l)).and.
     >              (pmax.gt.ak(k+1)+bk(k+1)*spval(l))) then
                  kcl0=k+1
                  goto 950
                endif
              enddo
 950          continue
            enddo
            if (kcl0.eq.2) then
              do l=1,nxy
                if (pmax.ge.ak(1)+bk(1)*spval(l)) then
                  kcl0=1
                  goto 952
                endif
              enddo
            endif
 952        continue
            do k=1,nz-1
              do l=1,nxy
                if ((pmin.lt.ak(k)+bk(k)*spval(l)).and.
     >              (pmin.ge.ak(k+1)+bk(k+1)*spval(l))) then
                  kcl1=k
                  goto 951
                endif
              enddo
 951          continue
            enddo
          endif
          nlev=kcl1-kcl0+1
          write(*,*)'number of levels in the clipped region: ',nlev
          write(*,*)'from indices ',kcl0,' to ',kcl1
        endif
      else if (vcflag.eq.2) then
        kcl0=nint(pmin)
        kcl1=nint(pmax)
        nlev=kcl1-kcl0+1
        write(*,*)'number of levels in the clipped region: ',nlev
        write(*,*)'from indices ',kcl0,' to ',kcl1
      endif

 960  continue

C     Serach grid points within clipped region

      ntra=0

      do k=1,nlev
      do l=1,nxy
        if (nlev.eq.1) then
          pval=pmin
        else
          if (grid) then
            pval=ak(k+kcl0-1)+bk(k+kcl0-1)*spval(l)
          else
            pval=pmin+(k-1)*delp
          endif
        endif

        if (vcflag.eq.1) then
          if ((levtyp.eq.1).and.(pval.le.spval(l))) then
            ntra=ntra+1
            xx(ntra)=xval(l)
            yy(ntra)=yval(l)
            pp(ntra)=pval
          else if (levtyp.eq.2) then
            ntra=ntra+1
            xx(ntra)=xval(l)
            yy(ntra)=yval(l)
            pp(ntra)=pval
          else if (levtyp.eq.3) then
            if (grid) then
              if ((pval.ge.pmin).and.(pval.le.pmax)) then
                ntra=ntra+1
                xx(ntra)=xval(l)
                yy(ntra)=yval(l)
                pp(ntra)=pval
              endif
            else
              if (pval.le.spval(l)) then
                ntra=ntra+1
                xx(ntra)=xval(l)
                yy(ntra)=yval(l)
                pp(ntra)=pval
              endif
            endif
          endif
        else if (vcflag.eq.2) then
          if ((levtyp.eq.1).and.(pval.le.spval(l))) then
            ntra=ntra+1
            xx(ntra)=xval(l)
            yy(ntra)=yval(l)
            pp(ntra)=pval
          else
            ntra=ntra+1
            xx(ntra)=xval(l)
            yy(ntra)=yval(l)
            pp(ntra)=pval
          endif
        endif
      enddo
      enddo
          
      write(*,*)'number of grid points in clipped region: ',ntra
 
C     Inform about eventually corrected boundary values
 
      if ((levtyp.eq.1).and.(nlev.gt.1)) then
        write(*,*)pmin,' corrected from ',pmin,' to --> ',ak(kcl0)
        write(*,*)pmax,' corrected from ',pmax,' to --> ',ak(kcl1)
        pmin=ak(kcl0)
        pmax=ak(kcl1)
      endif
 
      return
 
  980 ierr=1
      return
      end
      subroutine timediff(date1,date2,diff)
C     =====================================
 
C     New version with hour and minutes! (for hour and step [in hours]
C     use the routine oldtimediff!)
C
C     Calculates the time difference in hours (and minutes) for the two
C     dates specified by the two arrays date1 and date2.
C     They are expected to contain the following date information:
C     year      month   day     hour    minute.
C
C     date1     array specifying the first date
C     date2     array specifying the second date
C     diff      time differenc between date1 and date2 in hours
C
 
      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real      diff
      integer   ixday,imdiff,ihdiff,iddiff,j
      integer   yy,yy1,yy2
 
      idays(1)=31
      idays(2)=28
      idays(3)=31
      idays(4)=30
      idays(5)=31
      idays(6)=30
      idays(7)=31
      idays(8)=31
      idays(9)=30
      idays(10)=31
      idays(11)=30
      idays(12)=31
 
C     Check format of year (YYYY or YY - in case of YY assume 19YY)

      if (date1(1).lt.100) date1(1)=1900+date1(1)
      if (date2(1).lt.100) date2(1)=1900+date2(1)

C     Determine if the period between date1 and date2 contains a Feb.29
 
      ixday=0   ! extra day flag
 
      yy1=min(date1(1),date2(1))
      yy2=max(date1(1),date2(1))
      if (yy1.eq.yy2) then
        if (mod(yy1,4).eq.0) then
          idays(2)=29
        endif
      else
        if (mod(yy1,4).eq.0) then
          if (((yy1.eq.date1(1)).and.(date1(2).le.2)).or.
     >        ((yy1.eq.date2(1)).and.(date2(2).le.2))) then
            ixday=ixday+1
          endif
        endif
        if (mod(yy2,4).eq.0) then
          if (((yy2.eq.date1(1)).and.(date1(2).gt.2)).or.
     >        ((yy2.eq.date2(1)).and.(date2(2).gt.2))) then
            ixday=ixday+1
          endif
        endif
        if (yy2-yy1.gt.1) then
          do yy=yy1+1,yy2-1
            if (mod(yy,4).eq.0) then
              ixday=ixday+1
            endif
          enddo
        endif
      endif
 
      ihdiff=0  ! diff. in hours between date1/date2
      iddiff=0  ! diff. in days  between date1/date2
 
      if (date1(1).gt.date2(1)) then            ! compare years
        do j=date2(1),date1(1)-1
          iddiff=iddiff+365
        enddo
        iddiff=iddiff+ixday
      else if (date1(1).lt.date2(1)) then
        do j=date1(1),date2(1)-1
          iddiff=iddiff-365
        enddo
        iddiff=iddiff-ixday
      endif
 
      if (date1(2).gt.date2(2)) then            ! compare monthes
        do j=date2(2),date1(2)-1
          iddiff=iddiff+idays(j)
        enddo
      else if (date1(2).lt.date2(2)) then
        do j=date1(2),date2(2)-1
          iddiff=iddiff-idays(j)
        enddo
      endif
 
      iddiff=iddiff+date1(3)-date2(3)
      ihdiff=iddiff*24+date1(4)-date2(4)
      imdiff=ihdiff*60+date1(5)-date2(5)
 
      ihdiff=imdiff/60
      imdiff=mod(imdiff,60)
 
      diff=real(ihdiff)+real(imdiff)/60.
 
      return
      end
      subroutine oldtimediff(date1,date2,diff)
C     ========================================
 
C     Calculates the time difference in hours for the two dates specified
C     by the two arrays date1 and date2. They are expected to contain the
C     following date information:
C     year      month   day     time    step.
C
C     date1     array specifying the first date
C     date2     array specifying the second date
C     diff      time differenc between date1 and date2 in hours
C
C     Warning:  ihdiff is equal to 0 for date1/2 = 880203_12_00 and
C               880202_12_24 !!!
 
      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real      diff
      integer   ixday,ihdiff,iddiff,j
 
      data idays/31,28,31,30,31,30,31,31,30,31,30,31/
 
C     Determine if the period between date1 and date2 contains a Feb.29
 
      ixday=0   ! extra day flag
 
      if (((mod(date1(1),4).eq.0).and.(date1(1).ne.0)).or.
     >    ((mod(date2(1),4).eq.0).and.(date2(1).ne.0))) then
        if (((date1(2).le.2).and.(date2(2).gt.2)).or.
     >      ((date1(2).gt.2).and.(date2(2).le.2))) then
          ixday=1       ! set extra day flag to 1
          idays(2)=29   ! february has 29 days
        endif
      endif
 
      ihdiff=0  ! diff. in hours between date1/date2
      iddiff=0  ! diff. in days  between date1/date2
 
      if (date1(1).gt.date2(1)) then            ! compare years
        do j=date2(1),date1(1)-1
          iddiff=iddiff+365+ixday
        enddo
      else if (date1(1).lt.date2(1)) then
        do j=date1(1),date2(1)-1
          iddiff=iddiff-365-ixday
        enddo
      endif
 
      if (date1(2).gt.date2(2)) then            ! compare monthes
        do j=date2(2),date1(2)-1
          iddiff=iddiff+idays(j)
        enddo
      else if (date1(2).lt.date2(2)) then
        do j=date1(2),date2(2)-1
          iddiff=iddiff-idays(j)
        enddo
      endif
 
      iddiff=iddiff+date1(3)-date2(3)
      ihdiff=iddiff*24+date1(4)-date2(4)+date1(5)-date2(5)
 
      diff=real(ihdiff)
 
      return
      end
      subroutine newdate(date1,diff,date2)
C     ====================================
C
C     Routine calculates the new date when diff (in hours) is added to
C     date1.
C     date1	int	input	array contains a date in the form
C				year,month,day,hour,step
C     diff	real	input	timestep in hours to go from date1
C     date2	int	output	array contains new date in the same form

      integer   date1(5),date2(5)
      integer   idays(12)       ! array containing the days of the monthes
      real	diff
      logical	yearchange

      data idays/31,28,31,30,31,30,31,31,30,31,30,31/

      yearchange=.false.

      if ((mod(date1(1),4).eq.0).and.(date1(2).le.2)) idays(2)=29

      date2(1)=date1(1)
      date2(2)=date1(2)
      date2(3)=date1(3)
      date2(4)=date1(4)
      date2(5)=0
      date2(4)=date1(4)+int(diff)+date1(5)

      if (date2(4).ge.24) then
        date2(3)=date2(3)+int(date2(4)/24)
        date2(4)=date2(4)-int(date2(4)/24)*24
      endif
      if (date2(4).lt.0) then
        if (mod(date2(4),24).eq.0) then
          date2(3)=date2(3)-int(abs(date2(4))/24)
          date2(4)=date2(4)+int(abs(date2(4))/24)*24
        else
          date2(3)=date2(3)-(1+int(abs(date2(4))/24))
          date2(4)=date2(4)+(1+int(abs(date2(4))/24))*24
        endif
      endif

  100 if (date2(3).gt.idays(date2(2))) then
        if ((date2(2).eq.2).and.(mod(date2(1),4).eq.0)) idays(2)=29
        date2(3)=date2(3)-idays(date2(2))
        if (idays(2).eq.29) idays(2)=28
        date2(2)=date2(2)+1
        if (date2(2).gt.12) then
*         date2(1)=date2(1)+int(date2(2)/12)
*         date2(2)=date2(2)-int(date2(2)/12)*12
          date2(1)=date2(1)+1
          date2(2)=date2(2)-12
        endif
        if (date2(2).lt.1) then
          date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
          date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
        endif
        goto 100
      endif     
  200 if (date2(3).lt.1) then
        date2(2)=date2(2)-1
        if (date2(2).gt.12) then
          date2(1)=date2(1)+int(date2(2)/12)
          date2(2)=date2(2)-int(date2(2)/12)*12
        endif
        if (date2(2).lt.1) then
          date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
          date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
        endif
        if ((date2(2).eq.2).and.(mod(date2(1),4).eq.0)) idays(2)=29
        date2(3)=date2(3)+idays(date2(2))
        if (idays(2).eq.29) idays(2)=28
        goto 200
      endif

      if (date2(2).gt.12) then
        date2(1)=date2(1)+int(date2(2)/12)
        date2(2)=date2(2)-int(date2(2)/12)*12
      endif
      if (date2(2).lt.1) then
        date2(1)=date2(1)-(1+int(abs(date2(2)/12)))
        date2(2)=date2(2)+(1+int(abs(date2(2)/12)))*12
      endif

      if (date2(1).lt.1000) then
      if (date2(1).ge.100) date2(1)=date2(1)-100
      endif

      return
      end
      subroutine indipo(lonw,lone,lats,latn,dx,dy,nx,ny,nz,ak,bk,
     >           sp0,sp1,reltpos,levtyp,xpo,ypo,ppo,xind,yind,pind)
C     =============================================================
C
C     Calculates the position in the grid space from the physical
C     position.
C
      integer   nx,ny,nz,i
      real      dx,dy
      real      ak(nz),bk(nz)
      real      sp0(nx,ny),sp1(nx,ny),spval0,spval1,spval
      real      lonw,lone,lats,latn
      real      xpo,ypo,ppo,reltpos
      real      xind,yind,pind
      integer	levtyp
      real      int2d

C     levtyp=1 pressure levels
C     levtyp=2 theta levels
C     levtyp=3 model levels

C     Calculate the grid location to be interpolated to

      xind=(xpo-lonw)/dx+1.
      yind=(ypo-lats)/dy+1.

      if (nz.eq.1) then
        pind=1.
        return
      endif

      if (levtyp.eq.1) then
        do i=1,nz-1
          if ((ak(i).ge.ppo).and.(ak(i+1).le.ppo)) then
            pind=real(i)+(ak(i)-ppo)/(ak(i)-ak(i+1))
            return
          endif
        enddo
      else if (levtyp.eq.2) then
        do i=1,nz-1
          if ((ak(i).le.ppo).and.(ak(i+1).ge.ppo)) then
            pind=real(i)+(ak(i)-ppo)/(ak(i)-ak(i+1))
            return
          endif
        enddo
      else
        if (ppo.eq.1050.) then
          pind=1.
          return
        else
          if (reltpos.eq.0.) then
            spval=int2d(sp0,nx,ny,xind,yind)
          else if (reltpos.eq.1.) then
            spval=int2d(sp1,nx,ny,xind,yind)
          else
            spval0=int2d(sp0,nx,ny,xind,yind)
            spval1=int2d(sp1,nx,ny,xind,yind)
            spval=(1.-reltpos)*spval0+reltpos*spval1
          endif
          do i=1,nz-1
            if ((ak(i)+bk(i)*spval.ge.ppo).and.
     >          (ak(i+1)+bk(i+1)*spval.le.ppo)) then
              pind=real(i)+(ak(i)+bk(i)*spval-ppo)/
     >             (ak(i)+bk(i)*spval-ak(i+1)-bk(i+1)*spval)
              return
            endif
          enddo
        endif
      endif
      if (levtyp.lt.3) then
        write(*,*)'error in indipo: ppo ',ppo,' is outside the levels'
      else
        if (ppo.gt.(ak(1)+bk(1)*spval)) then
*         pind=1.
          pind=(spval-ppo)/(spval-(ak(1)+bk(1)*spval))
          return
        else
          write(*,*)'pos ',xpo,ypo,ppo,' sp ',ak(1)+bk(1)*spval
          write(*,*)'error in indipo: ppo is outside the levels'
        endif
      endif

      return
      end
      subroutine nindipo(lonw0,lats0,lonw1,lats1,dx,dy,nx0,ny0,
     >		 nx1,ny1,nz,ak,bk,sp0,sp1,reltpos,prelevs,
     >           xpo,ypo,ppo,xind0,yind0,xind1,yind1,pind)
C     ==============================================================
C
C     Calculates the position in the grid space from the physical
C     position.
C
      integer   nx0,ny0,nx1,ny1,nz,i
      real      dx,dy
      real      ak(nz),bk(nz)
      real      sp0(nx0,ny0),sp1(nx1,ny1),spval0,spval1,spval
      real      lonw0,lats0,lonw1,lats1
      real      xpo,ypo,ppo,reltpos
      real      xind0,yind0,xind1,yind1,pind
      logical   prelevs
      real      int3d
 
C     Calculate the grid location to be interpolated to
 
      xind0=(xpo-lonw0)/dx+1.
      yind0=(ypo-lats0)/dy+1.
      xind1=(xpo-lonw1)/dx+1.
      yind1=(ypo-lats1)/dy+1.
 
      if (prelevs) then
        do i=1,nz-1
          if ((ak(i).ge.ppo).and.(ak(i+1).le.ppo)) then
            pind=real(i)+(ak(i)-ppo)/(ak(i)-ak(i+1))
            return
          endif
        enddo
      else
        if (ppo.eq.1050.) then
          pind=1.
          return
        else
          if (reltpos.eq.0.) then
            spval=int3d(sp0,nx0,ny0,nz,xind0,yind0,1.)
          else if (reltpos.eq.1.) then
            spval=int3d(sp1,nx1,ny1,nz,xind1,yind1,1.)
          else
            spval0=int3d(sp0,nx0,ny0,nz,xind0,yind0,1.)
            spval1=int3d(sp1,nx1,ny1,nz,xind1,yind1,1.)
            spval=(1.-reltpos)*spval0+reltpos*spval1
          endif
          do i=1,nz-1
            if ((ak(i)+bk(i)*spval.ge.ppo).and.
     >          (ak(i+1)+bk(i+1)*spval.le.ppo)) then
              pind=real(i)+(ak(i)+bk(i)*spval-ppo)/
     >             (ak(i)+bk(i)*spval-ak(i+1)-bk(i+1)*spval)
              return
            endif
          enddo
        endif
      endif
      if (prelevs) then
        write(*,*)'error in indipo: ppo is outside the levels'
      else
        if (ppo.gt.(ak(1)+bk(1)*spval)) then
          pind=1.
          return
        else
          write(*,*)'pos ',xpo,ypo,ppo,' sp ',ak(1)+bk(1)*spval
          write(*,*)'error in indipo: ppo is outside the levels'
        endif
      endif
 
      return
      end
      subroutine indipo_mc2(lonw,lone,lats,latn,dx,dy,nx,ny,nz,
     >           bklev,bklay,bkt,zb,
     >           xpo,ypo,ppo,xind,yind,pindlev,pindlay)
C     =========================================================
C
C     Calculates the position in the grid space from the physical
C     position.
C
      implicit none

      integer   nx,ny,nz,i
      real      dx,dy
      real      bklev(nz),bklay(nz),bkt
      real      zb(nx,ny),spval0,spval1,spval,zbval
      real      lonw,lone,lats,latn
      real      xpo,ypo,ppo
      real      xind,yind,pindlev,pindlay
      real      int3d

C     Calculate the grid location to be interpolated to

      xind=(xpo-lonw)/dx+1.
      yind=(ypo-lats)/dy+1.

      if (ppo.eq.0.) then
        pindlev=1.                  
        pindlay=1.
        return
      endif
      zbval=int3d(zb,nx,ny,nz,xind,yind,1.)
      do i=1,nz-1
        if ((bklev(i)+(1.-bklev(i)/bkt)*zbval.le.ppo).and.
     >       bklev(i+1)+(1.-bklev(i+1)/bkt)*zbval.gt.ppo) then
          pindlev=real(i)+(ppo-bklev(i)-(1.-bklev(i)/bkt)*zbval)/
     >       ((1.-zbval/bkt)*(bklev(i+1)-bklev(i)))
          goto 40
        endif
      enddo
      if (ppo.lt.zbval) then
        write(*,*)'error in indipo_mc2: ppo is below the ground'
        call exit(1)
      else if (ppo.lt.bklev(1)+(1.-bklev(1)/bkt)*zbval) then
        pindlev=(ppo-zbval)/(bklev(1)-bklev(1)/bkt*zbval)
        return
      else
        write(*,*)'pos ',xpo,ypo,ppo,bklev(1)+(1.-bklev(1)/bkt)*zbval
        write(*,*)'error in indipo_mc2: ppo is outside the levels'
        call exit(1)
      endif
   40 do i=1,nz-1
        if ((bklay(i)+(1.-bklay(i)/bkt)*zbval.le.ppo).and.
     >       bklay(i+1)+(1.-bklay(i+1)/bkt)*zbval.gt.ppo) then
          pindlay=real(i)+(ppo-bklay(i)-(1.-bklay(i)/bkt)*zbval)/
     >       ((1.-zbval/bkt)*(bklay(i+1)-bklay(i)))
          return
        endif
      enddo
      if (ppo.lt.zbval) then
        write(*,*)'error in indipo_mc2: ppo is below the ground'
        call exit(1)
      else if (ppo.lt.bklay(1)+(1.-bklay(1)/bkt)*zbval) then
        pindlay=(ppo-zbval)/(bklay(1)-bklay(1)/bkt*zbval)
        return
      else
        write(*,*)'pos ',xpo,ypo,ppo,bklay(1)+(1.-bklay(1)/bkt)*zbval
        write(*,*)'error in indipo_mc2: ppo is outside the layers'
        call exit(1)
      endif

      end
      subroutine ipo(nx,ny,nz,field0,field1,reltpos,
     >		     xind,yind,pind,mdv,ipomode,value)
C     ================================================
C
C     Decides what kind of interpolation should be done.
C
      integer   nx,ny,nz
      real      field0(nx,ny,nz)
      real      field1(nx,ny,nz)
      real      reltpos
      real      xind,yind,pind
      real      value,mdv
      integer	ipomode

      if (ipomode.eq.1) then
        call linipo(nx,ny,nz,field0,field1,reltpos,
     >		    xind,yind,pind,mdv,value)
      else if (ipomode.eq.2) then
        call bicipo(nx,ny,nz,field0,field1,reltpos,
     >              xind,yind,pind,value)
      else if (ipomode.eq.3) then
        call vcuipo(nx,ny,nz,field0,field1,reltpos,
     >              xind,yind,pind,mdv,value)
      else if (ipomode.eq.11) then
        call linipolog(nx,ny,nz,field0,field1,reltpos,
     >                 xind,yind,pind,mdv,value)
      else
        write(*,*)'*** error in SR ipo: unknown ipomode'
        stop
      endif

      return
      end
      subroutine nipo(nx0,ny0,nx1,ny1,nz,field0,field1,reltpos,
     >           xind0,yind0,xind1,yind1,pind,mdv,ipomode,value)
C     ==========================================================
C
C     Decides what kind of interpolation should be done.
C
      integer   nx0,ny0,nx1,ny1,nz
      real      field0(nx0,ny0,nz)
      real      field1(nx1,ny1,nz)
      real      reltpos
      real      xind0,yind0,xind1,yind1,pind
      real      value,mdv
      integer   ipomode
 
      if (ipomode.eq.1) then
        call nlinipo(nx0,ny0,nx1,ny1,nz,field0,field1,reltpos,
     >               xind0,yind0,xind1,yind1,pind,mdv,value)
      else
        write(*,*)'*** error in SR ipo: unknown ipomode'
        stop
      endif
 
      return
      end
      subroutine linipo(nx,ny,nz,field0,field1,reltpos,
     >                  xind,yind,pind,mdv,value)
C     =================================================
C
C     Does linear interpolation both in time and space.
C
      integer   nx,ny,nz
      real 	field0(nx,ny,nz)
      real 	field1(nx,ny,nz)
      real 	reltpos
      real 	xind,yind,pind
      real	value,val0,val1,mdv
      real	int3dm,int2dm

      if (reltpos.eq.0.) then
        if (nz.eq.1) then
          value=int2dm(field0,nx,ny,xind,yind,mdv)
        else
          value=int3dm(field0,nx,ny,nz,xind,yind,pind,mdv)
        endif
      else if (reltpos.eq.1.) then
        if (nz.eq.1) then
          value=int2dm(field1,nx,ny,xind,yind,mdv)
        else
          value=int3dm(field1,nx,ny,nz,xind,yind,pind,mdv)
        endif
      else
        if (nz.eq.1) then
          val0=int2dm(field0,nx,ny,xind,yind,mdv)
          val1=int2dm(field1,nx,ny,xind,yind,mdv)
        else
          val0=int3dm(field0,nx,ny,nz,xind,yind,pind,mdv)
          val1=int3dm(field1,nx,ny,nz,xind,yind,pind,mdv)
        endif
        if ((val0.eq.mdv).or.(val1.eq.mdv)) then
          value=mdv
        else
          value=(1.-reltpos)*val0+reltpos*val1
        endif
      endif
      
      return
      end
      subroutine linipolog(nx,ny,nz,field0,field1,reltpos,
     >                     xind,yind,pind,mdv,value)
C     ====================================================
C
C     Does linear interpolation both in time and space.
C
      integer   nx,ny,nz
      real      field0(nx,ny,nz)
      real      field1(nx,ny,nz)
      real      reltpos
      real      xind,yind,pind
      real      value,val0,val1,mdv
      real      int3dmlog

      if (reltpos.eq.0.) then
        value=int3dmlog(field0,nx,ny,nz,xind,yind,pind,mdv)
      else if (reltpos.eq.1.) then
        value=int3dmlog(field1,nx,ny,nz,xind,yind,pind,mdv)
      else
        val0=int3dmlog(field0,nx,ny,nz,xind,yind,pind,mdv)
        val1=int3dmlog(field1,nx,ny,nz,xind,yind,pind,mdv)
        if ((val0.eq.mdv).or.(val1.eq.mdv)) then
          value=mdv
        else
          value=(1.-reltpos)*val0+reltpos*val1
        endif
      endif

      return
      end
      subroutine nlinipo(nx0,ny0,nx1,ny1,nz,field0,field1,reltpos,
     >                  xind0,yind0,xind1,yind1,pind,mdv,value)
C     ============================================================
C
C     Does linear interpolation both in time and space.
C
      integer   nx0,ny0,nx1,ny1,nz
      real      field0(nx0,ny0,nz)
      real      field1(nx1,ny1,nz)
      real      reltpos
      real      xind0,yind0,xind1,yind1,pind
      real      value,val0,val1,mdv
      real      int3dm
 
      if (reltpos.eq.0.) then
        value=int3dm(field0,nx0,ny0,nz,xind0,yind0,pind,mdv)
      else if (reltpos.eq.1.) then
        value=int3dm(field1,nx1,ny1,nz,xind1,yind1,pind,mdv)
      else
        val0=int3dm(field0,nx0,ny0,nz,xind0,yind0,pind,mdv)
        val1=int3dm(field1,nx1,ny1,nz,xind1,yind1,pind,mdv)
        value=(1.-reltpos)*val0+reltpos*val1
      endif
 
      return
      end
      subroutine vcuipo(nx,ny,nz,field0,field1,reltpos,
     >                  xind,yind,pind,mdv,value)
C     =================================================
C
C     Does linear interpolation both in time and space, except for a cubic
C     interpolation in the vertical.
C
      integer   nx,ny,nz
      real      field0(nx,ny,nz)
      real      field1(nx,ny,nz)
      real      reltpos
      real      xind,yind,pind
      real      value,val0,val1,mdv
      real      int3dmvc
 
      if (reltpos.eq.0.) then
        value=int3dmvc(field0,nx,ny,nz,xind,yind,pind,mdv)
      else if (reltpos.eq.1.) then
        value=int3dmvc(field1,nx,ny,nz,xind,yind,pind,mdv)
      else
        val0=int3dmvc(field0,nx,ny,nz,xind,yind,pind,mdv)
        val1=int3dmvc(field1,nx,ny,nz,xind,yind,pind,mdv)
        value=(1.-reltpos)*val0+reltpos*val1
      endif
 
      return
      end
      subroutine linipom(nx,ny,nz,field0,field1,reltpos,
     >                  mdv,xind,yind,pind,value)
C     =================================================
C
C     Does linear interpolation both in time and space.
C
      integer   nx,ny,nz
      real      field0(nx,ny,nz)
      real      field1(nx,ny,nz)
      real      reltpos,mdv
      real      xind,yind,pind
      real      value,val0,val1
      real      int3dm

      if (reltpos.eq.0.) then
        value=int3dm(field0,nx,ny,nz,xind,yind,pind,mdv)
      else if (reltpos.eq.1.) then
        value=int3dm(field1,nx,ny,nz,xind,yind,pind,mdv)
      else
        val0=int3dm(field0,nx,ny,nz,xind,yind,pind,mdv)
        val1=int3dm(field1,nx,ny,nz,xind,yind,pind,mdv)
        value=(1.-reltpos)*val0+reltpos*val1
      endif

      return
      end
      subroutine bicipo(nx,ny,nz,fld0,fld1,reltpos,
     >                  xind,yind,pind,value)
C     ===================================================
C
C     Does linear interpolation in time and bicubic interpolation
C     in space (horizontally)
C
      integer   nx,ny,nz
      real      fld0(nx,ny,nz),fld1(nx,ny,nz)
      real      reltpos,value
      real      xind,yind,pind

C     Declarations of local variables

      integer   k1,k2
      real      frac1k,frac2k
      real      val0,val1,val10,val20,val11,val21
      real	cint2d

      k1=int(pind)
      k2=k1+1

      frac1k=pind-aint(pind)
      frac2k=1.-frac1k

C     Interpolate value on lower level for first time

      if ((reltpos.ne.1.).and.(frac2k.ne.0.)) then
        val10=cint2d(fld0,nx,ny,nz,xind,yind,k1)
      endif

C     Interpolate value on upper level for first time

      if ((reltpos.ne.1.).and.(frac1k.ne.0.)) then
        val20=cint2d(fld0,nx,ny,nz,xind,yind,k2)
      endif

C     Interpolate value on lower level for second time

      if ((reltpos.ne.0.).and.(frac2k.ne.0.)) then
        val11=cint2d(fld1,nx,ny,nz,xind,yind,k1)
      endif

C     Interpolate value on upper level for second time

      if ((reltpos.ne.0.).and.(frac1k.ne.0.)) then
        val21=cint2d(fld1,nx,ny,nz,xind,yind,k2)
      endif

C     Do vertical interpolation and interpolation in time

      if (reltpos.eq.0.) then
        value=frac2k*val10+frac1k*val20
      else if (reltpos.eq.1.) then
        value=frac2k*val11+frac1k*val21
      else
        val0=frac2k*val10+frac1k*val20
        val1=frac2k*val11+frac1k*val21
        value=(1.-reltpos)*val0+reltpos*val1
      endif

      return
      end
