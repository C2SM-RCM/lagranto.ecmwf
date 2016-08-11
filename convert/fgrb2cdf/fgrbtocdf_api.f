      program grbtocdf
C
C     Purpose.
C     --------
C
C           Unpacks GRIB coded data (from the ECMWF) and
C           writes the data to a NetCDF-file.
C
C     Interface.
C     ----------
C
C           File of GRIB coded data attached as unit 11.
C	    NetCDF-filename and mdv read from unit 15.
C	    Varfile attached as unit 14.
C
C     Reference.
C     ----------
C
C           WMO Manual on Codes for GRIB definition.
C           WMO Publication No. 9, Volume B, for grid catalogue numbers.
C
C     Author.
C     -------
C
C	    H. Wernli		16.03.93
C	    Based on the programs grbtoxx.f by H.Wernli/D.Luethi and
C	    swmcdf.f by Ch.Schaer.
C
C     Modifications.
C     --------------
C
C           H. Wernli		6.06.97
C	    Implementation of options close and onepertime
C
C     -----------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'grib_api_f77.h'

C     Arrays are dimensioned to accommodate T213/N160 data volumes.

      integer   nxmax,nymax,nlevmax
      parameter(nxmax=2881,nymax=1441,nlevmax=300)

      integer	jfield
      PARAMETER(jfield=nxmax*nymax)
      
 
C     Maximum number of variables (in varfile)
      INTEGER   maxvar
      PARAMETER(maxvar=100)
 
C     Local variables declarations

      character*80 cdfnam,cdfname,varname,ifname,lastvar,var
      character*80 cfn,wrcstflag
      CHARACTER*10 vnshrt
      CHARACTER*2  string2,string3,string4
      CHARACTER*4  string1
      CHARACTER*10  datestr,timestr
      real	field(jfield),dummy(jfield),zsec4(jfield)
      real	pvarr(2*nlevmax+4),mdv,tstep
      real	varmin(4),varmax(4),stag(4),dx,dy,y0,y1,lat,xmin,xmax
      real	level(nlevmax),aklev(nlevmax),bklev(nlevmax),
     >		aklay(nlevmax),bklay(nlevmax)
      integer	idate(5),idate2(5),stdate(5),datar(14),vardim(4)
      integer	nd,ndim,nvar,nvars,dim,idefine
      integer	numerr,irec,idata,ipos,iendf,ilev,nlev,nloop
      integer	glev,nvals,nlev2
      integer	ierr,iword,idgrb,ifile,lenout
      INTEGER	i,j,k,ind,ind1,ind2,j0,n,ishift,gltyp,nn
      integer	nx,ny,npv
      integer	ifilen,ipo
      integer	cstid,cdfid,grbednr

      character*(15) vnam(maxvar)
      character*(10) snames(maxvar)
      character*(13) unit(maxvar)
      character*(15) varnam(maxvar)
      integer   lnum(maxvar),
     >          stg(maxvar),tdep(maxvar),p(maxvar),lval(maxvar)
      integer   arrcnt(maxvar),idef(maxvar),levty(maxvar)
      real      factor(maxvar),bias(maxvar)

      integer   varcnt,pos
      logical   infile,opnfile,prelev,periodic
      integer	lastlev,version,onepertime,rotate,noclose

      write(*,*)'*** start of program GRBTOCDF ***'

C     Initialize some variables
 
      numerr = 0	! clear error counter
      irec=0		! # of record on unit 11
      idata=0		! # of words (not yet used) in array inbuff
      ipos=0		! # of words (already used) in array inbuff
      iendf=0           ! flag to determine end of file on unit 11
      nlev=0		! # of levels on grib file
      nloop=0		! # of decoded fields
 
C     Lengths of INBUFF and record length on unit 11
 
      ifname = 'fort.11'

C     Read the NetCDF filename, constants filename and the missing
C     data flag from tape 15

      read(15,'(a)')cdfnam
      read(15,'(a)')wrcstflag
      read(15,'(a)')cfn
      read(15,*)mdv
      read(15,*)version
      read(15,*)rotate
      read(15,*)noclose
      read(15,*)onepertime

C     Read the varfile 

      call rdvarfile(vnam,snames,levty,unit,
     >              factor,bias,lnum,stg,tdep,p,lval,varcnt,ierr)
      if (ierr.ne.0) goto 994


C     ====================================================
C     Read the whole file in order to determine the levels
C     loop 50 --> goto 50
C     ====================================================

C     Initialize level array with dummy value
C     (otherwise ilev1=0 (surface data) will not be defined)
      DO i=1,nlevmax
        level(i)=-999.98999
      ENDDO
 
C     Open file with GRID coded data as direct access file

      CALL grib_check(grib_open_file(ifile, ifname, 'r'))
 
 50   CONTINUE
      nloop=nloop+1		! next field gets decoded

C     Get new record from GRIB file
      ierr = grib_new_from_file(ifile, idgrb)
 
      IF (idgrb.EQ.-1) THEN
        IF (ierr.NE.-1) THEN
          CALL grib_check(ierr)
          GOTO 9993
        ELSE
          GOTO 60     ! end of file is reached
        ENDIF
      ENDIF

C     Loop to get all levels for the NetCDF file

      CALL grib_check(grib_get_int(idgrb,'editionNumber',grbednr))
      CALL grib_check(grib_get_string(idgrb,'shortName',vnshrt))
      CALL grib_check(grib_get_int(idgrb,'level',glev))
      IF (grbednr.EQ.1) THEN
        CALL grib_check(grib_get_int(idgrb,'indicatorOfTypeOfLevel',
     &                        gltyp))
      ELSE
        ierr = grib_get_int(idgrb,'typeOfFirstFixedSurface',
     &                        gltyp)
        IF (gltyp.ne.105) gltyp=1
      ENDIF

C     Test if variable is on varfile
      IF (infile(vnshrt,gltyp,snames,levty,varcnt,pos)) THEN
C       Test if level is desired for 2d fields (quit loop if not)
C       Test if actual level (isec1(8)) is already defined (quit loop if yes)
        DO i=1,nlevmax
          IF (glev.EQ.level(i)) GOTO 24
        ENDDO
C       Define new level
        nlev=nlev+1
C       Special if to define the first level
        IF (level(1).EQ.-999.98999) THEN
          level(1)=glev
          GOTO 24
        ENDIF
C       Get position for new level in level-array such that levels
C       in the  array are monotonically decreasing
        DO i=1,nlevmax-1
          IF ((glev.LT.level(i)).AND.(glev.GT.level(i+1))) THEN
            DO j=nlevmax-1,i+1,-1
              level(j+1)=level(j)
            ENDDO
            level(i+1)=glev
            GOTO 24
          ENDIF
        ENDDO
        IF (glev.GT.level(1)) THEN
          DO j=nlevmax-1,1,-1
            level(j+1)=level(j)
          ENDDO
          level(1)=glev
          GOTO 24
        ENDIF
      ENDIF
   24 CONTINUE

      CALL grib_check(grib_release(idgrb))
   
      GOTO 50

   60 CONTINUE

      CALL grib_check(grib_close_file(ifile))

 
C     ======================================================
C     Start of loop to write specified fields on NetCDF file
C     loop 70 --> goto 70
C     ======================================================

C     Initialize some variables

      numerr = 0        ! clear error counter
      irec=1            ! # of record on unit 11
      idata=0           ! # of words (not yet used) in array inbuff
      ipos=0            ! # of words (already used) in array inbuff
      iendf=0           ! flag to determine end of file on unit 11
      nloop=0           ! # of decoded fields
      ilev=0            ! # of actual level in time step

      CALL grib_check(grib_open_file(ifile, ifname, 'r'))
      PRINT *,'opened file'

 65   ierr = grib_new_from_file(ifile, idgrb)
 
      IF (idgrb.EQ.-1) THEN
        IF (ierr.NE.-1) THEN
          CALL grib_check(ierr)
        ENDIF
        STOP 2        
      ENDIF

      CALL grib_check(grib_get_string(idgrb,'shortName',vnshrt))
      CALL grib_check(grib_get_int(idgrb,'level',glev))
      IF (grbednr.EQ.1) THEN
        CALL grib_check(grib_get_int(idgrb,'indicatorOfTypeOfLevel',
     &                        gltyp))
      ELSE
        ierr = grib_get_int(idgrb,'typeOfFirstFixedSurface',
     &                        gltyp)
        IF (gltyp.LT.100) gltyp=1
      ENDIF

      IF (.NOT. infile(vnshrt,gltyp,snames,levty,varcnt,pos)) THEN
        PRINT *,'get next record'
        CALL grib_check(grib_release(idgrb))
        GOTO 65
      ENDIF

C     Define horizontal dimensions

      CALL grib_check(grib_get_int(idgrb,'Ni',nx)) ! number of points along x
      CALL grib_check(grib_get_int(idgrb,'Nj',ny)) ! number of points along y

C     Decide if data is on pressure or model levels

      prelev=.true.
      IF (gltyp.NE.100) prelev=.false.

C     Define staggering coefficients

      stag(1)=0.
      stag(2)=0.
      stag(3)=-0.5

C     Define data region

      CALL grib_check(grib_get_real4(idgrb,
     &                'longitudeOfFirstGridPointInDegrees',varmin(1)))
      CALL grib_check(grib_get_real4(idgrb,
     &                'latitudeOfLastGridPointInDegrees',varmin(2)))
      CALL grib_check(grib_get_real4(idgrb,
     &                'longitudeOfLastGridPointInDegrees',varmax(1)))
      CALL grib_check(grib_get_real4(idgrb,
     &                'latitudeOfFirstGridPointInDegrees',varmax(2)))
      IF (varmin(1).GT.varmax(1)) THEN
        IF (varmin(1).GE.180.) THEN
          varmin(1) = varmin(1) - 360.
        ELSE
          varmax(1) = varmax(1) + 360.
        ENDIF
      ENDIF

      IF (prelev) THEN
        varmin(3)=level(1)
        varmax(3)=level(nlev)
      ELSE
        varmin(3)=1050.
        varmax(3)=0.
      ENDIF

C     Define dimensions

      vardim(1)=nx
      vardim(2)=ny

      IF ((lnum(pos).EQ.3).OR.(lnum(pos).EQ.4)) THEN
        vardim(3)=1
      ELSE
        vardim(3)=nlev
      ENDIF
      IF (tdep(pos).EQ.0) THEN
        vardim(4)=1
      ENDIF

C     Calculate grid increments

      dx=(varmax(1)-varmin(1))/real(vardim(1)-1)
      dy=(varmax(2)-varmin(2))/real(vardim(2)-1)

C     Check if the domain is (almost) periodic

C      print*,'test ',varmin(1),varmax(1),dx,varmax(1)-varmin(1),360.-dx
      IF (ABS((varmax(1)-varmin(1))-(360.-dx)).LT.0.1) THEN
        periodic=.TRUE.
        IF (noclose.EQ.0) THEN
          varmax(1)=varmin(1)+360.
          vardim(1)=vardim(1)+1
        ENDIF
      ENDIF
C      print*,'test ',varmin(1),varmax(1),dx,varmax(1)-varmin(1),360.-dx

      IF ((rotate.EQ.1).AND.(periodic)) THEN
        xmin=varmin(1)
        xmax=varmax(1)
        varmin(1)=0.
        varmax(1)=360.
      ENDIF

C     If one file per date then define the proper filename

C           Define actual date as starting date
      CALL grib_check(grib_get_string(idgrb,'dataDate',datestr))
      CALL grib_check(grib_get_string(idgrb,'dataTime',timestr))

      IF (onepertime.EQ.0) THEN
        cdfname=trim(cdfnam)
      ELSE
        READ (datestr(1:4),'(I4.4)') idate(1)
        READ (datestr(5:6),'(I2.2)') idate(2)
        READ (datestr(7:8),'(I2.2)') idate(3)
        READ (timestr(1:2),'(I2.2)') idate(4)

        write(string1,'(i4)')idate(1)
        if (idate(2).ge.10) then
          write(string2,'(i2)')idate(2)
        else
          write(string2,'(a1,i1)')'0',idate(2)
        endif
        if (idate(3).ge.10) then
          write(string3,'(i2)')idate(3)
        else
          write(string3,'(a1,i1)')'0',idate(3)
        endif
        if (idate(4).ge.10) then
          write(string4,'(i2)')idate(4)
        else
          write(string4,'(a1,i1)')'0',idate(4)
        endif

        cdfname=trim(cdfnam)//
     >          string1//string2//string3//'_'//string4

        cfn=TRIM(cdfname)//'_cst'
      ENDIF

C     Define the name of the constants file

      if (cfn.eq.'default') then
        cfn=cdfname//'_cst'
      endif

C     Open or create the NetCDF file, read or define the start date
      if (.not.opnfile) then
        call opncdf(cdfname,cdfid,varmin,varmax,nd,varnam,nvar,cfn,ierr)
        if (ierr.ne.0) then
          call crecdf(trim(cdfname),cdfid,varmin,varmax,
     >                3,trim(cfn),ierr)
          if (ierr.ne.0) goto 996
          write(*,*)'*** NetCDF file ',trim(cdfname),' created'

C         If requested write the constants file

          if (wrcstflag.eq.'yes') then

C           Define actual date as starting date
            CALL grib_check(grib_get_string(idgrb,'dataDate',datestr))
            CALL grib_check(grib_get_string(idgrb,'dataTime',timestr))
            IF (grbednr.EQ.1) THEN
              CALL grib_check(grib_get_int(idgrb,'startStep',
     &                                     stdate(5)))
            ELSE
              CALL grib_check(grib_get_int(idgrb,'forecastTime',
     &                                     stdate(5)))
            ENDIF

            READ (datestr(1:4),'(I4.4)') stdate(1)
            READ (datestr(5:6),'(I2.2)') stdate(2)
            READ (datestr(7:8),'(I2.2)') stdate(3)
            READ (timestr(1:2),'(I2.2)') stdate(4)

C           Define the datar array (used to define constants file)

            datar(1)=vardim(1)
            datar(2)=vardim(2)
            datar(3)=nint(1000.*varmax(2))      ! latitude north
            datar(4)=nint(1000.*varmin(1))      ! longitude west
            datar(5)=nint(1000.*varmin(2))      ! latitude south
            datar(6)=nint(1000.*varmax(1))      ! longitude east
            datar(7)=nint(1000.*dx)
            datar(8)=nint(1000.*dy)
            datar(9)=vardim(3)
            datar(10)=2                 ! dattyp = forecast
            datar(11)=0                 ! data version
            datar(12)=0                 ! cstfile version
            datar(13)=0                 ! longitude of pole
            datar(14)=90000             ! latitude of pole

C           Get the 4 levels/layers arrays

            IF (prelev) then
              call prelevs(nlev,level,aklev,bklev,aklay,bklay)
            ELSE
              if (nlev.eq.1) nlev=60
              nlev2=nlev
C             nlev2 is quick&dirty for only 36 levels from 60 levels
              if (nlev.eq.36) nlev2=60
              if (nlev.eq.46) nlev2=60
*             call modlevs(nlev2,aklev,bklev,aklay,bklay)
              
              CALL grib_check(grib_get_size(idgrb,'pv',npv))
              CALL grib_check(grib_get_real4_array(idgrb,'pv',
     &                                                pvarr,npv))
              DO i=1,nlev
                aklev(i) = pvarr(npv/2 - i)/100.
                bklev(i) = pvarr(npv - i)
              ENDDO
              DO i=1,nlev-1
                aklay(i+1) = 0.5*(aklev(i)+aklev(i+1))
                bklay(i+1) = 0.5*(bklev(i)+bklev(i+1))
              ENDDO
              aklay(1)=0.5*(0. + aklev(1))
              bklay(1)=0.5*(1. + bklev(1))
            ENDIF

C           Write the constants file

            call wricst(cfn,datar,aklev,bklev,aklay,bklay,stdate)
            write(*,*)'*** cst-file ',trim(cfn),' created'

          else          ! assume that the constants-file already exists

C           Inquire start time

            call cdfopn(cfn,cstid,ierr)
            call getstart(cstid,stdate,ierr)
            call clscdf(cstid,ierr)
            nlev2=nlev
C           nlev2 is quick&dirty for only 36 levels from 60 levels
            if (nlev.eq.36) nlev2=60
            if (nlev.eq.46) nlev2=60
          endif
          idef=0
        else
          write(*,*)'*** NetCDF file ',trim(cdfname),' opened'

C         Inquire start time

          call cdfopn(cfn,cstid,ierr)
          call getstart(cstid,stdate,ierr)
          call clscdf(cstid,ierr)

          nlev2=nlev
C         nlev2 is quick&dirty for only 36 levels from 60 levels
          if (nlev.eq.36) nlev2=60
          if (nlev.eq.46) nlev2=60

        endif
*       print*,'ilev nlev ',ilev,nlev
        opnfile=.true.
        idef=0
      endif
*     print*,idef

C     Read file with predetermined variables

      read(16,*)nvars
      if (nvars.gt.0) then
        do n=1,nvars
          read(16,*)var
          read(16,*)dim
          if (dim.eq.3) then
            vardim(3)=nlev
          else
            vardim(3)=1
          endif
          call putdef(cdfid,trim(var),
     >                4,mdv,vardim,varmin,varmax,stag,ierr)
          write(*,*)'*** variable ',trim(var),' created'
        enddo
        idefine=0
      else
        idefine=1
      endif

   70 continue

      nloop=nloop+1		! next field gets decoded

      
C     Read next field if the actual one is not on the varfile

      CALL grib_check(grib_get_string(idgrb,'shortName',vnshrt))
      IF (grbednr.EQ.1) THEN
        CALL grib_check(grib_get_int(idgrb,'indicatorOfTypeOfLevel',
     &                        gltyp))
      ELSE
        ierr = grib_get_int(idgrb,'typeOfFirstFixedSurface',
     &                        gltyp)
        IF (gltyp.LT.100) gltyp=1
      ENDIF
      if (.not.infile(vnshrt,gltyp,snames,levty,varcnt,pos))
     >  goto 80
      CALL grib_check(grib_get_int(idgrb,'level',glev))
 
C     Don't write field if level is not desired

      if ((lnum(pos).eq.3).or.(lnum(pos).eq.4)) then
        if (glev.ne.lval(pos)) goto 80
      endif

C     Don't write field if variable is not time-dependent and it
C     appears not the first time

      IF ((tdep(pos).EQ.0).AND.(arrcnt(pos).NE.0)) GOTO 80

C     Define the array containing the date information

      CALL grib_check(grib_get_string(idgrb,'dataDate',datestr))
      CALL grib_check(grib_get_string(idgrb,'dataTime',timestr))
      IF (grbednr.EQ.1) THEN
        CALL grib_check(grib_get_int(idgrb,'startStep',
     &                                     idate(5)))
      ELSE
        CALL grib_check(grib_get_int(idgrb,'forecastTime',
     &                                     idate(5)))
      ENDIF

      READ (datestr(1:4),'(I4.4)') idate(1)
      READ (datestr(5:6),'(I2.2)') idate(2)
      READ (datestr(7:8),'(I2.2)') idate(3)
      READ (timestr(1:2),'(I2.2)') idate(4)

C     Convert forecast step into normal date

*     print*,idate(1),idate(2),idate(3),idate(4),idate(5)
      call newdate(idate,0.,idate2)
      idate(1)=idate2(1)
      idate(2)=idate2(2)
      idate(3)=idate2(3)
      idate(4)=idate2(4)
      idate(5)=idate2(5)
*     print*,idate(1),idate(2),idate(3),idate(4),idate(5)

C     Get index of actual level

      if ((lnum(pos).eq.3).or.(lnum(pos).eq.4)) then
        ilev=1
      else
        do k=1,nlev
          if (glev.eq.level(k)) then
            ilev=k
            goto 22
          endif
        enddo
        goto 995        ! actual level glev is not in array level
  22    continue
      endif

C     Set variable name and dimension
      CALL grib_check(grib_get_size(idgrb,'values',nvals))
      CALL grib_check(grib_get_real4_array(idgrb,'values',zsec4,
     &                                                 nvals))

      print*,'nvals,nx*ny ',nvals,nx*ny

      varname=vnam(pos)
      ndim=4

C     Get exponential of LNSP

      if (varname.eq.'LNSP') then
        varname='PS'
        do i=1,nx
        do j=1,ny
          zsec4(i+(j-1)*nx)=exp(zsec4(i+(j-1)*nx))
        enddo
        enddo
      endif

C     Turn the array and scale/shift the data with factor/bias

      do i=1,nx
      do j=1,ny
        field(i+(j-1)*nx)=(zsec4(i+(ny-j)*nx)+bias(pos))*factor(pos)
      enddo
      enddo

C     Close periodic arrays
     
      print*,'periodi,noclose ',periodic,noclose

      if (periodic.and.(noclose.eq.0)) then

        print*,'Hallo'

        if (rotate.eq.1) then
          ishift=nint(-xmin/dx)
          do j=1,ny
            do i=1+ishift,nx
              ind1=i+(j-1)*nx
              ind2=i-ishift+(j-1)*nx
              dummy(ind2)=field(ind1)
            enddo
            do i=1,ishift
              ind1=i+(j-1)*nx
              ind2=nx-ishift+i+(j-1)*nx
              dummy(ind2)=field(ind1)
            enddo
          enddo

          do i=1,nx*ny
            field(i)=dummy(i)
          enddo
        endif 

        do j=1,ny
          do i=1,nx
            ind1=i+(j-1)*(nx+1)
            ind2=i+(j-1)*nx
            dummy(ind1)=field(ind2)
          enddo
          dummy((nx+1)+(j-1)*(nx+1))=field(1+(j-1)*nx)
        enddo
 
        do i=1,(nx+1)*ny
          field(i)=dummy(i)
        enddo

      endif
 
*     print*,cdfname

C     Calculate tstep for actual date

      if (onepertime.eq.1) then
        tstep=0.
      else
        call timediff(idate,stdate,tstep)
*       print*,idate,stdate,tstep
      endif

C     If necessary define new variable
      if ((onepertime.eq.1).or.
     >   ((arrcnt(pos).eq.0).and.(tdep(pos).eq.1))) then
        do i=1,maxvar
          if (onepertime.eq.0) then
            if (trim(varname).eq.trim(varnam(i))) goto 40
          endif
        enddo
*       print*,'ilev level ',ilev,level(ilev),nlev,prelev
*       print*,pos,idef(pos)
        if ((idefine.eq.1).and.(idef(pos).eq.0)) then
          call putdef(cdfid,trim(varname),
     >                ndim,mdv,vardim,varmin,varmax,stag,ierr)
          write(*,*)'*** variable ',trim(vnam(pos)),' created'
          if (vardim(3).gt.1) idef(pos)=1
          lastvar=varname
        endif
      endif
  40  continue
      arrcnt(pos)=arrcnt(pos)+1		! counter for actual variable

C     Put data on file

      if (onepertime.eq.0) then
        IF ((glev.LT.lastlev).AND.(.NOT.prelev)) THEN
C         if is quick&dirty for only 36 levels from 60 levels
          IF ((lastlev.NE.36).AND.(lastlev.NE.46).AND.
     >        (lastlev.NE.nlev)) THEN
            WRITE(*,*)'unusual bug exit'
            GOTO 9999
          ENDIF
        ENDIF
      ENDIF

      write(*,41)'*** var ',trim(vnam(pos)),
     >           ' at time ',tstep,' on level ',glev,
     >           ' written on cdffile'
  41  format(a,a,a,f6.0,a,i4,a)
      if (tdep(pos).eq.0) then
        call putcdf(cdfid,trim(varname),ndim,
     >              mdv,vardim,varmin,varmax,stag,field,ierr)
      else
         print*,cdfid
         print*,trim(varname)
         print*,tstep
         print*,ilev
         do j=1,10
            print*,j,field(j*nvals/10)
        enddo
        call putdat(cdfid,trim(varname),
     >              tstep,ilev,field,ierr)
      endif
      lastlev=glev	! bug correction q+d 17.Dec.93
      
C     Close NetCDF file
      if (prelev) then
        if ((onepertime.eq.1).and.(ilev.eq.nlev)) then
          call clscdf(cdfid,ierr)
          opnfile=.false.
        endif
      else
*       print*,ilev,nlev,level(ilev),nlev2
*       print*,ilev,varname,lastvar
        if ((onepertime.eq.1).and.(varname.eq.lastvar).and.
     >      ((level(ilev).eq.real(nlev2)).or.(level(ilev).eq.0.))) then
*         print*,'closing file'
          call clscdf(cdfid,ierr)
*         print*,'file closed'
          opnfile=.false.
        endif
      endif

   80 continue

C     Get new record from GRIB file
      ierr = grib_new_from_file(ifile, idgrb)
 
      IF (idgrb.EQ.-1) THEN
        IF (ierr.NE.-1) THEN
          CALL grib_check(ierr)
          GOTO 9993
        ELSE
          GOTO 9999  ! end of file is reached
        ENDIF
      ELSE 
        GOTO 70
      ENDIF

  994 stop 'error from subroutine rvarfile: not proper varfile'
  995 stop 'error: actual level not defined on constants file'
  996 stop 'error: could not open NetCDF file'
  999 stop 'date not found on tape 9'

 9993 stop 'error from subroutine grib_new_from_file'

 9999 continue

C     Close NetCDF file
      if (onepertime.eq.0) call clscdf(cdfid,ierr)
      if (ierr.ne.0) stop 'error when closing file'

      WRITE(*,*)'*** end of program FGRBTOCDF status: NORMAL ***'
      END

      LOGICAL FUNCTION infile (vnshrt,gltyp,snames,levtyp,varcnt,pos)
c     ================================================================
c     testet, ob die Kombination (ltyp,ivar) in den Arrays levtyp und
c     gribnr vorkommen. varcnt ist die Anzahl Eintraege in den Arrays.
c     pos gibt die Array-Position des gefundenen Eintrages zurueck.

      IMPLICIT NONE
      INTEGER   maxvar
      PARAMETER (maxvar=100)

      INTEGER       gltyp,levtyp(maxvar),varcnt,pos
      CHARACTER*10  vnshrt,snames(maxvar)
      INTEGER       i, ltyp, ivar
      i = 1
300   IF (i.LE.varcnt) THEN
C        IF ((TRIM(snames(i)).EQ.TRIM(vnshrt)).AND.
C     &                                 (levtyp(i).EQ.gltyp)) THEN
        IF (TRIM(snames(i)).EQ.TRIM(vnshrt)) THEN
          pos = i
          infile = .TRUE.
          RETURN
        ENDIF
        i = i+1
      ELSE
        GOTO 301
      ENDIF
      GOTO 300
301   infile = .FALSE.
      END

      SUBROUTINE rdvarfile(vnam,snames,levty,unit,factor,bias,
     >                    lnum,stg,tdep,p,lval,varcnt,ierr)
C     =======================================================
C     Variablen-File in Arrays einlesen

      IMPLICIT NONE
      INTEGER   maxvar
      PARAMETER(maxvar=100)

      CHARACTER*(15) vnam(maxvar)
      CHARACTER*(10) snames(maxvar)
      
      CHARACTER*(13) unit(maxvar)
      CHARACTER*(1)  flag
      INTEGER   levty(maxvar),lnum(maxvar),
     >          stg(maxvar),tdep(maxvar),p(maxvar),lval(maxvar)
      REAL      factor(maxvar),bias(maxvar)

      INTEGER   i,varcnt,ierr,nt

      nt=14             ! number of tape
      i=1               ! initialize var-counter

C     Read first character of row and decide if it is comment or not

  100 READ(nt,10,ERR=123,END=126) flag
      IF (flag.EQ."#") GOTO 100         ! don't bother about comments
      BACKSPACE nt
  121 READ(nt,122, ERR=123, END=126) vnam(i), snames(i), levty(i),
     & unit(i), factor(i), bias(i), lnum(i), stg(i), tdep(i), p(i),
     & lval(i)
      i=i+1
      GOTO 100


   10 FORMAT(a1)
  122 FORMAT(A14,A10,1X,I3,4X,A13,F7.5,2X,F7.2,I7,I4,I6,I3,I5)
  123 GOTO 121
  126 CONTINUE
      varcnt=i-1        ! # of variables in varfile_i

C     Check some things

      ierr=0            ! initialize error flag
      DO i=1,varcnt
        IF ((lnum(i).NE.1).AND.(lnum(i).NE.2).AND.(lnum(i).NE.3)
     >      .AND.(lnum(i).NE.4)) ierr=1
        IF ((stg(i).NE.0).AND.(stg(i).NE.1).AND.(stg(i).NE.10).AND.
     >      (stg(i).NE.11)) ierr=1
        IF ((tdep(i).NE.0).AND.(tdep(i).NE.1)) ierr=1
        IF ((p(i).NE.0).AND.(p(i).NE.1).AND.(p(i).NE.2)) ierr=1
        IF ((lval(i).LT.0).OR.(lval(i).GT.1050)) ierr=1
      ENDDO

      RETURN
      END
