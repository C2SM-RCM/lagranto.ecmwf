      program ptos

C     ******************************************************************
C
C     NAME:
C     p2s
C
C     Purpose
C     -------
C
C     	Calculates secondary data files from primary data files
C     	(based upon IVE-routines).
C
C     calling:
C     p2s [-m] file variable-list [-s] [-o]
C
C     example:
C     p2s P911201_00 TH PV CH QX
C
C     p2s -m #returns man-page
C
C     Authors
C     ------
C
C     	H. Wernli	April 96
C       D.N. Bresch     980311
C
C     Modifications
C     -------------
C       completely rewritten by D.N. Bresch 980311 for F90 
C       and many more variables added (nearly all available in IVE)
C       
C     ADD YOUR OWN VARIABLES AT ALL PLACES WITH (++)
C     for easy simple calculations, see VEL, M, B
C     for complicated calculations, see VORT, QX and PVR
C
C     Remarks
C     -------
C     ZB is only read once when needed, i.e. for first time...
C
C-    ******************************************************************

      integer,parameter :: ntmax=200,nzmax=200,nvarmax=100
      real	time(ntmax),time2(ntmax)
      REAL,ALLOCATABLE, DIMENSION (:,:) :: sp,cl,tl,f,zb,t2m,td2m,vip,
     >     u10m,v10m,oro,gradpv
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: var,th,pv,lpv,the,rh,dhr,
     >     tt,qq,uu,vv,ww,rho,alpha,zz,mm,zlay,ug,vg,fl,ipv
      character*80 cdfnam,cstnam,outfnam
      integer	cdfid,cdfid1,cstid,ierr,ndim,vardim(4),stat
      integer	cdfid2,vardim2(4)
      real	dx,dy,mdv,varmin(4),varmax(4),stag(4)
      real      aklev(nzmax),bklev(nzmax),aklay(nzmax),bklay(nzmax),
     >		ak(nzmax),bk(nzmax)
      integer	nx,ny,nz,ntimes,ntimes2,i,j,k,n
      integer	stdate(5)
      character*(80) qmode,arg,vnam(nvarmax)
      integer	mode,zdef
      real	rlat,rlon,lat
      real	pollon,pollat,yphys
      real	phstoph
      logical	prelev

      real,parameter ::  pi=3.141592654

      integer   PS_out,TH_out,TH_calc,RH_out,RH_calc
      integer   PV_out,PV_calc,THE_out,THE_calc,VIP_out,VIP_calc
      integer   GRADPV_out,GRADPV_calc
      integer	LPV_out,LPV_calc
      integer   CH_out,CH_calc,PVR_out,PVR_calc
      integer   THW_out,THW_calc,Z_outP,Z_out,Z_calc
      integer	DIVQU_out,DIVQU_calc
      integer   NSQ_out,NSQ_calc,RHO_out,RHO_calc,ALPHA_out,ALPHA_calc,VEL_out,VEL_calc
      integer   NSQM_out,NSQM_calc,W_out,W_calc,M_out,M_calc
      integer   VORT_out,VORT_calc,UG_out,UG_calc,VG_out,VG_calc
      integer	AVO_out,AVO_calc,CURVO_out,CURVO_calc
      integer	DTHDP_out,DTHDP_calc
      integer	COS_out
      integer   ZLAY_out,ZLAY_calc,UA_out,UA_calc,VA_out,VA_calc
      integer   P_out,P_calc,PLEV_out,PLEV_calc
      integer   QXF_out,QXF_calc,QYF_out,QYF_calc
      integer   QX_out,QX_calc,QY_out,QY_calc,PSRED_calc,PSRED_out
      integer	RI_calc,RI_out,BLH_calc,BLH_out
      integer   GRADTH_out,GRADTH_calc,B_out,B_calc !(++)
      logical   verbose
      logical   ZonP,TonP,PSonP,UonP,VonP,OMEGAonP,ZBonP,ZBneed
      logical	T2MonP,TD2MonP,U10MonP,V10MonP,PVonP,PV3onP
      character*(80) zbfile

c     set defaults:
      verbose=.false.
      ZonP=.false.
      TonP=.false.
      PSonP=.false.
      ZBonP=.false.
      UonP=.false.
      VonP=.false.
      OMEGAonP=.false.
      ZBonP=.false.
      T2MonP=.false.
      TD2MonP=.false.
      U10MonP=.false.
      V10MonP=.false.
      PVonP=.false.
      PV3onP=.false.	! Lukas
      ZBneed=.false.
      zbfile=''

      qmode='QNone'
      PS_out=1
      TH_out=0
      TH_calc=0
      RH_out=0
      RH_calc=0
      PV_out=0
      PV_calc=0
      LPV_out=0
      LPV_calc=0
      VIP_out=0
      VIP_calc=0
      THE_out=0
      THE_calc=0
      CH_out=0
      CH_calc=0
      PVR_out=0
      PVR_calc=0
      THW_out=0
      THW_calc=0
      DIVQU_out=0
      DIVQU_calc=0
      Z_calc=0
      Z_out=0
      Z_outP=0
      NSQ_out=0
      NSQ_calc=0
      DTHDP_out=0
      DTHDP_calc=0
      NSQM_out=0
      NSQM_calc=0
      RHO_out=0
      RHO_calc=0
      ALPHA_out=0
      ALPHA_calc=0
      VEL_out=0
      VEL_calc=0
      M_out=0
      M_calc=0
      B_out=0
      B_calc=0
      W_out=0
      W_calc=0
      VORT_out=0
      VORT_calc=0
      AVO_out=0
      AVO_calc=0
      CURVO_out=0
      CURVO_calc=0
      COS_out=0
      UG_out=0
      UG_calc=0
      VG_out=0
      VG_calc=0
      UA_out=0
      UA_calc=0
      VA_out=0
      VA_calc=0
      ZLAY_out=0
      ZLAY_calc=0
      P_out=0
      P_calc=0
      PLEV_out=0
      PLEV_calc=0
      FL_out=0
      FL_calc=0
      PSRED_out=0
      PSRED_calc=0
      RI_out=0
      RI_calc=0
      BLH_out=0
      BLH_calc=0
      QXF_out=0
      QXF_calc=0
      QYF_out=0
      QYF_calc=0
      QX_out=0
      QX_calc=0
      QY_out=0
      QY_calc=0
      GRADTH_out=0
      GRADTH_calc=0
      GRADPV_out=0
      GRADPV_calc=0             !(++)

c     get arguments:
c     get parameters from command-line:
c     COUNT THE ARGUMENTS:
      if (iargc() .lt. 1) then
         print*,'USAGE: p2s [-m] file variable-list [-s] ',
     >        '[-o] [-zb file]'
         STOP
      endif
      
c     REQUESTD INPUT:
c     ---------------
c     GET WITH getarg DIRECTLY FROM SHELL:
      call getarg(1,cdfnam)
      if (trim(cdfnam).eq.'-m') then
         print*,' '
         print*,'computes derived variables from primary ones on'
         print*,'the input-file'
         print*,'if the output-file is already present, it will be'
         print*,'updated'
         print*,' '
         print*,'check the source-code itself about details...'
         print*,' '
         print*,'file: netCDF file with basic variables on it'
         print*,'    requested are the variables needed to calculate'
         print*,'    the requested output (variable-list)'
         print*,'    If file is a P-file (starting with P), the output-'
         print*,'    file will be an S-file, otherwise the extension'
         print*,'    _out will be appended, unless -o is used'
         print*,'    if the S-file exists already, it is tried to '
         print*,'    append the new variable'
         print*,'[P]date means that you can either give PYYMMDD_HH'
         print*,'    or YYMMDD_HH alone'
         print*,'variable-list: a list of variables to be calculated'
         print*,'    and written to the S_file, available are:'
         print*,'    TH,PV,LPV,RH,THE,THW,CH,PVR,Z,ZonP,GRADTH,NSQ,NSQM'
         print*,'    M,B,W,RHO,VEL,VORT,AVO,CURVO,UG,VG,ZLAY,UA,VA,P'
         print*,'    PLEV,FL,QX,QY,QXF,QYF,PSRED,RI,BLH,GRADPV,DTHDP'
         print*,'    ALPHA'
         print*,'-s: only small S-file, i.e. TH and PV'
         print*,'-zb: file with ZB (for PSRED)'
         print*,'-o output: filename of the output netCDF file'
         STOP
      endif
C     check, if cdfnam is with or without P:
      if (cdfnam(1:1).eq.'P') then
         outfnam='S'//cdfnam(2:len_trim(cdfnam))
      else        
         outfnam=trim(cdfnam)//'_out'
      endif
      i=2
      do while (iargc().ge.i)
         call getarg(i,arg)
         i=i+1
         if (arg.eq.'TH') TH_out=1
         if (arg.eq.'THE') THE_out=1
         if (arg.eq.'THW') THW_out=1
         if (arg.eq.'DIVQU') DIVQU_out=1
         if (arg.eq.'PV') PV_out=1
         if (arg.eq.'LPV') LPV_out=1
         if (arg.eq.'VIP') VIP_out=1
         if (arg.eq.'RH') RH_out=1
         if (arg.eq.'CH') CH_out=1
         if (arg.eq.'PVR') PVR_out=1
         if (arg.eq.'Z') Z_out=1
         if (arg.eq.'ZonP') Z_outP=1
         if (arg.eq.'GRADTH') GRADTH_out=1
         if (arg.eq.'GRADPV') GRADPV_out=1
         if (arg.eq.'VEL') VEL_out=1
         if (arg.eq.'RHO') RHO_out=1
         if (arg.eq.'ALPHA') ALPHA_out=1
         if (arg.eq.'W') W_out=1
         if (arg.eq.'M') M_out=1
         if (arg.eq.'B') B_out=1
         if (arg.eq.'VORT') VORT_out=1
         if (arg.eq.'AVO') AVO_out=1
         if (arg.eq.'CURVO') CURVO_out=1
         if (arg.eq.'COS') COS_out=1
         if (arg.eq.'UG') UG_out=1
         if (arg.eq.'VG') VG_out=1         
         if (arg.eq.'UA') UA_out=1
         if (arg.eq.'VA') VA_out=1         
         if (arg.eq.'QX') QX_out=1         
         if (arg.eq.'QY') QY_out=1         
         if (arg.eq.'QXF') QXF_out=1         
         if (arg.eq.'QYF') QYF_out=1     
         if (arg.eq.'ZLAY') ZLAY_out=1
         if (arg.eq.'P') P_out=1
         if (arg.eq.'PLEV') PLEV_out=1
         if (arg.eq.'FL') FL_out=1
         if (arg.eq.'PSRED') PSRED_out=1
         if (arg.eq.'RI') RI_out=1
         if (arg.eq.'BLH') BLH_out=1
         if (arg.eq.'NSQ') NSQ_out=1
         if (arg.eq.'DTHDP') DTHDP_out=1
         if (arg.eq.'NSQM') NSQM_out=1 ! (++)      
         if (arg.eq.'-v') verbose=.true.
         if (arg.eq.'-s') then
            TH_out=1
            PV_out=1
         endif
         if (arg.eq.'-o') then
            if (iargc().ge.i) then          
               call getarg(i,outfnam)
               i=i+1
            else
               print*,'option -o requires filename'
               STOP
            endif
         endif
         if (arg.eq.'-zb') then
            if (iargc().ge.i) then          
               call getarg(i,zbfile)
               ZBneed=.true.
               i=i+1
            else
               print*,'option -zb requires filename'
               STOP
            endif
         endif
         if (arg.eq.'-nops') then
           PS_out=0
         endif
      enddo

C     force calculations of requested fields:
c     ---------------------------------------
      if (QX_out.eq.1) QX_calc=1
      if (QY_out.eq.1) QY_calc=1
      if (QXF_out.eq.1) QXF_calc=1
      if (QYF_out.eq.1) QYF_calc=1
      if (ZLAY_out.eq.1) ZLAY_calc=1
      if (P_out.eq.1) P_calc=1
      if (PLEV_out.eq.1) PLEV_calc=1
      if (FL_out.eq.1) FL_calc=1
      if (UG_out.eq.1) UG_calc=1
      if (VG_out.eq.1) VG_calc=1
      if (UA_out.eq.1) UA_calc=1
      if (VA_out.eq.1) VA_calc=1
      if (VORT_out.eq.1) VORT_calc=1
      if (AVO_out.eq.1) AVO_calc=1
      if (CURVO_out.eq.1) CURVO_calc=1
      if (TH_out.eq.1) TH_calc=1
      if (PV_out.eq.1) PV_calc=1
      if (LPV_out.eq.1) LPV_calc=1
      if (VIP_out.eq.1) VIP_calc=1
      if (RH_out.eq.1) RH_calc=1
      if (THE_out.eq.1) THE_calc=1
      if (THW_out.eq.1) THW_calc=1
      if (DIVQU_out.eq.1) DIVQU_calc=1
      if (CH_out.eq.1) CH_calc=1
      if (PVR_out.eq.1) PVR_calc=1
      if (Z_out.eq.1) Z_calc=1
      if (Z_outP.eq.1) Z_calc=1
      if (GRADTH_out.eq.1) GRADTH_calc=1
      if (GRADPV_out.eq.1) GRADPV_calc=1
      if (RHO_out.eq.1) RHO_calc=1
      if (ALPHA_out.eq.1) ALPHA_calc=1
      if (RI_out.eq.1) RI_calc=1
      if (BLH_out.eq.1) BLH_calc=1
      if (VEL_out.eq.1) VEL_calc=1
      if (M_out.eq.1) M_calc=1
      if (B_out.eq.1) B_calc=1
      if (NSQM_out.eq.1) NSQM_calc=1
      if (W_out.eq.1) W_calc=1
      if (DTHDP_out.eq.1) DTHDP_calc=1
      if (NSQ_out.eq.1) NSQ_calc=1 !(++)

C     make dependencies for variable calculations:
c     --------------------------------------------
      if (NSQ_calc.eq.1) TH_calc=1 !(++)
      if (DTHDP_calc.eq.1) TH_calc=1
      if (PSRED_out.eq.1) PSRED_calc=1
      if (PSRED_calc.eq.1) ZBneed=.true.
      if (QX_out.eq.1) UG_calc=1
      if (QY_out.eq.1) VG_calc=1
      if (UA_calc.eq.1) UG_calc=1
      if (VA_calc.eq.1) VG_calc=1
      if (UG_calc.eq.1) ZLAY_calc=1
      if (VG_calc.eq.1) ZLAY_calc=1
      if (ZLAY_calc.eq.1) Z_calc=1
      if (B_calc.eq.1) M_calc=1
      if (M_calc.eq.1) ZLAY_calc=1
      if (NSQ_calc.eq.1) RHO_calc=1
      if (NSQM_calc.eq.1) RHO_calc=1
      if (NSQM_calc.eq.1) THE_calc=1
      if (W_calc.eq.1) RHO_calc=1
      if (GRADTH_calc.eq.1) TH_calc=1
      if (THW_calc.eq.1) THE_calc=1
      if (THE_calc.eq.1) TH_calc=1
      if (VIP_calc.eq.1) PV_calc=1
      if (PV_calc.eq.1) TH_calc=1
      if (LPV_calc.eq.1) PV_calc=1
      if (PVR_calc.eq.1) CH_calc=1
      if (CH_calc.eq.1) TH_calc=1
      if (CH_calc.eq.1) RH_calc=1
      if (RI_calc.eq.1) TH_calc=1
      if (ALPHA_calc.eq.1) RHO_calc=1
      
      print*,'processing: ',trim(cdfnam)
      
C     Open files and get infos about data domain
c     ------------------------------------------
      
      if (Z_outP.eq.1) then
         call cdfwopn(trim(cdfnam),cdfid1,ierr)
      else
         call cdfopn(trim(cdfnam),cdfid1,ierr)
      endif
      if (ierr.ne.0) then
         print*,'ERROR opening input file, stopped'
         stop
      endif
      call getcfn(cdfid1,cstnam,ierr)
      if (ierr /= 0) then
         cstnam = cdfnam
         cstid  = cdfid1
      else
         call cdfopn(trim(cstnam),cstid,ierr)
      endif
     
 
C     Inquire the variables on the netCDF file:
c     -----------------------------------------
C     Inquire number of variables and variable names
      call getvars(cdfid1,nvars,vnam,ierr)
C      print*,nvars,vnam(1),vnam(2)
      do i=1,nvars
         if (trim(vnam(i)).eq.'Q') qmode='Q'
         if (trim(vnam(i)).eq.'QD') qmode='QD'
         if (trim(vnam(i)).eq.'Z') ZonP=.true.
         if (trim(vnam(i)).eq.'T') TonP=.true.
         if (trim(vnam(i)).eq.'PS') PSonP=.true.
         if (trim(vnam(i)).eq.'U') UonP=.true.
         if (trim(vnam(i)).eq.'V') VonP=.true.
         if (trim(vnam(i)).eq.'ZB') ZBonP=.true.
         if (trim(vnam(i)).eq.'T2M') T2MonP=.true.
         if (trim(vnam(i)).eq.'TD2M') TD2MonP=.true.
         if (trim(vnam(i)).eq.'U10M') U10MonP=.true.
         if (trim(vnam(i)).eq.'V10M') V10MonP=.true.
         if (trim(vnam(i)).eq.'OMEGA') OMEGAonP=.true.
         if (trim(vnam(i)).eq.'PV') PVonP=.true.
         if (trim(vnam(i)).eq.'PV3') PV3onP=.true.
      enddo

     
      if (TonP) then
        call getdef(cdfid1,'T',ndim,mdv,vardim,varmin,varmax,stag,ierr)
      else
        call getdef(cdfid1,vnam(2),ndim,mdv,
     >        vardim,varmin,varmax,stag,ierr)
      endif
      if (ierr.ne.0) goto 920
      mdv=-999.98999
      
C     Get the levels, pole, etc.
c     --------------------------
      
      nx=vardim(1)
      ny=vardim(2)
      nz=vardim(3)
      
C      print*,nx,ny,nz
      call getgrid(cstid,dx,dy,ierr)
      call getlevs(cstid,nz,aklev,bklev,aklay,bklay,ierr)
      call getpole(cstid,pollon,pollat,ierr)
      call getstart(cstid,stdate,ierr)

C     Allocate all arrays:
C     --------------------
      allocate(sp(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating sp(nx,ny)'
      allocate(oro(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating oro(nx,ny)'
      allocate(cl(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating cl(nx,ny)'
      allocate(tl(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating tl(nx,ny)'
      allocate(f(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating f(nx,ny)'

c     allocate memory for results
      allocate(var(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating var(nx,ny,nz)'

c     allocate memory for variables on P-file (allocate all)
      allocate(tt(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating tt(nx,ny,nz)'
      allocate(qq(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating qq(nx,ny,nz)'
      allocate(uu(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating uu(nx,ny,nz)'
      allocate(vv(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating vv(nx,ny,nz)'
      allocate(ww(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating ww(nx,ny,nz)'

      allocate(t2m(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating t2m(nx,ny)'
      allocate(td2m(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating td2m(nx,ny)'
      allocate(u10m(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating u10m(nx,ny)'
      allocate(v10m(nx,ny),STAT=stat)
      if (stat.ne.0) print*,'error allocating v10m(nx,ny)'
      allocate(ipv(nx,ny,nz),STAT=stat)
      if (stat.ne.0) print*,'error allocating ipv(nx,ny,nz)'
      
C     Determine if data is on pressure or model levels
      
      prelev=.true.
      do k=1,nz
         if (bklev(k).ne.0.) prelev=.false.
         if (bklay(k).ne.0.) prelev=.false.
      enddo
C      print*,'prelev ',prelev
      
C     Calculate cos(latitude) array and the coriolis parameter
      
      if ((abs(pollon).gt.0.001).or.(abs(pollat-90.).gt.0.001)) then
         do j=1,ny
            rlat=varmin(2)+(j-1)*dy
            do i=1,nx
               rlon=varmin(1)+(i-1)*dx
               yphys=phstoph(rlat,rlon,pollat,pollon)
C              if I use sind(lat in deg): troubles at the N-pole
               lat=2.*pi*yphys/360.
               cl(i,j)=cosd(rlat)
               tl(i,j)=tan(lat)
               f(i,j)=0.000145444*sin(lat)
            enddo
         enddo
      else
         do j=1,ny
            lat=varmin(2)+(j-1)*dy
            lat=2.*pi*lat/360.
            do i=1,nx
               cl(i,j)=cos(lat)
               f(i,j)=0.000145444*sin(lat)
            enddo
         enddo
      endif
      
C     Determine if data is on levels or layers
      
      if (stag(3).eq.-0.5) then
         ak=aklay
         bk=bklay
      else
         ak=aklev
         bk=bklev
      endif
      
C     Get all the fields
C     ------------------
      
      call gettimes(cdfid1,time,ntimes,ierr)
      
C     Loop over all times
C     -------------------

C      print*,'loop over all times'
C      print*,ntimes,vardim,nx,ny,nz
      
      do n=1,ntimes
         
         if (verbose) print*,'time=',time(n)
         
         if (.not.prelev) then
            if (PSonP) then
               call getdat(cdfid1,'PS',time(n),0,sp,ierr)
               if (ierr.ne.0) goto 921
            else
               if (verbose) print*,'PS not on P-file'
            endif
         endif

         if (ZBonP) then
            call getdat(cdfid1,'ZB',time(n),0,oro,ierr)
            if (ierr.ne.0) goto 930
         else
            if (verbose) print*,'ZB not on P-file'
         endif
         
         if (TonP) then
            call getdat(cdfid1,'T',time(n),0,tt,ierr)
            if (ierr.ne.0) goto 920
         else
            if (verbose) print*,'T not on P-file'
         endif

         if (qmode.eq.'Q') then             
            call getdat(cdfid1,'Q',time(n),0,qq,ierr)
            if (ierr.ne.0) goto 922         
         elseif (qmode.eq.'QD') then    
            call getdat(cdfid1,'QD',time(n),0,qq,ierr)
            if (ierr.ne.0) goto 922
         else
            if (verbose) print*,'neither Q nor QD on P-file'
         endif
         
         if (ZonP) then
            allocate(zz(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating zz'
            call getdat(cdfid1,'Z',time(n),0,zz,ierr)
            if (ierr.ne.0) then
               print*,'error reading Z'
               STOP
            else
               Z_calc=0         !indicate Z read
            endif
         endif

         if (ZBneed) then
            ZBneed=.false.
            allocate(zb(nx,ny),STAT=stat)
            if (stat.ne.0) print*,'error allocating zb'               
            if (ZBfile.ne.'') then
               call cdfopn(trim(ZBfile),cdfid2,ierr)
               if (ierr.ne.0) then
                  print*,'error opening ',trim(ZBfile)
                  ZBneed=.true.
               endif
               call gettimes(cdfid2,time2,ntimes2,ierr)
               call getdat(cdfid2,'ZB',time2(1),0,zb,ierr)
            else if (ZBonP) then               
               call getdat(cdfid1,'ZB',time(1),0,zb,ierr)
            else
               print*,'ZB needed (use option -zb filename)'
               ZBneed=.true.
            endif
            if (ierr.ne.0) then
               print*,'error reading ZB'
               ZBneed=.true.
            endif
         endif
         
         if (UonP) then
            call getdat(cdfid1,'U',time(n),0,uu,ierr)
            if (ierr.ne.0) goto 923
         else
            if (verbose) print*,'U not on P-file'
         endif
         
         if (VonP) then
            call getdat(cdfid1,'V',time(n),0,vv,ierr)
            if (ierr.ne.0) goto 924
         else
            if (verbose) print*,'V not on P-file'
         endif
         
         if (OMEGAonP) then
            call getdat(cdfid1,'OMEGA',time(n),0,ww,ierr)
            if (ierr.ne.0) goto 925
         else
            if (verbose) print*,'OMEGA not on P-file'
         endif

         if (T2MonP) then
            call getdat(cdfid1,'T2M',time(n),0,t2m,ierr)
            if (ierr.ne.0) goto 926
         else
            if (verbose) print*,'T2M not on P-file'
         endif

         if (TD2MonP) then
            call getdat(cdfid1,'TD2M',time(n),0,td2m,ierr)
            if (ierr.ne.0) goto 926
         else
            if (verbose) print*,'TD2M not on P-file'
         endif

         if (U10MonP) then
            call getdat(cdfid1,'U10M',time(n),0,u10m,ierr)
            if (ierr.ne.0) goto 926
         else
            if (verbose) print*,'U10M not on P-file'
         endif

         if (V10MonP) then
            call getdat(cdfid1,'V10M',time(n),0,v10m,ierr)
            if (ierr.ne.0) goto 926
         else
            if (verbose) print*,'V10M not on P-file'
         endif

         if (PVonP) then
            call getdat(cdfid1,'PV',time(n),0,ipv,ierr)
            if (ierr.ne.0) goto 926
         else
            if (verbose) print*,'PV not on P-file'
         endif

         if (PV3onP) then
            call getdat(cdfid1,'PV3',time(n),0,ipv,ierr)
            if (ierr.ne.0) goto 926
         else
            if (verbose) print*,'PV3 not on P-file'
         endif

C     Calculation of the geopotential
         
         if (Z_calc.eq.1) then                        
            allocate(zz(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating zz'       
            call zlev(zz,tt,qq,oro,sp,nx,ny,nz,aklev,bklev)            
         endif

         if (Z_outP.eq.1) then
            if (n.eq.1) then
               stag(3)=0.0
               call putdef(cdfid1,'Z',4,mdv,
     >              vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable Z created on P-file'
               stag(3)=-0.5
            endif
            call putdat(cdfid1,'Z',time(n),0,zz,ierr)           
         endif
         
C     Create the secondary data file
         
         if (n.eq.1) then
            call cdfwopn(trim(outfnam),cdfid,ierr)
            if (ierr.ne.0) then
               call crecdf(trim(outfnam),cdfid,varmin,varmax,
     >              3,trim(cstnam),ierr)
               write(*,*)'*** NetCDF file ',trim(outfnam),' created'
            else
               write(*,*)'*** NetCDF file ',trim(outfnam),
     >              ' used for appending'
               PS_out=0
            endif
         endif
         

C     Put surface pressure on S-file
         
         if ((.not.prelev).and.(PS_out.eq.1)) then
            vardim(3)=1
            if (n.eq.1) then
            call putdef(cdfid,'PS',4,mdv,vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable PS created on ',trim(outfnam)
            endif
            call putdat(cdfid,'PS',time(n),0,sp,ierr)
            vardim(3)=nz
         endif

C     Put geopotential on S-file
        
         if (Z_out.eq.1) then
            if (n.eq.1) then
               stag(3)=0.0
               call putdef(cdfid,'Z',4,mdv,
     >              vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable Z created on ',trim(outfnam)
               stag(3)=-0.5
            endif               
            call putdat(cdfid,'Z',time(n),0,zz,ierr)
         endif                 
         
C     Calculate the secondary data variables
C     --------------------------------------
         
C     Calculation of potential temperature
         
         if (TH_calc.eq.1) then
            allocate(th(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating th'
            call pottemp(th,tt,sp,nx,ny,nz,ak,bk)
            if (TH_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'TH',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable TH created on ',trim(outfnam)
               endif               
               call putdat(cdfid,'TH',time(n),0,th,ierr)
            endif
         endif
         
C     Calculation of relative humidity
         
         if (RH_calc.eq.1) then
            allocate(rh(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating rh'
            call relhum(rh,qq,tt,sp,nx,ny,nz,ak,bk)
            if (RH_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'RH',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable RH created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'RH',time(n),0,rh,ierr)
            endif
         endif
         
C     Calculation of potential vorticity (equals PV3 in IVE)
         
         if (PV_calc.eq.1) then
            allocate(pv(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating pv'
            call potvort(pv,uu,vv,th,sp,cl,f,nx,ny,nz,ak,bk,
     >           varmin,varmax)
            if (PV_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'PV',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable PV created on ',trim(outfnam)
               endif           
               call putdat(cdfid,'PV',time(n),0,pv,ierr)
            endif
         endif

C     Calculation of LAIT potential vorticity

         if (LPV_calc.eq.1) then
            allocate(lpv(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating lpv'
            call laitpv(lpv,pv,th,nx,ny,nz)
            if (LPV_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'LPV',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable LPV created on ',trim(outfnam)
               endif
               call putdat(cdfid,'LPV',time(n),0,lpv,ierr)
            endif
         endif
         
C     Calculation of equivalent potential temperature
         
         if (THE_calc.eq.1) then
            allocate(the(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating the'
            call equpot(the,tt,qq,sp,nx,ny,nz,ak,bk)
            if (THE_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'THE',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable THE created on ',trim(outfnam)
               endif               
               call putdat(cdfid,'THE',time(n),0,the,ierr)
            endif
         endif
         
C     Calculation of wet-bulb potential temperature
         
         if (THW_calc.eq.1) then
            call wetbpt(var,the,sp,nx,ny,nz,ak,bk)           
            if (THW_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'THW',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable THW created on ',trim(outfnam)
               endif               
               call putdat(cdfid,'THW',time(n),0,var,ierr)
            endif
         endif

C     Calculation of water flux divergence

         if (DIVQU_calc.eq.1) then
            call divqu(var,uu,vv,qq,cl,nx,ny,nz,varmin,varmax,mdv)
            if (DIVQU_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'DIVQU',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable DIVQU created on ',
     >		  trim(outfnam)
               endif
               call putdat(cdfid,'DIVQU',time(n),0,var,ierr)
            endif
         endif
         
C     Calculation of diabatic heating rate (also DHR in IVE)
         
         if (CH_calc.eq.1) then
            allocate(dhr(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating dhr'
            call diabheat(dhr,th,ww,rh,sp,nx,ny,nz,ak,bk)
            if (CH_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'CH',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable CH created on ',trim(outfnam)
               endif               
               call putdat(cdfid,'CH',time(n),0,dhr,ierr)
            endif
         endif
         
C     Calculation of diabatic PV rate (also DPVR in IVE)
         
         if (PVR_calc.eq.1) then
            call diabpvr(var,uu,vv,dhr,sp,cl,f,nx,ny,nz,ak,bk,
     >           varmin,varmax)
            if (PVR_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'PVR',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable PVR created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'PVR',time(n),0,var,ierr)
            endif
         endif
                  
C     Calculation of the theta gradient
         
         if (GRADTH_calc.eq.1) then    
         call gradth(var,th,sp,prelev,cl,nx,ny,nz,ak,bk,varmin,varmax)
            if (GRADTH_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'GRADTH',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable GRADTH created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'GRADTH',time(n),0,var,ierr)
            endif
         endif

C     Calculation of density RHO

         if (RHO_calc.eq.1) then  
            allocate(rho(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating rho'
            call calc_rho(rho,tt,sp,nx,ny,nz,ak,bk)         
            if (RHO_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'RHO',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable RHO created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'RHO',time(n),0,rho,ierr)
            endif
         endif

C     Calculation of specific volume ALPHA

         if (ALPHA_calc.eq.1) then  
	    print*,'calculate alpha',n 
            allocate(alpha(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating rho'
            alpha = 1. / rho   
            if (ALPHA_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'ALPHA',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable ALPHA created on ',
     >            trim(outfnam)
               endif            
               call putdat(cdfid,'ALPHA',time(n),0,alpha,ierr)
            endif
         endif

C     Calculation of Brunt-Vaisala frequ. N^2
         
         if (NSQ_calc.eq.1) then           
            call nsq(var,rho,th,sp,nx,ny,nz,ak,bk)
            if (NSQ_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'NSQ',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable NSQ created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'NSQ',time(n),0,var,ierr)
            endif
         endif

C     Calculation of Brunt-Vaisala frequ. N^2 with ThetaE

         if (NSQM_calc.eq.1) then           
            call nsq(var,rho,the,sp,nx,ny,nz,ak,bk)
            if (NSQM_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'NSQM',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable NSQM created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'NSQM',time(n),0,var,ierr)
            endif
         endif

C     Calculation of vertical gradient of potential temp. (-g d(th)/dp)

         if (DTHDP_calc.eq.1) then
            call dthetadp(var,th,sp,nx,ny,nz,ak,bk)
            if (DTHDP_out.eq.1) then
              if (n.eq.1) then
                call putdef(cdfid,'DTHDP',4,mdv,
     >               vardim,varmin,varmax,stag,ierr)
                write(*,*)'*** variable DTHDP created on ',trim(outfnam)
              endif
              call putdat(cdfid,'DTHDP',time(n),0,var,ierr)
            endif
         endif

C     Calculation of vertical velocity W
c     get the vertical wind velocity in Cart. coordinates. 
c     Omega is given in hPa/s.
         
         if (W_calc.eq.1) then 
            where (rho.ne.0.) var=-100.*ww/(9.80616*rho)
            if (W_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'W',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable W created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'W',time(n),0,var,ierr)
            endif
         endif
         
C     Calculation of velocity VEL
         
         if (VEL_calc.eq.1) then 
            var=sqrt(uu**2+vv**2)          
            if (VEL_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'VEL',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable VEL created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'VEL',time(n),0,var,ierr)
            endif
         endif
         
C     Calculation of ZLAY

         if (ZLAY_calc.eq.1) then
            allocate(zlay(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating zlay'
            call zlayer(zlay,tt,zz,sp,nx,ny,nz,aklev,bklev,aklay,bklay)
            if (ZLAY_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'ZLAY',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable ZLAY created on ',trim(outfnam)
               endif
               call putdat(cdfid,'ZLAY',time(n),0,zlay,ierr)
            endif
         endif

C     Calculation of Montgomery-Potential M
         
         if (M_calc.eq.1) then    
            allocate(mm(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating mm' 
            mm=1004.*(tt+273.15)+9.80616*1000.*zlay
            if (M_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'M',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable M created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'M',time(n),0,mm,ierr)
            endif
         endif

C     Calculation of Bernoulli-Function B
         
         if (B_calc.eq.1) then    
            var=mm+0.5*(uu**2+vv**2)
            if (B_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'B',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable B created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'B',time(n),0,var,ierr)
            endif
         endif

C     Calculation of relative vertical vorticity VORT

         if (VORT_calc.eq.1) then 
            call vort(var,uu,vv,sp,cl,nx,ny,nz,ak,bk,varmin,varmax)
            if (VORT_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'VORT',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable VORT created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'VORT',time(n),0,var,ierr)
            endif
         endif

C     Calculation of absolute vertical vorticity AVO

         if (AVO_calc.eq.1) then
            call avort(var,uu,vv,sp,cl,f,nx,ny,nz,ak,bk,varmin,varmax)
            if (AVO_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'AVO',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable AVO created on ',trim(outfnam)
               endif
               call putdat(cdfid,'AVO',time(n),0,var,ierr)
            endif
         endif

C     Calculation of vertical curvature vorticity CURVO

         if (CURVO_calc.eq.1) then
            call curvo(var,uu,vv,sp,cl,nx,ny,nz,ak,bk,varmin,varmax)
            if (CURVO_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'CURVO',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable CURVO created on ',trim(outfnam)
               endif
               call putdat(cdfid,'CURVO',time(n),0,var,ierr)
            endif
         endif

C     Output of cos(latitude) array
 
         if (COS_out.eq.1) then
            if (n.eq.1) then
               call putdef(cdfid,'COS',4,mdv,
     >              vardim,varmin,varmax,stag,ierr)
            write(*,*)'*** variable COS created on ',trim(outfnam)
            endif
            call putdat(cdfid,'COS',time(n),0,cl,ierr)
         endif

C     Output of vertically integrated PV (VIP)

         if (VIP_out.eq.1) then
            allocate(vip(nx,ny),STAT=stat)
            if (stat.ne.0) print*,'error allocating vip'
            call vint(vip,pv,sp,150.,500.,nx,ny,nz,ak,bk)
            if (n.eq.1) then
               vardim(3)=1
               call putdef(cdfid,'VIP',4,mdv,
     >              vardim,varmin,varmax,stag,ierr)
               vardim(3)=nz
            write(*,*)'*** variable VIP created on ',trim(outfnam)
            endif
            call putdat(cdfid,'VIP',time(n),0,vip,ierr)
         endif

C     Calculation of PV gradient (GRADPV)

         if (GRADPV_out.eq.1) then
C            write(*,*)'*** calling gradpv ***'
            allocate(gradpv(nx,ny),STAT=stat)
            if (stat.ne.0) print*,'error allocating gradpv'
            call calc_gradpv(gradpv,ipv,cl,nx,ny,nz,akd,bkd,
     >           varmin,varmax)
            if (n.eq.1) then
               vardim(3)=1
               call putdef(cdfid,'GRADPV',4,mdv,
     >              vardim,varmin,varmax,stag,ierr)
               vardim(3)=nz
            write(*,*)'*** variable GRADPV created on ',trim(outfnam)
            endif
            call putdat(cdfid,'GRADPV',time(n),0,gradpv,ierr)
         endif

C     Calculation of P

         if (P_calc.eq.1) then
            call pressure(var,sp,-0.5,nx,ny,nz,aklev,bklev,aklay,bklay)
            if (P_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'P',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable P created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'P',time(n),0,var,ierr)
            endif
         endif

C     Calculation of PLEV

         if (PLEV_calc.eq.1) then
            call pressure(var,sp,0.,nx,ny,nz,aklev,bklev,aklay,bklay)
            if (PLEV_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'PLEV',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable PLEV created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'PLEV',time(n),0,var,ierr)
            endif
         endif

C     Calculation of FL

         if (FL_calc.eq.1) then
            allocate(fl(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating fl'
            call pres(var,sp,nx,ny,nz,aklay,bklay)
            call flightlev(fl,var,nx,ny,nz)
            if (FL_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'FL',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable FL created on ',trim(outfnam)
               endif
               call putdat(cdfid,'FL',time(n),0,fl,ierr)
            endif
         endif

C     Calculation of geostrophic wind UG

         if (UG_calc.eq.1) then
            allocate(ug(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating ug' 
            call calc_ug(ug,zlay,sp,cl,f,nx,ny,nz,ak,bk,varmin,varmax)
            if (UG_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'UG',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable UG created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'UG',time(n),0,ug,ierr)
            endif
         endif

C     Calculation of geostrophic wind VG

         if (VG_calc.eq.1) then
            allocate(vg(nx,ny,nz),STAT=stat)
            if (stat.ne.0) print*,'error allocating vg' 
            call calc_vg(vg,zlay,sp,cl,f,nx,ny,nz,ak,bk,varmin,varmax)
            if (VG_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'VG',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable VG created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'VG',time(n),0,vg,ierr)
            endif
         endif

C     Calculation of ageostrophic geostrophic wind UA

         if (UA_calc.eq.1) then
            var=uu-ug
            if (UA_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'UA',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable UA created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'UA',time(n),0,var,ierr)
            endif
         endif

C     Calculation of ageostrophic geostrophic wind VA

         if (VA_calc.eq.1) then
            var=vv-vg
            if (VA_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'VA',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
               write(*,*)'*** variable VA created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'VA',time(n),0,var,ierr)
            endif
         endif


C     Calculation of x-component of the Q-vector (calc. with U)

         if (QXF_calc.eq.1) then
         call calc_qx(var,uu,vv,tt,sp,cl,nx,ny,nz,ak,bk,varmin,varmax)
            if (QXF_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'QXF',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable QXF created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'QXF',time(n),0,var,ierr)
            endif
         endif
         
C     Calculation of y-component of the Q-vector (calc. with V)
         
         if (QYF_calc.eq.1) then
      call calc_qy(var,uu,vv,tt,sp,cl,tl,nx,ny,nz,ak,bk,varmin,varmax)
            if (QYF_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'QYF',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable QYF created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'QYF',time(n),0,var,ierr)
            endif
         endif

C     Calculation of x-component of the Q-vector (calc. with UG)

         if (QX_calc.eq.1) then
         call calc_qx(var,ug,vg,tt,sp,cl,nx,ny,nz,ak,bk,varmin,varmax)
            if (QX_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'QX',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable QX created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'QX',time(n),0,var,ierr)
            endif
         endif
         
C     Calculation of y-component of the Q-vector (calc. with VG)
         
         if (QY_calc.eq.1) then
      call calc_qy(var,ug,vg,tt,sp,cl,tl,nx,ny,nz,ak,bk,varmin,varmax)
            if (QY_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'QY',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
                  write(*,*)'*** variable QY created on ',trim(outfnam)
               endif            
               call putdat(cdfid,'QY',time(n),0,var,ierr)
            endif
         endif

C     Calculation of reduced surface pressure PSRED

         if (PSRED_calc.eq.1) then
            if (.not.ZBneed) then
               call calc_psred(var,sp,tt,zb,nx,ny,nz,aklay,bklay)
               if (PSRED_out.eq.1) then
                  if (n.eq.1) then
                     vardim2=vardim
                     vardim2(3)=1 ! force one vertical-level
                     call putdef(cdfid,'PSRED',4,mdv,
     >                    vardim2,varmin,varmax,stag,ierr)
               write(*,*)'*** variable PSRED created on ',trim(outfnam)
                  endif            
                  call putdat(cdfid,'PSRED',time(n),0,var,ierr)
               endif
            else
               print*,'ZB needed to calculate PSRED'
            endif
         endif

C     Calculation of Richardson number

         if (RI_calc.eq.1) then
            call rich(var,sp,uu,vv,th,nx,ny,nz,ak,bk)
            if (RI_out.eq.1) then
               if (n.eq.1) then
                  call putdef(cdfid,'RI',4,mdv,
     >                 vardim,varmin,varmax,stag,ierr)
            write(*,*)'*** variable RI created on ',trim(outfnam)
               endif
               call putdat(cdfid,'RI',time(n),0,var,ierr)
            endif
         endif

C     Calculation of boundary layer height BLH
 
         if (BLH_calc.eq.1) then
            call calc_blh(var,sp,tt,qq,uu,vv,t2m,td2m,u10m,v10m,
     >			  nx,ny,nz,aklay,bklay)
            if (BLH_out.eq.1) then
               if (n.eq.1) then
                  vardim2=vardim
                  vardim2(3)=1 ! force one vertical-level
                  call putdef(cdfid,'BLH',4,mdv,
     >                 vardim2,varmin,varmax,stag,ierr)
            write(*,*)'*** variable BLH created on ',trim(outfnam)
               endif
               call putdat(cdfid,'BLH',time(n),0,var,ierr)
            endif
         endif
         
c     add calculations of your own variables here: (++)
         
      enddo
      
C     Close the NetCDF files
      
      call clscdf(cdfid,ierr)
 900  continue
      call clscdf(cdfid1,ierr)
      call clscdf(cstid,ierr)
      
      goto 999
      
 920  stop '*** error: variable T not found on P-file ***'
 921  stop '*** error: variable PS not found on P-file ***'
 922  stop '*** error: variable Q not found on P-file ***'
 923  stop '*** error: variable U not found on P-file ***'
 924  stop '*** error: variable V not found on P-file ***'
 925  stop '*** error: variable OMEGA not found on P-file ***'
 926  stop '*** error: variable T2M not found on P-file ***'
 927  stop '*** error: variable TD2M not found on P-file ***'
 928  stop '*** error: variable U10M not found on P-file ***'
 929  stop '*** error: variable V10M not found on P-file ***'
 930  stop '*** error: variable ZB not found on P-file ***'
      
 996  stop '*** error: could not create S-file ***'
 999  continue
      end
      
c-----------------------------------------------------------------------
c     following subroutines
c-----------------------------------------------------------------------

      subroutine calc_rho(var,tt,sp,ie,je,ke,ak,bk)
c     ============================================
c     computation of density
c
c     argument declaration
      integer	ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: tt(ie,je,ke),sp(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke)

      integer	i,j,k
      real,parameter ::	tzero=273.15
      
c     statement-functions for the computation of pressure
      real	pp,psrf
      integer	is
      pp(is)=ak(is)+bk(is)*psrf
      
c     computation of potential temperature
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke
c              distinction of temperature in K and deg C
               if (tt(i,j,k).lt.100.) then
                  var(i,j,k)=pp(k)/(2.87*(tt(i,j,k)+tzero))
               else
                  var(i,j,k)=pp(k)/(2.87*tt(i,j,k))                 
               endif
            enddo
         enddo
      enddo      
      end


      subroutine nsq(var,rho,th,sp,ie,je,ke,ak,bk)
C     ======================================================
      
c     argument declaration
      integer,intent(IN) :: ie,je,ke
      real,intent(OUT)   :: var(ie,je,ke)
      real,intent(IN)    :: rho(ie,je,ke),th(ie,je,ke),sp(ie,je)
      real,intent(IN)    :: ak(ke),bk(ke)
      
c     variable declaration
      integer   stat
      REAL,ALLOCATABLE, DIMENSION (:,:,:) ::  dthdp !3D array
      
      allocate(dthdp(ie,je,ke),STAT=stat)
      IF (stat.ne.0) PRINT*,'nsq: error allocating dthdp(ie,je,ke)'

      call ddp(th,dthdp,sp,ie,je,ke,ak,bk)

      where(th.ne.0.) var=-96.24*rho/th*dthdp

      IF (ALLOCATED(dthdp)) DEALLOCATE(dthdp)

      end

      subroutine dthetadp(var,th,sp,ie,je,ke,ak,bk)
C     ======================================================

c     argument declaration
      integer,intent(IN) :: ie,je,ke
      real,intent(OUT)   :: var(ie,je,ke)
      real,intent(IN)    :: th(ie,je,ke),sp(ie,je)
      real,intent(IN)    :: ak(ke),bk(ke)

c     variable declaration
      integer   stat
      REAL,ALLOCATABLE, DIMENSION (:,:,:) ::  dthdp !3D array

      allocate(dthdp(ie,je,ke),STAT=stat)
      IF (stat.ne.0) PRINT*,'nsq: error allocating dthdp(ie,je,ke)'

      call ddp(th,dthdp,sp,ie,je,ke,ak,bk)

c     factor 100 is because VORT is in units of f and so VORT*DTHDP
c     is in pvu
      var=-9.80616*100.*dthdp

      IF (ALLOCATED(dthdp)) DEALLOCATE(dthdp)

      end

      subroutine geopot(psi,q,t,oro,sp,ie,je,ke,ak,bk)
c     ================================================
      
c     argument declaration
      integer   ie,je,ke
      real      psi(ie,je,ke),t(ie,je,ke),q(ie,je,ke),oro(ie,je),
     >     sp(ie,je),ak(ke),bk(ke)
      
c     variable declaration
      integer   i,j,k
      real      r,c,g
      data	r,c,g /287.,0.608,9.80616/
      
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
      
c     integration of geopotential height(special for first layer)
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            psi(i,j,1)=1./g*(oro(i,j)
     >           +r*(t(i,j,1)+273.15)*(1.+c*q(i,j,1))*
     >           (psrf-pp(1))/(0.5*(psrf+pp(1))))
         enddo
      enddo
      do j=1,je
         do i=1,ie
            psrf=sp(i,j)
            do k=2,ke
               psi(i,j,k)=psi(i,j,k-1)+r/g*
     >              ((t(i,j,k-1)+273.15)*(1.+c*q(i,j,k-1))+
     >              (t(i,j,k)+273.15)*(1.+c*q(i,j,k)))*
     >              (pp(k-1)-pp(k))/(pp(k-1)+pp(k))
            enddo
         enddo
      enddo
      end
      

      subroutine getoro(or,dx,dy,starty,lonmin,latmin,ie,je,ierr)
c     ===========================================================
c     reads the orography for the actual data domain from the file
c     "/home/henry/ecmwf/cdf/orogr", which contains the orography
c     for the whole northern hemisphere.
      
c     argument declaration
      integer  ie,je
      real     or(ie,je)
      
c     variable declaration
      integer  i,j,cdfid,ierr,i0,j0,err,starty,nx
      real     gloro(480*240)
      real     dx,dy,dd,time,lonmin,latmin
      character*(80)filnam
      integer  vardim(3),ndim
      real     varmin(3),varmax(3),stag(3),mdv

      time=0.
      
c     open file with orography values and get them
      call ncpopt(NCVERBOS)
      
c     q&d bug fix (dx, dy were changed below when calculating i0)
      dd=dx
      if ((dx.eq.0.75).and.(dy.eq.0.75)) then
         if (starty.lt.96) then
            if (starty.le.92) then
               filnam="/home/henry/cdf/oro/9209orogr0_75"
            else if (starty.le.93) then
               filnam="/home/henry/cdf/oro/9309orogr0_75"
            else if (starty.le.94) then
               filnam="/home/henry/cdf/oro/9411orogr0_75"
            else if (starty.le.95) then
               filnam="/home/henry/cdf/oro/9509orogr0_75"
            else
               filnam="/home/henry/cdf/oro/orogr0_75"
            endif
         else
            filnam="/home/henry/cdf/oro/95orogr0_75"
         endif
      else if ((dx.eq.1.0).and.(dy.eq.1.0)) then
         if (starty.lt.95) then
            filnam="/home/henry/cdf/oro/orogr1_00"
         else
            filnam="/home/henry/cdf/oro/gloro1_00"
         endif
      else if ((dx.eq.0.5).and.(dy.eq.0.5)) then
         filnam="/home/henry/cdf/oro/95orogr0_50"
      else
         print*,'unappropriate data for the orography file'
         print*,'grid interval should be 0.5, 0.75 or 1.00
     >        and data only on the NH'
         ierr=2
         return
      endif
      print*,filnam
      call cdfopn(trim(filnam),cdfid,err)
      call getdat(cdfid,'ORO',time,0,gloro,err)
      call getdef(cdfid,'ORO',ndim,mdv,vardim,varmin,varmax,stag,err)
*     call clscdf(cdfid,err)
      
      i0=nint((lonmin-varmin(1))/dd)
      j0=nint((latmin-varmin(2))/dd)
      if (j0.lt.0) j0=-j0
      
      nx=nint(360./dd)
      do i=1,ie
      do j=1,je
        if (i0+i.le.nx) then
          or(i,j)=gloro(i0+i+(j0+j-1)*nx)
        else
          or(i,j)=gloro(i0+i-nx+(j0+j-1)*nx)
        endif
        if (or(i,j).gt.1.e20) print*,i,j,i0+i+(j0+j-1)*nx
      enddo
      enddo      
      end

      subroutine calc_qx(var,uu,vv,tt,sp,cl,ie,je,ke,ak,bk,vmin,vmax)
C     ===============================================================
c     calculate x-comp Q-vector: -9.8/273.*(D[U:X]*D[T:X]+D[V:X]*D[T:Y])

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),tt(ie,je,ke)
      real,intent(IN)  :: sp(ie,je),cl(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke),vmin(4),vmax(4)

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dudx,dtdx,dvdx,dtdy
      integer	stat
      integer	k
      real	mdv
      logical	prelev

c     test if pressure or model levels
      prelev=.true.
      do k=1,ke
        if (bk(k).ne.0.) prelev=.false.
      enddo
      mdv=-999.98999

      if (.not.prelev) then
        allocate(dspdx(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qx: error allocating dspdx'
        allocate(dspdy(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qx: error allocating dspdy'
      endif

      allocate(dudx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dudx'
      allocate(dtdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dtdx'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dvdx'
      allocate(dtdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dtdy'

      if (prelev) then
        do k=1,ke
          call ddh2m(uu(1,1,k),dudx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(vv(1,1,k),dvdx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
        enddo
      else
        call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
        call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
        call ddh3(uu,dudx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      endif

      var=-1.e9*9.80616/273.*(dudx*dtdx+dvdx*dtdy)
      if (prelev) then
        do i=1,ie
        do j=1,je
        do k=1,ke
          if ((dudx(i,j,k).eq.mdv).or.
     >        (dtdx(i,j,k).eq.mdv).or.
     >        (dvdx(i,j,k).eq.mdv).or.
     >        (dtdy(i,j,k).eq.mdv)) then
            var(i,j,k)=mdv
          endif
        enddo
        enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dudx)) DEALLOCATE(dudx)
      IF (ALLOCATED(dtdx)) DEALLOCATE(dtdx)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dtdy)) DEALLOCATE(dtdy)

      end


      subroutine calc_qy(var,uu,vv,tt,sp,cl,tl,ie,je,ke,ak,bk,vmin,vmax)
C     ===============================================================
c     calculate y-comp Q-vector: -9.8/273.*(D[U:Y]*D[T:X]+D[V:Y]*D[T:Y]
c                              -U/6.37E6*TAN*D[T:X]-V/6.37E6*TAN*D[T:Y])

c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),tt(ie,je,ke)
      real,intent(IN)  :: sp(ie,je),cl(ie,je),tl(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke),vmin(4),vmax(4)

c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dudy,dtdx,dvdy,dtdy,tl3d
      integer	k,stat
      real      mdv
      logical   prelev
 
c     test if pressure or model levels
      prelev=.true.
      do k=1,ke
        if (bk(k).ne.0.) prelev=.false.
      enddo
      mdv=-999.98999
      
      if (.not.prelev) then
        allocate(dspdx(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qy: error allocating dspdx'
        allocate(dspdy(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qy: error allocating dspdy'
      endif

      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dudy'
      allocate(dtdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dtdx'
      allocate(dvdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dvdy'
      allocate(dtdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dtdy'
      allocate(tl3d(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating tl3d'

      if (prelev) then
        do k=1,ke
          call ddh2m(uu(1,1,k),dudy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(vv(1,1,k),dvdy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
        enddo
      else
        call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
        call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
        call ddh3(uu,dudy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(vv,dvdy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      endif


      do k=1,ke
         tl3d(1:ie,1:je,k)=tl(1:ie,1:je)
      enddo

      var=-1.e9*9.80616/273.*(dudy*dtdx+dvdy*dtdy
     >    -uu/6.37E6*tl3d*dtdx-vv/6.37E6*tl3d*dtdy)
      if (prelev) then
        do i=1,ie
        do j=1,je
        do k=1,ke
          if ((dudy(i,j,k).eq.mdv).or.
     >        (dtdx(i,j,k).eq.mdv).or.
     >        (dvdy(i,j,k).eq.mdv).or.
     >        (dtdy(i,j,k).eq.mdv)) then
            var(i,j,k)=mdv
          endif
        enddo
        enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dtdx)) DEALLOCATE(dtdx)
      IF (ALLOCATED(dvdy)) DEALLOCATE(dvdy)
      IF (ALLOCATED(dtdy)) DEALLOCATE(dtdy)
      IF (ALLOCATED(tl3d)) DEALLOCATE(tl3d)

      end


      subroutine calc_ug(var,zz,sp,cl,f,ie,je,ke,ak,bk,vmin,vmax)
C     ===========================================================
C     calculate geostrophic wind
      
c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: zz(ie,je,ke),
     >     sp(ie,je),cl(ie,je),f(ie,je)
      real	ak(ke),bk(ke),vmin(4),vmax(4)
      
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dzzdy
      integer	k,stat
      
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'calc_ug: error allocating dspdy'
      allocate(dzzdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_ug: error allocating dzzdy'
      
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(zz,dzzdy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
         var(1:ie,1:je,k)=-dzzdy(1:ie,1:je,k)
     >        *9810.*(f(1:ie,1:je)/((ABS(f(1:ie,1:je))+1.e-12)**2))
      enddo  

      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dzzdy)) DEALLOCATE(dzzdy)

      end


      subroutine calc_vg(var,zz,sp,cl,f,ie,je,ke,ak,bk,vmin,vmax)
C     ===========================================================
C     calculate geostrophic wind
      
c     argument declaration
      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: zz(ie,je,ke),
     >     sp(ie,je),cl(ie,je),f(ie,je)
      real	ak(ke),bk(ke),vmin(4),vmax(4)
      
c     variable declaration
      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dzzdx
      integer	k,stat
      
      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'calc_vg: error allocating dspdx'
      allocate(dzzdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_vg: error allocating dzzdx'
      
      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh3(zz,dzzdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
         var(1:ie,1:je,k)=dzzdx(1:ie,1:je,k)
     >        *9810.*(f(1:ie,1:je)/((ABS(f(1:ie,1:je))+1.e-12)**2))
      enddo 

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dzzdx)) DEALLOCATE(dzzdx)
     
      end

      subroutine zlayer(ap,t,z,sp,ie,je,ke,aklev,bklev,aklay,bklay)
c     =============================================================
c     argument declaration
      integer  ie,je,ke
      real,intent(OUT) :: ap(ie,je,ke)
      real,intent(IN)  :: t(ie,je,ke),z(ie,je,ke),sp(ie,je)
      real,intent(IN)  :: aklev(ke),bklev(ke),aklay(ke),bklay(ke)
      
c     variable declaration
      integer  i,j,k
      real     psrf
      real,parameter :: rdg=29.271,tzero=273.15
 
c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf
 
c     computation of height of main levels
      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke-1
               ap(i,j,k) = (z(i,j,k) - rdg*(t(i,j,k)+tzero)*
     +              alog(prlay(k)/prlev(k)))/1000.
            enddo
            ap(i,j,ke) = (z(i,j,ke-1) + rdg*(t(i,j,ke)+tzero)*
     +           alog(prlev(ke-1)/prlay(ke)))/1000.
         enddo
      enddo
      end      

      
      subroutine calc_psred(psr,ps,t,zb,ie,je,ke,aklay,bklay)
c     =======================================================
c     calculate PSRED given T, ZB and ak,bk

c     argument declaration
      integer  ie,je,ke
      real,intent(OUT) :: psr(ie,je)
      real,intent(IN)  :: ps(ie,je),t(ie,je,ke),zb(ie,je)
      real,intent(IN)  :: aklay(ke),bklay(ke)

c     variable declaration
      integer  i,j
      real     psrf
      real     ztstar,zalpha,zt0
      real,parameter :: rdcp=0.286,tzero=273.15,r=287.05,g=9.80665
     
c     statement-functions for the computation of pressure
      real      prlay
      integer   is
      prlay(is)=aklay(is)+bklay(is)*psrf

      do i=1,ie
         do j=1,je
            psrf=ps(i,j)
            if (zb(i,j).lt.1.) then
               psr(i,j)=psrf
            else
               ztstar = (t(i,j,1)+tzero)*(1. + 0.0065*r/g*
     &              (ps(i,j)/prlay(1) - 1.0))
               zalpha = 0.0065*r
               zt0    = ztstar + 0.0065*zb(i,j)
               if (zt0.gt.290.5) then
                  if (ztstar.gt.290.5) then
                     zalpha = 0.0
                     ztstar = 0.5*(ztstar+290.5)
                  else
                     zalpha = r*(290.5-ztstar)/zb(i,j)
                  endif
               else if (ztstar.lt.255.) then
                  ztstar = 0.5*(255.0+ztstar)
               endif
               psr(i,j) = ps(i,j)* exp(g*zb(i,j)/(r*ztstar)*
     &              (1.0 - 0.5*(zalpha*zb(i,j)/(r*ztstar)) +
     &              0.333*(zalpha*zb(i,j)/(r*ztstar))**2))
            endif
         enddo
      enddo
      end

      subroutine calc_blh(blh,ps,t,q,u,v,t2m,td2m,u10m,v10m,
     >			  ie,je,ke,aklay,bklay)
c     ======================================================
c     calculate BLH with routines from Andi Stohl
 
c     argument declaration
      integer  ie,je,ke
      real,intent(OUT) :: blh(ie,je)
      real,intent(IN)  :: ps(ie,je),t2m(ie,je),td2m(ie,je),
     >			  u10m(ie,je),v10m(ie,je)
      real,intent(IN)  :: t(ie,je,ke),q(ie,je,ke),u(ie,je,ke),
     >		          v(ie,je,ke)
      real,intent(IN)  :: aklay(ke),bklay(ke)
 
c     variable declaration
      integer  i,j
      real     psrf
      real     ztstar,zalpha,zt0
      real,parameter :: rdcp=0.286,tzero=273.15,r=287.05,g=9.80665
    
c     statement-functions for the computation of pressure
      real      prlay
      integer   is
      prlay(is)=aklay(is)+bklay(is)*psrf
 
      do i=1,ie
         do j=1,je
           call richardson(ps,ust,t(i,j,1),q(i,j,1),u(i,j,1),v(i,j,1),
     >                     ke,aklay,bklay,hf,t2m,td2m,blh(i,j),wst)
         enddo
      enddo
      end

      subroutine richardson(psurf,ust,ttlev,qvlev,ulev,vlev,nuvz,
     +akz,bkz,hf,tt2,td2,h,wst)
*****************************************************************************
*                                                                           *
*     Calculation of mixing height based on the critical Richardson number. *
*     Calculation of convective time scale.                                 *
*     For unstable conditions, one iteration is performed. An excess        *
*     temperature (dependent on hf and wst) is calculated, added to the     *
*     temperature at the lowest model level. Then the procedure is repeated.*
*                                                                           *
*     Author: A. Stohl                                                      *
*                                                                           *
*     22 August 1996                                                        *
*                                                                           *
*     Literature:                                                           *
*     Vogelezang DHP and Holtslag AAM (1996): Evaluation and model impacts  *
*     of alternative boundary-layer height formulations. Boundary-Layer     *
*     Meteor. 81, 245-269.                                                  *
*                                                                           *
*     Update: 1999-02-01 by G. Wotawa                                       *
*                                                                           *
*     Two meter level (temperature, humidity) is taken as reference level   *
*     instead of first model level.                                         *
*     New input variables tt2, td2 introduced.                              *
*                                                                           *
*****************************************************************************
*                                                                           *
* Variables:                                                                *
* h                          mixing height [m]                              *
* hf                         sensible heat flux                             *
* psurf                      surface pressure at point (xt,yt) [Pa]         *
* tv                         virtual temperature                            *
* wst                        convective velocity scale                      *
*                                                                           *
* Constants:                                                                *
* ric                        critical Richardson number                     *
*                                                                           *
*****************************************************************************
 
*     include 'includepar'
 
      integer i,k,nuvz,itmax,iter
      real tv,tvold,zref,z,zold,pint,pold,theta,thetaref,wm,ri,riold
      real akz(nuvz),bkz(nuvz),ulev(nuvz),vlev(nuvz),hf,wst,tt2,td2,ew
      real psurf,ust,ttlev(nuvz),qvlev(nuvz),h,const,ric,b,excess,bs
      real thetaold,zl,ul,vl,thetal,ril
      parameter(r_air=287.05,ga=9.81)
      parameter(const=r_air/ga,ric=0.25,b=100.,bs=8.5,itmax=3)
 
 
      excess=0.0
      iter=0
 
C Compute virtual temperature and virtual potential temperature at
C reference level (2 m)
******************************************************************
 
30    iter=iter+1
 
      pold=psurf
      tvold=tt2*(1.+0.378*esat(td2)/psurf)
      zold=2.0
      zref=zold
 
      thetaref=tvold*(100000./pold)**(r_air/cpa)+excess
      riold=0.
 
 
C Integrate z up to one level above zt
**************************************
 
      do 10 k=3,nuvz
        pint=akz(k)+bkz(k)*psurf  ! pressure on model layers
        wm=qvlev(k)/(1.-qvlev(k))
        tv=ttlev(k)*(1.+0.378*wm/(wm+.622))
 
        if (abs(tv-tvold).gt.0.2) then
          z=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          z=zold+const*log(pold/pint)*tv
        endif
 
        theta=tv*(100000./pint)**(r_air/cpa)
 
 
Calculate Richardson number at each level
*****************************************
 
        ri=ga/thetaref*(theta-thetaref)*(z-zref)/
     +  max(((ulev(k)-ulev(2))**2+(vlev(k)-vlev(2))**2+b*ust**2),0.1)
 
        if (ri.gt.ric) goto 20
 
        riold=ri
        tvold=tv
        pold=pint
        thetaold=theta
10      zold=z
 
20    continue
 
C Determine Richardson number between the critical levels
*********************************************************
 
      do 15 i=1,20
        zl=zold+float(i)/20.*(z-zold)
        ul=ulev(k-1)+float(i)/20.*(ulev(k)-ulev(k-1))
        vl=vlev(k-1)+float(i)/20.*(vlev(k)-vlev(k-1))
        vl=vlev(k-1)+float(i)/20.*(vlev(k)-vlev(k-1))
        thetal=thetaold+float(i)/20.*(theta-thetaold)
        ril=ga/thetaref*(thetal-thetaref)*(zl-zref)/
     +  max(((ul-ulev(2))**2+(vl-vlev(2))**2+b*ust**2),0.1)
        if (ril.gt.ric) goto 25
15      continue
 
25    continue
      h=zl
 
 
C Calculate convective velocity scale
*************************************
 
      if (hf.lt.0.) then
        wst=(-h*ga/thetaref*hf/cpa)**0.333
        excess=-bs*hf/cpa/wst
        if (iter.lt.itmax) goto 30
      else
        wst=0.
      endif
 
      return
      end
