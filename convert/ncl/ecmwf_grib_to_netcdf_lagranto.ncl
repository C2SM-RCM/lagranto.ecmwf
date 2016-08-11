; ecmwf_grib_to_netcdf_lagranto.ncl
; 
; NCL script to read in ECMWF grib data and output in netCDF format 
; for lagranto.
;
; Written by Andrew Penny - November 26, 2013
;=========================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  ml_grb_DIR 		= "/Volumes/Quadra4/model_data_4/ECMWF_YOTC_GLOBAL/model_level/"
  sfc_grb_DIR 	= "/Volumes/Quadra4/model_data_4/ECMWF_YOTC_GLOBAL/surface/"
  out_DIR		= "/Volumes/Quadra4/model_data_4/ECMWF_YOTC_GLOBAL/netcdf/"
  basetime 		= "hours since 1979-01-01 00:00" 
  write_const_file 	= True
  
; START TIME
  s_yr 	= 2008
  s_mo 	= 8
  s_day	= 25
  s_hr		= 0
  s_mn		= 0
  s_sec	= 0

; END TIME  
  e_yr 	= 2008
  e_mo 	= 8
  e_day	= 30
  e_hr		= 12
  e_mn		= 0
  e_sec	= 0
  
; PARSE STARTING AND ENDING DATES
  start_date 	= ut_inv_calendar(s_yr, s_mo, s_day, s_hr, s_mn, s_sec, basetime, 0)
  start_date_utc = ut_calendar(start_date,0) 
  start_date_out = start_date_utc(:,0) + "-" + where(start_date_utc(:,1).lt.10,"0"+start_date_utc(:,1),start_date_utc(:,1)) + "-" + where(start_date_utc(:,2).lt.10,"0"+start_date_utc(:,2),start_date_utc(:,2))+ "-" + where(start_date_utc(:,3).lt.10,"0"+start_date_utc(:,3),start_date_utc(:,3)) + where(start_date_utc(:,4).lt.10,"0"+start_date_utc(:,4),start_date_utc(:,4))

  end_date 	= ut_inv_calendar(e_yr, e_mo, e_day, e_hr, e_mn, e_sec, basetime, 0)
  end_date_utc = ut_calendar(end_date,0) 
  end_date_out = end_date_utc(:,0) + "-" + where(end_date_utc(:,1).lt.10,"0"+end_date_utc(:,1),end_date_utc(:,1)) + "-" + where(end_date_utc(:,2).lt.10,"0"+end_date_utc(:,2),end_date_utc(:,2))+ "-" + where(end_date_utc(:,3).lt.10,"0"+end_date_utc(:,3),end_date_utc(:,3)) + where(end_date_utc(:,4).lt.10,"0"+end_date_utc(:,4),end_date_utc(:,4))      
  
  ml_FILES_ALL 	= systemfunc ("csh -c 'ls "+ ml_grb_DIR + "*.grib'")    
  ml_FILEnames 	= str_sub_str(systemfunc ("csh -c 'ls "+ ml_grb_DIR + "*.grib'"),ml_grb_DIR,"") 
    
  sfc_FILES_ALL 	= systemfunc ("csh -c 'ls "+ sfc_grb_DIR + "*.grib'")    
  sfc_FILEnames 	= str_sub_str(systemfunc ("csh -c 'ls "+ sfc_grb_DIR + "*.grib'"),sfc_grb_DIR,"")     
  
; PARSE FILE DATES
  ml_files_yr  	= stringtoint(str_get_cols(ml_FILEnames, 0,3))
  ml_files_mo  	= stringtoint(str_get_cols(ml_FILEnames, 4,5))
  ml_files_day  	= stringtoint(str_get_cols(ml_FILEnames, 6,7))
  ml_files_hr  	= stringtoint(str_get_cols(ml_FILEnames, 8,9))      
  ml_files_mn  	= ispan(0,dimsizes(ml_files_yr)-1,1)*0
  ml_files_sec  	= ispan(0,dimsizes(ml_files_yr)-1,1)*0
  ml_files_date  	= ut_inv_calendar(ml_files_yr, ml_files_mo, ml_files_day, ml_files_hr, ml_files_mn, ml_files_sec, basetime, 0)
  ml_files_date_utc  = ut_calendar(ml_files_date, 0)
  ml_files_date_out	= ml_files_date_utc(:,0) + "-" + where( ml_files_date_utc(:,1).lt.10, "0" + ml_files_date_utc(:,1), ml_files_date_utc(:,1)) + "-" + where( ml_files_date_utc(:,2).lt.10, "0" + ml_files_date_utc(:,2), ml_files_date_utc(:,2)) + " " + where( ml_files_date_utc(:,3).lt.10, "0" + ml_files_date_utc(:,3), ml_files_date_utc(:,3)) + where( ml_files_date_utc(:,4).lt.10, "0" + ml_files_date_utc(:,4), ml_files_date_utc(:,4))         

  sfc_files_yr  	= stringtoint(str_get_cols(sfc_FILEnames, 0,3))
  sfc_files_mo  	= stringtoint(str_get_cols(sfc_FILEnames, 4,5))
  sfc_files_day  	= stringtoint(str_get_cols(sfc_FILEnames, 6,7))
  sfc_files_hr  	= stringtoint(str_get_cols(sfc_FILEnames, 8,9))      
  sfc_files_mn  	= ispan(0,dimsizes(sfc_files_yr)-1,1)*0
  sfc_files_sec  	= ispan(0,dimsizes(sfc_files_yr)-1,1)*0
  sfc_files_date  	= ut_inv_calendar(sfc_files_yr, sfc_files_mo, sfc_files_day, sfc_files_hr, sfc_files_mn, sfc_files_sec, basetime, 0)
  sfc_files_date_utc	= ut_calendar(sfc_files_date, 0)
  sfc_files_date_out	= sfc_files_date_utc(:,0) + "-" + where( sfc_files_date_utc(:,1).lt.10, "0" + sfc_files_date_utc(:,1), sfc_files_date_utc(:,1)) + "-" + where( sfc_files_date_utc(:,2).lt.10, "0" + sfc_files_date_utc(:,2), sfc_files_date_utc(:,2)) + " " + where( sfc_files_date_utc(:,3).lt.10, "0" + sfc_files_date_utc(:,3), sfc_files_date_utc(:,3)) + where( sfc_files_date_utc(:,4).lt.10, "0" + sfc_files_date_utc(:,4), sfc_files_date_utc(:,4))         
  
; TAG FILES WITHINN START AND END DATES  
  ml_ind_files	= ind(ml_files_date.ge.start_date.and.ml_files_date.le.end_date)
  ml_FILES		= ml_FILES_ALL(ml_ind_files)
  ml_numFILES 	= dimsizes(ml_FILES)    

  sfc_ind_files	= ind(sfc_files_date.ge.start_date.and.sfc_files_date.le.end_date)
  sfc_FILES		= sfc_FILES_ALL(sfc_ind_files)
  sfc_numFILES 	= dimsizes(sfc_FILES)   
  
  if(ml_numFILES.ne.sfc_numFILES)then
    print("model level / surface file number mismatch")
    exit
  end if
  
; DEFINE SOME CONSTANTS
  R      		= 287.05
  g   			= 9.80665  
  
; OPEN ECMWF YOTC COEFFICIENTS FILE AND READ IN VALUES
  strs = asciiread("/Users/abpenny/code/ncl/ecmwf_coeffs",-1,"string")

  delim = ","                                              
  N  = stringtoint(str_get_field(strs, 1, delim))	
  a  = stringtofloat(str_get_field(strs, 2, delim))	
  b  = stringtofloat(str_get_field(strs, 3, delim))	  
  
  dims_N = dimsizes(N)  
  
; FLIP SO BOTTOM --> TOP
  aklev = a(dims_N-2:0)*0.01	; need to multiply "a" coefficients by 0.01 so that they match surface pressure in hPa
  bklev = b(dims_N-2:0)

  aklay = new((/(dims_N-1)/),float,"No_FillValue") 
  bklay = new((/(dims_N-1)/),float,"No_FillValue")   

  do k = 0, dims_N - 3  
    aklay(k+1) = 0.5*(aklev(k)+aklev(k+1))
    bklay(k+1) = 0.5*(bklev(k)+bklev(k+1))
  end do
    
  aklay(0)=0.5*(0. + aklev(0))
  bklay(0)=0.5*(1. + bklev(0))    
  
  aklev!0    		= "nz"    
  copy_VarCoords(aklev,bklev)        
  copy_VarCoords(aklev,aklay)        
  copy_VarCoords(aklev,bklay)            
      
;<><><><><><><><><><><><><><>
; --- LOOP THROUGH FILES ---

do ff = 0, ml_numFILES - 1

  print(" ")
  print("PROCESSING DATE: " + ml_files_date_out(ml_ind_files(ff)) ) 
     
  ml_grb_file 	= addfile(ml_FILES(ff),"r")
  sfc_grb_file 	= addfile(sfc_FILES(ff),"r")  
;  print(getfilevarnames(ml_grb_file)) 
;  print(getfilevarnames(sfc_grb_file)) 

; LATITUDE/LONGITUDE
  lat_ml 		= (/ml_grb_file->g0_lat_0(::-1)/)			; latitude [degrees_north] (Cylindrical Equidistant Projection Grid), lat
  lon_subset 		= (/ml_grb_file->g0_lon_1/)				; longitude [degrees_east] (Cylindrical Equidistant Projection Grid), lon
  dims_y		= dimsizes(lat_ml)
  dims_x		= dimsizes(lon_subset)+1
  dims_t		= 1
  dims_sfc		= 1
   
; MAKE LONGITUDE WRAP AROUND
  lon_ml		= new((/dims_x/),float,-999.99)   
  lon_ml(0:dims_x-2)	= lon_subset
  lon_ml(dims_x-1)	= 180.0  
  
; VERTICAL LEVELS
  lev_ml 		= (/ml_grb_file->lv_HYBL2/)				; Hybrid level (atmosphere hybrid sigma pressure coordinate)
  dims_z		= dimsizes(lev_ml)
  
; 3D VARIABLES (flip latitude so monotonically increasing)
  t_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)   
  u_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)   
  v_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)   
  q_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)   
  w_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)   

  t_ml(0,:,:,0:dims_x-2) 	= (/ml_grb_file->T_GDS0_HYBL(:,::-1,:)/)		; temperature [K], lv_HYBL2/lat/lon
  u_ml(0,:,:,0:dims_x-2) 	= (/ml_grb_file->U_GDS0_HYBL(:,::-1,:)/)		; u wind [m s**-1], lv_HYBL2/lat/lon
  v_ml(0,:,:,0:dims_x-2) 	= (/ml_grb_file->V_GDS0_HYBL(:,::-1,:)/)		; v wind [m s**-1], lv_HYBL2/lat/lon
  q_ml(0,:,:,0:dims_x-2) 	= (/ml_grb_file->Q_GDS0_HYBL(:,::-1,:)/)		; specific humidity [kg kg**-1], lv_HYBL2/lat/lon
  w_ml(0,:,:,0:dims_x-2) 	= (/ml_grb_file->W_GDS0_HYBL(:,::-1,:)/)		; vertical velocity (omega) [Pa s**-1], lv_HYBL2/lat/lon

; WRAP AROUND
  t_ml(0,:,:,dims_x-1) 	= (/ml_grb_file->T_GDS0_HYBL(:,::-1,0)/)		; temperature [K], lv_HYBL2/lat/lon
  u_ml(0,:,:,dims_x-1) 	= (/ml_grb_file->U_GDS0_HYBL(:,::-1,0)/)		; u wind [m s**-1], lv_HYBL2/lat/lon
  v_ml(0,:,:,dims_x-1) 	= (/ml_grb_file->V_GDS0_HYBL(:,::-1,0)/)		; v wind [m s**-1], lv_HYBL2/lat/lon
  q_ml(0,:,:,dims_x-1) 	= (/ml_grb_file->Q_GDS0_HYBL(:,::-1,0)/)		; specific humidity [kg kg**-1], lv_HYBL2/lat/lon
  w_ml(0,:,:,dims_x-1) 	= (/ml_grb_file->W_GDS0_HYBL(:,::-1,0)/)		; vertical velocity (omega) [Pa s**-1], lv_HYBL2/lat/lon

; SURFACE VARIABLES
  sst_sfc 		= new((/dims_t,dims_sfc,dims_y,dims_x/),float,-999.99) 
  z_sfc 		= new((/dims_t,dims_sfc,dims_y,dims_x/),float,-999.99)   
  ps_sfc 		= new((/dims_t,dims_sfc,dims_y,dims_x/),float,-999.99)   
  mslp_sfc 		= new((/dims_t,dims_sfc,dims_y,dims_x/),float,-999.99)   
  sst_sfc(0,0,:,0:dims_x-2)		= (/sfc_grb_file->SSTK_GDS0_SFC(::-1,:)/)		; sea surface temperature [K], lat/lon
  z_sfc(0,0,:,0:dims_x-2) 		= (/sfc_grb_file->Z_GDS0_SFC(::-1,:)/)		; geopotential [m**2 s**-2], lat/lon  
  ps_sfc(0,0,:,0:dims_x-2) 		= (/sfc_grb_file->SP_GDS0_SFC(::-1,:)/)		; surface pressure [Pa], lat/lon - keep in Pa for GHT calculation below
  mslp_sfc(0,0,:,0:dims_x-2) 	= (/sfc_grb_file->MSL_GDS0_SFC(::-1,:)/)		; mean sea-level pressure [Pa], lat/lon    

; WRAP AROUND
  sst_sfc(0,0,:,dims_x-1) 	= (/sfc_grb_file->SSTK_GDS0_SFC(::-1,0)/)		; sea surface temperature [K], lat/lon
  z_sfc(0,0,:,dims_x-1) 	= (/sfc_grb_file->Z_GDS0_SFC(::-1,0)/)		; geopotential [m**2 s**-2], lat/lon  
  ps_sfc(0,0,:,dims_x-1) 	= (/sfc_grb_file->SP_GDS0_SFC(::-1,0)/)		; surface pressure [Pa], lat/lon - keep in Pa for GHT calculation below
  mslp_sfc(0,0,:,dims_x-1) 	= (/sfc_grb_file->MSL_GDS0_SFC(::-1,0)/)		; mean sea-level pressure [Pa], lat/lon    
  
;<><><><><><><><><><><><><><><>  
; COMPUTE GEOPOTENTIAL HEIGHT
;<><><><><><><><><><><><><><><>  
  z_ml   	= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)   
  hgtprev   	= z_sfc(0,0,:,:)

  do j = 1, dims_N-1	 
  k = (dims_N-1) - j + 1  
  pend   = 0.5 * (a(k-1) + a(k)) + ps_sfc(0,0,:,:) * 0.5 * (b(k-1) + b(k))   
   
; IF USING vinth2p_ecmwf where we need to work with geopotential, not geopotential height
  if (k .eq. (dims_N-1) ) then               
    pstart = ps_sfc(0,0,:,:)
    z_ml(0,k-1,:,:) = hgtprev + R * t_ml(0,k-1,:,:) * (1. + 0.61 * q_ml(0,k-1,:,:)) * log(pstart/pend)     
  else
    pstart = 0.5 * (a(k) + a(k+1)) + ps_sfc(0,0,:,:) * 0.5 * (b(k) + b(k+1))     
    z_ml(0,k-1,:,:) = hgtprev + R * (t_ml(0,k,:,:) + t_ml(0,k-1,:,:))/2.*(1. + 0.61 * (q_ml(0,k,:,:)+q_ml(0,k-1,:,:))/2.) * log(pstart/pend) 
  end if      
    hgtprev = (/z_ml(0,k-1,:,:)/)    
  end do    
  
  z_ml = (/z_ml/g/)
  
; CONVERT SURFACE PRESSURE AND SLP TO hPA FROM Pa
  ps_sfc   = ps_sfc*0.01
  mslp_sfc = mslp_sfc*0.01
  
; COMPUTE PRESSURE
  prs_half 	= new((/dims_N,dims_y,dims_x/),float,"No_FillValue") 
  p_ml   	= new((/(dims_N-1),dims_y,dims_x/),float,"No_FillValue") 
  
  do k = 0, dims_N - 1
    prs_half(k,:,:) = a(k) + b(k)*ps_sfc(0,0,:,:)*100.0    
  end do
  
  do k = 0, dims_N - 2
    p_ml(k,:,:) = 0.5*(prs_half(k,:,:) + prs_half(k+1,:,:))	
  end do

; COMPUTE MIXING RATIO
  mixr_ml  		= q_ml/(1.0-q_ml)
  
; COMPUTE RELATIVE HUMIDITY  
  rh_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)     
  rh_ml(0,:,:,:)	= (/relhum(t_ml(0,:,:,:), mixr_ml(0,:,:,:), p_ml)/)
  rh_ml		= where(rh_ml.gt.100.0,100.0,rh_ml)

; COMPUTE THETA-E
  eth_ml  		= new((/dims_t,dims_z,dims_y,dims_x/),float,-999.99)     
  eth_ml(0,:,:,:)	= (/wrf_eth(q_ml(0,:,:,:), t_ml(0,:,:,:), p_ml)/)
  
; FLIP SO BOTTOM --> TOP
  eth_ml(0,:,:,:) 	= eth_ml(0,::-1,:,:)
  q_ml(0,:,:,:) 	= q_ml(0,::-1,:,:)
  rh_ml(0,:,:,:) 	= rh_ml(0,::-1,:,:)
  t_ml(0,:,:,:) 	= t_ml(0,::-1,:,:)  
  z_ml(0,:,:,:) 	= z_ml(0,::-1,:,:)    
  u_ml(0,:,:,:) 	= u_ml(0,::-1,:,:)    
  v_ml(0,:,:,:) 	= v_ml(0,::-1,:,:)    
  w_ml(0,:,:,:) 	= w_ml(0,::-1,:,:)          
  
; PARSE TIME AGAIN JUST FOR KICKS
  time = ml_grb_file->T_GDS0_HYBL@initial_time  
  mo  	= str_get_cols(time,0,1)
  day 	= str_get_cols(time,3,4)  
  year = str_get_cols(time,6,9) 
  hr 	= str_get_cols(time,12,13)
  mn 	= str_get_cols(time,15,16)
  sec	= 0
  
  date 	= ut_inv_calendar(stringtoint(year), stringtoint(mo), stringtoint(day), stringtoint(hr), stringtoint(mn), sec, basetime, 0)  
  date_utc 	= ut_calendar(date,0)   
  
; SETUP COORDINATES FOR 2D VARIABLES 
  ps_sfc!0    			= "time"  
  ps_sfc&time  			= tofloat(date)
  ps_sfc&time@units  		= basetime
  ps_sfc&time@calendar 		= "proleptic_gregorian"
  ps_sfc&time@standard_name 	= "time"

  ps_sfc!1    			= "dimz_PS"  

  ps_sfc!2    			= "dimy_Q"  
  ps_sfc&dimy_Q  			= lat_ml
  ps_sfc&dimy_Q@units  		= "degrees_north"
  ps_sfc&dimy_Q@standard_name 	= "latitude"
  ps_sfc&dimy_Q@_FillValue  	= -999.99  

  ps_sfc!3    			= "dimx_Q"  
  ps_sfc&dimx_Q  			= lon_ml
  ps_sfc&dimx_Q@units  		= "degrees_east"
  ps_sfc&dimx_Q@standard_name 	= "longitude"
  ps_sfc&dimx_Q@_FillValue  	= -999.99  

  ps_sfc@xmin 			= min(lon_ml)
  ps_sfc@xmax 			= max(lon_ml)  
  ps_sfc@xstag 			= 0.0

  ps_sfc@ymin 			= min(lat_ml)
  ps_sfc@ymax 			= max(lat_ml)  
  ps_sfc@ystag 			= 0.0

  ps_sfc@zmin 			= 1050.0
  ps_sfc@zmax 			= 0.0
  ps_sfc@zstag 			= -0.5

  ps_sfc@missing_data 		= ps_sfc@_FillValue

  copy_VarMeta(ps_sfc,z_sfc)
  copy_VarMeta(ps_sfc,sst_sfc)  
  copy_VarMeta(ps_sfc,mslp_sfc)    
  
; SETUP COORDINATES FOR 3D VARIABLES
  t_ml!0    			= "time"  
  t_ml&time  			= tofloat(date)
  t_ml&time@units  		= basetime
  t_ml&time@calendar 	= "proleptic_gregorian"
  t_ml&time@standard_name 	= "time"
  
  t_ml!1    			= "dimz_Q" 
  
  t_ml!2    			= "dimy_Q"  
  t_ml&dimy_Q  		= lat_ml
  t_ml&dimy_Q@units  	= "degrees_north"
  t_ml&dimy_Q@standard_name = "latitude"
  t_ml&dimy_Q@_FillValue  	= -999.99  

  t_ml!3    			= "dimx_Q"  
  t_ml&dimx_Q  		= lon_ml
  t_ml&dimx_Q@units  	= "degrees_east"
  t_ml&dimx_Q@standard_name = "longitude"
  t_ml&dimx_Q@_FillValue  	= -999.99  
  
  t_ml@xmin 			= min(lon_ml)
  t_ml@xmax 			= max(lon_ml)  
  t_ml@xstag 			= 0.0

  t_ml@ymin 			= min(lat_ml)
  t_ml@ymax 			= max(lat_ml)  
  t_ml@ystag 			= 0.0

  t_ml@zmin 			= 1050.0
  t_ml@zmax 			= 0.0
  t_ml@zstag 			= -0.5

  t_ml@missing_data 		= t_ml@_FillValue

  copy_VarMeta(t_ml,q_ml)
  copy_VarMeta(t_ml,u_ml)  
  copy_VarMeta(t_ml,v_ml)    
  copy_VarMeta(t_ml,w_ml)  
  copy_VarMeta(t_ml,z_ml)    
  copy_VarMeta(t_ml,rh_ml)      
  copy_VarMeta(t_ml,eth_ml)      
          
; PRIMARY "P" FILE     
  print("WRITING " + "P"+year+mo+day+"_"+hr)
  outfile_P = addfile(out_DIR + "P"+year+mo+day+"_"+hr+".nc","c")  
  filedimdef(outfile_P, (/"time","dimz_PS","dimz_Q","dimy_Q","dimx_Q"/), (/dims_t,dims_sfc,dims_z,dims_y,dims_x/), (/True,False,False,False,False/)) 
        
; define global attributes:
		outfile_P@domxmin = min(lon_ml)
		outfile_P@domxmax = max(lon_ml)
		outfile_P@domymin = min(lat_ml)
		outfile_P@domymax = max(lat_ml)
		outfile_P@domzmin = 1050.0
		outfile_P@domzmax = 0.0
		outfile_P@constants_file_name = "yotc_cst"
  outfile_P->$("U")$		= u_ml      
  outfile_P->$("V")$		= v_ml        
  outfile_P->$("OMEGA")$	= w_ml          
  outfile_P->$("PS")$	= ps_sfc   
  
  system("mv -f " + out_DIR + "P"+year+mo+day+"_"+hr+".nc" + " " + out_DIR + "P"+year+mo+day+"_"+hr)                 

; SECONDARY "S" FILE     
  print("WRITING " + "S"+year+mo+day+"_"+hr)
  outfile_S = addfile(out_DIR + "S"+year+mo+day+"_"+hr+".nc","c")  
  filedimdef(outfile_S, (/"time","dimz_PS","dimz_Q","dimy_Q","dimx_Q"/), (/dims_t,dims_sfc,dims_z,dims_y,dims_x/), (/True,False,False,False,False/)) 
        
; define global attributes:
		outfile_S@domxmin = min(lon_ml)
		outfile_S@domxmax = max(lon_ml)
		outfile_S@domymin = min(lat_ml)
		outfile_S@domymax = max(lat_ml)
		outfile_S@domzmin = 1050.0
		outfile_S@domzmax = 0.0
		outfile_S@constants_file_name = "yotc_cst"		   
  outfile_S->$("T")$		= t_ml
  outfile_S->$("Q")$		= q_ml
  outfile_S->$("Z")$		= z_ml
  outfile_S->$("RH")$	= rh_ml
  outfile_S->$("THE")$	= eth_ml
  outfile_S->$("ZB")$	= z_sfc
  outfile_S->$("SST")$	= sst_sfc
  outfile_S->$("SLP")$	= mslp_sfc
  outfile_S->$("PS")$	= ps_sfc
  
  system("mv -f " + out_DIR + "S"+year+mo+day+"_"+hr+".nc" + " " + out_DIR + "S"+year+mo+day+"_"+hr)
  
  end do
  
; <><><><><><><><><><><><><><>  
; - GENERATE CONSTANTS FILE -
; <><><><><><><><><><><><><><>  
  if(write_const_file.eq.True)then
  constants_file = "yotc_cst"
  print("CREATING CONSTANTS FILE: " + constants_file)
  outfile_cnst 			= addfile(out_DIR + constants_file + ".nc","c")  
  filedimdef(outfile_cnst, (/"nx","ny","nz"/), (/dims_x,dims_y,dims_z/), (/False,False,False/))   

  lonmin = new(1,"float","No_FillValue")
  lonmax = new(1,"float","No_FillValue")
  latmin = new(1,"float","No_FillValue")
  latmax = new(1,"float","No_FillValue")
  dellon = new(1,"float","No_FillValue")
  dellat = new(1,"float","No_FillValue")  
  pollon = new(1,"float","No_FillValue")  
  pollat = new(1,"float","No_FillValue")      
  starty = new(1,"integer","No_FillValue")    
  startm = new(1,"integer","No_FillValue")      
  startd = new(1,"integer","No_FillValue")      
  starth = new(1,"integer","No_FillValue")      
  starts = new(1,"integer","No_FillValue")        
  dattyp = new(1,"integer","No_FillValue")          
  datver = new(1,"integer","No_FillValue")            
  cstver = new(1,"integer","No_FillValue")  
  
  lonmin = (/min(lon_ml)/)
  lonmax = (/max(lon_ml)/)
  latmin = (/min(lat_ml)/)
  latmax = (/max(lat_ml)/)
  dellon = 0.25
  dellat = 0.25
  pollon = 0.0
  pollat = 90.0
  starty = (/toint(sfc_files_date_utc(ml_ind_files(0),0))/)
  startm = (/toint(sfc_files_date_utc(ml_ind_files(0),1))/)
  startd = (/toint(sfc_files_date_utc(ml_ind_files(0),2))/)
  starth = (/toint(sfc_files_date_utc(ml_ind_files(0),3))/)
  starts = (/toint(sfc_files_date_utc(ml_ind_files(0),4))/)
  dattyp = 2
  datver = 0
  cstver = 0
   
  lonmin!0 = "ncl_scalar"              
  copy_VarMeta(lonmin,lonmax)       
  copy_VarMeta(lonmin,latmin)       
  copy_VarMeta(lonmin,latmax)       
  copy_VarMeta(lonmin,dellon)       
  copy_VarMeta(lonmin,dellat)       
  copy_VarMeta(lonmin,pollon)       
  copy_VarMeta(lonmin,pollat)
         
  starty!0 = "ncl_scalar"                       
  copy_VarMeta(starty,startm)       
  copy_VarMeta(starty,startd)       
  copy_VarMeta(starty,starth)       
  copy_VarMeta(starty,starts)       
  copy_VarMeta(starty,dattyp)       
  copy_VarMeta(starty,datver)       
  copy_VarMeta(starty,cstver)       
    
  outfile_cnst->$("lonmin")$	= lonmin
  outfile_cnst->$("lonmax")$	= lonmax  
  outfile_cnst->$("latmin")$	= latmin
  outfile_cnst->$("latmax")$	= latmax
  outfile_cnst->$("dellon")$	= dellon
  outfile_cnst->$("dellat")$	= dellat
  outfile_cnst->$("pollon")$	= pollon
  outfile_cnst->$("pollat")$	= pollat
  outfile_cnst->$("starty")$	= starty
  outfile_cnst->$("startm")$	= startm
  outfile_cnst->$("startd")$	= startd
  outfile_cnst->$("starth")$	= starth
  outfile_cnst->$("starts")$	= starts
  outfile_cnst->$("dattyp")$	= dattyp
  outfile_cnst->$("datver")$	= datver
  outfile_cnst->$("cstver")$	= cstver
  outfile_cnst->$("aklev")$		= aklev
  outfile_cnst->$("bklev")$		= bklev
  outfile_cnst->$("aklay")$		= aklay
  outfile_cnst->$("bklay")$		= bklay
  
  system("mv -f " + out_DIR + constants_file + ".nc" + " " + out_DIR + constants_file)  
  end if  
  
  end