;-------------------------------------------------------------------------------------;
;		Trajectory Plot
; Assumes Lagranto lsl-file style output with five
; header lines
;-------------------------------------------------------------------------------------;

;-------------------------------------------------------------------------------------;
; Loading Libraries
;-------------------------------------------------------------------------------------;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;-------------------------------------------------------------------------------------;
; Read trajectory file 
;-------------------------------------------------------------------------------------;

; --------- Name of trajectory file --------------------------------------------------;
filename="trajectory.lsl"                          

;---------- Read file ----------------------------------------------------------------;
delim	= " "                                              ; Delimiting character

nhead   = 5                                                ; # header lines
head	= readAsciiHead(filename,nhead-1)		   ; Read header lines 
ncol    = str_fields_count(head(2),delim)                  ; # columns

data	= asciiread(filename,-1,"float")		   ; Read data in as float
n       = dimsizes(data);                                  ; Number of entries in data
ifirst  = 4                                                ; First data element
ifinal  = n-1                                              ; Last data element

timestp = data(ifirst+ncol) - data(ifirst)                 ; Time step of trajectories
period  = data(ifinal-ncol+1)                              ; Time period of trajectories
ntim    = floattoint( period/timestp + 1 )                 ; # times
ntra    = floattoint( (n-ifirst) / (ntim * ncol ) )        ; # trajectories

traj    = onedtond(data(ifirst:ifinal),(/ntra,ntim,ncol/)) ; 3d trajectory array 
vars    = str_split(head(2),delim)                         ; Field names (columns)

;-------------- Add meta information --------------------------------------------------;

; Add some attributes to the variable
traj@HOUR        = stringtoint( str_get_cols(head(0),22,23) )
traj@DAY         = stringtoint( str_get_cols(head(0),19,20) )
traj@MONTH       = stringtoint( str_get_cols(head(0),17,18) )
traj@YEAR        = stringtoint( str_get_cols(head(0),13,16) )
traj@FIELD_NAMES = vars
traj@MODEL       = "Lagranto"

; Set names of dimensions and assign values
traj!0    = "TRA"
traj!1    = "TIME"
traj!2    = "FIELD"
traj&TIME = fspan(0,period,ntim)

;-------------------------------------------------------------------------------------;
; Plot trajectories
;-------------------------------------------------------------------------------------;

; --------------- Open Workstation and Colormap --------------------------------------

wks  = gsn_open_wks("png","tra")		   ; open workstation
gsn_define_colormap(wks,"temp_19lev")		   ; use the BlueDarkRed18 colormap
cmap = gsn_retrieve_colormap(wks)	           ; colormap

; --------------- General resources for plot -----------------------------------------

mres		          = True	           ; map resource
mres@gsnDraw    	  = False                  ; don't draw
mres@gsnFrame   	  = False                  ; don't advance frame
mres@tiMainString         = "Trajectories"	   ; set the main title
mres@tiMainFont	          = "0_times_roman"        ; font of main title
mres@gsnMaximize	  = True		   ; Maximize plot


; ----------------- World map --------------------------------------------------------

mres@mpMaxLatF            = 70   		   ; choose subregion           
mres@mpMinLatF            = 0                      ; of world map
mres@mpMaxLonF            = 100
mres@mpMinLonF            = 0
mres@mpFillOn	          = False                  ; color land
mres@mpOutlineOn	  = True                   ; outline land/sea boarders

map = gsn_csm_map_ce(wks,mres)                     ; Draw world map 
draw(map)

; ----------------- Polyline (including coloring of trajectories ) --------------------


; --- Specify the coloring field
cnField                  = "p"                          ; Coloring according to pressure
cnIndex                  = 3                            ; Column index for coloring field
cnLevels 		 = fspan(200,1000,17)		; Coloring levels
cnLabels                 = flt2string(cnLevels)         ; Labels


; --- Resources for the polylines (except for color)
pres		  	 = True                         ; polyline resource
pres@tfDoNDCOverlay      = True 
pres@gsLineThicknessF    = 3.0                          ; line thickness    


; --- Resources for the colorbar/labelbar
lres                          = True
lres@vpWidthF                 = 0.60                    ; Width of Labelbar
lres@vpHeightF                = 0.10                    ; Height of Labelbar
lres@lbPerimOn                = False                   ; Turn off perimeter.
lres@lbOrientation            = "Horizontal"	        ; Default is "Vertical"
lres@amLabelBarOrthogonalPosF = -1.0
lres@lbLabelAlignment         = "InteriorEdges"         ; Default is "BoxCenters".
lres@lbFillColors             = cmap(2:,:)              ; Colors for boxes.
lres@lbMonoFillPattern        = True                    ; Fill them all solid.
lres@lbLabelFontHeightF       = 0.015                   ; Font Height
lres@lbLabelAutoStride        = True                    ; Auto correct labels
lres@vpYF                     = 0.15                    ; location of left edge
lres@vpXF                     = 0.22                    ; location of top edge


; --- Draw polylines 
do i = 0,ntra-1                                         ; Loop over all trajecories

  x = traj(i,0:(ntim-1),1)                              ; Longitudes for trajectory i
  y = traj(i,0:(ntim-1),2)                              ; Latitudes for trajectory i

  do k = 0,ntim-2                                       ; Loop over all times

      pres@gsLineColor = GetFillColor(cnLevels,cmap,traj(i,k,cnIndex))   ; color of polyline
      gsn_polyline  (wks,map,(/x(k),x(k+1)/),(/y(k),y(k+1)/),pres)       ; draw polyline
 
  end do	

end do


; --- Draw colorbar/labelbar
gsn_labelbar_ndc(wks,dimsizes(cnLevels)+1,cnLabels,0.2,0.2,lres)


; ----------------- Markers ----------------------------------------------------------------

; --- Resources for the markers
kres			 = True				; marker resource
kres@gsMarkerSizeF       = 5.0         			; marker size
kres@gsMarkerColor       = "black"			; marker color


; --- Draw markers
do i = 0,ntra-1                                         ; Loop over all trajecories

  x = traj(i,0:(ntim-1),1)                              ; Longitudes for trajectory i
  y = traj(i,0:(ntim-1),2)                              ; Latitudes for trajectory i

  do k = 0,ntim-2,5                                       ; Loop over all times

      gsn_polymarker(wks,map,x(k),y(k),kres)		; mark position
 
  end do	

end do


; ----------------- Advance frame ----------------------------------- --------------------

frame(wks) 

