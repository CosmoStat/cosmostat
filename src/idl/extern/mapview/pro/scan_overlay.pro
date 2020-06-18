;+
; NAME:
;     SCAN_OVERLAY
;
; PURPOSE:
;     Plot a single hour of WMAP's scan path based on quaternion data.
;     The user specifies the GMT time at which the hour starts, and must first
;     read in an IDL data structure from the appropriate TOD file.
;     
; CALLING SEQUENCE:
;     SCAN_OVERLAY, gmt_start, tod, [,duration=duration, PARFILE= ' ', $
;                  proj={'M','Z'},coord={'E','G','C'} ,color=num /ps, $
;                  image=image, side={'A','B'} ,da=' ', /Noannot
;
; INPUTS:
;     gmt_start - char string - The GMT time at which the scan trace begins.
;                         e.g., '20023070000' (year 2002 day 307 00H 00M)
; 
;     tod - IDL structure, read from a time-ordered data file, e.g. using
;            READ_FITS_TOD.   
; OUTPUTS:  
;     output consists of a plot, either to screen or postscript file.
;
; INPUT KEYWORDS:
;     duration - scalar float - Time interval over which to plot the
;                               scan pattern, specified in minutes
;                               Default = 60 minutes (1hr).
;
;     parfile - char string - Name of the programs.pars file to use for
;                          DA LOS positions and TS reference.
;                          Default = '$MAP_REF/programs.pars'
;
;     proj - char string - Single character specifying the projection type
;                          of the plot.  Either 'M' (Mollweide) or
;                          'Z' (Zenithal Equal Area) are allowed.  Case 
;                          insensitive.    Defaults to 'M'.
; 
;     coord - char string - Single character specifying the coordinate 
;                          system of the projection. Three systems are
;                          recognized: 'E' (Ecliptic J2000), 'G' (Galactic) or
;                          'C' (Celestial J2000).  Case insensitive.
;                          Defaults to 'G'.
;
;     color - byte -       The color (range 0-255) used to plot the 
;                          scan pattern. Defaults to 0.
;
;     ps -                 Set this keyword to 1 to write the output to a 
;                          postscript file with the name of'scan_overlay.ps'.  
;                          Alternatively, one can supply the name of the 
;                          postscript file, e.g. ps='scanKA.ps'.   Default
;                          is to plot to the screen  (windows device).
;
;     image - bytarr -     A byte-scaled image over which the scan pattern will
;                          be plotted.  If not specfied, the procedure will
;                          plot within the current window (screen) or plot
;                          a scan pattern with no image to the postscript 
;                          file (/ps).
;
;     side - char string - Single character specifying which side of
;                          the spacecraft to plot the scan pattern for: choose
;                          either 'A' or 'B'.  Only one side may be plotted
;                          at a time.  Defaults to 'A'. 
;
;     da - char string -   One of ten valid character strings specifying
;                          that the LOS for a specific differencing assembly
;                          should be used rather than the mean spacecraft LOS.
;                          A valid string is one of: 'K1','Ka1','Q1','Q2','V1',
;                          'V2','W1','W2','W3','W4'.  May be combined with
;                          the side keyword.
;     /Noannot -          Set this keyword to *not* annotate the overlay with
;                         the starting GMT 
;
; ROUTINES CALLED:
;     gmt2ts, timestamp_addtime, ts2gmt, gmt2jul, quat_to_sky_coords
;     q2m, coortrans, mollweide_xy, zea_xy, scrtv, load_map_params
;     
;
; EXAMPLE:
;     Overlay 1 hour of WMAP's Side A scan path beginning at 00 GMT on day 217 
;     of 2002 on a 1024 x 512 Mollweide WMAP K-band year 1 skymap.
;
;     IDL> fits_read_map,'map_k1_imap_yr1_v1.fits',t
;     IDL> reproj_healpix,t,kmoll,coord=1,proj=1,size=2
;     IDL> fits_read_tod, 'MAP_tod_20022162357_20022172357.fits', tod
;     IDL> scan_overlay,'20022171000',tod, image=bytscl(kmoll,0,300),da='K1'
;     
; COMMENTS:
;     The routine scales the overlay to fill the screen window. If the
;     user is trying to overlay an image on the screen which does not
;     fill the window, then the overlay will be scaled incorrectly.
;     The user need not worry about this when creating a postscript image,
;     as the image and overlay are scaled together automatically when
;     producing the .ps file.
;     
; MODIFICATION HISTORY:
;     initial version J. Weiland, 04 May 1999
;     Major rewrite to work with WMAP archival data  W. Landsman Feb. 2003
;     Use CONCAT_DIR() to concatenate directories, allow for non-X device
;        W. Landsman May 2003
;-
;======================================================================
;
pro scan_overlay,gmt_start,tod,duration=duration, file_prefix=file_prefix, $
                 parfile=parfile, NoAnnot= NoAnnot, $
                 proj=proj,coord=coord,color=color,ps=ps,image=image, $
                 side=side,da=da
;
;
 if (N_params() LT 2 ) then begin
    print,'Syntax: SCAN_OVERLAY, gmt_start, tod, [,optional keywords] '
    return
  endif
;
; Pad input gmt start time out to full 20 character length
;
 full_gmt_start = '00000000000000000000'
 strput,full_gmt_start,gmt_start,0        ;works regardless of input length

;
; default for Prefix for files (includes disk location)
;
 mref = strtrim(getenv('MAP_REF'),  2)
 if not keyword_set(parfile) then parfile = concat_dir( mref, 'programs.pars')
;
; Load in the programs.pars file of choice
;
 load_map_params,param_str,file=parfile
;
; Establish valid da list
;
 da_list = ['K1','KA1','Q1','Q2','V1','V2','W1','W2','W3','W4']
;
 if N_elements(proj) EQ 0 then proj = 'M'               ;default = Mollweide
 if N_elements(coord) EQ 0 then coord = 'G'             ;default = Galactic
;
; Color (defaults to 0 if not set)
;
if N_elements(color) eq 0 then color=0
;
; SC side defaults to 'A' if not specified
;
if N_elements(side) eq 0 then side='A'

;
; Save device name the user comes in with; restore at end of proc
;
user_device = !d.name

;
; Determine the end_time of the scan (start + duration)
;

 N = N_elements(tod.sci.day)
 ts_tod = transpose([[reform(tod.sci.day,n)],[reform(tod.sci.time,n)]])
 jd_tod = ts2jul(ts_tod)

 if N_elements(duration) EQ 0 then duration = 60.0
 jd1 = GMT2JUL(full_gmt_start)
 jd2 = jd1 + duration/60.0d/24.0d

 g = where( (jd_tod GE jd1) and (jd_tod LE jd2), Ngood)
 if Ngood EQ 0 then begin
     message,/INF,'No valid data within TOD structure'
     print,'Times must be between ',jul2gmt(jd_tod[0]), $
                            ' and ' ,jul2gmt(jd_tod[n-1])
     return
 endif
;
; Define pointing direction: S/C side and (mean SC LOS or DA LOS)
;
  quat = reform( tod[*].quaternions[*,1:30], 4, N)
  quat = quat[*,g]
  side = strupcase(strtrim(side,2))

 if not keyword_set(da) then begin

; Use mean optical axis

  if (side eq 'B') then $
     dir_SC = [0.d0,-0.94264149d0, -0.33380686d0]  else $
     dir_SC = [0.d0, 0.94264149d0, -0.33380686d0] 
 

 endif else begin
  my_da = strupcase(strtrim(da,2))
  da_index = where(da_list eq my_da, Nda)
  if (Nda EQ 0) then begin
     message,'Oops! That is not a valid DA!',/CON
     return
  endif
  if (side eq 'B') then $
     dir_SC = reform(param_str.Dir_B_LOS[*,da_index]) else $
     dir_SC = reform(param_str.Dir_A_LOS[*,da_index])

endelse

; Extract astronomical coordinates from input quaternions

case coord of

 'E':  QUAT_TO_SKY_COORDS, quat, da, ecl=out_uvec, side = side
 'G':  QUAT_TO_SKY_COORDS, quat, da, gal=out_uvec, side = side
 'C':  QUAT_TO_SKY_COORDS, quat, da, cel=out_uvec, side = side
 else: begin
        message,/CON,"ERROR - Coordinates must be 'G','E', or 'C'"
	return
       end
endcase

   lon = out_uvec[*,0]
   lat = out_uvec[*,1]

  ; plot to chosen device: postscript or window
  ;

  if N_elements(ps) GT 0 then begin
      if size(ps,/tname) EQ 'STRING' then psname = ps else $
          if keyword_set(ps) then psname = 'scan_overlay.ps'
  endif else psname = ''
  do_ps = strlen(psname) GT 0      
  if do_ps then begin
    set_plot,'ps'
    device,/land,/color,bits=8,/inches,xs=10,ys=5,yoff=10.5
    device,file=psname
    loadct,27
    message,/INF,'Writing image to a file ' + psname    
    if keyword_set(image) then tv,image,/inches,xs=10,ys=5
  endif else begin
  
    if !VERSION.OS_FAMILY EQ 'Windows' then  set_plot,'WIN' else set_plot,'X'

    if keyword_set(image) then begin
       dimen = size(image,/dimen)
       device,get_screen = screen
       if (dimen[0] gt screen[0] or dimen[1] gt screen[1]) then begin
          SCRTV,image 
       endif else begin
          window,xs= dimen[0],ys = dimen[1],/free  
          plot,indgen(10),indgen(5),/nodata,xstyle=5,ystyle=5, $
	       xmar=0.0,ymar=0.0
          tv,image
       endelse
    endif
 endelse
 
 case proj of       ;Convert to X,Y coordinates using specified projection

 'M':   mollweide_xy,lon,lat,x,y
 'Z':   zea_xy,lon,lat,x,y
 else: begin
     message,/CON, $
       'Only Mollweide (M) or Zenithal Equal Area (Z) projections allowed.'
     return
     end
 endcase
  ;

; Scale to normalized coordinates
;
  x = (x +2.)/4.
  y = (y +1.)/2.

; Plot the grid point by point and suppress the line when crossing the 
; equator (ZEA projection) or when crossing the 180 degree meridian 
; (Mollweide projection).

 if proj EQ 'M' then $
    plotcontin = abs(lon -shift(lon,1)) LT 180.0 else $
    plotcontin = lat*shift(lat,1) GT 0

    PLOTS, X[0], Y[0], /norm, Color=color
    FOR J = 1,n_elements(x)-1 DO BEGIN
        IF plotcontin[j] THEN $
          PLOTS, X[J], Y[J], /norm, /Continue, Color=color $ 
        ELSE $
          PLOTS, X[J], Y[J], /norm, Color=color 
    ENDFOR
 if not keyword_set(Noannot) then begin 
; Make sure annotation color differs from the background
     acol = !D.TABLE_SIZE-1 
     if keyword_set(image) then $
          if image[0,0] EQ acol then acol = 0
     xyouts,0.7,0.01,'GMT ' + full_GMT_start,charsize=1.4, $
                      charthick=2,/norm,color=acol
 endif 

 if do_ps then device,/close
;
; Reset device to user value prior to calling this routine
;
 set_plot, user_device
;

return
end

