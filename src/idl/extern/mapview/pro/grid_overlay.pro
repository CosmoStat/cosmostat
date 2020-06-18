;+
; NAME:
;     GRID_OVERLAY
;
; PURPOSE:
;     Plot coordinate grid lines on a map projection.
;     
; CALLING SEQUENCE:
;     grid_overlay, [ PROJ={'M','Z'}, IMCOORD={'G','E','C'}, 
;                    GRCOORD={'G','E','C", COLOR =, /PS, IMAGE=,
;                    LON= , LAT=, /ZBUFF ]
;
; INPUTS:
;     none required, if defaults are satisfactory.
;
; OUTPUTS:  
;     output consists of a plot, either to screen or postscript file.
;
; OPTIONAL INPUT KEYWORDS:
;     proj - char string - Single character specifying the projection type
;                          of the plot.  Either 'M' (Mollweide) or
;                          'Z' (Zenithal Equal Area) are allowed.  Case 
;                          insensitive.   Defaults to 'M'.
; 
;     imcoord - char string - Single character specifying the coordinate 
;                          system of the image projection. Three systems are
;                          recognized: 'E' (Ecliptic J2000), 'G' (Galactic) or
;                          'C' (Celestial J2000).  Case insensitive.
;                          Defaults to 'G'.
;
;     grcoord - char string - Single character specifying the coordinate 
;                          system of the grid overlay. Three systems are
;                          recognized: 'E' (Ecliptic J2000), 'G' (Galactic) or
;                          'C' (Celestial J2000).  Case insensitive.
;                          Defaults to 'G'.
;
;     color - byte -       The color (range 0-255) used to plot the 
;                          grid pattern. Defaults to 0.
;
;     /ps - value=0 or 1 -  Set this keyword to direct the plot to a postscript
;                          file.  The output file is named 'grid_overlay.ps'.  
;                          If not set, plot will be directed to the screen 
;                          (windows device).
;
;     image - bytarr -     A byte-scaled image over which the grid will
;                          be plotted.  If not specfied, the procedure will
;                          plot within the current window (screen) or plot
;                          a grid with no image to the postscript file (/ps).
;
;     lon - fltarr -       lines of constant longitude, in degrees.
;                          Default is every 15 degrees centered on 0.
;                          Specifying lon=-1000. results in NO longitude lines.
;
;     lat - fltarr -       lines of constant latitude, in degrees.
;                          Default is every 10 degrees centered on 0.
;                          Specifying lat=-1000. results in NO latitude lines.
;
;     /zbuff - value=0 or 1 - Set this keyword if plotting to the z buffer
;                          is desired. 
;
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     coortrans, mollweide_xy, zea_xy, scrtv
;     
;
; EXAMPLE:
;     grid_overlay,image=bytscl(kmoll,0,300),/ps
;     
; COMMENTS:
;     The routine scales the overlay to fill the screen window. If the
;     user is trying to overlay an image on the screen which does not
;     fill the window, then the overlay will be scaled incorrectly.
;     The user need not worry about this when creating a postscript image,
;     as the image and overlay are scaled together automatically when
;     producing the .ps file.
;
;     The Zenithal equal area plot has a hard time plotting the lat=0 edge 
;     border.  This can be remedied by plotting lat=[-0.00000001,0.0000000001]
;     in its place (and this is what the code does by default).
;
;     based on G. Hinshaw's grid_overplot.pro & grid_overlay.pro
;
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 04 May 1999
;     added zbuff keyword, JW, 15 Dec 1999
;     changed default lon & lat grid. JP, 09 Dec 2001.
;     Allow for non-X device  W. Landsman   May 2003
;
;-
;======================================================================
;
pro grid_overlay,proj=proj,imcoord=imcoord,grcoord=grcoord, $
                 lon=lon, lat=lat, color=color,ps=ps,image=image, $
                 zbuff=zbuff
;
;
; Projection Type
;
if keyword_set(proj) then begin
   proj = strupcase(strtrim(proj,2))
   if ((proj ne 'M') and (proj ne 'Z')) then begin
     print, 'Only Mollweide (M) or Zenithal Equal-Area (Z) projections allowed.'
     return
   endif
endif else begin
   proj = 'M'                      ;default = Mollweide
endelse
;
; Coordinate Type for image 
;
if keyword_set(imcoord) then begin
   imcoord = strupcase(strtrim(imcoord,2))
   if ((imcoord ne 'E') and (imcoord ne 'G') and (imcoord ne 'C')) then begin
     message, 'Ecliptic(E), Celestial(C) or Galactic(G) coordinates only ',/CON
     return
   endif
endif else begin
   imcoord = 'G'                      ;default = Galactic
endelse
;
; Coordinate Type for grid overlay 
;
if keyword_set(grcoord) then begin
   grcoord = strupcase(strtrim(grcoord,2))
   if ((grcoord ne 'E') and (grcoord ne 'G') and (grcoord ne 'C')) then begin
     message, 'Ecliptic(E), Celestial(C) or Galactic(G) coordinates only ',/CON
     return
   endif
endif else begin
   grcoord = 'G'                      ;default = Galactic
endelse

;
; define lon,lat grid lines in the GRID coordsys, units= degrees
; Note special case allowing NO plot of a lon or lat array if -1000 used
;
if keyword_set(lon) then begin
   lon = float(lon)            ;ensure integer input carries through
   nlon = n_elements(lon)
endif else begin
   ; old default: lon = [-180., -120., -60., 0., 60., 120., 180.]
   lon  = findgen(25)*15. - 180.    ;every 15 degrees
   nlon = n_elements(lon)
endelse
if keyword_set(lat) then begin
   lat = float(lat)            ;ensure integer input carries through
   nlat = n_elements(lat)
endif else begin
   ; old default: lat = [-60., -30., -.0000000001,0.000000001, 30., 60.]
   lat  = [ [findgen(8)*10. -80.], -.0000000001,0.000000001,[findgen(8)*10.+10.] ]
   nlat = n_elements(lat)
endelse
;
; Determine the coordinate transform code for use in coortrans, IF ANY.
;
if (imcoord eq grcoord) then begin
    code = ''
endif else begin
    code = strlowcase(grcoord)+'2'+strlowcase(imcoord)
endelse

;
;
;
; Save device name the user comes in with; restore at end of proc
;
user_device = !d.name
;
;
; plot to chosen device: postscript or window
;
if keyword_set(ps) then begin
  set_plot,'ps'
  device,/land,/color,bits=8,/inches,xs=10,ys=5,yoff=10.5
  device,file='grid_overlay.ps'
  loadct,27
  if keyword_set(image) then tv,image,/inches,xs=10,ys=5
endif else if keyword_set(zbuff) then begin
  set_plot,'z'
  if keyword_set(image) then begin
     size_info = size(image)
     device,z_buff=0,set_resolution=[size_info[1],size_info[2]]
     tv,image
  endif else begin
     device,z_buff=0,set_resolution=[!d.x_size,!d.y_size]
  endelse
endif else begin
  if !VERSION.OS_FAMILY EQ 'Windows' then  set_plot,'WIN' else set_plot,'X'
  if keyword_set(image) then begin
     sz_screen = get_screen_size()
     dimen = size(image,/dimen)
     if N_elements(dimen) LT 2 then message, $
        'ERROR - Image keyword must be a 2-dimension array'
     if (dimen[0] gt sz_screen[0] or dimen[1] gt  sz_screen[1]) then begin
        scrtv,image 
     endif else begin
        window,xs= dimen[0],ys = dimen[1],/free
        tv,image
     endelse
  endif
endelse   

; Color (defaults to !P.COLOR if not set)
;
if N_elements(color) eq 0 then color=!P.COLOR


case proj of

'M': begin
  ;
  ; plot the lines of constant latitude
  ;
  if (lat[0] ne -1000) then begin
  for i = 0,nlat-1 do begin
    phi   = findgen(361)-180.
    theta = replicate(lat[i],361)
    ;
    ; transform coordinates, if necessary
    ;
    if (code ne '') then begin
        coortrans,[[phi],[theta]],ll,code,/lonlat
        phi = reform(ll[*,0])
        sel = where(phi gt 180.)
        if (sel[0] ne -1) then phi[sel] = phi[sel] -360.
        theta = reform(ll[*,1])
    endif
    ;
    ; get Mollweide x,y from lon, lat in degrees
    ;
    Mollweide_XY, phi, theta, X, Y
    ;
    ; Scale to normalized coordinates
    ;
    x = (x +2.)/4.
    y = (y +1.)/2.
    ;
    ; Plot the grid point by point and suppress the line when crossing the 
    ; 180 degree meridian
    PLOTS, X[0], Y[0], /norm, Color=color
    FOR J = 1,n_elements(x)-1 DO BEGIN
        IF( ABS(phi[J]-phi[J-1]) LT 180. )THEN BEGIN
          PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
        ENDIF ELSE BEGIN
          PLOTS, X[J], Y[J], /norm, Color=color 
        ENDELSE
    ENDFOR
  endfor
  endif

  ;
  ; plot the lines of constant longitude
  ;
  if (lon[0] ne -1000) then begin
  for i = 0,nlon-1 do begin
    phi   = replicate(lon[i],181)
    theta = findgen(181) - 90.
    ;
    ; transform coordinates, if necessary
    ;
    if (code ne '') then begin
        coortrans,[[phi],[theta]],ll,code,/lonlat
        phi = reform(ll[*,0])
        sel = where(phi gt 180.)
        if (sel[0] ne -1) then phi[sel] = phi[sel] -360.
        theta = reform(ll[*,1])
    endif
    ;
    ; get Mollweide x,y
    ;
    Mollweide_XY, phi, theta, X, Y
    ;
    ; Scale to normalized coordinates
    ;
    x = (x +2.)/4.
    y = (y +1.)/2.
    ;
    ; plot to chosen device, which has already been established
    ;
    PLOTS, X[0], Y[0], /norm, Color=color
    FOR J = 1,n_elements(x)-1 DO BEGIN
        IF( ABS(phi[J]-phi[J-1]) LT 180. )THEN BEGIN
          PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
        ENDIF ELSE BEGIN
          PLOTS, X[J], Y[J], /norm, Color=color 
        ENDELSE
    ENDFOR
  endfor
  endif
  end   ; end Mollweide case


'Z': begin
  ;
  ; plot the lines of constant latitude
  ;
  if (lat[0] ne -1000) then begin
  for i = 0,nlat-1 do begin
    phi   = findgen(361)-180.
    theta = replicate(lat[i],361)
    ;
    ; transform coordinates, if necessary
    ;
    if (code ne '') then begin
        coortrans,[[phi],[theta]],ll,code,/lonlat
        phi = reform(ll[*,0])
        sel = where(phi gt 180.)
        if (sel[0] ne -1) then phi[sel] = phi[sel] -360.
        theta = reform(ll[*,1])
    endif
    ;
    ; given lon, lat in degrees, get x in [-2,2] and y in [-1,1]
    ;
     zea_xy,phi,theta,x,y
    ;
    ; Scale to normalized coordinates
    ;
    x = (x +2.)/4.
    y = (y +1.)/2.
    ;
    ; Plot the grid point by point and suppress the line when crossing the 
    ; equator
    PLOTS, X[0], Y[0], /norm, Color=color
    FOR J = 1,n_elements(x)-1 DO BEGIN
        IF(theta[J]*theta[J-1] gt 0) THEN BEGIN
          PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
        ENDIF ELSE BEGIN
          PLOTS, X[J], Y[J], /norm, Color=color 
        ENDELSE
    ENDFOR
  endfor
  endif
  ;
  ; plot the lines of constant longitude
  ;
  if (lon[0] ne -1000) then begin
  for i = 0,nlon-1 do begin
    phi   = replicate(lon[i],181)
    theta = findgen(181) - 90.
    ;
    ; transform coordinates, if necessary
    ;
    if (code ne '') then begin
        coortrans,[[phi],[theta]],ll,code,/lonlat
        phi = reform(ll[*,0])
        sel = where(phi gt 180.)
        if (sel[0] ne -1) then phi[sel] = phi[sel] -360.
        theta = reform(ll[*,1])
    endif
    ;
    ; given lon, lat in degrees, get x in [-2,2] and y in [-1,1]
    ;
    zea_xy,phi,theta,x,y
    ;
    ; Scale to normalized coordinates
    ;
    x = (x +2.)/4.
    y = (y +1.)/2.
    ;
    ;
    ; Plot the grid point by point and suppress the line when crossing the 
    ; equator
    PLOTS, X[0], Y[0], /norm, Color=color
    FOR J = 1,n_elements(x)-1 DO BEGIN
        IF(theta[J]*theta[J-1] gt 0) THEN BEGIN
          PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
        ENDIF ELSE BEGIN
          PLOTS, X[J], Y[J], /norm, Color=color 
        ENDELSE
    ENDFOR
  endfor
  endif
  end   ; end zenithal equal area case

endcase

if keyword_set(ps) then device,/close
;
; Reset device to user value prior to calling this routine
;
set_plot, user_device
;
return
end

