;;+
; NAME:
;     PLANET_OVERLAY
;
; PURPOSE:
;     Plot a cross at the position of the chosen planet over an existing
;     image within a screen window.
;     
; CALLING SEQUENCE:
;    planet_overlay,gmt=gmt,jd=jd [,planet=planet] [,coord={'G','E','C'}] 
;                    [,proj={'M','Z'}] [,color=color]
;
; INPUTS:
;    At minimum, use must specify a time, using either the GMT or JD
;    keywords (but not both).     
;
; OUTPUTS:  
;     output consists of a plot to screen.
;
; KEYWORDS:
;     gmt - char string  - MAP GMT time for which the planet position is
;                          to be computed, e.g., '2000345000000'. Valid
;                          GMTS lie between 1 Nov 2000 and 20 April 2003.
;                          Do not specify GMT if you are already using the
;                          JD keyword.
;
;     jd - scalar or vector - Julian date for which the planet position is
;                          to be computed.  Only JDs between 2451849.5 and
;                          2452749.5 are accepted.  It is acceptable to
;                          input a MAP reduced Julian date with values
;                          between 1849.4 and 2749.5. 
;                          Do not specify JD if you are already using the
;                          GMT keyword.
;
;     planet - char string - Name of the planet.  Only one planet may
;                           be specified.  Valid planet names are
;                          'mars', 'jupiter','saturn','uranus','neptune'.
;                          Default planet = 'jupiter'.
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
;
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     is_ieee_big(), gmt2jul, coortrans, mollweide_xy, zea_xy
;     
; EXTERNAL FILES REFERENCED:
;     The FITS file JPLEPH.405 $MAP_REF/planet_overlay/ provides the JPL
;     DE405 ephemeris (as Chebyshev polynomials). 
;
; EXAMPLE:
;     To overplot Jupiter's positions for 500 days on a Galactic Mollweide
;     projection:
;       IDL> jd0= 1850.0
;       IDL> planet_overlay,jd=jd0+dindgen(500)
;     
; COMMENTS:
;     The routine scales the overlay to fill the screen window. If the
;     user is trying to overlay an image on the screen which does not
;     fill the window, then the overlay will be scaled incorrectly.
;     
; PROCEDURES USED:
;     CONCAT_DIR(), EULER,  MOLLWEIDE_XY, PLANET_COORDS, ZEA_XY
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 07 May 1999
;     conversion for SGI, JW, 26 Apr 2000 
;     use getenv to extract environment variables.  MRG, SSAI, 20 August 2002.
;     Let JD accept vector values  WBL  December 2002
;     Use CONCAT_DIR() to concatanate arrays  WBL    May 2003
;-
;======================================================================
pro planet_overlay,gmt=gmt,jd=jd,planet=planet,coord=coord,proj=proj, $
                   color=color
;
; Plots a cross where the requested planet is on the requested date.
;
; Image on screen MUST fill the plot window in order for this to work.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; Get the time
;
if keyword_set(gmt) then begin
   full_gmt = '00000000000000000000'
   strput,full_gmt,gmt,0        ;works regardless of input length
   jd = 2450000.D0 + gmt2jul(gmt)
endif else if keyword_set(jd) then begin
   jd = double(jd)
   red = where(jd LT 2450000., Nred)
   if Nred GT 0 then jd[red] =2450000.D0 + jd[red]       ;handle reduced JD 
endif else begin
   print,"Syntax - PLANET_OVERLAY, [gmt=,jd= ,planet= ,coord={'G','E','C'}]" 
   print,'Either a gmt or JD must be specified'
   return
endelse
;
; Valid times lie between 1 Nov 2000 and 20 Apr 2003
;
minjd = min(jd,max= maxjd)
if (minjd lt 2451543.5d0  OR maxjd gt 2455196.5d0) then $
   message,/warn,'Supplied dates outside of the range 2000-2010'
   
;
; planet id
;
if keyword_set(planet) then $
   planet = strupcase(strtrim(planet,2)) else $
   planet = 'JUPITER'                      ;default = jupiter

;
; Coordinate Type and coortrans code defn
;
if keyword_set(coord) then begin
   coord = strlowcase(strtrim(coord,2))
   if ((coord ne 'e') and (coord ne 'g') and (coord ne 'c')) then begin
     print, 'Ecliptic(E), Celestial(C) or Galactic(G) coordinates only '
     return
   endif
endif else begin
   coord = 'g'                      ;default = Galactic
endelse

if coord ne 'c' then transcode = 'c2'+coord
if coord eq 'c' then transcode = 'u2ll'

;
; Projection Type
;
if keyword_set(proj) then begin
   proj = strupcase(strtrim(proj,2))
   if ((proj ne 'M') and (proj ne 'Z')) then begin
     print, 'Only Mollweide (M) or Zenithal Equal Area (Z) projections allowed.'
     return
   endif
endif else begin
   proj = 'M'                      ;default = Mollweide
endelse

;
; Color (defaults to 0 if not set)
;
if keyword_set(color) eq 0 then color=0

jplfile = concat_dir('$MAP_REF', concat_dir('planet_overlay','JPLEPH.405'))
planet_coords,jd,ra,dec,planet=planet,jpl=jplfile,/jd
;
case coord of 
'g': euler,ra,dec,lon,lat,1
'e': euler,ra,dec,lon,lat,3
else: begin
     lon = ra & lat = dec
      end
endcase


case proj of

'M': begin
  ;
  mollweide_xy,lon,lat,x,y
  ;
  ; Scale to normalized coordinates
  ;
  x = (x +2.)/4.
  y = (y +1.)/2.
  ;
  ; Place a cross at the central position
  ; x:y = 2:1
  ;
  for i=0,N_elements(x)-1 do begin
  plots,[x[i]-.005,x[i]+.005],[y[i],y[i]],/norm,color=color
  plots,[x[i],x[i]],[y[i]-.01,y[i]+.01],/norm,color=color
  endfor
  
  end

'Z': begin
  ;
  ; given lon, lat in degrees, get x in [-2,2] and y in [-1,1]
  ;
  zea_xy,lon,lat,x,y
  ;
  ; Scale to normalized coordinates
  ;
  x = (x +2.)/4.
  y = (y +1.)/2.
  ;
  ; Place a cross at the central position
  ; x:y = 2:1
  ;
  for i=0,N_elements(x)-1 do begin
  plots,[x[i]-.005,x[i]+.005],[y[i],y[i]],/norm,color=color
  plots,[x[i],x[i]],[y[i]-.01,y[i]+.01],/norm,color=color
  endfor

  end

endcase

;
end


