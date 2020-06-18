pro glactc,ra,dec,year,gl,gb,j, degree=degree
;+
; NAME:  
;       GLACTC
; PURPOSE:
;        Convert between celestial and Galactic coordinates.
; EXPLANATION:
;       Program to convert right ascension (ra) and declination (dec) to
;       Galactic longitude (gl) and latitude (gb) (j=1) or vice versa (j=2).
;
; CALLING SEQUENCE: 
;       glactc, ra, dec, year, gl, gb, j, [ /DEGREE ]
;
; INPUT PARAMETERS: 
;       year     equinox of ra and dec, scalar       (input)
;       j        direction of conversion     (input)
;               1:  ra,dec --> gl,gb
;               2:  gl,gb  --> ra,dec
;
; INPUTS OR OUTPUT PARAMETERS: ( depending on argument J )
;       ra       Right ascension, hours (or degrees if /DEGREES is set), 
;                         scalar or vector
;       dec      Declination, degrees,scalar or vector
;       gl       Galactic longitude, degrees, scalar or vector
;       gb       Galactic latitude, degrees, scalar or vector
;
;       All results forced double precision floating.
;
; KEYWORD PARAMETER:
;       /DEGREE - If set, then the RA parameter (both input and output) is 
;                given in degrees rather than hours. 
; COMMON BLOCKS:   
;      gal      See Side Effects.    
;
; SIDE EFFECTS:
;       Year and galaxy orientation saved in common to make repeated     
;       computations more efficient.
;
; EXAMPLES:
;       Find the Galactic coordinates of Altair (RA (2000): 19,50,47 
;       Dec (2000): 08 52 06)
;
;       IDL> glactc, ten(19,50,47),ten(8,52,6),2000,gl,gb,1
;       ==> gl = 47.74, gb = -8.91
;
; HISTORY: 
;       FORTRAN subroutine by T. A. Nagy, 21-MAR-78.
;       Conversion to IDL, R. S. Hill, STX, 19-OCT-87.
;       Modified to handle vector input, E. P. Smith, GSFC, 14-OCT-94
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Added DEGREE keyword, C. Markwardt, Nov 1999
;-
common gal,oldyr,rapol,decpol,dlon,sdp,cdp,radhrs
if N_params() lt 6 then begin
     print,'Syntax -  glactc, ra, dec, year, gl, gb, j
     print,'j = 1: ra,dec --> gl,gb   j = 2:  gl,gb -->ra,dec
     return
endif
sz = size(oldyr) 
first = sz[1] eq 0 
new = 0
if not first then new=(double(year) ne oldyr)
radeg = 180.0d/!DPI
if first or new then begin
   oldyr = double(year)
   ;
   ; Galactic pole at ra 12 hrs 49 mins, dec 27.4 deg, equinox 1950.0
   ; position angle from galactic center to equatorial pole = 123 degs.
   rapol = 12.0d0 + 49.0d0/60.0d0 + 8.13d-4 * (year-1950.0d0)
   decpol = 27.4d0 - 5.44d-3 * (year-1950.0d0)
   dlon = 123.0d0 - 1.33d-3 * (year-1950.0d0)
   sdp = sin(decpol/radeg)
   cdp = sqrt(1.0d0-sdp*sdp)
   radhrs=radeg/15.0d0
endif
;
; Branch to required type of conversion.
case j of                   
    1:  begin
        sdec = sin(dec/radeg)
        cdec = sqrt(1.0d0-sdec*sdec)
        if keyword_set(degree) then      ras = ra/radeg - rapol/radhrs $
        else                             ras = (ra-rapol)/radhrs
        sgb = sdec*sdp + cdec*cdp*cos(ras)
        gb = radeg * asin(sgb)
        cgb = sqrt(1.0d0-sgb*sgb)
        sine = cdec * sin(ras) / cgb
        cose = (sdec-sdp*sgb) / (cdp*cgb)
        gl = dlon - radeg*atan(sine,cose)
        ltzero=where(gl lt 0.0, Nltzero)
        if Nltzero ge 1 then gl[ltzero]=gl[ltzero]+360.0d0
        return
        end
    2:  begin
        sgb = sin(gb/radeg)
        cgb = sqrt(1.0d0-sgb*sgb)
        sdec = sgb*sdp + cgb*cdp*cos((dlon-gl)/radeg)
        dec = radeg * asin(sdec)
        cdec = sqrt(1.0d0-sdec*sdec)
        sinf = cgb * sin((dlon-gl)/radeg) / cdec
        cosf = (sgb-sdp*sdec) / (cdp*cdec)
        ra = rapol + radhrs*atan(sinf,cosf)
        gt24 = where(ra gt 24.0, Ngt24)
        if Ngt24 ge 1 then ra[gt24] = ra[gt24] - 24.0d0
        if keyword_set(degree) then      ra = ra * 15.0D0
        return
        end
endcase
end
