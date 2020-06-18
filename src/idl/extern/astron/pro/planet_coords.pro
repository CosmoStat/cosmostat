pro planet_coords, date, ra, dec, planet=planet, jd = jd
;+
; NAME:
;    PLANET_COORDS
; PURPOSE:  
;    Find low-precision RA and DEC for the planets given a date
;
; EXPLANATION:
;    This routine uses HELIO to get the heliocentric ecliptic coordinates
;    of the planets at the given date, then converts these to geocentric
;    ecliptic coordinates ala "Astronomical Algorithms" by Jean Meeus
;    (1991, p 209). These are then converted to RA and Dec using EULER.
;    The accuracy between the years 1800 and 2050 is better than 1 arcminute for 
;    the terrestial planets, but reaches 10 arcminutes for Saturn.    Before
;    1850 or after 2050 the accuracy can get much worse.   
;
; CALLING SEQUENCE:
;    PLANET_COORDS, DATE, RA, DEC, [ PLANET = , /JD]
;
; INPUTS:
;       DATE - If /JD is not set, then date is a 3-6 element vector containing
;              year,month (1-12), day, and optionally hour, minute, & second.
;              If /JD is set then DATE is a Julian date.   An advantage of the
;              /JD option is that it allows the use of vector dates.
; OUTPUTS:
;       RA - right ascension of planet(s), J2000 degrees
;       DEC - declination of   planet(s), J2000 degrees
;
; OPTIONAL INPUT KEYWORD:
;       PLANET - scalar string giving name of a planet. Default is coords for 
;               all of them except Earth.
;       /JD - If set, then the date parameter should be supplied as Julian date
; EXAMPLES:
;    (1)  Find the RA, Dec of Venus on 1992 Dec 20
;          IDL> planet_coords, [1992,12,20], ra,dec    ;Compute for all planets
;          IDL> print,adstring(ra[1],dec[1],1)         ;Venus is second planet
;     ====> RA = 21 05  2.66  Dec = -18 51 45.7
;    This position is 40" from the full DE2000 ephemeris position of
;          RA = 21 05  5.38        -18 51 35.6
;
;    (2) Return the current RA and Dec of all 8 planets
;          IDL> get_juldate, jd                 ;Get current Julian Date
;          IDL> planet_coords,jd,ra,dec,/jd     ;Find positions of all planets
;          IDL> forprint,adstring(ra,dec,0)     ;Display positions   
;
;    (3) Plot the declination of Mars for every day in the year 2001
;          IDL> jdcnv,2001,1,1,0,jd      ;Get Julian date of midnight on Jan 1 
;               Now get Mars RA,Dec for 365 consecutive days
;          IDL> planet_coords,jd+indgen(365),ra,dec,/jd, planet = 'mars'     
;          IDL> plot,indgen(365)+1,dec
; NOTES:
;       (1) HELIO is based on the two-body problem and neglects interactions 
;           between the planets.   This is why the worst results are for
;           Saturn.  See http://ssd.jpl.nasa.gov/cgi-bin/eph for a more 
;           accurate ephemeris generator online.
; PROCEDURES USED:
;        EULER, JULDATE, HELIO  
;
; REVISION HISTORY:
;        Written P.Plait & W. Landsman     August 2000
;        Fixed Julian date conversion   W. Landsman August 2000
;-   
 On_error,2
 if N_params() LT 1 then begin
     print,'Syntax - planet_coords, date, ra,dec, [PLANET =, /JD ]'
     print,'      date - either 3-6 element date or Julian date (if /JD is set)'
     print,'      ra,dec - output ra and dec in degrees'
     print,'      PLANET - name of planet (optional)'
     return
 endif

;convert input date to real JD
  
  if keyword_set(jd) then begin 
       jj = date
       if N_elements(jj) GT 0 then if N_elements(planet) GT 1 then $
       message,'ERROR - A planet name must be supplied for vector dates'
  endif else begin 
           juldate,date,jj
            jj = jj + 2400000.0d
  endelse 

;make output arrays to include each planet
; note that we need Earth to convert from heliocentric
; ecliptic coordinates to geocentric and then to RA and DEC

   if keyword_set(planet) then begin
       planetlist = ['MERCURY','VENUS','MARS', $
        'JUPITER','SATURN','URANUS','NEPTUNE','PLUTO']
        index = 1+ where(planetlist eq strupcase(strtrim(planet,2)), Nfound) 
        if index[0] GE 3 then index = index + 1
       if Nfound EQ 0 then message,'Unrecognized planet of ' + planet
   endif else index = [1,2,4,5,6,7,8,9]

   helio,jj,index,rad,lon,lat,/radian

; extract Earth's info

   helio,jj,3,rade,lone,late,/radian

;get rectangular coords of planets

   x = rad * cos(lat) * cos(lon) - rade * cos(late) * cos(lone)
   y = rad * cos(lat) * sin(lon) - rade * cos(late) * sin(lone)
   z = rad * sin(lat)            - rade * sin(late)

;get geocentric longitude lambda and geo latitude, beta

   lambda = atan(y,x) * !radeg
   beta   = z/sqrt(x*x + y*y) * !radeg

;convert to Ra and Dec

   euler, lambda, beta, ra, dec, 4

   return
   end
