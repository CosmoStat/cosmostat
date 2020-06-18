pro xyad, hdr, x, y, a, d, PRINT = print        ;X, Y to Ra, Dec
;+
; NAME:
;       XYAD
; PURPOSE:
;       Use a FITS header to convert pixel (X,Y) to celestial coordinates
; EXPLANATION: 
;       Use astrometry in a FITS image header to compute R.A. and Dec in decimal
;       degrees from X and Y.  
;
; CALLING SEQUENCE:
;       XYAD, HDR               ;Prompt for X and Y positions
;       XYAD, HDR, X, Y, A, D, [ /PRINT]
;
; INPUTS:
;       HDR - FITS Image header containing astrometry info
;
; OPTIONAL INPUTS:
;       X     - row position in pixels, scalar or vector
;       Y     - column position in pixels, scalar or vector
;
;       X and Y should be in IDL convention, (first pixel is (0,0)).
;
; OPTIONAL OUTPUT:
;       A - Right ascension in decimal DEGREES, same number of elements as
;               X and Y
;       D - Declination in decimal DEGREES
;
; OPTIONAL KEYWORD INPUT:
;       /PRINT - If this keyword is set and non-zero, then results are displayed
;               at the terminal.
;
; OPERATIONAL NOTES:
;       If less than 5 parameters are supplied, or if the /PRINT keyword is
;       set, then then the X and Y positions are displayed at the terminal.
;
;       If this procedure is to be used repeatedly with the same header,
;       then it would be faster to use XY2AD.
;
; PROCEDURES CALLED
;       ADSTRING(), EXTAST, GSSSXYAD, XY2AD
; REVISION HISTORY:
;       W. Landsman                 STX          Jan, 1988
;       Use astrometry structure  W. Landsman    Jan, 1994
;       Recognize GSSS header  W. Landsman       June, 1994
;       Changed ADSTRING output format   W. Landsman    September 1995
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use vector call to ADSTRING() W. Landsman February 2000
;-
 On_error,2

 npar = N_params()
 if ( npar EQ 0 ) then begin
        print,'Syntax -  xyad, hdr, [x, y, a, d, /PRINT]'
        print,'HDR - FITS header (string array) containing astrometry'
        print,'X,Y - Input X and Y positions (scalar or vector)'
        print,'A,D - Ouput RA and Dec in decimal degrees'
        return
 endif                                                         

  extast, hdr, astr, noparams              ;Extract astrometry structure

  if ( noparams LT 0 ) then $ 
        message,'ERROR - No astrometry info in supplied FITS header'

  if ( npar lt 3 ) then read,'XYAD: Enter X and Y positions: ',x,y

  
  case strmid(astr.ctype[0],5,3)  of 
        'GSS': gsssxyad, astr, x, y, a, d
         else: xy2ad, x, y, astr, a, d
  endcase

  if (npar lt 5) or keyword_set(PRINT) then begin
        npts = N_elements(X)
        fmt = '(1x ,2F8.2,2x,2F9.4,A)'
        print,'    X       Y         RA       DEC       RA         DEC'
        str = adstring(a,d,1)
        for i=0l, npts-1 do $
        print,FORMAT=fmt, float(x[i]), float(y[i]), a[i], d[i], str[i]
   endif
   
   return
   end
