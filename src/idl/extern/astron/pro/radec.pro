pro radec,ra,dec,ihr,imin,xsec,ideg,imn,xsc ;Convert from decimal to hrs,min
;+
; NAME:
;	RADEC
; PURPOSE:
;	To convert RA and Dec  from decimal to sexigesimal units.
; EXPLANATION: 
;	The conversion is to sexigesimal hours for RA,  and sexigesimal 
;	degrees for declination.
;
; CALLING SEQUENCE:
;	radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc
;
; INPUTS:
;	ra   - right ascension in decimal DEGREES, scalar or vector
;	dec  - declination in decimal DEGREES, scalar or vector, same number
;		of elements as RA
;
; OUTPUTS:
;	ihr  - right ascension hours   (INTEGER*2)
;	imin - right ascension minutes (INTEGER*2)
;	xsec - right ascension seconds  (REAL*4 or REAL*8)
;	ideg - declination degrees (INTEGER*2)
;	imn  - declination minutes (INTEGER*2)
;	xsc  - declination seconds (REAL*4 or REAL*8)
;
; RESTRICTIONS:
;	RADEC does minimal parameter checking.
;
; REVISON HISTORY:
;	Written by B. Pfarr, STX, 4/24/87
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
  On_error,2

  if (N_params() LT 2 ) then begin
    print,'Syntax - radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc'
    return
  endif
  
;    Compute RA

  ra = ra mod 360.          ;Make sure between 0 and 24 hours
  ra = ra + 360*(ra lt 0)
  ihr = fix(ra/15.)
  xmin =abs(ra*4.0-ihr*60.0)
  imin = fix(xmin)
  xsec = (xmin-imin)*60.0

;    Compute Dec

  ideg = fix(dec)
  xmn = abs(dec-ideg)*60.0
  imn = fix(xmn)
  xsc = (xmn-imn)*60.0

; Now test for the special case of zero degrees

  zero_deg = ( ideg EQ 0 ) and (dec LT 0)
  imn = imn - 2*imn*fix( zero_deg*(imn NE 0) )
  xsc = xsc - 2*xsc*zero_deg*(imn EQ 0)

  return
  end
