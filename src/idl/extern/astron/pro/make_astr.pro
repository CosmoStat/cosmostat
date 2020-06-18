pro make_astr,astr, CD=cd, DELTA = cdelt, CRPIX = crpix, CRVAL = crval, $
                    CTYPE = ctype, LONGPOLE = longpole, PROJP1 = projp1
	    	    PROJP2 = projp2
;+
; NAME:
;       MAKE_ASTR
; PURPOSE:
;       Build an astrometry structure from input parameter values
; EXPLANATION:
;       This structure can be subsequently placed in a FITS header with 
;       PUTAST
;
; CALLING SEQUENCE:
;       MAKE_ASTR, astr, CD = , DELT =, CRPIX =, CRVAL =, CTYPE =,
;               LONGPOLE =, PROJP1 =, PROJP2 =    
;
; OUTPUT PARAMETER:
;       ASTR - Anonymous structure containing astrometry info.  See the 
;              documentation for EXTAST for descriptions of the individual
;              tags
;
; REQUIRED INPUT KEYWORDS
;       CRPIX - 2 element vector giving X and Y coordinates of reference pixel
;               (def = NAXIS/2)
;       CRVAL - 2 element double precision vector giving R.A. and DEC of 
;               reference pixel in DEGREES
; OPTIONAL INPUT KEYWORDS
;       CD -  2 x 2 array containing the astrometry parameters CD1_1 CD1_2
;              in DEGREES/PIXEL                                CD2_1 CD2_2
;       DELT - 2 element vector giving physical increment at reference pixel
;              CDELT default = [1.0D, 1.0D].
;       CTYPE - 2 element string vector giving projection types, default
;              ['RA---TAN','DEC--TAN']
;       LONGPOLE - scalar longitude of north pole, default = 180
;       PROJP1 - Scalar parameter needed in some projections, default = -1.
;       PROJP2 - Scalar parameter needed in some projections, default = -2.
;
; NOTES:
;       (1) An anonymous structure is created to avoid structure definition
;               conflicts.    This is needed because some projection systems
;               require additional dimensions (i.e. spherical cube
;               projections require a specification of the cube face).
;       (2) The name of the keyword for the CDELT parameter is DELT because
;               the IDL keyword CDELT would conflict with the CD keyword
; REVISION HISTORY:
;       Written by   W. Landsman              Mar. 1994
;       Converted to IDL V5.0                 Jun  1998
;-
 On_error,2

 if ( N_params() LT 1 ) then begin
	print,'Syntax - MAKE_ASTR, astr, CD = , DELT =, CRPIX =, CRVAL =, '
        print,'			   CTYPE =, LONGPOLE =, PROJP1 =, PROJP2 = ]'
	return
 endif


 if N_elements( cd ) EQ 0 then cd = [ [1.,0.], [0.,1.] ]

 if N_elements( crpix) EQ 0 then message, $
	'ERROR - CRPIX is a required keyword for a new astrometry structure'
 
 if N_elements( crval) EQ 0 then message, $
	'ERROR - CRVAL is a required keyword for a new astrometry structure'

 if N_elements( ctype)  EQ 0 then ctype = ['RA---TAN','DEC--TAN']

 if N_elements( cdelt) EQ 0 then cdelt = [1.0D, 1.0D]

 if N_elements(longpole) EQ 0 then longpole = 180.0D

 if N_elements(projp1) EQ 0 then projp1 = -1.

 if N_elements(projp2) EQ 0 then projp2 = -2.

 ASTR = {CD: double(cd), CDELT: double(cdelt), $
		CRPIX: float(crpix), CRVAL:double(crval), $
		CTYPE: string(ctype), LONGPOLE: float( longpole[0]),  $
		PROJP1: float(projp1[0]), PROJP2: float(projp2[0])}

  return
  end
