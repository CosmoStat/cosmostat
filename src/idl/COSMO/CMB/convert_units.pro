FUNCTION CONVERT_UNITS, data_in, units_out, units_in=units_in, nu_in=nu_in,$
                header=header, new_header=new_header

; NAME:
;      CONVERT_UNITS
; PURPOSE:
;      Convert values in selected units
;
; CALLING SEQUENCE:
;      result = CONVERT_UNITS( data_in, units_out,
;               [units_in=units_in, nu_in=nu_in, header=header, new_header=new_header])   
;
; INPUTS:
;      header =  FITS header array, (e.g. as returned by READFITS) 
;             string array, each element should have a length of 80 characters
;      data_in = initial value(s) (can be an array)
;
; INPUT KEYWORDS
;      units_in and nu_in, or header must be specified
;
; OPTIONAL INPUT KEYWORDS:
;      units_in = input unit ('MJY/SR' or 'DT/T_CMB' or 'T_RJ' or 'Y_SZ')
;      nu_in = frequency (in Hz)
;      header =  FITS header array, (e.g. as returned by READFITS) 
;             string array, each element should have a length of 80 characters
;
; OUTPUTS:
;      result = value expressed in new units;
; OPTIONAL OUTPUT KEYWORD:
;      new_header = FITS header array in new units.
;      Can be used if header is defined
;
; EXAMPLES:
;
; RESTRICTIONS:
;      Not completely tested yet.
;
; PROCEDURES CALLED:
;      SXPAR(),
;      SXADDPAR,
;      PLANCK_DBNU_DT()
; 
; MODIFICATION HISTORY:
;      Written, Guillaume Patanchon July 2000
;      Now, units_in can be equal to units_out, G.P. October 2000
;      Header updated, interactive help added,
;      explicit call to !const removed, JD, april 2001
;      unit 'Y_SZ' added, JBM, may 2001
;      now new_header is an output keyword, not an optionnal output, JD, july 2001
;      jd_dbnu_dt replaced by planck_dbnu_dt, JBM, Feb 2006 
;

IF ( N_PARAMS() EQ 0 ) THEN BEGIN
print, 'CALL IS : result = CONVERT_UNITS(data_in, units_out, [new_header])'
print, "Keywords : units_in = ('MJY/SR' or 'DT/T_CMB' or 'T_RJ' or 'Y_SZ'),"
print, "nu_in=nu_in, header=header"
GOTO, closing
ENDIF

const= {c:299792458d0, h:6.6260693d-34,$
hbar:1.05457168d-34,k:1.3806505d-23,g:6.6742d-11,msol:1.98844d30,tcmb:2.725}
defsysv, '!const', const



IF DEFINED(!const) THEN BEGIN
	Tcmb = !const.tcmb
	c_light = !const.c
ENDIF ELSE BEGIN
	Tcmb = 2.726			; COBE FIRAS by default
	c_light = 2.9979246e+08	        ; speed of light
ENDELSE

IF KEYWORD_SET(header) THEN BEGIN
	units_in = SXPAR(header,'BUNIT')
	length_in = SXPAR(header,'CRVAL3')
	nu_in = c_light/length_in
	new_header = header
	SXADDPAR, new_header, 'BUNIT', units_out
ENDIF
uni_in=STRCOMPRESS(units_in,/REMOVE_ALL)
uni_out=STRCOMPRESS(units_out,/REMOVE_ALL)
nu_in=DOUBLE(nu_in)

x = !const.h*nu_in/!const.k/Tcmb

CASE uni_in OF

uni_out : BEGIN
	    data_out=data_in
          END
'MJY/SR': BEGIN
            IF (uni_out EQ 'DT/T_CMB') THEN BEGIN
              data_out=1d-20/PLANCK_DBNU_DT(nu_in)*data_in/Tcmb
            ENDIF
            IF (uni_out EQ 'T_RJ') THEN BEGIN
              data_out=1d-20*(!const.c)^2/2/!const.k/(nu_in)^2*data_in
            ENDIF
            IF (uni_out EQ 'Y_SZ') THEN BEGIN
              data_out=1d-20/PLANCK_DBNU_DT(nu_in)*data_in/Tcmb/ $
                       (x*(EXP(x)+1.)/(EXP(x)-1.)-4.)
            ENDIF
            END
'DT/T_CMB': BEGIN
             IF (uni_out EQ 'MJY/SR') THEN BEGIN
               data_out=1d20*PLANCK_DBNU_DT(nu_in)*data_in*Tcmb
             ENDIF
             IF (uni_out EQ 'T_RJ') THEN BEGIN
               data_out=PLANCK_DBNU_DT(nu_in)*(!const.c)^2/2/!const.k/(nu_in)^2* $
                        data_in*Tcmb
             ENDIF
             IF (uni_out EQ 'Y_SZ') THEN BEGIN
               data_out=1./(x*(EXP(x)+1.)/(EXP(x)-1.)-4.)*data_in
             ENDIF
            END
'T_RJ': BEGIN
         IF (uni_out EQ 'MJY/SR') THEN BEGIN
           data_out=1d20*2*!const.k*(nu_in)^2/(!const.c)^2*data_in
         ENDIF
         IF (uni_out EQ 'DT/T_CMB') THEN BEGIN
           data_out=(2d*(!const.k)/(!const.c)^2)*(nu_in)^2/PLANCK_DBNU_DT(nu_in)* $
                    data_in/Tcmb
         ENDIF
         IF (uni_out EQ 'Y_SZ') THEN BEGIN
            data_out=1./x^2/EXP(x)*(EXP(x)-1)^2/ $
                     (x*(EXP(x)+1.)/(EXP(x)-1.)-4.)*data_in/Tcmb
         ENDIF
        END
'Y_SZ': BEGIN
         IF (uni_out EQ 'MJY/SR') THEN BEGIN
            data_out=1d20*PLANCK_DBNU_DT(nu_in)*(x*(EXP(x)+1.)/(EXP(x)-1.)-4.)* $
                     data_in*Tcmb
         ENDIF
         IF (uni_out EQ 'DT/T_CMB') THEN BEGIN
            data_out=(x*(EXP(x)+1.)/(EXP(x)-1.)-4.)*data_in
         ENDIF
         IF (uni_out EQ 'T_RJ') THEN BEGIN
            data_out=x^2*EXP(x)/(EXP(x)-1)^2* $
                     (x*(EXP(x)+1.)/(EXP(x)-1.)-4.)*data_in*Tcmb
         ENDIF
        END
ENDCASE

RETURN, data_out

closing:

END
