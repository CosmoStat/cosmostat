pro i_whitenoise, MAP, NSIDE, SIGMAPIX,$
                    HITMAP = hitmap,seed=seed

;+
; NAME:
;	I_WHITENOISE
;
; PURPOSE:
;       Generates a white noise I map.
;
; CALLING SEQUENCE:
; 	I_WHITENOISE, MAP, NSIDE, SIGMAPIX
;                        
;
; INPUTS:
;      NSIDE       INTEGER    Healpix NSIDE parameter.
;      SIGMAPIX    REAL       Pixel RMS noise.      
;
; OPTIONAL INPUTS:
;      HITMAP      ARRAY      Healpix hit map.
;      SEED        INTEGER    Seed for random number generator.
;      
; OPTIONAL INPUT KEYWORDS:
;
; NOTES
;
; SIDE EFFECTS
;
; EXAMPLES
; i_whitenoise, whitenoisemap,2048,1
;
; COMMONS USED : 
;
; PROCEDURES USED: 
;
; MODIFICATION HISTORY:
; 	May 2006, Samuel Leach, SISSA
;-
bad_data= -1.63750000e+30


npix=nside2npix(nside)
MAP = make_array(npix,1,value=bad_data,/float)

if(defined(hitmap)) then begin
;Generate non-stationary white noise map
   print,'Generating non-stationary white noise map'
   obs    = WHERE( hitmap[*] GT 0 )
   nobs=n_elements(obs)
   MAP[obs]=RANDOMN(seed, nobs)
   MAP[obs]=MAP[obs]*1/sqrt(hitmap[obs])*sigmapix
endif else begin
   ;Generate pure white noise map
   print,'Generating stationary white noise map'
   MAP[*]=RANDOMN(seed, npix)
   MAP[*]=MAP[*]*sigmapix
endelse

END
