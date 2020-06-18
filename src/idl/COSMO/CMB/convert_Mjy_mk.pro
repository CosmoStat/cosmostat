pro convert_Mjy_mk, lam, sig_Mjy, sig_mk, CMB_TH=CMB_TH, RJ=RJ

;------------------------------------------------
; Convert MJy/sr in mK, either themodynamic or
; Rayleigh-Jeans

; Input: Lam     = wavelength in micron
;        sig_MjY = Signal in MJy/sr 
;
; Written Guilaine Lagache, 2001
;------------------------------------------------

IF KEYWORD_SET(RJ) THEN BEGIN
 ; mK Rayleigh-Jeans
 sig_mk=sig_Mjy/(1.e-3*2*1.38e-23*1e20/(lam*1.e-6)^2.) 
ENDIF

IF KEYWORD_SET(CMB_TH) THEN BEGIN
  ; Thermodynamique
  PLANCK_VALUES=  PLANCK (2.726, lam, dBdT_VALUES,/MJY)
  sig_mk=sig_Mjy/(dbdt_values*1.e-3)
ENDIF


END
