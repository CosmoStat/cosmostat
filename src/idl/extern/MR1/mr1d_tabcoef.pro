FUNCTION MR1D_TABCOEF, Type
;+ 
; NAME: 
;     M1RD_TABCOEF
;
; PURPOSE: 
;     Return the pre-computed table for noise behaviour in the
;      multiresolution transform.
;
; CALLING SEQUENCE: 
;   output=MR1D_TABCOEF(TypeTransform)
;
; INPUTS: 
;   Type -- string: type of multiresolution
;                      Type = "atrou","pyrmed" or "pavemed"
;
;
; OUTPUTS: 
;    output: 1D table: pre-computed table
;                     The ouput is a float array which
;                     gives the standard deviation of the noise at each scale,
;                     when we apply a multiresolution transform to a signal
;                     following a Gaussian distribution  with a standard
;                     deviation equal to 1.
;
; MODIFICATION HISTORY: 
;    8-Jan-1996 JL Starck 
;-
 
  
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=MR1D_TABCOEF(Type)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
TabLevel = [0.721618,0.286209,0.178281,0.121010,0.0831040,0.0588543, 0.0446171, 0.0289142, 0.0177307]

;TabLevel = [0.726336,0.283629,0.177784,0.1261,0.1,0.1, 0.1, 0.1, 0.1]

TabLevelMed = [0.915152,0.357905,0.248778,0.209270,0.162813,0.111927, 0.0786654, 0.0594114, 0.0643011]
TabPaveMed = [0.938905, 0.333192, 0.244986, 0.207760, 0.136711]

return_Tab=-1
case strupcase(Type) of
   "ATROU": return_Tab = TabLevel
   "PYRMED": return_Tab = TabLevelMed
   "PAVEMED":return_Tab = TabPaveMed
   else: print, "Unknown Type in mr1d_tabcoef: ", Type
endcase

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, return_Tab
 
 END
