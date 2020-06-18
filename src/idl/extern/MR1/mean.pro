FUNCTION MEAN, data
;+ 
; NAME: 
;       MEAN
;
; PURPOSE: 
;       mean routine of an IDL n-array
;                output = total(data) / N_ELEMENTS(data)
;
; CALLING SEQUENCE: 
;   output=MEAN(data)
;
; INPUTS: 
;   data
;
; OUTPUTS: 
;    output
;
; MODIFICATION HISTORY: 
;    30-Mar-1995 JL Starck written with template_gen 
;-
 
;------------------------------------------------------------
; initialization
;------------------------------------------------------------
 
 output=''
 
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=MEAN(data)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
 output = total(data) / N_ELEMENTS(data)
 
;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
 
  RETURN, output
 
 END
