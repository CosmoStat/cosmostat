FUNCTION sigma_clip, Data, sigma_clip=sigma_clip, mean=mean, Niter=Niter
;+ 
; NAME: 
;       sigma_clip
;
; PURPOSE: 
;       return the sigma obtained by k-sigma. Default sigma_clip value is 3. 
;       if mean is set, the mean (taking into account outsiders) is returned.
;
; CALLING SEQUENCE: 
;   output=sigma_clip(Data, sigma_clip=sigma_clip, mean=mean)
;
; INPUTS: 
;   Data -- IDL array: data
;
; OPTIONAL INPUT PARAMETERS: 
;   none
;
; KEYED INPUTS: 
;   sigma_clip -- float : sigma_clip value 
;
; KEYED OUTPUTS: 
;   mean -- float : mean value 
;
; OUTPUTS: 
;    output
;
; EXAMPLE: 
;    output_sigma = sigma_clip(Image, sigma_clip=2.5)
;
; MODIFICATION HISTORY: 
;    25-Jan-1995 JL Starck written with template_gen 
;-
 
 
 output=''
 
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'output=sigma_clip(Data, sigma_clip=sigma_clip)'
   GOTO, CLOSING
 ENDIF
 
vsize = size(Data)
if vsize[0] EQ  0 then BEGIN
                  PRINT, 'Error: Data must be an array ... '
                  GOTO, CLOSING
                  ENDIF
    
;------------------------------------------------------------
; function body
;------------------------------------------------------------

if keyword_set(sigma_clip) then k = sigma_clip else k = 3.
if not keyword_set(Niter) then Ni = 2. $
else Ni = Niter - 1

Sig = 0.
Buff = Data

m = TOTAL(Buff)/ N_ELEMENTS(Buff)
Sig = sigma(Buff)
index = where((abs(Buff-M) LT k*Sig), count)

for i = 1,Ni do begin
   if (count Gt 0) then begin
      m = TOTAL(Buff[index]) / N_ELEMENTS(Buff(index))
      Sig = sigma(Buff[index])
      index = where((abs(Buff-M) LT k*Sig), count)
   end else Sig = 0.
   endfor
output = Sig

 mean = m

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, output
 
 END
