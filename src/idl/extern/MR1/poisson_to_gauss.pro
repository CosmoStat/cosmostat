FUNCTION POISSON_TO_GAUSS, Data, poisson=poisson, inv=inv
;+ 
; NAME: 
;     POISSON_TO_GAUSS
;
; PURPOSE: 
;     Transform data with poisson noise in data with gaussian noise.
;
;        T = 2/gain * sqrt(gain*Data + 3./8.*gain^2 + sigma^2 - gain*mean)
;         where
;           gain = gain of the CCD
;           sigma = standard deviation of the read out noise
;           mean = mean of the read out noise
;
;     For a pure poisson noise, gain=1, mean=0, and sigma=0
;
;     if inv is set, then the inverse transfom is applied:
;        T = Data^2 /4. * gain - (3./8.*gain^2 + sigma^2 - gain*mean)/gain
;
;     by default, poisson=[1.,0.,0.].
;
;     more details about this transform ca be found in 
;     J.L. Starck, A. Bijaoui, and F. Murtagh, 
;     "Multiresolution Support Applied to Image Filtering and Deconvolution",
;      in CVIP: Graphical Models and Image Processing, Vol. 57, 5, 
;      pp 420-431, Sept. 1995.
;
; CALLING SEQUENCE: 
;   output=POISSON_TO_GAUSS(Data, poisson=poisson, inv=inv)
;
; INPUTS: 
;   Data -- IDL array: data
;
; OPTIONAL INPUT PARAMETERS: 
;   none
;
; KEYED INPUTS: 
;   poisson -- fltarr(3):  poisson(0) = gain
;                          poisson(1) = sigma read out noise
;                          poisson(2) = mean read out noise
;            or fltarr(2): poisson(0) = gain 
;                          and poisson(1) = sigma read out noise
;                          and mean is set to 0
;            or fltarr(1) : poisson(0) = gain
;                          and mean and sigma are set to 0
;                          
;   inv -- scalar: if set, the inverse transfomr is applied
;
; OUTPUTS: 
;    output -- IDL array: transformed data
;
; MODIFICATION HISTORY: 
;    8-Jan-1996 JL Starck written with template_gen 
;-

;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=POISSON_TO_GAUSS(Data, poisson=poisson, inv=inv)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
if keyword_set(poisson) then gain = poisson(0) else gain = 1
if keyword_set(poisson) and (size(poisson))(1) GE 2 then $
                                sigmap = poisson(1) else sigmap = 0
if keyword_set(poisson) and (size(poisson))(1) GE 3 then $
                                meanp = poisson(2) else meanp = 0

b = 3./8.*gain^2 +sigmap^2 - gain*meanp

if not keyword_set(inv) then trans = 2./gain*sqrt(Data*gain + b) $
else trans = Data^2 / 4. * gain - b / gain

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, trans
 
 END
