;+
; NAME:
;        mrs_GetWTLocalVariance
;
; PURPOSE:
;       Compute the local variance for each position and for each wavelet scale. of a given map.
;       if the keyword Cl is set (Cl is assumed to a power spectrum), then the variance relative to CL
;       at a given wavelet scale is added to output.
;
; CALLING:
;     WTLV = mrs_GetWTLocalVariance(InputMap, NbrScale=NbrScale, BS=BS, Cl=Cl, WTTrans=WTTrans)
;
; INPUTS:
;     InputMap -- Healpix Map = Input Map    
;    
; OUTPUTS:
;     WTLM -- 2D IDL array = Variance per spatial position and per scale. WTLV
;
; INPUT KEYWORD: 
;     NbrScale: int = Number of scales to use in the Wavelet Transform. Default is 5.
;     BS: int = block size at the finest scale to compute the local variance. Default is 16.
;     Cl: IDL 1D array = Power Spectrum to be added to the calculated variance. 
;
; OUTPUT KEYWORD: 
;    WTTrans: wavelet transform of the input map
;
; EXTERNAL CALL:
;  mrs_almtrans         in mrs_almtrans.pro
;  mrs_bandpass_sigma   in mrs_init.pro

; EXAMPLE:
;       Compute the local RMS at each scale of a given pyramidal wavelet transform: 
;               mrs_pwttrans, Imag, Trans
;               mrs_wt_local_rms, Trans, RMSTrans
;         
; HISTORY:
;       Written: Jerome Bobin & Jean-Luc Starck, 2012
;--------------------------------------------------------------------------------------------------------

 
function mrs_GetWTLocalVariance, Input, NbrScale=NbrScale, BS=BS, Cl=Cl, WTTrans=WTTrans

if not keyword_set(NbrScale) then NbrScale = 5
; if keyword_set(Cl) then print,'Add CMB variance'
if not keyword_set(BS) then BS = 16

nside = gnside(input)

mrs_wttrans,Input, WTTrans,NbrScale=NbrScale
cn = WTTrans.coef
T = WTTrans.TabPsi
TPhi = WTTrans.TABPHI
out = 0

varmap = 0.*cn

Tall = 0.*T

BSl = BS

for s=0,NbrScale-1 do begin

	tempn = cn[*,s]	
	
	CMBVar = 0
	if keyword_set(Cl) then begin
	  if s NE NbrScale-1 then CMBVar = mrs_bandpass_sigma(Cl, T[*,s])^2. $
	  else CMBVar = mrs_bandpass_sigma(Cl, TPhi[*,NbrScale-2])^2.
	end
	
	for facenum=0,11 do begin ;--scan all faces
	
		fn = get_one_face(tempn,facenum)
		gn = fn
		
		nb_bin = double(nside)/double(BSl)
		for bx = 0,nb_bin-1. do begin ;--- scan bx
			for by = 0,nb_bin-1. do begin
			
				pn = fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]	
				gn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = sigma(pn)^2. + CMBVar ;--- Attention aux normalisations !!
				
			endfor
		endfor ;--- patch loops
	
		put_one_face,tempn,gn,facenum

	endfor ;--- face loop

	varmap[*,s] = tempn
		
	if (nside GE BSl*2) then  BSl = 2.*BSl

endfor

return, varmap

end
