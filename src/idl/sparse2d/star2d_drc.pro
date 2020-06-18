;+
; NAME:
;        STAR2D_DRC
;
; PURPOSE:
;       Dynamic Range Compression using the 2D starlet transform of an image . 
;       The output is a 2D image.  
;
; CALLING:
;      DRC = STAR2D_DRC(Imag, Nscale=Nscale, Eps=Eps)
;       
; INPUTS:
;     Imag -- 2D IDL array: image we want to transform
;
; OUTPUTS:
;     DRC -- 2D IDL array: Image
;  
; KEYWORDS:
;      Nscale -- int: Number of scales. Default is 4.
; 
; EXAMPLE:
;       DRC_Image = STAR2D_DRC(I)
;
; REFERENCES:
;    J.-L. Starck, M. Elad, and D.L. Donoho,
;    "Redundant Multiscale Transforms and their Application for Morphological
;    Component Analysis", Advances in Imaging and Electron Physics, 132, 2004.
;
; EXTERNAL CALL:  
;     star2d.pro
;
; AUTHOR:
;    J.L. Starck
;    Service d'Astrophysique, Centre d'Etudes de SACLAY,
;    Orme des Merisiers, 91191 GIF-Sur-YVETTE CEDEX, France 
;    Email: jstarck@cea.fr        Tel:  +33 (0) 169 085 764
;    http://jstarck.free.fr       Fax:  +33 (0) 169 086 577
;-
;-------------------------------------------------------------------------------

function mad, data
  n = n_elements(data)
  d = fltarr(n)
  d[*] = abs(data-median(data))
  ind = sort(d)
  d = d[ind]
  return, d[n/2] / 0.6747
end


function star2d_drc, Ima, Nscale=Nscale, Eps=Eps, NoSign=NoSign  

        DRC = -1

	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: WT = star2d_drc(Ima, Nscale=Nscale, Eps=Eps)'
        goto, DONE
	end

	if not keyword_set(Nscale) then Nscale =4
	if not keyword_set(Eps) then Eps = 0.1

	if Nscale LT 2 then  Nscale =4
	
	vsize = size(Ima)
	if vsize[0] NE 2 then begin
        print, 'Error: First parameter is not an image ...'
        print, 'CALLING SEQUENCE:  WT = star2d_drc(Ima, Nscale=Nscale, Eps=Eps)'
        goto, DONE
	end
	
	;-------------------------------------------------------------------------------
	TabNormPaveB3Spline = [0.889434,  0.200105, 0.0857724, 0.0413447, $
	                       0.0202689, 0.00995628, 0.00513504, 0., 0., 0., 0, $
			       0,0,0,0,0,0,0,0,0,0]
			    
	TabSigmaNoise = fltarr(Nscale)	       
	w = star2d(Ima, Nscale=Nscale)
	Noise = mad( w[*,*,0] )
	for j=0,Nscale-1 do begin
	   TabSigmaNoise[j] = Noise * TabNormPaveB3Spline[j]
	   if TabSigmaNoise[j] EQ 0 then TabSigmaNoise[j] = TabSigmaNoise[j-1]
	   W[*,*,j] = W[*,*,j] / TabSigmaNoise[j] 
	end
	
 	SignW = w
	SignW[*] = 1
	ind = where(w LT 0, c)
	if c GT 0 then SignW[ind] = -1
	info,  SignW
	if not keyword_set(NoSign) then WL = SignW * alog(ABS(W) + Eps) $
	else WL =  alog(ABS(W) + Eps)
	DRC = total(WL, 3)
	
DONE:
    return, DRC
    
end	
	
	
