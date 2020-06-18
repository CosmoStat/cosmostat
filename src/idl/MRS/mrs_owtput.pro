;+
; NAME:
;        mrs_owtput
;
; PURPOSE:
;		Put a scale, 2D [*,*] array of coefficients, in the wavelet transform obtained by the command mrs_owttrans.pro, i.e. (bi-) orthogonal wavelet transform on the sphere
;
; CALLING:
;
;		mrs_owtput, W, Scale, ScaleNumber, DirectNumber=DirectNumber, BorderSize=BorderSize  
;       
; INPUTS:
;		W -- IDL structure: wavelet coefficients of an image, see mrs_owttrans.pro for more details.
;		Scale -- IDL 2D array: wavelet scale coefficients we want use in the decomposition. 
;		ScaleNumber -- int: Scale number, The scale number must be 
;							between 0 and W.NbrScale-1
;
; KEYWORDS:
;		DirectNumber -- int: Select the direction of the band, 0 for horizontal band, 1 for vertical band and 2 for diagonal band. Default value is 0.
;		BorderSize -- int: Parameter used for ignoring borders of bands, default value 0.
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Olivier Fourt, 2009
;	October, 2009 File creation
;
;--------------------------------------------------------------------------------------------------------------------------------

pro mrs_owtput, W, Scale, ScaleNumber, DirectNumber=DirectNumber, BorderSize=BorderSize

if not keyword_set(BorderSize) then BorderSize=0
if not keyword_set(DirectNumber) then DirectNumber =0

  EndX = w.nx / 2^ScaleNumber
  EndY = w.ny / 2^ScaleNumber
  HalfX = EndX / 2
  HalfY = EndY / 2
  
	size_x1 = HalfX+BorderSize
    if size_x1 LT 0 then size_x1 = 0
    size_x2 = HalfX-1-BorderSize
    if size_x2 LT 0 then size_x2 = 0
    size_x3 = EndX-1-BorderSize
    if size_x3 LT 0 then size_x3 = 0
    size_y1 = HalfY+BorderSize
    if size_y1 LT 0 then size_y1 = 0
    size_y2 = HalfY-1-BorderSize
    if size_y2 LT 0 then size_y2 = 0
    size_y3 = EndY-1-BorderSize
    if size_y3 LT 0 then size_y3 = 0

if DirectNumber EQ 0 then w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(BorderSize,size_y2) : max(BorderSize,size_y2), * ] = Scale; w.coef[HalfX+BorderSize:EndX-1-BorderSize, BorderSize:HalfY-1-BorderSize, *]

if DirectNumber EQ 1 then w.coef[ min(BorderSize,size_x2) : max(BorderSize,size_x2), min(size_y1,size_y3) : max(size_y1,size_y3), * ] = Scale; w.coef[BorderSize:HalfX-1-BorderSize,HalfY+BorderSize:EndY-1-BorderSize, *]

if DirectNumber EQ 2 then w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(size_y1,size_y3) : max(size_y1,size_y3), * ] = Scale; w.coef[HalfX+BorderSize:EndX-1-BorderSize, HalfY+BorderSize:EndY-1-BorderSize, *]
  
end

;====================================================================================================
