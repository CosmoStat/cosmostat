;+
; NAME:
;        mrs_owtget
;
; PURPOSE:
;		Return in a 2D [*,*] array a scale coefficients of the wavelet transform obtained by the command mrs_owttrans.pro, i.e. (bi-) orthogonal wavelet transform on the sphere
;
; CALLING:
;		Scale = mrs_owtget( W, ScaleNumber, DirectNumber=DirectNumber, BorderSize=BorderSize, Isotrop=Isotrop )
;
; INPUTS:
;		W -- IDL structure: wavelet coefficients of an image, see mrs_owttrans.pro for more details.
;		ScaleNumber -- int: Scale number, The scale number must be 
;							between 0 and W.NbrScale-1
;    
;
; KEYWORDS:
;		DirectNumber -- int: Select the direction of the band, 0 for horizontal band, 1 for vertical band and 2 for diagonal band. Default value is 0.
;		BorderSize -- int: Parameter used for ignoring borders of bands, default value 0.
;		Isotropic -- scalar: If set, ignore the keyword DirectNumber and gets the three bands in a 3D [*,*,3] array
;
; EXAMPLE:
;       Compute the orthogonal wavelet transform of an image with five scales and then extract the second scale
;               mrs_owttrans, Imag, WT, NbrScale=5
;               scale_2 = mrs_owtget( WT, 1 )
;         
; HISTORY:
;	Written:  Olivier Fourt, 2009
;	October, 2009 File creation
;-
;-----------------------------------------------------------------

;====================================================================================================

function mrs_owtget, W, ScaleNumber, DirectNumber=DirectNumber, BorderSize=BorderSize, Isotropic=Isotropic

if not keyword_set(BorderSize) then BorderSize=0
if not keyword_set(DirectNumber) then DirectNumber =0

  EndX = w.nx
  EndY = w.ny
  
  EndX = EndX / 2^ScaleNumber
  EndY = EndY / 2^ScaleNumber
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
     
  if keyword_set(isotropic) then begin
		BandHorizontal = w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(BorderSize,size_y2) : max(BorderSize,size_y2), * ]; w.coef[HalfX+BorderSize:EndX-1-BorderSize, BorderSize:HalfY-1-BorderSize, *]
		
		BandVertical = w.coef[ min(BorderSize,size_x2) : max(BorderSize,size_x2), min(size_y1,size_y3) : max(size_y1,size_y3), * ]; w.coef[BorderSize:HalfX-1-BorderSize,HalfY+BorderSize:EndY-1-BorderSize, *]
		
		BandDiagonal = w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(size_y1,size_y3) : max(size_y1,size_y3), * ]; w.coef[HalfX+BorderSize:EndX-1-BorderSize, HalfY+BorderSize:EndY-1-BorderSize, *]
		
		Scale = [ BandHorizontal, BandVertical, BandDiagonal ]
  end else begin
      if DirectNumber EQ 0 then Scale = w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(BorderSize,size_y2) : max(BorderSize,size_y2), * ]
      									;w.coef[HalfX+BorderSize:EndX-1-BorderSize, BorderSize:HalfY-1-BorderSize, *]
      									
      if DirectNumber EQ 1 then Scale = w.coef[ min(BorderSize,size_x2) : max(BorderSize,size_x2), min(size_y1,size_y3) : max(size_y1,size_y3), * ]
      									;w.coef[BorderSize:HalfX-1-BorderSize,HalfY+BorderSize:EndY-1-BorderSize, *]
      									
      if DirectNumber EQ 2 then Scale = w.coef[ min(size_x1,size_x3) : max(size_x1,size_x3), min(size_y1,size_y3) : max(size_y1,size_y3), * ]
      									;w.coef[HalfX+BorderSize:EndX-1-BorderSize, HalfY+BorderSize:EndY-1-BorderSize, *]
   end

return, Scale
end

