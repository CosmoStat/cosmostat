;+
; NAME:
;        STAR2D_ADJOINT
;
; PURPOSE:
;       Computes the adjoint of the starlet transform of an image (i.e. undecimated isotropic wavelet transform). 
;       The output is a 1D IDL array. If the keyword Gen2 is set, then it is
;       the 2nd generation starlet transform which is computed:  i.e. g = Id - h*h 
;       instead of g = Id - h.
;       
; CALLING:
;      Imag = STAR2D_ADJOINT(DataTrans, Gen2=Gen2)
;       
; INPUTS:
;      DataTransf -- 3D IDL array: Wavelet Transform
;                                 DataTransf(*,*,i) = ith band of the 
;                                                    wavelet transform
;
; OUTPUTS:   
;     Imag -- 2D IDL array: image we want to transform
;
; 
; 
; EXAMPLE:
;       Compute the starlet transform of an image I with default options
;       (i.e. a trou algorithm with 4 scales).  
;               Output = STAR2D(I)
;       Reconstruction using the adjoint operator:
;               RecI = STAR2D_ADJOINT(Output)
;
; REFERENCES:
;    [1] J.L. Starck and F. Murtagh, 
;    "Image Restoration with Noise Suppression 
;    Using the Wavelet Transform",
;    Astronomy and Astrophysics, 288, pp-343-348, 1994.
;    
;
;    For the modified STARLET transform:
;    [2] J.-L. Starck, J. Fadili and F. Murtagh, "The Undecimated Wavelet Decomposition 
;        and its Reconstruction", IEEE Transaction on Image Processing,  16,  2, pp 297--309, 2007.
; AUTHOR:
;    J.L. Starck
;    Service d'Astrophysique, Centre d'Etudes de SACLAY,
;    Orme des Merisiers, 91191 GIF-Sur-YVETTE CEDEX, France 
;    Email: jstarck@cea.fr        Tel:  +33 (0) 169 085 764
;    http://jstarck.free.fr       Fax:  +33 (0) 169 086 577
;-
;-------------------------------------------------------------------------------

pro b3spline_star2d, Im_in, Im_Out,  Step
   
   Nl = (size(Im_in))(1)
   Nc = (size(Im_in))(2)
   Im_Out = dblarr(Nl,Nc)
   
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Buffs = dblarr(Nl,Nc)
   Im_Out= dblarr(Nl,Nc)
   for i = 0, Nl-1 do begin
    for j = 0, Nc-1 do begin
      jm = test_ind(j-Step, Nc)
      jp = test_ind(j+Step, Nc)
      jm2 = test_ind(j-2*Step, Nc)
      jp2 = test_ind(j+2*Step, Nc)
      Buffs (i,j) =$ 
		C3 * Im_in(i,j) + C2 * (Im_in(i,jm) + Im_in(i,jp))$
		+ C1 * (Im_in(i,jm2) + Im_in(i,jp2)) 
	endfor
   endfor

   for i = 0, Nl-1 do begin
    for j = 0, Nc-1 do begin
      im = test_ind(i-Step, Nl)
      ip = test_ind(i+Step, Nl)
      im2 = test_ind(i-2*Step, Nl)
      ip2 = test_ind(i+2*Step, Nl)
      Im_Out (i,j) =$
			C3 * Buffs(i,j) + C2 * (Buffs(im,j) + Buffs(ip,j))$
			+ C1 * (Buffs(im2,j) + Buffs(ip2,j))
    endfor
   endfor
end


pro fast_b3spline_star2d, Im_in, Im_Out,  Step
; /edge_wrap
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   KSize = 4*Step+1
   KS2 = KSize/2
   Kernel = fltarr(KSize)
   Kernel[0] = C1
   Kernel[KSize-1] = C1
   Kernel[KS2+Step] = C2
   Kernel[KS2-Step] = C2
   Kernel[KS2]=C3
   ; /EDGE_MIRROR does not work with old IDL version
   ; z = convol(Im_in,Kernel, /EDGE_MIRROR, /CENTER)  
   z = convol(Im_in,Kernel, /EDGE_TRUNCATE, /CENTER) 
   kernelY = transpose(Kernel)
   ; Im_Out = convol(z, kernelY, /EDGE_MIRROR, /CENTER)
   Im_Out = convol(z, kernelY, /EDGE_TRUNCATE, /CENTER)
end


function star2d_adjoint,  WT, Gen2=Gen2
RecIma = -1
Fast=1
	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: Rec = star2d_adjoint(WT, Gen2=Gen2)'
        goto, DONE
	end

	vsize = size(WT)
	if vsize[0] NE 3 then begin
        print, 'Error: First parameter is not an cube ...'
        print, 'CALLING SEQUENCE:  WT = star2d_adjoint(WT, Gen2=Gen2)'
        goto, DONE
	end
	;-------------------------------------------------------------------------------
		
vs = size(WT)
Nscale = vs[3]
NStep=Nscale-1
Step_trou = 1
for i=0, NStep-2 do Step_trou = Step_trou * 2
	
RecIma = WT[*,*,Nscale-1]
Im_Out = fltarr(vs[1],vs[2])
if keyword_set(Gen2) then  Im_Out2 = fltarr(vs[1],vs[2])

for j=Nscale-2,0,-1 do begin
   ; print, j, Step_trou
   ; fast_b3spline_star2d,recIma,dum,Step_trou
   if keyword_set(fast) then fast_b3spline_star2d, recIma,dum,Step_trou  $
   else b3spline_star2d, recIma,dum,Step_trou 
   recIma = dum

  ; fast_b3spline_star2d, WT[*,*,j], Im_Out, Step_trou
   if keyword_set(fast) then fast_b3spline_star2d, WT[*,*,j], Im_Out, Step_trou  $
   else b3spline_star2d, WT[*,*,j], Im_Out, Step_trou 

  if keyword_set(Gen2) then  begin    
      ; fast_b3spline_star2d, Im_Out, Im_Out2, Step_trou
      if keyword_set(fast) then  fast_b3spline_star2d, Im_Out, Im_Out2, Step_trou  $
      else  b3spline_star2d, Im_Out, Im_Out2, Step_trou
      RecIma = RecIma + WT[*,*,j] - Im_Out2
  end else RecIma = RecIma + WT[*,*,j] - Im_Out
  Step_trou = Step_trou / 2
end

DONE:
    return,  RecIma
end

;-------------------------------------------------------------------------------
