;+
; NAME:
;        ISTAR2D
;
; PURPOSE:
;       Reconstruct an image from its starlet transform. 
;       The output is a 2D IDL array. If the keyword Gen2 is set, then it is
;       the 2nd generation starlet transform which was is used computed:  i.e. g = Id - h*h 
;       instead of g = Id - h. 
;       If  the keyword Gen2 is not set, the reconstruction is done by simple co-addition 
;       of all frames: Rec = total(output, 3)
;       
; CALLING:
;      Rec = ISTAR2D(Trans, Gen2=Gen2)
;       
; INPUTS:
;     Trans -- 3D IDL array: Starlet transform of an image
;
; OUTPUTS:
;     Rec -- 2D IDL array: Reconstructed image
;  
; KEYWORDS:
;     Gen2 -- scalar: if set, it is the seconde generation of the starlet transform which is used.
; 
; EXAMPLE:
;       Compute the multiresolution of the image I with default options
;       (i.e. a trou algorithm with 4 scales).  
;               Output = STAR2D(I)
;       Reconstruction can be done by simple co-addition of all frames:
;               Rec = total(output, 3)
;
; REFERENCES:
;    [1] J.L. Starck and F. Murtagh, 
;    "Image Restoration with Noise Suppression 
;    Using the Wavelet Transform",
;    Astronomy and Astrophysics, 288, pp-343-348, 1994.
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
;-------------------------------------------------------------------------------

pro b3spline_istar2d, Im_in, Im_Out,  Step
; /edge_wra

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
   z = convol(Im_in,Kernel, /EDGE_TRUNCATE) 
   kernelY = transpose(Kernel)
   Im_Out = convol(z, kernelY, /EDGE_TRUNCATE)
end


function istar2d,  WT, Gen2=Gen2
RecIma = -1

	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: Rec = istar2d(WT, Gen2=Gen2)'
        goto, DONE
	end

	vsize = size(WT)
	if vsize[0] NE 3 then begin
        print, 'Error: First parameter is not an cube ...'
        print, 'CALLING SEQUENCE:  WT = istar2d(WT, Gen2=Gen2)'
        goto, DONE
	end
	;-------------------------------------------------------------------------------
		
vs = size(WT)
Nscale = vs[3]
NStep=Nscale-1
Step_trou = 1
for i=0, NStep-2 do Step_trou = Step_trou * 2
	
RecIma = WT[*,*,Nscale-1]
if keyword_set(Gen2) then Im_Out = fltarr(vs[1],vs[2])

for j=Nscale-2,0,-1 do begin
  ; print, j, Step_trou
  if keyword_set(Gen2) then  begin
      b3spline_istar2d, RecIma, Im_Out, Step_trou
      RecIma = Im_Out + WT[*,*,j]
  end else RecIma = RecIma + WT[*,*,j]
  Step_trou = Step_trou / 2
end

DONE:
    return,  RecIma
end

;-------------------------------------------------------------------------------


