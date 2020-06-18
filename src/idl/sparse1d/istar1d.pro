;+
;
; NAME:
;        ISTAR1D
;
; PURPOSE:
;       Reconstruct a signal from its starlet transform. 
;       The output is a 1D IDL array. If the keyword Gen2 is set, then it is
;       the 2nd generation starlet transform which was is used computed:  i.e. g = Id - h*h 
;       instead of g = Id - h. 
;       If  the keyword Gen2 is not set, the reconstruction is done by simple co-addition 
;       of all scales: Rec = total(output, 2)
;       
; CALLING:
;      Rec = ISTAR1D(Trans, Gen2=Gen2)
;       
; INPUTS:
;     Trans -- 2D IDL array: Starlet transform of an image
;
; OUTPUTS:
;     Rec -- 1D IDL array: Reconstructed image
;  
; KEYWORDS:
;     Gen2 -- scalar: if set, it is the seconde generation of the starlet transform which is used.
; 
; EXAMPLE:
;       Compute the starlet transform of a signal I with default options
;       (i.e. a trou algorithm with 4 scales).  
;               Output = STAR2D(I)
;       Reconstruction can be done by simple co-addition of all scales:
;               Rec = total(output, 2)
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

function test_ind_1D, ind, N
	;this function implements mirror like limit conditions on the edges of the input signal
	;ATTENTION : the output may still be out of range ie not in [0, N-1]
	;refinements may be necessary, although not meaningful in practice
	
	ret = ind
	if ind LT 0 then ret = -ind $
	else if ind GE N then ret = 2*N-ind-2
	
	return, ret
	
end

;-------------------------------------------------------------------------------

pro b3spline_istar1d, Sig_in, Sig_out,  Step
   N = double( (size(Sig_in))(1) )
   
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Sig_out  = Sig_in

   for i = 0., N-1. do begin
      im = test_ind_1D(i-Step, N)
      ip = test_ind_1D(i+Step, N)
      im2 = test_ind_1D(i-2*Step, N)
      ip2 = test_ind_1D(i+2*Step, N)
      Sig_out(i) = C3 * Sig_in(i) + C2 * (Sig_in(im) + Sig_in(ip)) + C1 * (Sig_in(im2) + Sig_in(ip2)) 
   endfor
end


pro b3spline_1D, Sig_in, Sig_out,  Step
; help, Sig_in
   N = double( (size(Sig_in))(1) )
   
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Sig_out  = Sig_in

   for i = 0., N-1. do begin
      im = test_ind_1D(i-Step, N)
      ip = test_ind_1D(i+Step, N)
      im2 = test_ind_1D(i-2*Step, N)
      ip2 = test_ind_1D(i+2*Step, N)
      Sig_out(i) = C3 * Sig_in(i) + C2 * (Sig_in(im) + Sig_in(ip)) + C1 * (Sig_in(im2) + Sig_in(ip2)) 
   endfor
end

;-------------------------------------------------------------------------------

function istar1d,  WT, Gen2=Gen2
RecIma = -1

	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: Rec = istar1d(WT, Gen2=Gen2)'
        goto, DONE
	end

	vsize = size(WT)
	if vsize[0] NE 2 then begin
        print, 'Error: First parameter is not a 2D array ...'
        print, 'CALLING SEQUENCE:  WT = istar1d(WT, Gen2=Gen2)'
        goto, DONE
	end
	;-------------------------------------------------------------------------------
		
vs = size(WT)
Nscale = vs[2]
NStep=Nscale-1
Step_trou = 1
for i=0, NStep-2 do Step_trou = Step_trou * 2
	
RecIma = WT[*,Nscale-1]
if keyword_set(Gen2) then Im_Out = fltarr(vs[1])

for j=Nscale-2,0,-1 do begin
  ; print, j, Step_trou
  if keyword_set(Gen2) then  begin
      b3spline_istar1d, RecIma, Im_Out, Step_trou
      ; b3spline_istar1d, RecIma, Im_Out, Step_trou
      RecIma = Im_Out + WT[*,j]
  end else RecIma = RecIma + WT[*,j]
  Step_trou = Step_trou / 2
end

DONE:
    return,  RecIma
end

;-------------------------------------------------------------------------------


