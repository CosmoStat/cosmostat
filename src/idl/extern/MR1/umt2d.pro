;+
; NAME:
;        UMT2D
;
; PURPOSE:
;       Computes the undecimated median transform of an image. 
;       The output is a 3D IDL array.  
;       
; CALLING:
;      DataTransf = UMT2D(Imag, Nscale=Nscale)
;       
; INPUTS:
;     Imag -- 2D IDL array: image we want to transform
;
; OUTPUTS:
;     DataTransf -- 3D IDL array: Undecimated Median Transform
;                                 DataTransf(*,*,i) = ith band of the transform
;  
; KEYWORDS:
;      Nscale -- int: Number of scales. Default is 4.
;                     DataTransf is a cube of Nscale scales: 
;                                      DataTransf(*,*,0: Nscale-1)
;			There is no test of the validity of the number of scales with respect
;			to the size of the input images. 
; 
; EXAMPLE:
;       Compute the undecimated median transform of an image I with default options
;               Output = UMT2D(I)
;       Reconstruction can be done by simple co-addition of all frames:
;               Rec = total(output, 3)
;
; REFERENCES:
;    J.L. Starck, F. Murtagh, B. Pirenne, and M. Albrecht, 
;   "Astronomical Image Compression based on Noise Suppression",  
;   Publications of the Astronomical Society of Pacific, 108, pp 446--455, 1996.
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

function umt2d, Ima, Nscale=Nscale, wt=wt

        Wtrans = -1
	FirstWindowSize = 2	

	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: WT = umt2d(Ima, Nscale=Nscale)'
        goto, DONE
	end

	if not keyword_set(Nscale) then Nscale =4

	if Nscale LT 2 then  Nscale =4
	
	vsize = size(Ima)
	if vsize[0] NE 2 then begin
        print, 'Error: First parameter is not an image ...'
        print, 'CALLING SEQUENCE:  WT = umt2d(Ima, Nscale=Nscale)'
        goto, DONE
	end
	
	;-------------------------------------------------------------------------------
	;recovering a few parameter values, initializing the transform loop ...
	Nx = vsize[1]
	Ny = vsize[2]
	Nz = Nscale
	NStep=Nscale-1
  	Wtrans = fltarr(Nx,Ny,Nz)
	Im_in = Ima
	Im_Out = fltarr(Nx,Ny)
	WinM = FirstWindowSize
	for i=0, NStep-1 do begin
 		WindowSize = 2*WinM+1
		Im_Out = median(Im_in, WindowSize)
 		Wtrans[*,*,i] = Im_in - Im_Out
		if keyword_set(WT) then begin
 		   Band = Wtrans[*,*,i]
		   ind = where(ABS(BAND) GT 5*mad(band), c)
		   if c GT 0 then Band[ind] = 0
		   Im1 = Im_Out + Band
 		   W1 = star2d(Im1, Nscale=i+2)
		   Im_Out  = W1[*,*,i+1]
		   Wtrans[*,*,i] = Im_in - Im_Out
		end
  		Im_in = Im_Out
 		WinM = WinM * 2
        endfor
	; Smooth array
 	Wtrans[*,*, NStep] = Im_Out
	
DONE:
    return, Wtrans
    
end

;-------------------------------------------------------------------------------


