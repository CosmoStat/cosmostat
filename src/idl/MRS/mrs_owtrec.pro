;+
; NAME: 
;       MRS_OWTREC
;
; PURPOSE: 
;        Reconstruct an image on the Sphere from its (bi-) orthogoanl wavelet 
;        transform (see mrs_owttrans).   
;
; CALLING:
;       MRS_OWTREC, WT_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Wavelet transform structure (see MRS_OWTTRANS) 
;          
; OUTPUTS:
;      Result:  1D array of an Healpix image (NESTED format)
;             
; EXTERNAL CALLS
;            bwt01_lift written by Olivier Forni and Nabila Aghanim
;
; EXAMPLE:
;   mrs_owttrans, Imag, WT, NbrScale=5
;   mrs_owtrec, WT, RecIma
;       Wavelet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

pro mrs_owtrec, Trans, HealpixIma 

COMMON MR1ENV


if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mrs_owtrec, WT_Struct, result'
        goto, DONE
        end
        
CubeFace = fltarr(Trans.Nx, Trans.Ny, 12)
NameSig = 'xx_signal.fits'

for f=0,11 do begin						;begin loop over 12 healpix faces

  if Trans.bin EQ 1 then begin				;if package mre is available and was used in the direct transform
		
		Coef = Trans.Coef[*,*,f]
		mr_recons, Coef, Rec, Header=Trans.HeadTrans
		CubeFace[*,*,f] = Rec 
    
  endif else begin							;use mrs for reconstruction as it was used for the direct transform
  
		Coef = Trans.Coef[*,*,f]
		;rec = BWT01_INVERSE( Coef, Trans.NbrScale)
		rec = bwt01_lift( Coef, Trans.NbrScale-1, /inverse)
  		CubeFace[*,*,f] = Rec 

  endelse
end										;end loop over 12 healpix faces


put_all_faces, CubeFace, HealpixIma

delete, NameSig
DONE:
end


