;+
; NAME: 
;       MRS_MODPHASE_DWT_REC
;
; PURPOSE: 
;        Reconstruct a vector field image on the Sphere from its decimated wavelet 
;        transform (see mrs_modphase_dwt_trans).   
;
; CALLING:
;       MRS_MODPHASE_DWT_REC, WT_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Wavelet transform structure (see MRS_MODPHASE_DWT_TRANS) 
;          
; OUTPUTS:
;      Result:  2D array of an Healpix vector field (nested format)
;             
; EXTERNAL CALLS
;            bwt01_lift written by Olivier Forni and Nabila Aghanim
;
; EXAMPLE:
;   mrs_modphase_dwt_trans, Imag, WT, NbrScale=5
;   mrs_modphase_dwt_rec, WT, RecIma
;       Wavelet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck, May 2008.
;-
;-----------------------------------------------------------------

pro mrs_modphase_dwt_rec, Trans, RecIma 

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mrs_modphase_dwt_rec, WT_Struct, result'
        goto, DONE
        end
        
RecIma = fltarr(Trans.Npix, 2)
CubeFace1 = fltarr(Trans.nside, Trans.nside, 12)
CubeFace2 = fltarr(Trans.nside, Trans.nside, 12)
for f=0,11 do begin   ;begin loop over 12 healpix faces
        WT = {NBRSCALE: Trans.NbrScale, MODCOEFF: Trans.Coef[*,*,f, 0], ANGCOEFF: Trans.Coef[*,*,f, 1]}
        mr_modphase_dwt_rec, WT, Rec
        CubeFace1[*,*,f] = Rec[*,*,0] 
        CubeFace2[*,*,f]  = Rec[*,*,1]
 end    ;end loop over 12 healpix faces


put_all_faces, CubeFace1, Ima1
put_all_faces, CubeFace2, Ima2

RecIma[*,0] = Ima1
RecIma[*,1] = Ima2

DONE:

end


 



