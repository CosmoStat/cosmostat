;+
; NAME: 
;       MRS_MODPHASE_UWT_REC
;
; PURPOSE: 
;        Reconstruct a vector field image on the Sphere from its undecimated wavelet 
;        transform (see mrs_modphase_uwt_trans).   
;
; CALLING:
;       MRS_MODPHASE_UWT_REC, WT_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Wavelet transform structure (see MRS_MODPHASE_UWT_TRANS) 
;          
; OUTPUTS:
;      Result:  2D array of an Healpix vector field (nested format)
;             
; EXTERNAL CALLS
;            mr_modphase_dwt_rec
;
; EXAMPLE:
;   mrs_modphase_uwt_trans, Imag, WT, NbrScale=5
;   mrs_modphase_uwt_rec, WT, RecIma
;       Wavelet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck, May 2008.
;-
;-----------------------------------------------------------------

pro mrs_modphase_uwt_rec, Trans, RecIma 

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mrs_modphase_uwt_rec, WT_Struct, result'
        goto, DONE
        end
        
RecIma = fltarr(Trans.Npix, 2)

CubeFace1 = fltarr(Trans.nside, Trans.nside, 12)
CubeFace2 = fltarr(Trans.nside, Trans.nside, 12)
for f=0,11 do begin   ;begin loop over 12 healpix faces
        WT = {NBRSCALE: Trans.NbrScale, MODCOEFF: reform(Trans.Coef[*,*,*,f, 0]), ANGCOEFF: reform(Trans.Coef[*,*,*,f, 1])}
        ; hs, wt
        mr_modphase_uwt_rec, WT, Rec
        CubeFace1[*,*,f] = Rec[*,*,0] 
        CubeFace2[*,*,f]  = Rec[*,*,1]
 end    ;end loop over 12 healpix faces

put_all_faces, CubeFace1, Ima1
put_all_faces, CubeFace2, Ima2

RecIma[*,0] = Ima1
RecIma[*,1] = Ima2

DONE:

end


 



