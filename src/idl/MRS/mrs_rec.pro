;+
; NAME:
;        MRS_REC
;
; PURPOSE:
;	Reconstruct an image from its decomposition.
;
; CALLING:
;      MRS_REC, Trans, Rec
;
; INPUTS:
;       Trans : IDL structure; Transform structure (see MRS_TRANS) ;     
;    
; OUTPUTS:
;      Rec -- IDL array of a  healpix  image : Output image be reconstructed 
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of an image,  with five scales and reconstrcution
;               mrs_trans, Imag, WT,  /UWT, NbrScale=5
;               mrs_rec, WT, RecIma
;               tvscl, RecIma  ; plot the  reconstructed image   
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2008
;-
;-----------------------------------------------------------------

pro mrs_rec, Trans, Rec

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mrs_rec, WT_Struct, result'
        goto, DONE
        end

npix =Trans.npix
nside = Trans.nside

; Trans = { NbrScale : NbrScale, nside : nside, npix:npix, Dec: Dec,    $
;           lmax:lmax, MeyerWave:MeyerWave, $
;          DifInSH:DifInSH, pyrtrans:pyrtrans, TabNorm:TabNorm, TransChoice:TransChoice, TransTypeName:TransTypeName} 
;
; 
;TransTypeName = ['Alm','Bi-Orthogonal WT', 'A_Trous WT', 'Pyramidal WT', 'Undecimated WT', 'Ridgelet Transform', 'Curvelet', 'DCT']
;TabCodeTransform = ['T_ALM', 'T_OWT', 'T_AT', 'T_PyrWT', 'T_UWT', 'T_Ridgelet', 'T_CUR',  'T_DCT']

TransChoice = Trans.TransChoice
  
Imag = fltarr(npix, 3)

if  TransChoice EQ 'T_ALM' then BEGIN
  mrs_almrec, Trans.Dec, rec
END

; BI-Orthogonal WT
if TransChoice EQ 'T_OWT' then BEGIN
   mrs_owtrec, Trans.Dec, rec
   END

; Pyramidal WT
if TransChoice EQ 'T_PyrWT'  then BEGIN
   mrs_pwtrec,  Trans.Dec , rec
 END

; Undecimated WT
if TransChoice EQ 'T_UWT' then BEGIN
   mrs_wtrec,  Trans.Dec, rec
 END

; A trous wavelet Transform
if TransChoice EQ 'T_AT' then BEGIN
   mrs_atrec,  Trans.Dec, rec
 END

; Ridgelet Transform
if TransChoice EQ 'T_Ridgelet' then BEGIN
   mrs_ridrec,  Trans.Dec, rec
  END

; Curvelet
if TransChoice EQ 'T_CUR' then BEGIN
   mrs_currec,  Trans.Dec, rec 
 END

; Curvelet
if TransChoice EQ 'T_DCT' then BEGIN
   mrs_dctrec,  Trans.Dec,  rec 
 END

if TransChoice EQ 'T_WT1D' then BEGIN
  mrs_wt1d1drec, Trans.Dec,  rec 
 END

DONE:

END


;================================================

pro test_rec, Trans, RecIma

n = randomn(seed, 32L^2*12)*100.

mrs_trans, n, t, /alm
mrs_rec, t, r
print, "ALM REC", sigma(n-r)

mrs_trans, n, t, /owt
mrs_rec, t, r
print, "OWT REC", sigma(n-r)

mrs_trans, n, t, /pyrwt
mrs_rec, t, r
print, "PyrWT REC", sigma(n-r)

mrs_trans, n, t, /uwt
mrs_rec, t, r
print, "UWT REC", sigma(n-r)

mrs_trans, n, t, /at
mrs_rec, t, r
print, "AT REC", sigma(n-r)

mrs_trans, n, t, /Rid
mrs_rec, t, r
print, "RID REC", sigma(n-r)

mrs_trans, n, t, /DCT
mrs_rec, t, r
print, "DCT REC", sigma(n-r)

n = randomn(seed, 64L^2*12)*100.
mrs_trans, n, t, /Cur
mrs_rec, t, r
print, "CUR REC", sigma(n-r)

end


