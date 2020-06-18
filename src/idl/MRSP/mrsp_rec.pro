;+
; NAME:
;        MRSP_REC
;
; PURPOSE:
;	Reconstruct a polarized image from its wavelet transform.
;
; CALLING:
;      MRSP_REC, Trans, Rec
;
; INPUTS:
;       Trans : IDL structure; Wavelet transform structure (see MRSP_TRANS) ;     
;    
; OUTPUTS:
;      Imag -- IDL array of a polarized image fltarr[*,*,3] : Output image be reconstructed 
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of a vector field I with five scales and reconstrcution
;               mrsp_trans, Imag, /UWT, NbrScale=5
;               mrsp_rec, WT, RecIma
;               tvscl, RecIma.coef[*,*,0] ; plot the first reconstructed component   
;               tvscl, RecIma.coef[*,*,1] ; plot the second reconstructed component   
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2008
;-
;-----------------------------------------------------------------

pro mrsp_rec, Trans, RecIma

if N_PARAMS() LT 2 then begin 
        print, 'CALL SEQUENCE: mrsp_rec, WT_Struct, result'
        goto, DONE
        end

npix =Trans.npix
nside = Trans.nside

; Trans = { NbrScale : NbrScale, nside : nside, npix:npix, Dec1: Dec1, Dec2:Dec2, Dec3:Dec3,  $
;          ebdec:ebdec, lmax:lmax, MeyerWave:MeyerWave, $
;          DifInSH:DifInSH, pyrtrans:pyrtrans, TabNorm:TabNorm, TransChoice:TransChoice, TransTypeName:TransTypeName} 
;
;TransTypeName = ['EBDEC','Bi-Orthogonal WT', 'Pyramidal WT', 'Undecimated WT', 'Module-Phase Decimated Transform', $
;                  'Module-Phase Undecimated Transform', 'Curvelet']

TransChoice = Trans.TransChoice
  
Imag = fltarr(npix, 3)

if  TransChoice EQ 'T_EBDEC' then BEGIN
  Imag[*,0] = Trans.Dec1
  Imag[*,1] = Trans.Dec2
  Imag[*,2] = Trans.Dec3
END

; BI-Orthogonal WT
if TransChoice EQ 'T_OWT' then BEGIN
   mrs_owtrec, Trans.Dec1, rec
   Imag[*,0] = rec
   mrs_owtrec,  Trans.Dec2 , rec
   Imag[*,1] = rec
   mrs_owtrec,  Trans.Dec3 , rec
   Imag[*,2] = rec
END

; Pyramidal WT
if TransChoice EQ 'T_PyrWT'  then BEGIN
   mrs_pwtrec,  Trans.Dec1 , rec
   Imag[*,0] = rec
   mrs_pwtrec,  Trans.Dec2 , rec
   Imag[*,1] = rec
   mrs_pwtrec,  Trans.Dec3 , rec
   Imag[*,2] = rec
END

; Undecimated WT
if TransChoice EQ 'T_UWT' then BEGIN
   mrs_wtrec,  Trans.Dec1, rec
   Imag[*,0] = rec
   mrs_wtrec,  Trans.Dec2, rec 
   Imag[*,1] = rec
   mrs_wtrec, Trans.Dec3, rec 
   Imag[*,2] = rec
END

; Module-Phase Decimated Transform
if TransChoice EQ 'T_MPDWT' then BEGIN
   mrs_owtrec,  Trans.Dec1 , rec
   Imag[*,0] = rec
   mrs_modphase_dwt_rec, Trans.Dec2, rec
   Imag[*,1] = rec[*,0]
   Imag[*,2] = rec[*,1]
END

; Module-Phase Undecimated Transform
if TransChoice EQ 'T_MPUWT' then BEGIN
   mrs_atrec,  Trans.Dec1, rec
   Imag[*,0] = rec
   mrs_modphase_uwt_rec, Trans.Dec2 , rec
   Imag[*,1] = rec[*,0]
   Imag[*,2] = rec[*,1]
 END

; Curvelet
if TransChoice EQ 'T_CUR' then BEGIN
   mrs_currec,  Trans.Dec1, rec 
   Imag[*,0] = rec
   mrs_currec,  Trans.Dec2, rec
   Imag[*,1] = rec
   mrs_currec, Trans.Dec3 , rec
   Imag[*,2] = rec
END

if  Trans.EBDEC EQ 1 then mrsp_teb2tqu, Imag,  RecIma $
else RecIma  = Imag

DONE:

END


;================================================

pro test_rec, Trans, RecIma

n = randomn(seed, 32L^2*12, 3)*100.

  mrsp_tqu2teb, n, neb
  mrsp_teb2tqu, neb,  r
  mrsp_tqu2teb, r, neb
  mrsp_teb2tqu, neb,  rn
print, "EB REC", sigma(n-r)


mrsp_trans, n, t, /ebdec
mrsp_rec, t, r
print, "EB REC", sigma(n-r)

mrsp_trans, n, t, /owt
mrsp_rec, t, r
print, "OWT REC", sigma(n-r)

mrsp_trans, n, t, /pyrwt
mrsp_rec, t, r
print, "PyrWT REC", sigma(n-r)

mrsp_trans, n, t, /uwt
mrsp_rec, t, r
print, "UWT REC", sigma(n-r)

mrsp_trans, n, t, /MPDWT
mrsp_rec, t, r
print, "MPDWT REC", sigma(n-r)

mrsp_trans, n, t, /MPUWT
mrsp_rec, t, r
print, "MPUWT REC", sigma(n-r)

n = randomn(seed, 64L^2*12, 3)*100.
mrsp_trans, n, t, /Cur
mrsp_rec, t, r
print, "CUR REC", sigma(n-r)

end


