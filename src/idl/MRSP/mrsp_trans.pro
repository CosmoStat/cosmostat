; NAME:
;        MRSP_TRANS
;
; PURPOSE:
;	Compute a transform on polarized spherical data.
;   The transoform can be:
;       0) A E/B decomposition using the spin_2 transform.
;       1) An orthogonal wavelet on each T,Q,U component
;       2) A pyramidal isotropic wavelet on each T,Q,U component 
;       3) An undecimated wavelet transform on each component
;       4) A decimated modulus-Phase wavelet transform 
;       5) A undecimated modulus-Phase wavelet transform 
;       6) A curvelet transform
;       Input image must be in the healPix pixel representation (nested data representation). 
;       The output is a IDL structure.
;
; CALLING:
;     mrsp_trans, InImag, Trans, ebdec=ebdec, NbrScale=NbrScale, Cur=Cur, UWT=UWT, PyrWT=PyrWT, OWT=OWT, $ 
;                   MPDWT=MPDWT, MPUWT=MPUWT, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, 
;                   NeedletWave=NeedletWave, B_NeedletParam=B_NeedletParam, Overlap=Overlap, FirstBlockSize=FirstBlockSize
;
; INPUTS:
;     Imag -- IDL array of a healpix polarized image  fltarr[*,3] : Input image be transformed 
;    
; OUTPUTS: 
;     Trans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;             Nside -- Nside value
;              Npix -- Number of pixels 
;             DEC1 -- IDL structure: First component transformation (depends on the chosen transform)
;             DEC2 -- IDL structure: Second component transformation (depends on the chosen transform)
;             DEC3 -- IDL structure: Third component transformation (depends on the chosen transform)
;            ebdec -- int: set to 1 if an EB decomposiiton has been applied
;            lmax  -- int: maximum used spherical harmonic (for isoptropic wavelet transform only)
;			pyrtrans -- int set to 1 if a pyramidal decomposition has been applied
;	    MeyerWave : int = 1 if the keyword MeyerWave used, otherwise 0
;           NeedletWave: int = 1 if the  keyword  NeedletWave  is used, otherwise 0
;           B_NeedletParam: int = B_NeedletParam (default is 2)
;	    DifInSH : int = 1 if the keyword DifInSH used, otherwise 0
;           TransChoice -- string: Code of the chosen transform.
;           TransTypeName -- string array: array of transform names
;	    TransName -- string transform's name
;	    TabCodeTransform -- string array: array of transform code
;
; INPUT KEYWORDS:
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;         EBDEC -- int: if set, an E/B decomposition is applied before the chosen multiscale decomposition
;         CUR -- int: if set, perform a curvelet transform
;         UWT  -- int: if set, perform a undecimated isotropic transform
;         PyrWT  -- int: if set, perform a pyramidal isotropic transform
;         OWT  -- int: if set, perform a bi-orthogonal wavelet transform on each face
;         MPDWT  -- int: if set, perform a decimated module-phase wavelet transform
;         MPUWT -- int: if set, perform a undecimated module-phase wavelet transform
;		  Overlap -- LONG: is equal to 1 if blocks are overlapping, only used with curvelet transform
;		  FirstBlockSize -- INT: Block size in the ridgelet transform at the finest scale (default is 16), only used with curvelet transform
;		  lmax  -- int: maximum used spherical harmonic (for isoptropic wavelet transform only)
;      DifInSH   : If set, compute the wavelet coefficients as the
;					difference between two resolution in the spherical harmonics representation.
;					Otherwise, the wavelet coefficients are computed as the difference between two resolutions
;					in the initial representation. Only used with /UWT or /PyrWT.
;	   MeyerWave : If set, use Meyer wavelets and set the keyword DifInSH. Only used with /UWT or /PyrWT.
;          NeedletWave: int = 1 if the  keyword  NeedletWave  is used, otherwise 0
;          B_NeedletParam: int = B_NeedletParam (default is 2)
;
; EXTERNAL CALLS:
;         mrsp_tqu2teb
;         mrs_attrans
;         mr_modphase_uwt_trans
;         mrs_owttrans
;         mr_modphase_dwt_trans
;         mrs_pwttrans
;         mrs_wttrans
;         mrs_curtrans
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of a vector field I with five scales
;       The result is stored in WT
;               mrsp_trans, Imag, WT, NbrScale=5, /UWT
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2008
;-

pro mrsp_trans, InImag, Trans, ebdec=ebdec, NbrScale=NbrScale, Cur=Cur, UWT=UWT, PyrWT=PyrWT, OWT=OWT, $ 
               MPDWT=MPDWT, MPUWT=MPUWT, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, $
               Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrsp_trans,InImag, Trans, ebdec=ebdec, NbrScale=NbrScale, Cur=Cur, UWT=UWT, PyrWT=PyrWT, OWT=OWT, MPDWT=MPDWT,'
        print, '                  MPUWT=MPUWT, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, Overlap=Overlap, FirstBlockSize=FirstBlockSize, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam'
        goto, DONE
end
	 
if  not keyword_set(NbrScale) then NbrScale = 4
if NbrScale le 1 or NbrScale ge 20 then begin 
     print,'Error: Number of scales should be between 2 and 20'
     goto, DONE
end

NbrBand = NbrScale
nside=0
pyrtrans=0
npix = (size(Inimag))[1]
nside = npix2nside(npix)
if not keyword_set(ebdec)  then ebdec = 0
if not keyword_set(lmax)  then lmax = nside *3
if not keyword_set(DifInSH)  then DifInSH = 0
if not keyword_set(MeyerWave) then MeyerWave = 0
if keyword_set(PyrWT)  then pyrtrans = 1

TransTypeName = ['EBDEC','Bi-Orthogonal WT', 'Pyramidal WT', 'Undecimated WT', 'Module-Phase Decimated Transform', $
                  'Module-Phase Undecimated Transform', 'Curvelet']
 
TabCodeTransform = ['T_EBDEC', 'T_OWT', 'T_PyrWT', 'T_UWT', 'T_MPDWT', 'T_MPUWT', 'T_CUR']

TransChoice = TabCodeTransform[0]
if keyword_set(OWT) then TransChoice= TabCodeTransform[1]                  
if keyword_set(PyrWT) then TransChoice= TabCodeTransform[2]             
if keyword_set(UWT) then TransChoice= TabCodeTransform[3]                
if keyword_set(MPDWT) then TransChoice= TabCodeTransform[4]                  
if keyword_set(MPUWT) then TransChoice= TabCodeTransform[5]                  
if keyword_set(Cur) then TransChoice= TabCodeTransform[6]                
Coef=0
 
if keyword_set(EBDEC) then mrsp_tqu2teb,InImag, Imag $
else Imag = InImag

if TransChoice EQ 'T_EBDEC'  then begin
  ebdec=1
  mrsp_tqu2teb,InImag, Imag
  Dec1 = Imag[*,0]
  Dec2 = Imag[*,1]
  Dec3 = Imag[*,2]
  NbrScale=1
end

TabNbrBandPerScale = intarr(NbrScale)
TabNbrBandPerScale[*] = 1

; BI-Orthogonal WT
if TransChoice EQ 'T_OWT' then BEGIN
   mrs_owttrans, Imag[*,0], Dec1,  NbrScale=NbrScale
   mrs_owttrans, Imag[*,1], Dec2,  NbrScale=NbrScale
   mrs_owttrans, Imag[*,2], Dec3,  NbrScale=NbrScale
  TabNbrBandPerScale[0:NbrScale-2] = 3
END

; Pyramidal WT
if TransChoice EQ 'T_PyrWT'  then BEGIN
   mrs_pwttrans, Imag[*,0], Dec1,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave
   mrs_pwttrans, Imag[*,1], Dec2,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave
   mrs_pwttrans, Imag[*,2], Dec3,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave
END

; Undecimated WT
if TransChoice EQ 'T_UWT' then BEGIN
   mrs_wttrans, Imag[*,0], Dec1,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam
   mrs_wttrans, Imag[*,1], Dec2,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam
   mrs_wttrans, Imag[*,2], Dec3,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam
END

; Module-Phase Decimated Transform
if TransChoice EQ 'T_MPDWT' then BEGIN
   mrs_owttrans, Imag[*,0], Dec1,  NbrScale=NbrScale
   mrs_modphase_dwt_trans, Imag[*,1:2], Dec2,  NbrScale=NbrScale 
   Dec3 = 0
END

; Module-Phase Undecimated Transform
if TransChoice EQ 'T_MPUWT' then BEGIN
   mrs_attrans, Imag[*,0], Dec1,  NbrScale=NbrScale, modif=modif
   mrs_modphase_uwt_trans, Imag[*,1:2], Dec2,  NbrScale=NbrScale 
   ; Dec3 = Coef[*,*,*, *, 1]
   Dec3 = 0
END

; Curvelet

if TransChoice EQ 'T_CUR' then BEGIN
   mrs_curtrans, Imag[*,0], Dec1, lmax=lmax, NbrScale=NbrScale, Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize, /Silent
   mrs_curtrans, Imag[*,1], Dec2, lmax=lmax, NbrScale=NbrScale, Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize, /Silent
   mrs_curtrans, Imag[*,2], Dec3, lmax=lmax, NbrScale=NbrScale, Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize, /Silent
   TabNbrBandPerScale[0:NbrScale-2] = Dec1.TABNBRSCALERID[0:NbrScale-2]
END

TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 

Ind = where ( TabCodeTransform EQ TransChoice, c)
IndChoiceTrans = Ind[0]

FirstComponentStructField = 7
Trans = { FirstComponentStructField : FirstComponentStructField, $
          TransChoice:TransChoice, $
          TransName: TransTypeName[IndChoiceTrans], $
          IndexTrans: IndChoiceTrans, $
          NbrScale : NbrScale, $
          nside : nside, $
          npix:npix, $
          Dec1: Dec1,  $    ; Field Position of the first decomposition
          Dec2:Dec2, Dec3:Dec3,  $ 
          ebdec:ebdec, lmax:long(lmax), MeyerWave:MeyerWave, $
          DifInSH:DifInSH, pyrtrans:pyrtrans, TabNorm:TabNorm, TransTypeName:TransTypeName, $
          TabCodeTransform:TabCodeTransform, $
          TabNbrBandPerScale: TabNbrBandPerScale }

DONE:

END




