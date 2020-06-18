; NAME:
;        MRS_TRANS
;
; PURPOSE:
;	Compute a transform on polarized spherical data.
;   The transoform can be:
;       0) Spherical Harmonic Transform.
;       1) Orthogonal wavelet  (per face)
;       2) A trou undecimated wavelet transform (per face)
;       2) Pyramidal isotropic wavelet  on the sphere
;       4) Undecimated wavelet transform  on the sphere
;       5) Ridgelet transform  on the sphere
;       6) Curvelet transform  on the sphere
;       7) Discrete Cosine Transform (per face).
;       Input image must be in the healPix pixel representation (nested data representation). 
;       The output is a IDL structure.
;
; CALLING:
;     mrs_trans, InImag, Trans,  NbrScale=NbrScale, Cur=Cur, AT=AT, UWT=UWT, PyrWT=PyrWT, OWT=OWT, $ 
;                   Rid= Rid, DCT= DCT, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave, B_NeedletParam=B_NeedletParam, 
;                   Overlap=Overlap, FirstBlockSize=FirstBlockSize
;
; INPUTS:
;     Imag -- IDL array of a healpix  image : Input image be transformed 
;    
; OUTPUTS: 
;     Trans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;             Nside -- Nside value
;              Npix -- Number of pixels 
;             DEC1 -- IDL structure: Transformation (depends on the chosen transform)
;            lmax  -- int: maximum used spherical harmonic (for isoptropic wavelet transform only)
;			pyrtrans -- int set to 1 if a pyramidal decomposition has been applied
;	     	MeyerWave : int = 1 if the keyword MeyerWave used, otherwise 0
;           NeedletWave: int = 1 if the  keyword  NeedletWave  is used, otherwise 0
;           B_NeedletParam: int = B_NeedletParam (default is 2)
;		DifInSH : int = 1 if the keyword DifInSH used, otherwise 0
;      TransChoice -- string: Code of the chosen transform.
;      TransTypeName -- string array: array of transform names
;		TransName -- string transform's name
;		TabCodeTransform -- string array: array of transform code
;
; INPUT KEYWORDS:
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;         CUR -- int: if set, perform a curvelet transform
;         UWT  -- int: if set, perform a undecimated isotropic transform
;         PyrWT  -- int: if set, perform a pyramidal isotropic transform
;         OWT  -- int: if set, perform a bi-orthogonal wavelet transform on each face
;         WT1D - int: if set, a 1d-1d undecimated wavelet transform is applied.
;		  Overlap -- LONG: is equal to 1 if blocks are overlapping, only used with curvelet transform
;		  FirstBlockSize -- INT: Block size in the ridgelet transform at the finest scale (default is 16), only used with curvelet transform
;		  lmax  -- int: maximum used spherical harmonic (for isoptropic wavelet transform only)
;      DifInSH   : If set, compute the wavelet coefficients as the
;					difference between two resolution in the spherical harmonics representation.
;					Otherwise, the wavelet coefficients are computed as the difference between two resolutions
;					in the initial representation. Only used with /UWT or /PyrWT.
;	   MeyerWave : If set, use Meyer wavelets and set the keyword DifInSH. Only used with /UWT or /PyrWT.
;      NeedletWave:  If set, use  Needlet wavelet instead of Cubic spline
;      B_NeedletParam: float: needlet parameter.
;
; EXTERNAL CALLS:
;         mrs_attrans
;         mrs_owttrans
;         mrs_pwttrans
;         mrs_wttrans
;         mrs_ridtrans
;         mrs_curtrans
;         mrs_dcttrans
;         mrs_wt1d1d
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of an image
;       The result is stored in WT
;               mrs_trans, Imag, WT, NbrScale=5, /UWT
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2008
;-

pro mrs_trans, Imag, Trans,  NbrScale=NbrScale, Cur=Cur, Alm=Alm, AT=AT, Rid=Rid, UWT=UWT, PyrWT=PyrWT, OWT=OWT, DCT=DCT, $ 
               lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam, $
               Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize,wt1d=wt1d,  gen2=gen2

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_trans,InImag, Trans, Alm=Alm, NbrScale=NbrScale, Cur=Cur, UWT=UWT, PyrWT=PyrWT, OWT=OWT, AT=AT,'
        print, '                  Rid=Rid, DCT=DCT, WT1D=WT1D, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam, Overlap=Overlap, FirstBlockSize=FirstBlockSize,  gen2=gen2'
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
npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(Alm)  then Alm = 0
if not keyword_set(lmax)  then lmax = nside *3
if not keyword_set(DifInSH)  then DifInSH = 0
if not keyword_set(MeyerWave) then MeyerWave = 0
if not keyword_set(NeedletWave) then NeedletWave = 0
if not keyword_set(B_NeedletParam) then B_NeedletParam = 2

if keyword_set(PyrWT)  then pyrtrans = 1

TransTypeName = ['Alm','Bi-Orthogonal WT', 'A_Trous WT', 'Pyramidal WT', 'Undecimated WT', 'Ridgelet Transform', $
                  'Curvelet', 'DCT','WT-1D1D']
 
TabCodeTransform = ['T_ALM', 'T_OWT', 'T_AT', 'T_PyrWT', 'T_UWT', 'T_Ridgelet', 'T_CUR',  'T_DCT','T_WT1D']

TransChoice = TabCodeTransform[0]
if keyword_set(OWT) then TransChoice= TabCodeTransform[1]       
if keyword_set(AT) then TransChoice= TabCodeTransform[2]             
if keyword_set(PyrWT) then TransChoice= TabCodeTransform[3]             
if keyword_set(UWT) then TransChoice= TabCodeTransform[4]                
if keyword_set(Rid) then TransChoice= TabCodeTransform[5]                  
if keyword_set(Cur) then TransChoice= TabCodeTransform[6]                  
if keyword_set(DCT) then TransChoice= TabCodeTransform[7]    
if keyword_set(WT1D) then TransChoice= TabCodeTransform[8]                
            
Coef=0

if TransChoice EQ 'T_ALM'  then begin
   mrs_almtrans, Imag, Dec1, lmax=lmax
   NbrScale=1
end

TabNbrBandPerScale = intarr(NbrScale)
TabNbrBandPerScale[*] = 1

; BI-Orthogonal WT
if TransChoice EQ 'T_OWT' then BEGIN
   mrs_owttrans, Imag, Dec1,  NbrScale=NbrScale
   TabNbrBandPerScale[0:NbrScale-2] = 3
END

; Pyramidal WT
if TransChoice EQ 'T_PyrWT'  then BEGIN
   mrs_pwttrans, Imag, Dec1,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave
END

; Undecimated WT
if TransChoice EQ 'T_UWT' then BEGIN
   mrs_wttrans, Imag, Dec1,  NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam
  END

; A trous Transform
if TransChoice EQ 'T_AT' then BEGIN
   mrs_attrans, Imag, Dec1,  NbrScale=NbrScale, Opt=Opt, modif=modif, healpix=healpix
  END

; Module-Phase Undecimated Transform
if TransChoice EQ 'T_Ridgelet' then BEGIN
   mrs_ridtrans, Imag, Dec1, NbrScale=NbrScale, overlap=overlap, blocksize= FirstBlockSize, Opt=Opt  
 END

; Curvelet
if TransChoice EQ 'T_CUR' then BEGIN
   mrs_curtrans, Imag, Dec1, lmax=lmax, NbrScale=NbrScale, Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize, /Silent
   TabNbrBandPerScale[0:NbrScale-2] = Dec1.TABNBRSCALERID[0:NbrScale-2]
END

; DCT
if TransChoice EQ 'T_DCT' then BEGIN
  mrs_dcttrans, Imag, Dec1,  blockSize= FirstBlockSize
  NbrScale=1
END

if TransChoice EQ 'T_WT1D' then BEGIN
  mrs_wt1d1dtrans, Imag, Dec1,  gen2=gen2
  NbrScale=1
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
          Dec: Dec1,  $    ; Field Position of the first decomposition
          lmax:long(lmax), MeyerWave:MeyerWave, NeedletWave:NeedletWave, B_NeedletParam:B_NeedletParam,  $
          DifInSH:DifInSH, pyrtrans:pyrtrans, TabNorm:TabNorm, TransTypeName:TransTypeName, $
          TabCodeTransform:TabCodeTransform, $
          TabNbrBandPerScale: TabNbrBandPerScale }

DONE:

END




