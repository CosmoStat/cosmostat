;+
; NAME:
;        mrs_wttrans
;
; PURPOSE:
;	Computes the undecimated isotropic wavelet transform on the sphere, 
;   using the HEALPix representation (NESTED DATA REPRESENTATION) or Glesp Data representation.
;	The wavelet function is zonal and its spherical harmonics coefficients a_l0 follow 
;   a cubic box-spline profile. If DifInSH is set, wavelet coefficients are derived in the Spherical
;   Harmonic Space, otherwise (default) they are derived in the direct space.
;
;   If the keyword MeyerWave then the Meyer wavelet function is used instead of the Cubic spline.
;   If the keyword NeedletWave then the Needlet wavelet function is used instead of the Cubic spline, and 
;   the B-needlet parameter of the wavelet functuion can be modified using the keyword B_NeedletParam.
;
;
; CALLING:
;
;     mrs_wttrans, Imag, Trans, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, NeedletWave=NeedletWave, Healpix_with_Glesp=Healpix_with_Glesp, B_NeedletParam=B_NeedletParam
;       
;
; INPUTS:
;     Imag -- IDL array of healpix map or Glesp image IDL structure: Input image to be transformed
;
; INPUT/OUTPUT:
;	  lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter (0 for a Glesp image)
;                      lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;                      npix : long = Number of pixels of the input image (12*nside*nside)
;                      Coef : fltarr[npix,NbrScale] = wavelet transform of the data
;                             Coef[*,0] = wavelet coefficients of the finest scale (highest frequencies).
;                             Coef[*,NbrScale-1] = coarsest scale (lowest frequencies). 
;                      lmax : int = lmax parameter at the first scale
;		           Healpix_with_Glesp : int = 1 if the keyword Healpix_with_Glesp used, otherwise 0
;				   UseGLESP : int = 1 if the input image was in Glesp format, otherwise 0
;				  MeyerWave : int = 1 if the keyword MeyerWave is used, otherwise 0
;                 NeedletWave: int = 1 if the  keyword  NeedletWave  is used, otherwise 0
;                 B_NeedletParam: int = B_NeedletParam (default is 2)
;   		    tabphi -- IDL array[0:lmax, Nscale-1]  = Scaling function  at resolution level
;   		    tabpsi -- IDL array[0:lmax, Nscale-1]   = Wavelet function  at resolution level
;   		   TabFilterH --IDL array[0:lmax, Nscale-1] = filter H allowing to go from a resolution to the next one
;    		   TabFilterG -- IDL array[0:lmax, Nscale-1] = filter G allowing to compute the wavelet coeff from the previous resolution level
;					DifInSH : int = 1 if the keyword DifInSH is used, otherwise 0
;						 nx : int = number of rings (Glesp parameter), otherwise, 0
;						 np : int = max number of pixel on a ring (Glesp parameter), otherwise, 0
;					  x_sky : fltarray with COS( THETA ) for each ring (Glesp parameter), otherwise, 0
;					  y_sky : long 1D array number of pixels/ring (Glesp parameter), otherwise, 0
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4). If it is set to -1, then the number scales is:  log(lmax) / log(2)  - 2
;      Lmax      : Number of spherical harmonics computed in the decomposition
;					(default is 3*nside, should be between 2*nside and 4*nside)
;      DifInSH   : If set, compute the wavelet coefficients as the
;					difference between two resolution in the spherical harmonics representation.
;					Otherwise, the wavelet coefficients are computed as the difference between two resolutions
;					in the initial representation.
;	   MeyerWave : If set, use Meyer wavelets and set the keyword DifInSH
;	   Healpix_with_Glesp : If set, a copy of Imag is done in Glesp format in order to compute the wavelet transform
;      NeedletWave:  If set, use  Needlet wavelet instead of Cubic spline
;      B_NeedletParam: float: needlet parameter. Default is 2.
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;   	synfast (healpix software)
;   	alm_product2 (idl)
;   	compute_g (idl)
;   	compute_h (idl)
;
; EXAMPLE:
;
;       Compute the multiresolution of an image I with default options
;        The result is stored in Output
;               mrs_wttrans, Imag, Output, NbrScale=5
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2004
;	December, 2004 File creation
;---------------------------------------------------------------------------------------------------------------------------------------------

pro mrs_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, Healpix_with_Glesp=Healpix_with_Glesp, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam

COMMON C_PLANCK

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, Healpix_with_Glesp=Healpix_with_Glesp, NeedletWave=NeedletWave,  B_NeedletParam=B_NeedletParam'
        goto, DONE
        end
	 
				    
if type_code(Imag) EQ 8 then begin
   GLESP=1  
   Healpix_with_Glesp = 0
end else GLESP=0

nside=0
nx = 0
np = 0
x_sky = 0
y_sky = 0

if keyword_set(Healpix_with_Glesp) then  begin
   NpixOrig = (size(imag))[1]
   NsideOrig = npix2nside(npixorig)
   OrigImag=Imag
   Imag = healpix2glesp(Imag)
   GLESP =  1
   DifInSH=0
end else Healpix_with_Glesp = 0

if GLESP EQ 0 then begin
  npix = (size(imag))[1]
  nside = npix2nside(npix)
  if not keyword_set(lmax)  then lmax = nside *3
  if lmax GT P_LMAX then lmax = P_LMAX
end else begin
   npix = (size(imag.t_sky))[1]
   nx = imag.nx
   np = imag.np
   x_sky = imag.x_sky
   y_sky = imag.y_sky
   if not keyword_set(lmax)  then lmax = min([(nx-1)/2,np/4])
end

if not keyword_set(DifInSH)  then DifInSH = 0
if not keyword_set(MeyerWave) then MeyerWave = 0 else DifInSH = 1
if not keyword_set(NeedletWave) then NeedletWave = 0 else DifInSH = 1
if not keyword_set(B_NeedletParam) then B_NeedletParam = 2
if not keyword_set(NbrScale) then NbrScale = 4
if NbrScale  EQ -1 then NbrScale  = ceil(alog(Lmax)/alog(2.)) - 2.

if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end

ech = 1.
Hscale = imag
if keyword_set(Healpix_with_Glesp) then  H_Hscale = OrigImag 

mrs_almtrans, imag, ALM, lmax=lmax 
ALM_HighResolImag = ALM.alm
index=ALM.index
ech =1.
nlmax2 =lmax
							
if keyword_set(Healpix_with_Glesp) then   TabWavelet = dblarr(NpixOrig, NbrScale) $
else TabWavelet = dblarr(npix, NbrScale)

if  keyword_set(NeedletWave) then begin
  pyrtrans = 0
   
     NedF = mrs_needlet_filters(B=B_NeedletParam,  Lmax, tabphi=tabphi, tabpsi=tabpsi, TabFilterH=TabFilterH, TabFilterG= TabFilterG, Nscale=NbrScale)
     A = ALM
   	for j=0, NbrScale-1 do begin
         g = NedF[*,j]
         alm_product2, ALM_HighResolImag, g, Almg
         A.alm = Almg
        mrs_almrec, A,  WScale
        if GLESP EQ 0 then TabWavelet[*,j] = WScale $
	    else TabWavelet[*,j] = WScale.t_sky
     end
end else begin 
for j=0,NbrScale-2 do begin
  if keyword_set(MeyerWave) then begin
     hgmey, nlmax2, j, h, g, dif=dif 
  end else begin
     compute_g, lmax, ech, g
     compute_h, lmax, ech, h
  end
  
  if (j EQ 0) then begin
    TabFilterH = fltarr( N_ELEMENTS(h), NbrScale-1)
    TabFilterG = fltarr( N_ELEMENTS(h), NbrScale-1)
    TabPhi = fltarr( N_ELEMENTS(h), NbrScale-1)
    TabPsi = fltarr( N_ELEMENTS(h), NbrScale-1)
    TabFilterH[*,j] = h
    TabFilterG[*,j] = 1. - h
    TabPhi[*,j] = h
    TabPsi[*,j] =  1. - h
  end else begin
    TabFilterH[*,j] = h
     TabFilterG[*,j] =  TabFilterH[*,j-1] - h
     TabPhi[*,j] =  TabPhi[*,j-1] * h
     TabPsi[*,j]  =  (1. -h) * TabPhi[*,j-1]
 end
 
  alm_product2, ALM_HighResolImag ,h, alm_h  
  if keyword_set(DifInSH) then begin
        alm_product2, ALM_HighResolImag, g, alm_g  
        ALM.alm = alm_g
	mrs_almrec, ALM, WScale
 	if GLESP EQ 0 then TabWavelet[*,j] = WScale $
	else TabWavelet[*,j] = WScale.t_sky
  end else begin
        ALM.alm = alm_h
	mrs_almrec, ALM, LScale
	if GLESP EQ 0 then TabWavelet[*,j] = double(Hscale) - double(LScale) $
	else if keyword_set(Healpix_with_Glesp) then begin
	   H_LScale = glesp2healpix(LScale, nside=NsideOrig)
	   TabWavelet[*,j] =  double(H_Hscale) - double(H_LScale)
	   H_Hscale = H_LScale
	end else TabWavelet[*,j] =  Hscale.t_sky - LScale.t_sky
	Hscale = LScale
  end
  ALM_HighResolImag = alm_h
  ech = ech*2
endfor

j=NbrScale-1
 
if keyword_set(DifInSH) then begin
   ALM.alm = alm_h
   mrs_almrec, ALM, LScale
   if GLESP EQ 0 then TabWavelet[*,j] = LScale $
   else  TabWavelet[*,j] = LScale.t_sky
end else begin 
   if GLESP EQ 0 then TabWavelet[*,j] = LScale $
   else if not keyword_set(Healpix_with_Glesp) then TabWavelet[*,j] = LScale.t_sky $
   else begin
      H_LScale = glesp2healpix(LScale, nside=NsideOrig)
      TabWavelet[*,j] = H_LScale
  end
end
 pyrtrans = 0
if not keyword_set(MeyerWave) then MeyerWave = 0 else MeyerWave = 1
endelse ; needlet filters

if keyword_set(Healpix_with_Glesp) then  begin
   Imag = OrigImag
   GLESP = 0
   nside = NsideOrig
end

ll = findgen(lmax+1)
ll =2*ll + 1
TabNorm = fltarr(NbrScale)
for j=0, NbrScale -2 do  TabNorm [j] = sqrt((total(tabpsi[*, j]^2.*ll)) / double(npix))
TabNorm [NbrScale -1] =  sqrt((total(tabphi[*, NbrScale -2]^2.*ll)) / double(npix))

; TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 
out = {UseGLESP: GLESP, NbrScale : NbrScale, nside : nside, nx: nx, np:np, npix:npix, Coef : TabWavelet, lmax:long(lmax), MeyerWave:MeyerWave, $
       TabFilterH:TabFilterH, TabFilterG:TabFilterG, TabPhi:TabPhi, TabPsi:TabPsi, $
       DifInSH:DifInSH, pyrtrans:pyrtrans, x_sky:x_sky,  y_sky :y_sky, TabNorm:TabNorm, Healpix_with_Glesp: Healpix_with_Glesp, $
       NeedletWave:NeedletWave,  B_NeedletParam:B_NeedletParam}
DONE:

END

