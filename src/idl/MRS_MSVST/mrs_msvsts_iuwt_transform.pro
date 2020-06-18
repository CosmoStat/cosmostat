;+
; NAME:
;        mrs_msvsts_iuwt_transform
;
; PURPOSE:
;	Computes the multi-scale variance stabilising transform
;   on the sphere with undecimated isotropic wavelet transform, 
;   using the HEALPix representation (nested data
;   representation). The wavelet function  is zonal and 
;   its spherical harmonics coefficients a_l0 follow 
;   a cubic box-spline profile. If DifInSH is set, wavelet coefficients are derived in the Spherical
;   Harmonic Space, otherwise (default) they are derived in the direct space.
;
;
; CALLING:
;
;      mrs_msvsts_iuwt_transform, Imag, Trans, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH
;       
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image to be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter
;                      lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;                      npix : int = Number of pixels of the input image (12*nside*nside)
;                      Coef : fltarr[npix,NbrScale] = stabilised wavelet transform of the data
;                             Coef[*,0] = stabilised wavelet coefficients of the finest scale (highest frequencies).
;                             Coef[*,NbrScale-1] = coarsest scale (lowest frequencies). 
;                      lmax: int= lmax parameter at the first scale
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4)
;      Lmax      : Number of spherical harmonics computed in the decomposition
;					(default is 3*nside, should be between 2*nside and 4*nside)
;      DifInSH   : If set, compute the wavelet coefficients as the
;					difference between two resolution in the spherical harmonics representation.
;					Otherwise, the wavelet coefficients are computed as the difference between two resolutions
;					in the initial representation.
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
;               mrs_msvsts_iuwt_transform, Imag, Output, NbrScale=5
;         
; HISTORY:
;	Written: Jérémy Schmitt, 2010
;	February, 2010 File creation
;--------------------------------------------------------------------------------------------------------

function vst,x,b,c
y=b*sgn(x+c)*sqrt(abs(x+c))
return,y
end

pro mrs_msvsts_iuwt_transform, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, Healpix_with_Glesp=Healpix_with_Glesp 

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH'
        goto, DONE
        end
	 
if  not keyword_set(NbrScale) then NbrScale = 4
if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
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


mrs_msvsts_iuwt_param_computing,nbrscale,c,b,hh,tau1,tau2,tau3,sigma

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


ech = 1.
Hscale = imag
if keyword_set(Healpix_with_Glesp) then  H_Hscale = OrigImag 

mrs_almtrans, imag, ALM, lmax=lmax 
ALM_HighResolImag = ALM.alm
index=ALM.index
ech =1.
nlmax2 =lmax
							
if keyword_set(Healpix_with_Glesp) then   TabWavelet = fltarr(NpixOrig, NbrScale) $
else TabWavelet = fltarr(npix, NbrScale)

for j=0,NbrScale-2 do begin
  if keyword_set(MeyerWave) then begin
     hgmey, nlmax2, j, h, g, dif=dif 
  end else begin
     compute_g, lmax, ech, g
     compute_h, lmax, ech, h
  end
  alm_product2, ALM_HighResolImag ,h, alm_h  
  if keyword_set(DifInSH) then begin
        alm_product2, ALM_HighResolImag, g, alm_g  
        ALM.alm = alm_g
	mrs_almrec, ALM, WScale
 	if GLESP EQ 0 then TabWavelet[*,j] = vst(WScale,b[j],c[j]) $
	else TabWavelet[*,j] = vst(WScale.t_sky,b[j],c[j])
  end else begin
        ALM.alm = alm_h
	mrs_almrec, ALM, LScale
	if GLESP EQ 0 then TabWavelet[*,j] = vst(Hscale,b[j],c[j]) - vst(LScale,b[j+1],c[j+1]) $
	;if GLESP EQ 0 then TabWavelet[*,j] = Hscale - LScale $
	else if keyword_set(Healpix_with_Glesp) then begin
	   H_LScale = glesp2healpix(LScale, nside=NsideOrig)
	   TabWavelet[*,j] =  vst(H_Hscale,b[j],c[j]) - vst(H_LScale,b[j+1],c[j+1])
	   H_Hscale = H_LScale
	end else TabWavelet[*,j] =  vst(Hscale.t_sky,b[j],c[j]) - vst(LScale.t_sky,b[j+1],c[j+1])
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


if keyword_set(Healpix_with_Glesp) then  begin
   Imag = OrigImag
   GLESP = 0
   nside = NsideOrig
end

TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 
out = {UseGLESP: GLESP, NbrScale : NbrScale, nside : nside, nx: nx, np:np, npix:npix, Coef : TabWavelet, lmax:lmax, MeyerWave:MeyerWave, $
       DifInSH:DifInSH, pyrtrans:pyrtrans, x_sky:x_sky,  y_sky :y_sky, TabNorm:TabNorm, Healpix_with_Glesp: Healpix_with_Glesp}
DONE:

END

