;+
; NAME:
;        mrs_mmt
;
; PURPOSE:
;	Computes the multiscale median transform of an Healpix spherical image. 
;
; CALLING:
;
;     mrs_mmt, Imag, Trans, NbrScale=NbrScale
;
; INPUTS:
;     Imag -- IDL array of healpix image: Input image to be transformed
;
; INPUT/OUTPUT:
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter (0 for a Glesp image)
;                      npix : long = Number of pixels of the input image (12*nside*nside)
;                      Coef : fltarr[npix,NbrScale] = wavelet transform of the data
;                             Coef[*,0] = wavelet coefficients of the finest scale (highest frequencies).
;                             Coef[*,NbrScale-1] = coarsest scale (lowest frequencies). 
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4). If it is set to -1, then the number scales is:  log(lmax) / log(2)  - 2
;
; EXAMPLE:
;
;       Compute the multiscale median transform of an image I  
;        The result is stored in Output
;               mrs_mmt, Imag, Output, NbrScale=5
;         
; HISTORY:
;	Written:   Jean-Luc Starck, 2011
;---------------------------------------------------------------------------------------------------------------------------------------------

pro mrs_mmt, Imag, out, NbrScale=NbrScale, WindowSize=WindowSize

COMMON C_PLANCK

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_mmt, Imag, out, NbrScale=NbrScale'
        goto, DONE
        end
	    
GLESP=0
nside=0
nx = 0
np = 0
x_sky = 0
y_sky = 0
Healpix_with_Glesp = 0
pyrtrans = 0

npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(NbrScale) then NbrScale = 4
if NbrScale  EQ -1 then NbrScale  = ceil(alog(Lmax)/alog(2.)) - 2.

if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end

if not keyword_set(WindowSize) then Win = 3 else Win = WindowSize
Hscale = imag 							
TabWavelet = dblarr(npix, NbrScale)
 TabWavelet[*,0] = Hscale
for j=0,NbrScale-2 do begin
   print, ' Scale ', j+1
    TabWavelet[*,j+1] = mrs_median( TabWavelet[*,j], windowsize=Win)
    TabWavelet[*,j] = TabWavelet[*,j] - TabWavelet[*,j+1]
    Win = Win*2
endfor
 
TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 
out = {UseGLESP: GLESP, NbrScale : NbrScale, nside : nside, nx: nx, np:np, npix:npix, Coef : TabWavelet, lmax:0, MeyerWave: 0, $
       TabFilterH: 0, TabFilterG: 0, TabPhi: 0, TabPsi: 0, $
       DifInSH: 0, pyrtrans: pyrtrans, x_sky: 0,  y_sky : 0, TabNorm:TabNorm, Healpix_with_Glesp: 0, $
       NeedletWave: 0,  B_NeedletParam: 0}
DONE:

END

