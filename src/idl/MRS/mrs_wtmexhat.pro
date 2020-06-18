;+
; NAME:
;        mrs_wtmexhat
;
; PURPOSE:
;	Convolution with the mexican hat wavelet function
;
; CALLING:
;
;     Scale = mrs_wtmexhat( Image, ScaleParameter )    
;    
; INPUT:
;     Imag -- IDL array of healpix map in RING format: Input image to be convolved with the wavelet function
;     ScaleParameter -- float:  Scale parameter in arc minute. 
;
; OUTPUT:
;     map_out  -- IDL array of healpix map: wavelet coefficients  
;
; EXAMPLE:
;       Compute a  wavelet scale (here, the mexican hat wavelet will have an angular spread of about pi/6) 
;               coef = mrs_wtmexhat(Image, sqrt(3.)/3.)    
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck & Patricio Vielva
;	September, 2005 File creation
;-----------------------------------------------------------------


pro alm2al0,alm,al0

taille = size (alm)
taille = taille[1]
index2lm2,taille,l,m
al0 = fltarr(l)
for i=0l,l-1 do begin
   lm2index2,i,0,index
   al0[i] = real_part( alm[index] )				;attention au melange entre type reel et type complex.
endfor											;al0 est suppose etre un type complexe lors de l'appel a almproduct2
												;toutefois, idl ne bute pas sur ce probleme de typage
end

;==============================================================

function mrs_wtmexhat, Map, ScaleParameter 
out=-1
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: coef = mrs_wtmexhat(Data,ScaleParameter)'
        goto, DONE
        end

scale  = ScaleParameter / 60. * !dtor
 
npix = (size(map))[1]
nside = npix2nside(npix)
pi = !dpi
scale = double(scale)


;calcul de l'ondelette continu
tab = fltarr(npix)

;normalization factor
n = scale*sqrt(1.+scale^2./2.+scale^4./4.)


; building a theta phi map the same siz as the input image in  ring order
i = lindgen(npix)
pix2ang_ring, nside, i, theta, phi

; building a map of the wavelet function
r = 2. * tan (theta/2.)
tab(i) = 1./sqrt(2.0*pi) / n *(2.0-(r/scale)^2)*exp( - r^2 /(2.*scale^2) )*(1+(r/2.0)^2)^2

;convolution in the spherical harmonics domain.  
nlmax = 3*nside
alm_cxxtrans, Map, Almimg, nside=nside, nlmax=nlmax, index=index
alm_cxxtrans, tab, AlmWave, /ring, nside=nside, nlmax=nlmax, index=index
alm2al0, almwave, al0
alm_product2,almimg,al0,alm_out
alm_cxxitrans, alm_out, Out, nside=nside, nlmax=nlmax, index=index

DONE:

return, out
end



