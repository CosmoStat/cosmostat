;+
; NAME:
;        mrs_wt_coeff
;
; PURPOSE:
;  Computes the wavelet decomposition of a Healpix Map.
;
; CALLING:
;     mrs_wt_coeff, Imag, out, bkj

;
; INPUTS:
;     Imag       -- Healpix map = Input CMB image.
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;                   Cl.
;    
; OUTPUTS:
;     out        -- IDL structures with the following fields:  
;                         NScale : int = number of wavelet
;                                  scales/filters.
;                          nside : int = Healpix nside parameter.
;                           npix : int = Healpix npix parameter.
;                           Coef : fltarr[npix,NScale] = wavelet
;                                  transform of the data Coef[*,0] =
;                                  wavelet coefficients of the
;                                  coarsest scale (lowest
;                                  frequencies).
;                           lmax : int = maximum l value in the
;                                  Spherical Harmonic Space (Healpix).
;                           npix : long = Number of pixels of the
;                                  input image (12*nside*nside).
;                         TabPsi : IDL array[0:lmax, Nscale-1] =
;                                  wavelet function at resolution
;                                  level.
;
; EXAMPLE:
;      Wavelet transform of a CMB map:
;      mrs_wt_coeff, cmb, needt, bkj
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;--------------------------------------------------------------------------------------------------------

;================================ mrs_wt_coeff ================================

pro mrs_wt_coeff, Imag, out, bkj

COMMON C_PLANCK

out = -1
if N_PARAMS() LT 3  then begin 
   print, 'CALLING SEQUENCE: mrs_wt_coeff, Imag, out, bkj'
   goto, DONE
endif

Nscale = (size(bkj))[2]
lmax = (size(bkj))[1]-1
npix = (size(Imag))[1]
nside = npix2nside(npix)

;-----------------------------------------------------------------------
mrs_almtrans, imag, ALM, lmax=lmax
ALM_HighResolImag = ALM.alm
TabWavelet = dblarr(npix, NScale)

Hscale = imag
TabPsi = fltarr( lmax+1, NScale)

;-------------------- Compute wt coefficients --------------------------
for j=0, Nscale-1 do begin
   h = bkj[*,j]
   alm_product2, ALM_HighResolImag, h, alm_h
   ALM.alm = alm_h
   mrs_almrec, ALM, Hscale
   TabPsi[*,j] = h
   TabWavelet[*,j] = double(Hscale)
endfor

ll = findgen(lmax+1)
ll =2*ll + 1
TabNorm = fltarr(NScale)
for j=0, NScale -1 do  TabNorm [j] = sqrt( (total(tabpsi[*, j]^2.*ll)) / double(npix)  / (4.*!dpi) )

out = { NScale: NScale, nside: nside, npix:npix, Coef: TabWavelet, lmax:long(lmax), $
       TabPsi:TabPsi, TabNorm:TabNorm }

DONE:

end

;================================================================================================
