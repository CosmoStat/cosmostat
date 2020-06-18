;+
; NAME:
;        toucan_wavelet
;
; PURPOSE:
;  Computes the wavelet decomposition of a Healpix Map.
;
; CALLING:
;     toucan_wavelet, Imag, out, Nscale=Nscale, Firstl=Firstl,
;                     lmax=lmax, filtWidth=l, lin=lin, filters=bkj

;
; INPUTS:
;     Imag       -- Healpix map = Input CMB image.
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
;                      cutofFreq : limit l values of the filters.
;                        modFreq : mod of the wavelet filters.
;                       TabWidth : widths of wavelet filters in the
;                                  multipole domain.
;
; INPUT KEYWORD: 
;     Nscale      : int = number of scales/filters used in the wavelet
;                   decomposition = number of compressed measurements.
;     Firstl      : int = initial multipole value to be reconstructed in Cl.
;     lmax        : int = final multipole value to be reconstructed in Cl.
;     filtWidth   : int = minimum width of the wavelet filters in the
;                   multipole domain.
;     lin         : bool = if set, the width and spacing of the filters in the
;                   multipole domain is linear.
;     filters     : array of size lmax+1 by Nscale = wavelet filters in the
;                   multipole domain to be used to compute the compressed
;                   measurements.
;
; EXTERNAL CALLS:
;     newspline (idl)
;
; EXAMPLE:
;      Wavelet transform of a CMB map:
;      toucan_wavelet, cmb, needt, Nscale=Nscale,lmax=lmax
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;--------------------------------------------------------------------------------------------------------

;================================ newspline ================================
; compute a spline between firstl and l with cutting frequencies lminj
; and lmaxj and minimum width minwidth. The output is tab of length
; size.

pro newspline, size, firstl, l, lminj, lmaxj, tab, minwidth=minwidth

if not keyword_set(minwidth) then minwidth = lmaxj - lminj

width = min([lmaxj,l]) - max([lminj,firstl])
leftm = max([0.,lminj - firstl])
rightm = max([0.,l-lmaxj])
delt = max([0.,(minwidth - width)/2.0])
delt = delt + max([0., delt - leftm]) + max([0., delt - rightm])
lmaxj = lmaxj + delt
lminj = lminj - delt

res = l * dindgen(size+1)/size
res = 2.0* (2.0 * res - lmaxj-lminj)  / (lmaxj-lminj)
tab =  (3.0/2.0)*1.0 /12.0 * (( abs(res-2))^3 - 4.0* (abs(res-1))^3 + 6 *(abs(res))^3 - 4.0 *( abs(res+1))^3+(abs(res+2))^3)

end

;================================ toucan_wavelet ================================

pro ntousi_needlet, Imag, out, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtWidth=l, lin=lin, filters=bkj

COMMON C_PLANCK

out = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE: toucan_wavelet, Imag, out, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtWidth=l, lin=lin, filters=bkj'
   goto, DONE
endif

if keyword_set(bkj) then begin
   Nscale = (size(bkj))[2]
   lmax = (size(bkj))[1]-1
endif
npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(lmax)  then lmax = nside *3l
if lmax GT P_LMAX then lmax = P_LMAX
if not keyword_set(NScale) then NScale = fix(lmax/2.)
if keyword_set(bkj) then Nscale = (size(bkj))[2]
if not keyword_set(l) then l = 70
if not keyword_set(Firstl) then Firstl=0
if not keyword_set(lin) then log=1
lmarg = l/4

;------------- make logarithmically spaced intervals -----------------
if keyword_set(log) then begin
   log1 = alog(4*(l+Firstl))
   log2 = alog(long(lmax)-Firstl+4*(l+Firstl))
   logs = (log2 - log1)/(Nscale-1)
   Llims = exp( log1 + logs*(findgen(Nscale))) + Firstl - 4*(l+Firstl)

;---------------- make linearly spaced intervals ---------------------
endif else begin ;if lin eq 1 then begin

   Llims = Firstl + (lmax-Firstl)*(findgen(Nscale))/(Nscale-1)

endelse


;---------------------------------------------------------------------
mrs_almtrans, imag, ALM, lmax=lmax
ALM_HighResolImag = ALM.alm
TabWavelet = dblarr(npix, NScale)

Hscale = imag
TabPsi = fltarr( lmax+1, NScale)
modfreq = fltarr(Nscale)
TabWidth = fltarr(Nscale)

;-------------------- Compute needlet coefficients --------------------------
if not keyword_set(bkj) then begin
   if log eq 1 then begin
      ;;; for the first scale do:
      j = 0
      lmaxj = Llims[1]
      lminj = 2*Llims[0] - lmaxj
      
      newspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
      if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
      
      alm_product2, ALM_HighResolImag, h, alm_h
      ALM.alm = alm_h
      mrs_almrec, ALM, Hscale
      TabPsi[*,j] = h
      TabWavelet[*,j] = double(Hscale)
      modfreq[j] = (lmaxj+lminj)/2
      Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
      
   ;;; for intermediate scales do:
      for j=1, Nscale-2 do begin
         lmaxj = Llims[j+1]
         lminj = Llims[j-1]
         
         newspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
         if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
         
         alm_product2, ALM_HighResolImag, h, alm_h
         ALM.alm = alm_h
         mrs_almrec, ALM, Hscale
         TabPsi[*,j] = h
         TabWavelet[*,j] = double(Hscale)   
         modfreq[j] = (lmaxj+lminj)/2
         Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
      endfor
      
   ;;; for the last scale do:
      j = Nscale-1
      lminj = Llims[Nscale-2]
      lmaxj = 2*Llims[Nscale-1] - lminj
      
      newspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
      if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
      
      alm_product2, ALM_HighResolImag, h, alm_h
      ALM.alm = alm_h
      mrs_almrec, ALM, Hscale
      TabPsi[*,j] = h
      TabWavelet[*,j] = double(Hscale)   
      modfreq[j] = (lmaxj+lminj)/2
      Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
      
   endif
endif else begin
   for j=0, Nscale-1 do begin
      h = bkj[*,j]
      if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
      
      alm_product2, ALM_HighResolImag, h, alm_h
      ALM.alm = alm_h
      mrs_almrec, ALM, Hscale
      TabPsi[*,j] = h
      TabWavelet[*,j] = double(Hscale)   
      vmax = max(bkj[*,j],vind)
      modfreq[j] = vind
      lind = where(abs(bkj[0:vind,*]) lt 1.e-2, c)
      if c gt 0 then lminj = max(lind) else lminj = 0
      lind = where (abs(bkj[vind :*,*]) lt 1.e-5, c)
      if c gt 0 then lmaxj = min(lind) else lmaxj = lmax
      Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
   endfor
endelse

ll = findgen(lmax+1)
ll =2*ll + 1
TabNorm = fltarr(NScale)
for j=0, NScale -1 do  TabNorm [j] = sqrt( (total(tabpsi[*, j]^2.*ll)) / double(npix)  / (4.*!dpi) )

out = { NScale: NScale, nside: nside, npix:npix, Coef: TabWavelet, lmax:long(lmax), $
       TabPsi:TabPsi, TabNorm:TabNorm, cutofFreq: Llims, modFreq: modFreq, TabWidth:TabWidth }

DONE:

end

;================================================================================================

;.compile ntousi_needlet
; Nscale = 6
;ntousi_needlet, cmb2, needt, Nscale=Nscale,lmax=lmax
;window,/free & for i=0, Nscale-1 do if i eq 0 then plot, phi[*,i] else oplot, phi[*,i]

; window,/free & for i=0, Nscale-1 do if i eq 0 then plot, out.tabpsi[*,i] else oplot, out.tabpsi[*,i]

