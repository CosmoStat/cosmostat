;+
; NAME:
;        mrs_wt_powspec
;
; PURPOSE:
;        Computes the wavelet binned power spectrum, corresponding to 
;        the power spectrum integrated over filters of a wavelet decomposition.
;
; CALLING:
;     wtps  =  mrs_wt_powspec(Map,  filters=filters, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtWidth= filtWidth)
;
; INPUTS:
;     Map        -- Healpix map: Input CMB image to compute the
;                   compressed measurements from.
;    
; OUTPUTS:
;     wtps     -- IDL dblarr(Nscale, 3) = binned Cl 
;                   									wtps[i,0] = beginning of the band
;                                                       wtps[i,1] = end of the band
;                                                       wtps[i,2] = binned Cl value in the bad
;
; INPUT/OUTPUT KEYWORD: 
;     filters     : array of size lmax+1 by Nscale = wavelet filters
;                   in the multipole domain to be used to compute the
;                   compressed measurements.
;     Mask        : Healpix map = mask applied to the input Map.
;     Nscale      : int = number of scales used in the wavelet
;                   decomposition = number of binned measurements.
;     Firstl      : int = initial multipole value to be reconstructed in Cl.
;     filtWidth   : int = minimum width of the wavelet filters in the
;                   multipole domain.
;     lmax        : int = final multipole value to be reconstructed in Cl.
;
; EXAMPLE:
;      Estimate the binned power spectrum from a CMB map  :
;      wtps  =  mrs_wt_powspec(Map)
;
; HISTORY:
;       Written:  Aurele Balavoine & J.L. Starck , 2012
;----------------------------------------------------------------------------


pro touspline, size, firstl, l, lminj, lmaxj, tab, minwidth=minwidth

if not keyword_set(minwidth) then minwidth = lmaxj - lminj

width = min([lmaxj,l]) - max([lminj,firstl])
leftm = max([0.,lminj - firstl])
rightm = max([0.,l-lmaxj])
delt = max([0.,(minwidth - width)/2.0])
;; delt = delt + max([0., delt - leftm]) + max([0., delt - rightm])
lmaxj = lmaxj + delt
lminj = lminj - delt

res = l * dindgen(size+1)/size
res = 2.0* (2.0 * res - lmaxj-lminj)  / (lmaxj-lminj)
tab =  (3.0/2.0)*1.0 /12.0 * (( abs(res-2))^3 - 4.0* (abs(res-1))^3 + 6 *(abs(res))^3 - 4.0 *( abs(res+1))^3+(abs(res+2))^3)

end

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

;================================ mrs_get_ps_wt_filters ================================

function mrs_get_ps_wt_filters, Nscale, lmax, Firstl=Firstl, Filtwidth=l, lin=lin, Bmat=Bmat, cutofFreq=Llims, modFreq=modFreq, Tabwidth=Tabwidth, TabBeginEndBand=TabBeginEndBand

bkj = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE:  res = mrs_get_ps_wt_filters(Nscale, lmax, Firstl=Firstl, Filtwidth=l, lin=lin, Bmat=Bmat, cutofFreq=Llims, modFreq=modFreq)'
   goto, DONE
endif

if not keyword_set(Firstl) then Firstl = 0
if not keyword_set(Filtwidth) then Filtwidth = 100

;---------------- make linearly spaced intervals ---------------------
if keyword_set(lin) then begin
   Llims = Firstl + l/4. + (lmax-2.-l/2.-Firstl)*(findgen(Nscale))/(Nscale-1)
endif else begin
;------------- make logarithmically spaced intervals -----------------
   log1 = alog(4*(l+Firstl))
   log2 = alog(long(lmax)-2.-l/2.-Firstl+4*(l+Firstl))
   logs = (log2 - log1)/(Nscale-1)
   Llims = exp( log1 + logs*(findgen(Nscale))) + Firstl + l/4. - 4*(l+Firstl)
endelse

bkj = fltarr( lmax+1, NScale)
modfreq = fltarr(Nscale)
TabWidth = fltarr(Nscale)
TabBeginEndBand = fltarr(Nscale, 2)
;; for the first scale do:
j = 0
lmaxj = Llims[1]
lminj = 2*Llims[0] - lmaxj

touspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]

bkj[*,j] = h
modfreq[j] = (lmaxj+lminj)/2
Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
TabBeginEndBand[j,0] = min([lmaxj,lmax])
TabBeginEndBand[j,1] = max([lminj,Firstl])

;; for intermediate scales do:
for j=1, Nscale-2 do begin
   lmaxj = Llims[j+1]
   lminj = Llims[j-1]
   TabBeginEndBand[j,0] = lminj
   TabBeginEndBand[j,1] = lmaxj
   ; print, j+1, ' [ ', lminj, ', ', lmaxj, ']'
   touspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
   if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
   
   bkj[*,j] = h
   modfreq[j] = (lmaxj+lminj)/2
   Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
endfor

;; for the last scale do:
j = Nscale-1
lminj = Llims[Nscale-2]
lmaxj = 2*Llims[Nscale-1] - lminj

touspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]

bkj[*,j] = h
modfreq[j] = (lmaxj+lminj)/2
Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])
TabBeginEndBand[j,0] = min([lmaxj,lmax])
TabBeginEndBand[j,1] = max([lminj,Firstl])

Bmat = dblarr(lmax+1,Nscale)
lcl = findgen(lmax+1)*2+1
for i=0, Nscale-1 do begin
   Bmat(*,i) = bkj[*,i]^2*lcl
endfor
Bmat = Bmat / (4. * !dpi)

DONE:
return, bkj

end

;================================ toucan_wavelet ================================

pro wavelet_transform, Imag, out, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtWidth=l, lin=lin, filters=bkj

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

;-------------------- Compute wavelet coefficients --------------------------

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

;;==========================  binned powspec =============================

function mrs_wt_powspec, Map,  filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, filtwidth=filtwidth, scaledMask=scaledMask, Bmat=Bmat, rescaleMask=rescaleMask, EdgeMask=maskEdge, cutoffreq=cutoffreq, Trans=Trans, lin=lin

ang_ps = -1
if N_PARAMS() LT 1  then begin 
   print, 'CALLING SEQUENCE: Res =  mrs_wt_powspec(Map,  filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, scaledMask=scaledMask, Bmat=Bmat, niter=niter,  cutoffreq=cutoffreq)'
   goto, DONE
endif

npix = N_elements(Map)
nside = npix2nside(npix)

;;---------------------------------------- Set default parameters ----------------------------------------
if not keyword_set(lmax) then lmax =2* nside

if not keyword_set(Filtwidth) then Filtwidth = 70
if not keyword_set(Firstl) then Firstl=0
if not keyword_set(Nscale) then Nscale = round(4*alog(lmax))

bkj = mrs_get_ps_wt_filters(Nscale, lmax, Firstl=Firstl, Filtwidth=Filtwidth, lin=lin, Bmat=Bmat, modFreq=modFreq, Tabwidth=Tabwidth, TabBeginEndBand=TabBeginEndBand)
TabReducePowSpec = fltarr(Nscale, 3)

ll = dindgen(lmax+1)
l2 =  2*ll+1
TN = dblarr(Nscale)
for i=0, Nscale-1 do   TN[i] = 4. * !DPI    / total(bkj[*,i]^2 * l2)
 for i=0, Nscale-1 do bkj[*,i]= bkj[*,i] *  sqrt( ll * (ll+1) * TN[i] )

;; ---------------- Needlet transform  
if not keyword_set(Trans) then begin
   print, 'Computing Wavelet  decomposition over', Nscale, ' scales...'
   wavelet_transform, Map, Trans, NScale=Nscale, Firstl=Firstl, lmax=lmax, filters=bkj, filtwidth=filtwidth
endif

;; ---------------- Parameters of the needlet transform 
npix = Trans.Npix
Nscale = Trans.NScale
nside = Trans.Nside
lmax = Trans.lmax ; lmax = 3*nside
ell = findgen(lmax+1)*2+1
bkj = Trans.tabpsi
cutoffreq = Trans.cutoffreq

ang_ps = dblarr(Nscale)

;; ;------ If a mask is present
if keyword_set(Mask) then begin
      nMask = double( total(Mask))
      for i=0, Nscale-1 do  ang_ps[i] = total(Mask*Trans.coef[*,i]^2)/nMask   
endif else begin
      for i=0, Nscale-1 do   ang_ps[i] = total(Trans.coef[*,i]^2)/ double(npix)
endelse

if arg_present(Bmat) then begin
;   print, 'computing Bmat'
   Bmat = dblarr(lmax+1,Nscale)
   lcl = findgen(lmax+1)
   lcl = 2*lcl+1  
   for i=0, Nscale-1 do  Bmat(*,i) = bkj[*,i]^2*lcl
   Bmat = Bmat / (4. * !dpi)
endif

lcl = findgen(lmax+1)*2+1
print, Nscale
for i=0, Nscale-1 do begin
   TabReducePowSpec[i,0] = TabBeginEndBand[i,0]
   TabReducePowSpec[i,1] =  TabBeginEndBand[i,1]
   l =  (TabBeginEndBand[i,0] + TabBeginEndBand[i,1]) / 2.
   TabReducePowSpec[i,2] =  ang_ps[i]   ; * 4. * !DPI  * l* (l+1) /  total( lcl * bkj[*,i]^2 )
  ;  if i EQ Nscale-1 then  print, ang_ps[i], ' ' , lcl[i], ' ' , total( lcl * bkj[*,i]^2 ), ' ', ang_ps[i]   * 4. * !DPI    
end
  
plotclr, TabReducePowSpec
DONE:
return, TabReducePowSpec
end
