;+
; NAME:
;        toucan_filters
;
; PURPOSE:
;  Computes the filters that will be used to compute the integrated
;  power spectrum.
;
; CALLING:
;     toucan_filters(Nscale, lmax=lmax, Firstl=Firstl, Filtwidth=l,
;     lin=lin, Bmat=Bmat, cutofFreq=Llims)
;
; INPUTS:
;     Nscale     -- int = number of scales used in the wavelet
;                   decomposition = number of compressed measurements.
;    
; OUTPUTS:
;     filters    -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;                   Cl.
;
; INPUT KEYWORDS: 
;     lmax        : int = final multipole value to be reconstructed in
;                   Cl.
;     Firstl      : int = initial multipole value to be reconstructed.
;     filtWidth   : int = minimum width of the wavelet filters in the
;                   multipole domain.
;     lin         : bool = if set, the width and spacing of the filters in the
;                   multipole domain is linear.
;
; OUTPUT KEYWORD:
;            Bmat : dblarr(lmax+1,Nscale) = matrix formed from the
;                   filters. When multiplying Cl by Bmat, one obtains
;                   the integrated Cl.
;       cutofFreq : flt = limit values of the filters in the multipole
;                   domain.
;
; EXTERNAL CALLS:
;     newspline (idl)
;
; EXAMPLE:
;      Creation on toucan filters:
;      bkj = toucan_filters(Nscale, lmax=2*nside, Firstl=2, Filtwidth=140)
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
;; delt = delt + max([0., delt - leftm]) + max([0., delt - rightm])
lmaxj = lmaxj + delt
lminj = lminj - delt

res = l * dindgen(size+1)/size
res = 2.0* (2.0 * res - lmaxj-lminj)  / (lmaxj-lminj)
tab =  (3.0/2.0)*1.0 /12.0 * (( abs(res-2))^3 - 4.0* (abs(res-1))^3 + 6 *(abs(res))^3 - 4.0 *( abs(res+1))^3+(abs(res+2))^3)

end

;================================ toucan_filters ================================

function toucan_filters, Nscale, Firstl=Firstl, lmax=lmax, Filtwidth=l, lin=lin, Bmat=Bmat, cutoffreq=Llims
  
;;   if keyword_set(lin) then begin
;;      Llims = Firstl + (lmax-Firstl)*(findgen(Nscale))/(Nscale-1)
;;   endif else begin
;;      print, '--- Compute Toucan filters...'
;;      log1 = alog(4*(l+Firstl))
;;      log2 = alog(long(lmax)-Firstl+4*(l+Firstl))
;;      logs = (log2 - log1)/(Nscale-1)
;;      Llims = exp( log1 + logs*(findgen(Nscale))) + Firstl - 4*(l+Firstl)
;;   endelse

  if keyword_set(lin) then begin
     Llims = Firstl + l/4. + (lmax-2.-l/2.-Firstl)*(findgen(Nscale))/(Nscale-1)
  endif else begin
     print, '--- Compute Toucan filters...'
     log1 = alog(4*(l+Firstl))
     log2 = alog(long(lmax)-2.-l/2.-Firstl+4*(l+Firstl))
     logs = (log2 - log1)/(Nscale-1)
     Llims = exp( log1 + logs*(findgen(Nscale))) + Firstl + l/4. - 4*(l+Firstl)
  endelse

  bkj = fltarr( lmax+1, NScale)
  
  j = 0
  lmaxj = Llims[1]
  lminj = 2*Llims[0] - lmaxj
  
  newspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
  if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
  
  bkj[*,j] = h
  
   ;;; for intermediate scales do:
  for j=1, Nscale-2 do begin
     lmaxj = Llims[j+1]
     lminj = Llims[j-1]
     
     newspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
     if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
     
     bkj[*,j] = h
  endfor
  
   ;;; for the last scale do:
  j = Nscale-1
  lminj = Llims[Nscale-2]
  lmaxj = 2*Llims[Nscale-1] - lminj
  
  newspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
  if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
  
  bkj[*,j] = h
  
  Bmat = dblarr(lmax+1,Nscale)
  lcl = findgen(lmax+1)*2+1
  for i=0, Nscale-1 do begin
     Bmat(*,i) = bkj[*,i]^2*lcl
  endfor
  Bmat = Bmat / (4. * !dpi)
    
  return, bkj
  
end
