;+
; NAME:
;        mrs_toucan_wt_filters
;
; PURPOSE:
;  Computes the filters that will be used to compute the integrated
;  power spectrum.
;
; CALLING:
;     mrs_toucan_wt_filters(Nscale, lmax, Firstl=Firstl, Filtwidth=l,
;     lin=lin, Bmat=Bmat, cutofFreq=Llims, modFreq=modFreq)
;
; INPUTS:
;     Nscale     -- int = number of scales used in the wavelet
;                   decomposition = number of compressed measurements.
;     lmax       -- int = final multipole value to be reconstructed in
;                   Cl.
;    
; OUTPUTS:
;     filters    -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;                   Cl.
;
; INPUT KEYWORDS: 
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
;         modFreq : flt = mod of the filters.
;
; EXTERNAL CALLS:
;     touspline (idl)
;
; EXAMPLE:
;      Creation on toucan filters:
;      bkj = mrs_toucan_wt_filters(Nscale, lmax, Firstl=2, Filtwidth=140)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;--------------------------------------------------------------------------------------------------------

;================================ touspline ================================
; compute a spline between firstl and l with cutting frequencies lminj
; and lmaxj and minimum width minwidth. The output is tab of length
; size.

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

;================================ mrs_toucan_wt_filters ================================

function mrs_toucan_wt_filters, Nscale, lmax, Firstl=Firstl, Filtwidth=l, lin=lin, Bmat=Bmat, cutofFreq=Llims, modFreq=modFreq

bkj = -1
if N_PARAMS() LT 2  then begin 
   print, 'CALLING SEQUENCE:  res = mrs_toucan_wt_filters(Nscale, lmax, Firstl=Firstl, Filtwidth=l, lin=lin, Bmat=Bmat, cutofFreq=Llims, modFreq=modFreq)'
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

print, '--- Compute Toucan filters...'
bkj = fltarr( lmax+1, NScale)
modfreq = fltarr(Nscale)
TabWidth = fltarr(Nscale)

;; for the first scale do:
j = 0
lmaxj = Llims[1]
lminj = 2*Llims[0] - lmaxj

touspline, lmax, Firstl, lmax, lminj, lmaxj, h, minwidth=l
if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]

bkj[*,j] = h
modfreq[j] = (lmaxj+lminj)/2
Tabwidth[j] = min([lmaxj,lmax])-max([lminj,Firstl])

;; for intermediate scales do:
for j=1, Nscale-2 do begin
   lmaxj = Llims[j+1]
   lminj = Llims[j-1]
   
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

Bmat = dblarr(lmax+1,Nscale)
lcl = findgen(lmax+1)*2+1
for i=0, Nscale-1 do begin
   Bmat(*,i) = bkj[*,i]^2*lcl
endfor
Bmat = Bmat / (4. * !dpi)

DONE:
return, bkj

end
