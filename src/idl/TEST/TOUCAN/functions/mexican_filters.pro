;================================ mexicanWF ================================
; compute a mexican transform between firstl and lmax with cutting
; middle frequencies midl. The output is tab of length size.

pro mexicanWF, size, midl, Firstl, lmax, tab
  
  res = ( (lmax * dindgen(size+1)/size) - Firstl )/(midl - Firstl)
  res = res^2.
  tab =  res * exp( - res )
  tab = tab/exp(-1.)
  lstart = max([0,floor(Firstl*size/lmax)])
  tab[0:lstart]=0.

end

;================================ mexican_filters ================================

function mexican_filters, Nscale, Firstl=Firstl, lmax=lmax, Filtwidth=l, lin=lin, Bmat=Bmat
  
;;   if keyword_set(lin) then begin
;;      Llims = Firstl + (lmax-Firstl)*(findgen(Nscale+2))/(Nscale+1)
;;   endif else begin
;;      print, '--- Compute Mexican filters...'
;;      log1 = alog(4*(l+Firstl))
;;      log2 = alog(long(lmax)-Firstl+4*(l+Firstl))
;;      logs = (log2 - log1)/(Nscale+1)
;;      Llims = exp( log1 + logs*(findgen(Nscale+2))) + Firstl - 4*(l+Firstl)
;;   endelse
  
;;   Llims = Llims[1:Nscale]

if keyword_set(lin) then begin
    Llims = Firstl + l/4. + (lmax-2.-l/2.-Firstl)*(findgen(Nscale))/(Nscale-1)
endif else begin
    print, '--- Compute Mexican filters...'
    log1 = alog(4*(l+Firstl))
    log2 = alog(long(lmax)-2.-l/2.-Firstl+4*(l+Firstl))
    logs = (log2 - log1)/(Nscale-1)
    Llims = exp( log1 + logs*(findgen(Nscale))) + Firstl + l/4. - 4*(l+Firstl)
endelse
  
  bkj = fltarr( lmax+1, NScale)
  
  for j=0, Nscale-1 do begin
     midl = Llims[j]
     
     mexicanWF, lmax, midl, (midl-l/2.), lmax, h
     if Firstl gt 0 then h = [fltarr(Firstl),h[Firstl :*]]
     
     bkj[*,j] = h
  endfor
  
  Bmat = dblarr(lmax+1,Nscale)
  lcl = findgen(lmax+1)*2+1
  for i=0, Nscale-1 do begin
     Bmat(*,i) = bkj[*,i]^2*lcl
  endfor
  Bmat = Bmat / (4. * !dpi)
    
  return, bkj
  
end
