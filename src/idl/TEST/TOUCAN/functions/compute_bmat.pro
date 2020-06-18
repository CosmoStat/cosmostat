function compute_Bmat, filters=bkj, Nscale=Nscale, Firstl=Firstl, lmax=lmax, Filtwidth=l

if keyword_set(bkj) then begin
   Nscale = (size(bkj))[2]
   lmax = (size(bkj))[1]-1
endif

   print, '--- Compute sensing matrice Bmat...'

if not keyword_set(bkj) then begin

   log1 = alog(4*(l+Firstl))
   log2 = alog(long(lmax)-Firstl+4*(l+Firstl))
   logs = (log2 - log1)/(Nscale-1)
   Llims = exp( log1 + logs*(findgen(Nscale))) + Firstl - 4*(l+Firstl)

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

endif

Bmat = dblarr(lmax+1,Nscale)
lcl = findgen(lmax+1)*2+1
for i=0, Nscale-1 do begin
   Bmat(*,i) = bkj[*,i]^2*lcl
endfor
Bmat = Bmat / (4. * !dpi)

print, '--- Bmat computed.'
   
return, Bmat

end
