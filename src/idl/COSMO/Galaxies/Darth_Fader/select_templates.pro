pro select_templates, ntemp, zmax=zmax, wigglez=wigglez, test = test, $
                      indir = indir

;; this code randomly selects ntemp training set spectra up to a
;; limiting redshift zmax

if not keyword_set(zmax) then zmax = 0.1
if not keyword_set(indir) then indir = '/export/zen/aleonard/CMC_templates/'

loadct,5
if not keyword_set(wigglez) then $
   readcol,indir+'all_z_sorted.txt',z $
else $
   readcol,indir+'WiggleZfinalV04_99.dat',$
           f='x,x,x,f,f,a,x,x,x,f,x,x,x,x,x,x,f', z, zerr, q, gmag, $
           gerr, skipline = 5

if keyword_set(test) then zmax = max(z)

if not keyword_set(wigglez) then $
   idx = where(z le zmax,count) $
else idx = where(z gt 0 and zerr lt 0.0005 and q gt 3,count)

while n_elements(ind) lt ntemp do begin
   regen:   indt = round(randomu(seed)*count)
   if n_elements(ind) eq 0 then begin
      ind = indt 
      goto, regen
   endif else begin
      idx2 = where(ind eq indt, ct)
      if ct gt 0 then goto, regen
   endelse
   ind = [ind,indt]
endwhile
ind = ind[sort(ind)]
ind = idx[ind]
ztemps = z[ind]
print,min(ztemps),max(ztemps)

if not keyword_set(wigglez) then begin
   for i = 0, 9 do begin
      incl = where(ind ge long(i)*9026L and ind lt long(i+1)*9026L, count)
      if count gt 0 then begin
         inds = ind[incl] mod 9026
         print,"i range", long(i)*9026L, long(i+1)*9026L - 1
         print,n_elements(ind[incl])
         print,count
      endif
      if count gt 0 then begin
         rdfits_struct,indir+'specs_'+strtrim(i,2)+$
                       '.fits',inp 
         if n_elements(temps) eq 0 then temps=inp.tab1[*,inds] else $
            temps = [[temps],[inp.tab1[*,inds]]]
         if n_elements(lamtemps) eq 0 then lamtemps=inp.im0[*,inds] else $
            lamtemps = [[lamtemps],[inp.im0[*,inds]]]
      endif
   endfor
endif else begin
   temps = dblarr(5100,ntemp)
   lamtemps = dblarr(5100,ntemp)
   noise = dblarr(5100,ntemp)
   for i = 0, ntemp - 1 do begin
      if ind[i]+1 gt 0 and ind[i]+1 lt 10 then prefix = '00000'
      if ind[i]+1 ge 10 and ind[i]+1 lt 100 then prefix = '0000'
      if ind[i]+1 ge 100 and ind[i]+1 lt 1000 then prefix = '000'
      if ind[i]+1 ge 1000 and ind[i]+1 lt 10000 then prefix = '00'
      if ind[i]+1 ge 10000 and ind[i]+1 lt 100000 then prefix = '0'
      if ind[i]+1 ge 100000 then prefix = ''

      filename = indir+'spectra/wig'+prefix+strtrim(ind[i]+1,2)+'.fits'
      spec = read_wigglez_data(filename)
      temps[0:n_elements(spec.spectrum)-1,i] = spec.spectrum
      lamtemps[0:n_elements(spec.spectrum)-1,i] = spec.lambda
      noise[0:n_elements(spec.spectrum)-1,i] = spec.variance
   endfor
endelse

;; now blueshift the spectra

for i = 0, n_elements(ind) - 1 do $
   lamtemps[*,i] = lamtemps[*,i]/(1.+ztemps[i])

if not keyword_set(wigglez) then begin
   writefits,indir+'training_set.fits',temps
   writefits,indir+'training_set.fits',lamtemps,$
             /append
endif else begin
   writefits,indir+'training_set.fits',temps
   writefits,indir+'training_set.fits',lamtemps,$
             /append
   writefits,indir+'training_set.fits',noise,/append
endelse

return
end
