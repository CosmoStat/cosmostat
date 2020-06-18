function df_rebin_spectra, inlam, inspec, outlam, innoise=innoise, $
                           outnoise=outnoise, lanczos_a=lanczos_a

;; innoise and outnoise are *variances*

if not keyword_set(lanczos_a) then lanczos_a=4
nx = n_elements(inspec)

lam0 = inlam[0]
lam1 = inlam[nx-1]
in_pix_lstep = inlam[1]-inlam[0]

npix = n_elements(outlam)
outspec = dblarr(npix)
if keyword_set(outnoise) then outnoise=dblarr(npix)

if keyword_set(innoise) then begin
   weight = 1./innoise
   idx=where(finite(innoise) eq 0,count)
   if count gt 0 then weight[idx]=0.
endif
for i = 0, npix - 1 do begin

   dum = (outlam[i]-inlam)/in_pix_lstep
      
   idx=where(abs(dum) le lanczos_a,count)
   if count eq 0 then outspec[i] = 0. else begin
      outspec[i] = total(inspec[idx]*lanczos(dum[idx],lanczos_a)) 
      if n_elements(outnoise) gt 0 then outnoise[i] = 1./(total((lanczos(dum[idx],lanczos_a))^2/weight[idx]))
   endelse
         
endfor

if n_elements(outnoise) gt 0 then outnoise = 1./(outnoise)


return,outspec
end


