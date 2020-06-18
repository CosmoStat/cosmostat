pro templates_from_data, ntemp, sample=sample, indir = indir

if not keyword_set(indir) then indir = '/export/zen/aleonard/wigglez/'

if keyword_set(sample) then begin
   print,'Reading input redshift catalogue'
   readcol,indir+'WiggleZfinalV04_99.dat',$
           f='x,x,x,d,d,i,d,d,d,d,d,d,d,d,d,d,d,d',$
           z, zerr, q, fuv, nuv, umag, gmag, rmag, imag, zmag, $
           fuverr, nuverr, uerr, gerr, rerr, ierr, zerr, skipline=5

;; select indices matching certain criteria
   print,'Selecting subsample'
   idx = where(rmag lt 20,count)
   
;; now randomly sample from this subsample
   
   print,'Generating random sample'
   while n_elements(ind) lt ntemp do begin
      regen: indt = round(randomu(seed)*count)
      if n_elements(ind) eq 0 then begin
         ind = indt
         goto,regen
      endif else begin
         idx2 = where(ind eq indt,ct)
         if ct gt 0 then goto, regen
      endelse
      ind = [ind,indt]
   endwhile
   
   ind = ind[sort(ind)]
   ind = idx[ind]
   ztemps = z[ind]
   
   temps = dblarr(5100,ntemp)
   lamtemps = dblarr(5100,ntemp)
   noise = dblarr(5100,ntemp)
   print,'Reading in selected spectra'
   for i = 0, ntemp - 1 do begin
      if ind[i]+1 gt 0 and ind[i]+1 lt 10 then prefix = '00000'
      if ind[i]+1 ge 10 and ind[i]+1 lt 100 then prefix = '0000'
      if ind[i]+1 ge 100 and ind[i]+1 lt 1000 then prefix = '000'
      if ind[i]+1 ge 1000 and ind[i]+1 lt 10000 then prefix = '00'
      if ind[i]+1 ge 10000 and ind[i]+1 lt 100000 then prefix = '0'
      if ind[i]+1 ge 100000 then prefix = ''
      
      filename = '/export/zen/aleonard/wigglez/spectra/wig'+prefix+strtrim(ind[i]+1,2)+'.fits'
      
      spec = read_wigglez_data(filename)
      df_mask_lines, spec.lambda, [6362.,5577.,5893.,6300.], mask, tol=4

      temps[0:n_elements(spec.spectrum)-1,i] = spec.spectrum*(1-mask)
      lamtemps[0:n_elements(spec.spectrum)-1,i] = spec.lambda
      noise[0:n_elements(spec.spectrum)-1,i] = spec.variance
      
   endfor

   print,'Blueshifting selected spectra'
   for i = 0, ntemp -1  do $
      lamtemps[*,i] = lamtemps[*,i]/(1.+ztemps[i])
   
   writefits,indir+'training_set.fits',temps
   writefits, indir+'training_set.fits',lamtemps,$
             /append
   writefits,indir+'training_set.fits',noise,/append
   
endif else begin
   temps = mrdfits(indir+'training_set.fits',0)
   lamtemps = mrdfits(indir+'training_set.fits',1)
   noise = mrdfits(indir+'training_set.fits',2)

endelse

;; now rebin the templates

input = df_input_params_wigglez()
print,'Generating pixellation scheme'
outlam_temp = df_gen_log_pix(input.training_lmin, $
                             input.training_lmax, $
                             delta=input.lstep)

training = dblarr(n_elements(outlam_temp),ntemp)
outnoise = dblarr(n_elements(outlam_temp),ntemp)

print,'Rebinning spectra and noise'
for i = 0, ntemp - 1 do begin
   idx = where(lamtemps[*,i] gt 0)
   dum = dblarr(n_elements(outlam_temp))
   training[*,i] = df_rebin_spectra(lamtemps[idx,i], $
                                    temps[idx,i], $
                                    outlam_temp, $
                                    lanczos_a = input.lanczos_a, $
                                    innoise = noise[idx,i],$
                                    outnoise=dum)
   outnoise[*,i] = sqrt(dum)
   
   mask = dblarr(n_elements(lamtemps[*,i]))
   blah = where(finite(noise[*,i]) eq 0)
   mask[blah] = 1
   
   mask2 = df_rebin_spectra(lamtemps[idx,i],mask[idx],$
                            outlam_temp, $
                            lanczos_a = input.lanczos_a)
   
   dd = where(mask2 ne 0, count, complement = comp)
   if count gt 0 then begin
      training[dd,i] = 0.
      outnoise[dd,i] = max(outnoise[comp,i])
   endif 
   dd = where(finite(dum) eq 0, count,complement=comp)
   if count gt 0 then begin
      training[dd,i] = 0.
      outnoise[dd,i] = max(outnoise[comp,i])
   endif
endfor

;; now compute the components of each training spectrum
print,'Getting components of training set'
dectraining = df_get_spectra_components(training, rms = outnoise, $
                                        nscale = input.nscale, $
                                        nsigma = input.fdr, $
;                                        nsigma = 0.1, $
                                        regulparam = input.regulparam, $
                                        gen2starlet = input.gen2, $
                                        regulfilter = input.regul)

print,'Computing eigentemplates on dirty training spectra'
dirty_pca = df_get_eigenvector(dectraining.data_nobaseline, $
                               alleigen=alleigen, $
                               energpercent = input.energpercent, $
                               ntemp = input.ntemplates)

writefits,indir+'outs/dirty_pca2.fits',dirty_pca

return
end
