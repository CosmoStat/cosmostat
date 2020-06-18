pro wigglez_pipeline
loadct,5
close,101
input = df_input_params_wigglez()
ndata = input.ndata;225415L

;; read in redshift catalogue

readcol,input.indir+'WiggleZfinalV04_99.dat',f='x,x,x,d,d,i',z, zerr, q,skipline=5

;; generate the log pixel scheme required for the data and the templates

if input.rebin_spectra eq 1 or input.rebin_training eq 1 then begin
   print,'generating pixels'
   outlam_temp = df_gen_log_pix(input.training_lmin, input.training_lmax, $
                                delta=input.lstep)
   
   idx = where(outlam_temp ge input.data_lmin and outlam_temp le $
               input.data_lmax)
   outlam_data = outlam_temp[idx]
   input.shiftpar = round((alog10(double(max(outlam_temp))) - $
                          alog10(double(max(outlam_data))))/input.lstep)
   input.training_lmin = min(outlam_temp)
   input.training_lmax = max(outlam_temp)
   input.data_lmin = min(outlam_data)
   input.data_lmax = max(outlam_data)
endif

;; Read in and rebin training set data, if required
if input.rebin_training eq 1 then begin
   
   trainingspec = mrdfits(input.indir+'training_set.fits',0)
   traininglam = mrdfits(input.indir+'training_set.fits',1)
   ;; uncomment the line below if using noisy training data

;   trainingvar = mrdfits(input.indir+'training_set.fits',2)
 
   sz_training = size(trainingspec)
   ntrain = sz_training[2]
   
   training = dblarr(n_elements(outlam_temp),ntrain)
   for i = 0, ntrain - 1 do begin
      print, 'rebinning spectrum', i

      ;; the commented lines below should be used if noisy data are
      ;; selected as training data. 

 ;     outnoise = dblarr(n_elements(outlam_temp))
 ;     mask = dblarr(n_elements(sp.spectrum))
 ;     blah = where(finite(trainingvar) eq 0)
 ;     mask[blah] = 1
      idx = where(traininglam[*,i] gt 0)
      training[*,i] = df_rebin_spectra(traininglam[idx,i], $
                                       trainingspec[idx,i], $
                                       outlam_temp, $
                                       lanczos_a=input.lanczos_a);, $
                                       ;innoise=trainingvar[idx,i],$
                                       ;outnoise=outnoise)
      ;mask2 = df_rebin_spectra(traininglam[idx,i], mask[idx], outlam_temp, $
      ;                         lanczos_a = input.lanczos_a)
      
     ; maskout = where(mask2 ne 0, count,complement=comp)
     ; if count gt 0 then begin
     ;    training[maskout,i] = 0.
     ;    outnoise[maskout,i] = max(outnoise[comp])
     ; endif
     ; rms = sqrt(outnoise)
     ; maskout2 = where(finite(rms) eq 0,count,complement=maskin)
     ; if count gt 0 then begin
     ;    rms[maskout2] = max(rms[maskin])
     ;    spec[maskout2] = 0.
      ; endif
      
      
   endfor
   
   ;; write these to a file so that this stage can be skipped if
   ;; desired
   
   if input.outputtraining eq 1 then begin
      writefits, input.outdir+input.trainingcat, $
                 training
      writefits, input.outdir+input.trainingcat, outlam_temp, /append

      ;; uncomment line below if using noisy training data

     ; writefits, input.outdir+input.trainingcat, rms, /append
   endif
   
endif else begin
   ;; read in the training set from a file if rebinning is not needed
   training = mrdfits(input.outdir + input.trainingcat,0)
   outlam_temp = mrdfits(input.outdir + input.trainingcat,1)
   
   ;; uncomment line below if using noisy training data

   ;rms = mrdfits(input.outdir + input.trainingcat,2)
endelse

;; compute eigentemplates, if required
if input.computeeigentemplates eq 1 then begin
   print,'removing baseline from training catalogue'
   training_nobaseline = df_del_baseline(training, nscale = input.nscale)

   ;; uncomment line below if using noisy training data

   ;training_comp = df_get_spectra_components(training, rms = rms, $
   ;                                          nsigma = input.fdr, $
   ;                                          gen2starlet = input.gen2, $
   ;                                          regulfilter=input.regul, $
   ;                                          regulparam = input.regulparam)
   ; training_nobaseline = training_comp.data_nobaseline
   


   print,'computing eigentemplates'
   template = df_get_eigenvector(training_nobaseline, alleigen = alleigen, $
                                 energpercent = input.energpercent, $
                                 ntemp = input.ntemplates)
   
   if input.outputtemplates eq 1 then writefits, input.outdir+$
      input.tempfilename, template
   if input.outputalleigen eq 1 then writefits, input.outdir + $
      input.alleigenname, alleigen
endif else begin
   ;; or read them in if already generated
   if strmatch(input.templatecat,'*.fits') or  $
   strmatch(input.templatecat,'*.fits.gz') then $
      template = readfits(input.outdir + input.templatecat) else $
         message, 'No templates specified!'

   ;; you may need to uncomment the line below if you are using noisy
   ;; data as your templates

   ;template = template[*,0:19]

endelse

skip:
;; read in line information for line matching/flagging of spectra
readcol,input.linefile, f='d,a', linelam, linename
linepix = round(double(alog10(linelam) - alog10(input.data_lmin))/double(input.lstep))
nlinestomatch = n_elements(input.lines)
if strmatch(input.matchtype,'OR') then numcrit = 1
if strmatch(input.matchtype,'AND') then numcrit = nlinestomatch

;; initialise output arrays and open file to write
estred = dblarr(ndata)
estred2 = estred
flag = intarr(ndata)
openw,101,input.outdir+input.RedshiftFileName
snr = dblarr(ndata)
for i = 0L, long(ndata) - 1 do begin
;; loop over spectra
   if (i+1) lt 10 then prefix = '00000'
   if (i+1) ge 10 and (i+1) lt 100 then prefix = '0000'
   if (i+1) ge 100 and (i+1) lt 1000 then prefix = '000'
   if (i+1) ge 1000 and (i+1) lt 10000 then prefix = '00'
   if (i+1) ge 10000 and (i+1) lt 100000 then prefix = '0'
   if (i+1) ge 100000 then prefix = ''
   
   if input.rebin_spectra eq 1 then begin

      filename = input.indir + 'spectra/wig'+prefix+strtrim(i+1,2)+'.fits'

      print,'Rebinning spectrum ','wig'+prefix+strtrim(i+1,2)

      ;; read in spectrum & noise
      sp = read_wigglez_data(filename)
      
      ;; rebin spectrum and error curve
      outnoise = dblarr(n_elements(outlam_data))
      df_mask_lines, sp.lambda, [6362.,5577.,5893.,6300.], mask, tol=4
      idx = where(finite(sp.variance) eq 0,count)
      if count gt 0 then mask[idx] = 1
      spec = df_rebin_spectra(sp.lambda, sp.spectrum*(1-mask), outlam_data, $
                              innoise = sp.variance, outnoise = outnoise, $
                              lanczos_a = input.lanczos_a)
      mask2 = df_rebin_spectra(sp.lambda, mask, outlam_data, $
                               lanczos_a = input.lanczos_a)
      if input.OutputRebinnedSpectra eq 1 then begin
         writefits,input.outdir+'rebinned_spectra/'+input.RebSpecName+prefix+$
                   strtrim(i+1,2)+'.fits',spec
         writefits,input.outdir+'rebinned_spectra/'+input.RebSpecName+prefix+$
                   strtrim(i+1,2)+'.fits',outnoise,/append
         writefits,input.outdir+'rebinned_spectra/'+input.RebSpecName+prefix+$
                strtrim(i+1,2)+'.fits',outlam_data,/append
      endif 

      maskout = where(mask2 ne 0, count,complement=comp)
      if count gt 0 then begin
         spec[maskout] = 0.
         outnoise[maskout] = max(outnoise[comp])
      endif
      
 
   endif else begin
readspec:      
      filename = input.outdir + 'rebinned_spectra/'+input.RebSpecName+prefix+$
                 strtrim(i+1,2)+'.fits'
      print,'Reading spectrum ', input.RebSpecName+prefix+strtrim(i+1,2)+'.fits'
      outlam_data = mrdfits(filename,2)
      spec = mrdfits(filename,0)
;      spec = {spectrum:dum, lambda:outlam_data}
      outnoise = mrdfits(filename,1)
   endelse
         
   rms = sqrt(outnoise)
   maskout2 = where(finite(rms) eq 0,count,complement=maskin)
   if count gt 0 then begin
      rms[maskout2] = max(rms[maskin])
      spec[maskout2] = 0.
   endif
   idx = where(outlam_data ge 5600 and outlam_data le 6760)
   snr[i] = mean(spec[idx])/mean(rms[idx])
   goto,skip2

   ;; compute components of spectrum
   if input.computecomponents eq 1 then begin
      print,'computing components'
      
      decspectra = df_get_spectra_components(spec, rms = rms, $
                                             nsigma = input.fdr, $
                                             gen2starlet = input.gen2, $
                                             regulfilter=input.regul, $
                                             regulparam = input.regulparam)

      ;; estimate redshift
      
      if input.outputcomponents eq 1 then begin
         
         writefits,input.outdir+'rebinned_spectra/'+$
                   input.compfilename+prefix+$
                   strtrim(i+1,2)+'.fits', decspectra.baseline
         writefits,input.outdir+'rebinned_spectra/'+$
                   input.compfilename+prefix+$
                   strtrim(i+1,2)+'.fits', decspectra.posline,/append
         writefits,input.outdir+'rebinned_spectra/'+$
                   input.compfilename+prefix+$
                   strtrim(i+1,2)+'.fits', decspectra.negline,/append
         writefits,input.outdir+'rebinned_spectra/'+$
                   input.compfilename+prefix+$
                   strtrim(i+1,2)+'.fits', decspectra.noise,/append
         writefits,input.outdir+'rebinned_spectra/'+$
                   input.compfilename+prefix+$
                   strtrim(i+1,2)+'.fits', decspectra.data_nobaseline,/append
         writefits,input.outdir+'rebinned_spectra/'+$
                   input.compfilename+prefix+$
                   strtrim(i+1,2)+'.fits', decspectra.data_clean,/append
      endif
   endif else begin
      baseline = mrdfits(input.outdir+'rebinned_spectra/'+$
                         input.compfilename+prefix+strtrim(i+1,2)+$
                         '.fits',0)
      posline = mrdfits(input.outdir+'rebinned_spectra/'+$
                        input.compfilename+prefix+strtrim(i+1,2)+$
                        '.fits',1)
      negline = mrdfits(input.outdir+'rebinned_spectra/'+$
                        input.compfilename+prefix+strtrim(i+1,2)+$
                        '.fits',2)
      noise = mrdfits(input.outdir+'rebinned_spectra/'+$
                      input.compfilename+prefix+strtrim(i+1,2)+$
                      '.fits',3)
      data_nobaseline = mrdfits(input.outdir+'rebinned_spectra/'+$
                                input.compfilename+prefix+strtrim(i+1,2)+$
                                '.fits',4)
      data_clean = mrdfits(input.outdir+'rebinned_spectra/'+$
                           input.compfilename+prefix+strtrim(i+1,2)+$
                           '.fits',5)
      
      decspectra = {baseline:baseline, posline:posline, negline:negline,$
                    noise:noise, data_nobaseline:data_nobaseline, $
                    data_clean:data_clean}
   endelse
   
   print,'estimating redshift'
   chisquare = dblarr(n_elements(outlam_temp))
   if input.estimredshiftfromnoisydata eq 1 then $
      estred[i] = df_get_redshift(decspectra.data_nobaseline,template,$
                                  lstep=input.lstep, chisquare=chisquare,$
                                  shiftpar = input.shiftpar) $
   else $
      estred[i] = df_get_redshift(decspectra.data_clean, template, $
                                  chisquare = chisquare, $
                                  lstep = input.lstep, shiftpar = input.shiftpar)

   ;; decide whether to keep spectrum/flag for cleaning
   
   print,'Performing line matching'
   delpix = round(double(alog10(1.+(estred[i]))/input.lstep))
   peaks = df_get_peaks(decspectra.negline, decspectra.posline, noise=rms, $
                        nsigma = input.nsigmapeak, loc = loc)
   print,'Number of features:', (peaks.allpeaks)
   matches = df_linematch(loc, linepix, linename, delpix, $
                          thresh = input.matchthresh)
   nametest = (matches.name)[0]
   nmatches = intarr(nlinestomatch)
   if strmatch(strtrim(nametest,2),'0') then nmatches[*] = 0 else begin
      for match = 0, n_elements(matches.distance) - 1 do $
         for line = 0, nlinestomatch - 1 do begin
         if strmatch(matches.name[match],input.lines[line]) then $
            nmatches[line] = 1
      endfor
   endelse

   if total(nmatches) ge numcrit then flag[i] = 1 else flag[i] = 0

   ;; output to catalogue

   printf,101,format = '(f8.6,1x, i1)', estred[i], flag[i]

   skip2:
endfor

close,101
plothist,snr,bin=1
read,idum


;; plot DF redshift vs human redshift

window,0
plot,z,estred,psym=4, xr=[0,2],yr=[0,2]
oplot,[0,2],[0,2],color=100

idx=where(abs(estred-z) le 0.05*(1+z),ct)
print,'Success rate (noisy temps) before cleaning = ', double(ct)/double(ndata)*100.

id2 = where(flag eq 1,count)
oplot,z[id2],estred[id2],psym=4,color=100
dum = where(abs(estred[id2]-z[id2]) le 0.05*(1+z[id2]),ct2)
print,'Success rate (noisy temps) after cleaning = ', double(ct2)/double(count)*100.
print,'Retention rate (noisy temps) after cleaning = ', double(count)/double(ndata)*100.

stop

end
