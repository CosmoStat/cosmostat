;; ;; ======= Plot test results ======

;; ;; nside = 512
;; ;; nscale = 60
;; ;; width = 140
;; ;; savename = 'test_toucan_decomp_'+strcompress(Nside,/remove_all)+'_'+strcompress(Nscale,/REMOVE_ALL)+'_'+strcompress(width,/REMOVE_ALL)+'.save'
;; ;; restore, savename, /verbose

;; ======= ***** PRINT OPTION ***** ====== ;;
print = 0
;; ======= ***** PRINT OPTION ***** ====== ;;

npix = nside2npix(nside)
lmax = 2l*nside
ell = findgen(lmax+1)*2+1

extra = ''
if (lin eq 1) then extra = extra + 'lin_'
if (mex eq 1) then extra = extra + 'mex_'

loadct, 39

;; ;; --------------------------- Average estimated integrate cl - no mask -----------------------------------------

;; intcl_ave_pn = fltarr(Nscale)
;; intcl_ave_pn_MC = fltarr(Nscale)
;; intcl_ave_block = fltarr(Nscale)
;; intcl_ave_MC = fltarr(Nscale)

;; for i=0, Nscale-1 do begin
;;    intcl_ave_pn(i) = mean(intcl_pn_ns[i,*])
;;    intcl_ave_pn_MC(i) = mean(intcl_pn_MC_ns[i,*])
;;    intcl_ave_block(i) = mean(intcl_block_ns[i,*])
;;    intcl_ave_MC(i) = mean(intcl_MC_ns[i,*])
;; endfor

;; imagename = '~/images/toucan_intcl_ns_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
;; if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
;; else  window,/free
;; tit = 'Compressed measurements - no mask'
;; xtit = 'wavelet filter'
;; ytit = 'averaged compressed measurements'
;; plot, var_needt_ideal, title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, yrange=[min([ var_needt_ideal,intcl_ave_pn,intcl_ave_pn_MC,intcl_ave_block,intcl_ave_MC]),max([ var_needt_ideal,intcl_ave_pn,intcl_ave_pn_MC,intcl_ave_block,intcl_ave_MC])] 
;; oplot, intcl_ave_pn, color=250, line=2
;; oplot, intcl_ave_pn_MC, color=50, line=2
;; oplot, intcl_ave_block, color=200, linestyle=2
;; oplot, intcl_ave_MC, color=100, linestyle=2
;; legend, ['True','est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[0,2,2,2,2],colors=[0,250,50,200,100],/right, charsize=1.3
;; if (print eq 1) then cps

;; ;; --------------------------- Average estimated integrate cl - with mask -----------------------------------------

;; intcl_ave_pn = fltarr(Nscale)
;; intcl_ave_pn_MC = fltarr(Nscale)
;; intcl_ave_block = fltarr(Nscale)
;; intcl_ave_MC = fltarr(Nscale)

;; for i=0, Nscale-1 do begin
;;    intcl_ave_pn(i) = mean(intcl_pn_mask[i,*])
;;    intcl_ave_pn_MC(i) = mean(intcl_pn_MC_mask[i,*])
;;    intcl_ave_block(i) = mean(intcl_block_mask[i,*])
;;    intcl_ave_MC(i) = mean(intcl_MC_mask[i,*])
;; endfor

;; imagename = '~/images/toucan_intcl_mask_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
;; if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
;; else  window,/free
;; tit = 'Compressed measurements - with mask'
;; xtit = 'wavelet filter'
;; ytit = 'averaged compressed measurements'
;; plot, var_needt_ideal, title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, yrange=[min([ var_needt_ideal,intcl_ave_pn,intcl_ave_pn_MC,intcl_ave_block,intcl_ave_MC]),max([ var_needt_ideal,intcl_ave_pn,intcl_ave_pn_MC,intcl_ave_block,intcl_ave_MC])] 
;; oplot, intcl_ave_pn, color=250, line=2
;; oplot, intcl_ave_pn_MC, color=50, line=2
;; oplot, intcl_ave_block, color=200, linestyle=2
;; oplot, intcl_ave_MC, color=100, linestyle=2
;; legend, ['True','est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[0,2,2,2,2],colors=[0,250,50,200,100],/right, charsize=1.3
;; if (print eq 1) then cps

;; ;; ----------------------------------- Log plot of normalized MSE - no mask --------------------------------------------

intcl_pn_mean_nmse = fltarr(Nscale)
intcl_pn_var_nmse = fltarr(Nscale)
intcl_pn_MC_mean_nmse = fltarr(Nscale)
intcl_pn_MC_var_nmse = fltarr(Nscale)
intcl_block_mean_nmse = fltarr(Nscale)
intcl_block_var_nmse = fltarr(Nscale)
intcl_MC_mean_nmse = fltarr(Nscale)
intcl_MC_var_nmse = fltarr(Nscale)

for i=0, Nscale-1 do begin
   err_i = ( ( var_needt[i,*]-intcl_pn_ns[i,*] )/var_needt[i,*] )^2.
   intcl_pn_mean_nmse[i] = mean(err_i)
   intcl_pn_var_nmse[i] = sqrt(mean( (err_i - intcl_pn_mean_nmse(i) )^2 ))

   err_i = ( ( var_needt[i,*]-intcl_pn_MC_ns[i,*] )/var_needt[i,*] )^2.
   intcl_pn_MC_mean_nmse[i] = mean(err_i)
   intcl_pn_MC_var_nmse[i] = sqrt(mean( (err_i - intcl_pn_MC_mean_nmse(i) )^2 ))

   err_i = ( ( var_needt[i,*]-intcl_block_ns[i,*] )/var_needt[i,*] )^2.
   intcl_block_mean_nmse[i] = mean(err_i)
   intcl_block_var_nmse[i] = sqrt(mean( (err_i - intcl_block_mean_nmse(i) )^2 ))

   err_i = ( ( var_needt[i,*]-intcl_MC_ns[i,*] )/var_needt[i,*] )^2.
   intcl_MC_mean_nmse[i] = mean(err_i)
   intcl_MC_var_nmse[i] = sqrt(mean( (err_i - intcl_MC_mean_nmse(i) )^2 ))
endfor

;; imagename = '~/images/toucan_intcl_ns_nmse_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
;; if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
;; else  window,/free
;; tit = 'Normalized MSE of compressed measurements - no mask'
;; xtit = 'wavelet filter'
;; ytit = 'nMSE'
;; ploterror, intcl_pn_mean_nmse, intcl_pn_var_nmse, /ylog, yrange=[min([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse]) , max([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse])],/nohat, errstyle=2, title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, /nodata ;, xrange=[0,1500], yrange=[0,300]
;; oploterror, intcl_pn_mean_nmse, intcl_pn_var_nmse, /nohat, errstyle=2, errcolor=250, color=250
;; oploterror, intcl_pn_MC_mean_nmse, intcl_pn_MC_var_nmse, /nohat, errstyle=2, errcolor=50, color=50
;; oploterror, intcl_block_mean_nmse, intcl_block_var_nmse, /nohat, errstyle=2, errcolor=200, color=200
;; oploterror, intcl_MC_mean_nmse, intcl_MC_var_nmse, /nohat, errstyle=2, errcolor=100, color=100
;; legend, ['est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[2,2,2,2],colors=[250,50,200,100],/right, /bottom, charsize=1.3
;; if (print eq 1) then cps

imagename = '~/images/toucan_intcl_ns_nmse_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Normalized MSE of compressed measurements - no mask'
xtit = 'wavelet filter'
ytit = 'nMSE'
plot, intcl_pn_mean_nmse, /ylog, yrange=[min([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse]) , max([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse])], title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, /nodata ;, xrange=[0,1500], yrange=[0,300]
oplot, intcl_pn_mean_nmse, color=250
oplot, intcl_pn_MC_mean_nmse, color=50
oplot, intcl_block_mean_nmse, color=200, line=2
oplot, intcl_MC_mean_nmse, color=100, line=2
legend, ['est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[0,0,2,2],colors=[250,50,200,100],/left, charsize=1.3
if (print eq 1) then cps

;;-------------------------------- Gain 

gain_block = fltarr(Nscale)
var_gain_block = fltarr(Nscale)
gain_MC = fltarr(Nscale)
var_gain_MC = fltarr(Nscale)

for i=0, Nscale-1 do begin
   err_pn_i = ( ( var_needt[i,*]-intcl_pn_ns[i,*] )/var_needt[i,*] )^2.
   err_i = ( ( var_needt[i,*]-intcl_block_ns[i,*] )/var_needt[i,*] )^2.
   gain_i = ( err_pn_i - err_i )/intcl_block_var_nmse[i]
   gain_block(i) = mean(gain_i)
   var_gain_block(i) = sqrt(mean( (gain_i - gain_block(i) )^2 ))

   err_pn_i = ( ( var_needt[i,*]-intcl_pn_MC_ns[i,*] )/var_needt[i,*] )^2.
   err_i = ( ( var_needt[i,*]-intcl_MC_ns[i,*] )/var_needt[i,*] )^2.
   gain_i = ( err_pn_i - err_i )/intcl_MC_var_nmse[i]
   gain_MC(i) = mean(gain_i)
   var_gain_MC(i) = sqrt(mean( (gain_i - gain_MC(i) )^2 ))
endfor

imagename = '~/images/toucan_intcl_ns_gain_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Gain on measurements estimation - no mask'
xtit = 'wavelet filter'
ytit = 'gain on compressed measurements' ; * cut_of_freq'
plot, gain_block, title = tit, xtitle=xtit, ytitle=ytit, yrange=[min([min(gain_block),min(gain_MC)]),max([max(gain_block),max(gain_MC)])], /nodata, charsize=1.3 ;, BACKGROUND = 255, COLOR = 0
oplot, gain_block, color=200, line=2
oplot, gain_MC, color=100, linestyle=2
oplot, fltarr(Nscale), color=250, line=2
legend, ['est. by patch','est. by MC'],linestyle=[2,2],colors=[200,100],/left, charsize=1.3
if (print eq 1) then cps

;; imagename = '~/images/toucan_intcl_ns_gain_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
;; if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
;; else  window,/free
;; tit = 'Gain on measurements estimation - no mask'
;; xtit = 'wavelet filter'
;; ytit = 'gain on compressed measurements' ; * cut_of_freq'
;; ploterror, gain_block, var_gain_block, xtitle=xtit, ytitle=ytit, yrange=[min([min(gain_block),min(gain_MC)]),max([max(gain_block),max(gain_MC)])], /nohat, /nodata, charsize=1.3 ;, BACKGROUND = 255, COLOR = 0
;; oploterror, gain_block, var_gain_block, /nohat, errstyle=2, errcolor=200, color=200, line=2
;; oploterror, gain_MC, var_gain_block, /nohat, errstyle=2, errcolor=100, color=100, linestyle=2
;; legend, ['est. by patch','est. by MC'],linestyle=[2,2],colors=[200,100],/left, charsize=1.3
;; if (print eq 1) then cps

;; ;; ----------------------------------- Log plot of normalized MSE - with mask --------------------------------------------

intcl_pn_mean_nmse = fltarr(Nscale)
intcl_pn_var_nmse = fltarr(Nscale)
intcl_pn_MC_mean_nmse = fltarr(Nscale)
intcl_pn_MC_var_nmse = fltarr(Nscale)
intcl_block_mean_nmse = fltarr(Nscale)
intcl_block_var_nmse = fltarr(Nscale)
intcl_MC_mean_nmse = fltarr(Nscale)
intcl_MC_var_nmse = fltarr(Nscale)

for i=0, Nscale-1 do begin
   err_i = ( ( var_needt[i,*]-intcl_pn_mask[i,*] )/var_needt[i,*] )^2.
   intcl_pn_mean_nmse[i] = mean(err_i)
   intcl_pn_var_nmse[i] = sqrt(mean( (err_i - intcl_pn_mean_nmse(i) )^2 ))

   err_i = ( ( var_needt[i,*]-intcl_pn_MC_mask[i,*] )/var_needt[i,*] )^2.
   intcl_pn_MC_mean_nmse[i] = mean(err_i)
   intcl_pn_MC_var_nmse[i] = sqrt(mean( (err_i - intcl_pn_MC_mean_nmse(i) )^2 ))

   err_i = ( ( var_needt[i,*]-intcl_block_mask[i,*] )/var_needt[i,*] )^2.
   intcl_block_mean_nmse[i] = mean(err_i)
   intcl_block_var_nmse[i] = sqrt(mean( (err_i - intcl_block_mean_nmse(i) )^2 ))

   err_i = ( ( var_needt[i,*]-intcl_MC_mask[i,*] )/var_needt[i,*] )^2.
   intcl_MC_mean_nmse[i] = mean(err_i)
   intcl_MC_var_nmse[i] = sqrt(mean( (err_i - intcl_MC_mean_nmse(i) )^2 ))
endfor

;; imagename = '~/images/toucan_intcl_mask_nmse_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
;; if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
;; else  window,/free
;; tit = 'Normalized MSE of compressed measurements - with mask'
;; xtit = 'wavelet filter'
;; ytit = 'nMSE'
;; ploterror, intcl_pn_mean_nmse, intcl_pn_var_nmse, /ylog, yrange=[min([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse]) , max([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse])],/nohat, errstyle=2, title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, /nodata ;, xrange=[0,1500], yrange=[0,300]
;; oploterror, intcl_pn_mean_nmse, intcl_pn_var_nmse, /nohat, errstyle=2, errcolor=250, color=250
;; oploterror, intcl_pn_MC_mean_nmse, intcl_pn_MC_var_nmse, /nohat, errstyle=2, errcolor=50, color=50
;; oploterror, intcl_block_mean_nmse, intcl_block_var_nmse, /nohat, errstyle=2, errcolor=200, color=200
;; oploterror, intcl_MC_mean_nmse, intcl_MC_var_nmse, /nohat, errstyle=2, errcolor=100, color=100
;; legend, ['est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[2,2,2,2],colors=[250,50,200,100],/right, charsize=1.3
;; if (print eq 1) then cps

imagename = '~/images/toucan_intcl_mask_nmse_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Normalized MSE of compressed measurements - with mask'
xtit = 'wavelet filter'
ytit = 'nMSE'
plot, intcl_pn_mean_nmse, /ylog, yrange=[min([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse]) , max([intcl_pn_mean_nmse,intcl_pn_MC_mean_nmse,intcl_block_mean_nmse,intcl_MC_mean_nmse])], title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, /nodata ;, xrange=[0,1500], yrange=[0,300]
oplot, intcl_pn_mean_nmse, color=250
oplot, intcl_pn_MC_mean_nmse, color=50
oplot, intcl_block_mean_nmse, color=200, line=2
oplot, intcl_MC_mean_nmse, color=100, line=2
legend, ['est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[0,0,2,2],colors=[250,50,200,100],/left, charsize=1.3
if (print eq 1) then cps

;;-------------------------------- Gain 

gain_block = fltarr(Nscale)
var_gain_block = fltarr(Nscale)
gain_MC = fltarr(Nscale)
var_gain_MC = fltarr(Nscale)

for i=0, Nscale-1 do begin
   err_pn_i = ( ( var_needt[i,*]-intcl_pn_mask[i,*] )/var_needt[i,*] )^2.
   err_i = ( ( var_needt[i,*]-intcl_block_mask[i,*] )/var_needt[i,*] )^2.
   gain_i = ( err_pn_i - err_i )/intcl_block_var_nmse[i]
   gain_block(i) = mean(gain_i)
   var_gain_block(i) = sqrt(mean( (gain_i - gain_block(i) )^2 ))

   err_pn_i = ( ( var_needt[i,*]-intcl_pn_MC_mask[i,*] )/var_needt[i,*] )^2.
   err_i = ( ( var_needt[i,*]-intcl_MC_mask[i,*] )/var_needt[i,*] )^2.
   gain_i = ( err_pn_i - err_i )/intcl_MC_var_nmse[i]
   gain_MC(i) = mean(gain_i)
   var_gain_MC(i) = sqrt(mean( (gain_i - gain_MC(i) )^2 ))
endfor

imagename = '~/images/toucan_intcl_mask_gain_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Gain on measurements estimation - with mask'
xtit = 'wavelet filter'
ytit = 'gain on compressed measurements' ; * cut_of_freq'
plot, gain_block, title = tit, xtitle=xtit, ytitle=ytit, yrange=[min([min(gain_block),min(gain_MC)]),max([max(gain_block),max(gain_MC)])], /nodata, charsize=1.3 ;, BACKGROUND = 255, COLOR = 0
oplot, gain_block, color=200, line=2
oplot, gain_MC, color=100, linestyle=2
oplot, fltarr(Nscale), color=250, line=2
legend, ['est. by patch','est. by MC'],linestyle=[2,2],colors=[200,100],/left, charsize=1.3
if (print eq 1) then cps

;; imagename = '~/images/toucan_intcl_mask_gain_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
;; if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
;; else  window,/free
;; tit = 'Gain on measurements estimation - with mask'
;; xtit = 'wavelet filter'
;; ytit = 'gain on compressed measurements' ; * cut_of_freq'
;; ploterror, gain_block, var_gain_block, xtitle=xtit, ytitle=ytit, yrange=[min([min(gain_block),min(gain_MC)]),max([max(gain_block),max(gain_MC)])], /nohat, /nodata, charsize=1.3 ;, BACKGROUND = 255, COLOR = 0
;; oploterror, gain_block, var_gain_block, /nohat, errstyle=2, errcolor=200, color=200, line=2
;; oploterror, gain_MC, var_gain_block, /nohat, errstyle=2, errcolor=100, color=100, linestyle=2
;; legend, ['est. by patch','est. by MC'],linestyle=[2,2],colors=[200,100],/left, charsize=1.3
;; if (print eq 1) then cps

end
