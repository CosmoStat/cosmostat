;; ;; ======= Plot test results ======

;; ;; nside = 512
;; ;; nscale = 60
;; ;; width = 140
;; ;; savename = 'test_toucan_decomp_'+strcompress(Nside,/remove_all)+'_'+strcompress(Nscale,/REMOVE_ALL)+'_'+strcompress(width,/REMOVE_ALL)+'.save'
;; ;; restore, savename, /verbose

;; ======= ***** PRINT OPTION ***** ====== ;;
print = 1
;; ======= ***** PRINT OPTION ***** ====== ;;

npix = nside2npix(nside)
lmax = 2l*nside
ell = findgen(lmax+1)*2+1

extra = ''
if (lin eq 1) then extra = extra + 'lin_'
if (mex eq 1) then extra = extra + 'mex_'

filename = '~/images/toucan_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'_'+strcompress(niter,/remove_all)

;; loadct, 15
loadct, 39

;; ;; ----------------------------------- Average estimated cl - no mask --------------------------------------------

tousi_cl_ave = fltarr(lmax+1)
tousi_cl_MC_ave = fltarr(lmax+1)
toucan_cl_block_ave = fltarr(lmax+1)
toucan_cl_MC_ave = fltarr(lmax+1)

for i=2, lmax do begin
    tousi_cl_ave(i) = mean(tousi_cl_pn_ns[i,*])
    tousi_cl_MC_ave(i) = mean(tousi_cl_pn_MC_ns[i,*])
    toucan_cl_block_ave(i) = mean(toucan_cl_block_ns[i,*])
    toucan_cl_MC_ave(i) = mean(toucan_cl_MC_ns[i,*])
endfor

imagename = filename + '/toucan_cl_ns_' + extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width, /remove_all)+'.eps'
if (print eq 1) then ops, file=imagename, /COLOR,/landscape $
else window,/free
plotcl, cl
oplotcl, tousi_cl_ave, color=250, line=2
oplotcl, tousi_cl_MC_ave, color=50, line=2
oplotcl, toucan_cl_block_ave, color=200, line=2
oplotcl, toucan_cl_MC_ave, color=100, line=2
legend, ['True Cl', 'est. from tousi', 'est. from MC tousi', 'est. by patch', 'est. by MC'], linestyle=[0,2,2,2,2], colors=[0,250,50,200,100],/right,charsize=1.3

;; ;; ----------------------------------- Average estimated cl - with mask --------------------------------------------

tousi_cl_ave = fltarr(lmax+1)
tousi_cl_MC_ave = fltarr(lmax+1)
toucan_cl_block_ave = fltarr(lmax+1)
toucan_cl_MC_ave = fltarr(lmax+1)

for i =2, lmax do begin
    tousi_cl_ave(i) = mean(tousi_cl_pn_mask[i,*])
    tousi_cl_MC_ave(i) = mean(tousi_cl_pn_MC_mask[i,*])
    toucan_cl_block_ave(i) = mean(toucan_cl_block_mask[i,*])
    toucan_cl_MC_ave(i) = mean(toucan_cl_MC_mask[i,*])
endfor

imagename = filename + '/toucan_cl_mask_' + extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops, file=imagename, /COLOR,/landscape $
else window,/free
plotcl, cl
oplotcl, tousi_cl_ave, color=250, line=2
oplotcl, tousi_cl_MC_ave, color=50, line=2
oplotcl, toucan_cl_block_ave, color=200, line=2
oplotcl, toucan_cl_MC_ave, color=100, line=2
legend, ['True Cl', 'est. from tousi', 'est. from MC tousi', 'est. by patch', 'est. by MC'], linestyle=[0,2,2,2,2], colors=[0, 250, 50, 200, 100],/right,charsize=1.3

;; ;; ----------------------------------- Log plot of normalized MSE - no mask --------------------------------------------

tousi_cl_mean_nmse = fltarr(lmax+1)
tousi_cl_var_nmse = fltarr(lmax+1)
tousi_cl_MC_mean_nmse = fltarr(lmax+1)
tousi_cl_MC_var_nmse = fltarr(lmax+1)
toucan_cl_block_mean_nmse = fltarr(lmax+1)
toucan_cl_block_var_nmse = fltarr(lmax+1)
toucan_cl_MC_mean_nmse = fltarr(lmax+1)
toucan_cl_MC_var_nmse = fltarr(lmax+1)

for i=2, lmax do begin
   err_i = ( ( cl[i]-tousi_cl_pn_ns[i,*] )/cl[i] )^2.
   tousi_cl_mean_nmse[i] = mean(err_i)
   tousi_cl_var_nmse[i] = sqrt(mean( (err_i - tousi_cl_mean_nmse(i) )^2 ))

   err_i = ( ( cl[i]-tousi_cl_pn_MC_ns[i,*] )/cl[i] )^2.
   tousi_cl_MC_mean_nmse[i] = mean(err_i)
   tousi_cl_MC_var_nmse[i] = sqrt(mean( (err_i - tousi_cl_MC_mean_nmse(i) )^2 ))

   err_i = ( ( cl[i]-toucan_cl_block_ns[i,*] )/cl[i] )^2.
   toucan_cl_block_mean_nmse[i] = mean(err_i)
   toucan_cl_block_var_nmse[i] = sqrt(mean( (err_i - toucan_cl_block_mean_nmse(i) )^2 ))

   err_i = ( ( cl[i]-toucan_cl_MC_ns[i,*] )/cl[i] )^2.
   toucan_cl_MC_mean_nmse[i] = mean(err_i)
   toucan_cl_MC_var_nmse[i] = sqrt(mean( (err_i - toucan_cl_MC_mean_nmse(i) )^2 ))
endfor

imagename = filename+'/toucan_cl_ns_nmse_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Normalized MSE of reconstructed Cl - no mask'
xtit = 'multipole l'
ytit = 'nMSE'
plot, tousi_cl_mean_nmse, /ylog, yrange=[min([tousi_cl_mean_nmse[2:*],tousi_cl_MC_mean_nmse[2:*],toucan_cl_block_mean_nmse[2:*],toucan_cl_MC_mean_nmse[2:*]]) , max([tousi_cl_mean_nmse[2:*],tousi_cl_MC_mean_nmse[2:*],toucan_cl_block_mean_nmse[2:*],toucan_cl_MC_mean_nmse[2:*]])], title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, /nodata ;, xrange=[0,1500], yrange=[0,300]
oplot, tousi_cl_mean_nmse, color=250 ;210
oplot, tousi_cl_MC_mean_nmse, color=50
oplot, toucan_cl_block_mean_nmse, color=200, line=2
oplot, toucan_cl_MC_mean_nmse, color=100, line=2
legend, ['est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[0,0,2,2],colors=[250,50,200,100],/left, charsize=1.3
if (print eq 1) then cps

;;-------------------------------- Gain 

gain_block = fltarr(lmax+1)
var_gain_block = fltarr(lmax+1)
gain_MC = fltarr(lmax+1)
var_gain_MC = fltarr(lmax+1)

for i=2, lmax do begin
   err_pn_i = ( ( cl[i]-tousi_cl_pn_ns[i,*] )/cl[i] )^2.
   err_i = ( ( cl[i]-toucan_cl_block_ns[i,*] )/cl[i] )^2.
   gain_i = ( err_pn_i - err_i )/toucan_cl_block_var_nmse[i]
   gain_block(i) = mean(gain_i)
   var_gain_block(i) = sqrt(mean( (gain_i - gain_block(i) )^2 ))

   err_pn_i = ( ( cl[i]-tousi_cl_pn_MC_ns[i,*] )/cl[i] )^2.
   err_i = ( ( cl[i]-toucan_cl_MC_ns[i,*] )/cl[i] )^2.
   gain_i = ( err_pn_i - err_i )/toucan_cl_MC_var_nmse[i]
   gain_MC(i) = mean(gain_i)
   var_gain_MC(i) = sqrt(mean( (gain_i - gain_MC(i) )^2 ))
endfor

imagename = filename+'/toucan_cl_ns_gain_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Gain on Cl estimation - no mask'
xtit = 'multipole l'
ytit = 'gain on Cl' ; * cut_of_freq'
plot, gain_block, title = tit, xtitle=xtit, ytitle=ytit, yrange=[min([gain_block, gain_MC]),max([gain_block, gain_MC])], /nodata, charsize=1.3 ;, BACKGROUND = 255, COLOR = 0
oplot, gain_block, color=200, line=2
oplot, gain_MC, color=100, linestyle=2
oplot, fltarr(lmax+1), color=250, line=2
legend, ['est. by patch','est. by MC'],linestyle=[2,2],colors=[200,100],/left, charsize=1.3
if (print eq 1) then cps

;; ;; ----------------------------------- Log plot of normalized MSE - no mask --------------------------------------------

tousi_cl_mean_nmse = fltarr(lmax+1)
tousi_cl_var_nmse = fltarr(lmax+1)
tousi_cl_MC_mean_nmse = fltarr(lmax+1)
tousi_cl_MC_var_nmse = fltarr(lmax+1)
toucan_cl_block_mean_nmse = fltarr(lmax+1)
toucan_cl_block_var_nmse = fltarr(lmax+1)
toucan_cl_MC_mean_nmse = fltarr(lmax+1)
toucan_cl_MC_var_nmse = fltarr(lmax+1)

for i=2, lmax do begin
   err_i = ( ( cl[i]-tousi_cl_pn_mask[i,*] )/cl[i] )^2.
   tousi_cl_mean_nmse[i] = mean(err_i)
   tousi_cl_var_nmse[i] = sqrt(mean( (err_i - tousi_cl_mean_nmse(i) )^2 ))

   err_i = ( ( cl[i]-tousi_cl_pn_MC_mask[i,*] )/cl[i] )^2.
   tousi_cl_MC_mean_nmse[i] = mean(err_i)
   tousi_cl_MC_var_nmse[i] = sqrt(mean( (err_i - tousi_cl_MC_mean_nmse(i) )^2 ))

   err_i = ( ( cl[i]-toucan_cl_block_mask[i,*] )/cl[i] )^2.
   toucan_cl_block_mean_nmse[i] = mean(err_i)
   toucan_cl_block_var_nmse[i] = sqrt(mean( (err_i - toucan_cl_block_mean_nmse(i) )^2 ))

   err_i = ( ( cl[i]-toucan_cl_MC_mask[i,*] )/cl[i] )^2.
   toucan_cl_MC_mean_nmse[i] = mean(err_i)
   toucan_cl_MC_var_nmse[i] = sqrt(mean( (err_i - toucan_cl_MC_mean_nmse(i) )^2 ))
endfor

imagename = filename+'/toucan_cl_mask_nmse_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Normalized MSE of reconstructed Cl - with mask'
xtit = 'multipole l'
ytit = 'nMSE'
plot, tousi_cl_mean_nmse, /ylog, yrange=[min([tousi_cl_mean_nmse[2:*],tousi_cl_MC_mean_nmse[2:*],toucan_cl_block_mean_nmse[2:*],toucan_cl_MC_mean_nmse[2:*]]) , max([tousi_cl_mean_nmse[2:*],tousi_cl_MC_mean_nmse[2:*],toucan_cl_block_mean_nmse[2:*],toucan_cl_MC_mean_nmse[2:*]])], title=tit, xtitle=xtit, ytitle=ytit, charsize=1.3, /nodata ;, xrange=[0,1500], yrange=[0,300]
oplot, tousi_cl_mean_nmse, color=250
oplot, tousi_cl_MC_mean_nmse, color=50
oplot, toucan_cl_block_mean_nmse, color=200, line=2
oplot, toucan_cl_MC_mean_nmse, color=100, line=2
legend, ['est. from noise powspec','est. from MC noise powspec','est. by patch', 'est. by MC'],linestyle=[0,0,2,2],colors=[250,50,200,100],/left, charsize=1.3
if (print eq 1) then cps

;;-------------------------------- Gain 

gain_block = fltarr(lmax+1)
var_gain_block = fltarr(lmax+1)
gain_MC = fltarr(lmax+1)
var_gain_MC = fltarr(lmax+1)

for i=2, lmax do begin
   err_pn_i = ( ( cl[i]-tousi_cl_pn_mask[i,*] )/cl[i] )^2.
   err_i = ( ( cl[i]-toucan_cl_block_mask[i,*] )/cl[i] )^2.
   gain_i = ( err_pn_i - err_i )/toucan_cl_block_var_nmse[i]
   gain_block(i) = mean(gain_i)
   var_gain_block(i) = sqrt(mean( (gain_i - gain_block(i) )^2 ))

;;    err_pn_i = ( ( cl[i]-tousi_cl_pn_MC_mask[i,*] )/cl[i] )^2.
   err_i = ( ( cl[i]-toucan_cl_MC_mask[i,*] )/cl[i] )^2.
   gain_i = ( err_pn_i - err_i )/toucan_cl_MC_var_nmse[i]
   gain_MC(i) = mean(gain_i)
   var_gain_MC(i) = sqrt(mean( (gain_i - gain_MC(i) )^2 ))
endfor

imagename = filename+'/toucan_cl_mask_gain_nobars_'+extra+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'_'+strcompress(width,/remove_all)+'.eps'
if (print eq 1) then ops,file=imagename, /COLOR,/landscape $
else  window,/free
tit = 'Gain on Cl estimation - with mask'
xtit = 'multipole l'
ytit = 'gain on Cl' ; * cut_of_freq'
plot, gain_block, title = tit, xtitle=xtit, ytitle=ytit, yrange=[min([gain_block, gain_MC]),max([gain_block, gain_MC])], /nodata, charsize=1.3 ;, BACKGROUND = 255, COLOR = 0
oplot, gain_block, color=200, line=2
oplot, gain_MC, color=100, linestyle=2
oplot, fltarr(lmax+1), color=250, line=2
legend, ['est. by patch','est. by MC'],linestyle=[2,2],colors=[200,100],/left, charsize=1.3
if (print eq 1) then cps

end
