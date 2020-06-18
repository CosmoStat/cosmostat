;; .compile mrs_get_cl_theo_powspec
;; .compile /functions/master_anafast_h.pro

;; goto, PLOTS
;; goto, Here
;; goto, Begining
;; goto, MSE
;; goto, bestit

;; ============ Number of Simulations =============
Ntest = 100

;; ============ Parameters of Toucnan =============
width = 140
nside = 512 ; 512; 256
Nscale = 60
Firstl = 2 ;800
lmax = 2l*nside
Niter = 180
lin = 1
mex = 0
tousil = 150

;; ============ Parameters of cmb =============
npix = nside2npix(nside)
ell = findgen(lmax+1)*2.+1.
fsky = 1

;; ============ savename of experiments =============
savename = 'test_toucan_decomp_'
if (lin eq 1) then savename = savename + 'lin_'
if (mex eq 1) then savename = savename + 'mex_'
if (lin eq 1) then savename = savename + strcompress(Nside,/remove_all)+'_'+strcompress(Nscale,/REMOVE_ALL)+'_'+strcompress(width,/REMOVE_ALL)+'.save' else savename = 'test_toucan_decomp_'+strcompress(Nside,/remove_all)+'_'+strcompress(Nscale,/REMOVE_ALL)+'_'+strcompress(width,/REMOVE_ALL)+'.save'

stop
;; ============ Generate the wavelet filters =============
bkj = mrs_toucan_wt_filters(Nscale, lmax, Firstl=Firstl, filtwidth=width, Bmat=Bmat, lin=lin) 
print, '--- wavelet filters generated.'

;;  ================== Put in 15 arcmin  ==================
h = getbeam(FWHM=15, lmax=lmax)

;; ============ Non-stationary noise =============
sigmaNoise = rims('images/WMAP7_rms_map_94GHz.fits')
if nside ne npix2nside(N_elements(sigmaNoise)) then sigmaNoise = mrs_resize(sigmaNoise, nside=nside)
sigfact = 533 ;sigma(cmb)/sigma(sigmaNoise)/25.
sigmaNoise = sigfact*sigmaNoise
;; pn = mrs_powspec(randomn(seed, npix)*sigmaNoise) *** -> do a new realization each time instead
print, '--- non-stationary noise generated.'

;; ================== Mask ====================
if not keyword_set(mask) then mask = rims('images/mask64.fits')
if nside ne npix2nside(N_elements(mask)) then  mask = mrs_resize(mask,nside=nside)
mask = mask[*]
nmask = total(total(mask))
fsky = nmask/npix
print, '--- mask restored.'

;; ================== restore estimated noise variances  ==================
;; ------- method by patches *** -> do a new realization each time instead
;; if (lin eq 1) then blockname = 'block_varnoise_lin_'+strcompress(nside,/remove_all)+'_'+strcompress(Nscale,/remove_all)+'_'+strcompress(width,/REMOVE_ALL)+'.save' else blockname = 'block_varnoise_'+strcompress(nside,/remove_all)+'_'+strcompress(Nscale,/remove_all)+'_'+strcompress(width,/REMOVE_ALL)+'.save'
;; filepresent = file_search(blockname,count=c)
;; if c eq 0 then begin 
;;    sigmaNoise2 = rims('images/WMAP7_rms_map_94GHz.fits')
;;    if npix ne N_elements(sigmaNoise2) then sigmaNoise2 = mrs_resize(sigmaNoise2, nside=nside)
;;    n = randomn(seed,npix)
;;    NoiseReal = n*SigmaNoise2
;;    var_noise =  mrs_var_noise_patch(NoiseReal, bkj)
;;    save, var_noise, filename=blockname 
;; endif else restore, blockname

;; var_noise_block = var_noise
;; var_noiseBlock = var_noise_block*(sigfact)^2. 
;; print, '--- noise variance estimated by patch method.'

;; -------- method by Monte-Carlo
nMC = 20
if (lin eq 1) then MCname = 'MCvarnoise_lin_'+strcompress(nside,/remove_all)+'_'+strcompress(Nscale,/remove_all)+'_'+strcompress(width,/REMOVE_ALL)+'_'+strcompress(nMC,/remove_all)+'.save' else MCname = 'MCvarnoise_'+strcompress(nside,/remove_all)+'_'+strcompress(Nscale,/remove_all)+'_'+strcompress(width,/REMOVE_ALL)+'_'+strcompress(nMC,/remove_all)+'.save'
filepresent = file_search(MCname,count=c)
if c eq 0 then begin 
   sigmaNoise2 = rims('images/WMAP7_rms_map_94GHz.fits')
   if npix ne N_elements(sigmaNoise2) then sigmaNoise2 = mrs_resize(sigmaNoise2, nside=nside)
   var_noise =  mrs_var_noise_mc(sigmaNoise2, bkj, nMC=nMC)
   save, var_noise, filename=MCname 
endif else restore, MCname

var_noise_MC = var_noise
var_noiseMC = var_noise_MC*(sigfact)^2.
print, '--- noise variance estimated by Monte-Carlo method.'

;; --------- Monte-Carlo estimation of the noise power spectrum
pnMC = 0. ; ?????? change to use same noise realization as MCvarnoise ?
for i=0, nMC-1 do begin
   pnMC = pnMC + mrs_powspec(randomn(seed, npix)*sigmaNoise)
endfor
pnMC = pnMC/nMC
print, '--- noise power spectrum estimated by Monte-Carlo method.'

;; ;;;;----------------- View cl's ------------------
;; cln = mrs_powspec(cmb2)
;; window,/free
;; plotcl, cl_inf[0:lmax]
;; oplotcl, cl[0:lmax], color=20
;; oplotcl, cln[0:lmax], color=100, line=1
;; stop

Here:

;;----- integrated Cl variables
;; true integrated Cl
var_needt = fltarr(Nscale, Ntest)
var_needt_ideal = fltarr(Nscale)
;; without mask
intcl_pn_ns = fltarr(Nscale, Ntest)
intcl_pn_MC_ns = fltarr(Nscale, Ntest)
intcl_block_ns = fltarr(Nscale, Ntest)
intcl_MC_ns = fltarr(Nscale, Ntest)
;; with mask
intcl_pn_mask = fltarr(Nscale, Ntest)
intcl_pn_MC_mask = fltarr(Nscale, Ntest)
intcl_block_mask = fltarr(Nscale, Ntest)
intcl_MC_mask = fltarr(Nscale, Ntest)

;;----- reconstructed Cl variables
;; reconstructed with tousi without mask
tousi_cl_pn_ns = fltarr(lmax+1, Ntest)
tousi_cl_pn_MC_ns = fltarr(lmax+1, Ntest)
;; reconstructed with toucan without mask with patch method
toucan_cl_block_ns = fltarr(lmax+1, Ntest)
toucan_ucl_block_ns = fltarr(lmax+1, Ntest)
toucan_part_nMSE_block_ns = fltarr(Niter, Ntest)
toucan_true_resi_block_ns = fltarr(Niter, Ntest)
;; reconstructed with toucan without mask with Monte-Carlo method
toucan_cl_MC_ns = fltarr(lmax+1, Ntest)
toucan_ucl_MC_ns = fltarr(lmax+1, Ntest)
toucan_part_nMSE_MC_ns = fltarr(Niter, Ntest)
toucan_true_resi_MC_ns = fltarr(Niter, Ntest)

;; reconstructed with tousi with mask
tousi_cl_pn_mask = fltarr(lmax+1, Ntest)
tousi_cl_pn_MC_mask = fltarr(lmax+1, Ntest)
;; reconstructed with toucan with mask with patch method
toucan_cl_block_mask = fltarr(lmax+1, Ntest)
toucan_ucl_block_mask = fltarr(lmax+1, Ntest)
toucan_part_nMSE_block_mask = fltarr(Niter, Ntest)
toucan_true_resi_block_mask = fltarr(Niter, Ntest)
;; reconstructed with toucan with mask with Monte-Carlo method
toucan_cl_MC_mask = fltarr(lmax+1, Ntest)
toucan_ucl_MC_mask = fltarr(lmax+1, Ntest)
toucan_part_nMSE_MC_mask = fltarr(Niter, Ntest)
toucan_true_resi_MC_mask = fltarr(Niter, Ntest)

Begining:
Zerotime = systime(1)
currenttime = STRMID(SYSTIME(0), 11, 5)
testtime = 60.*Ntest
for k=0, Ntest-1 do begin
   
   print, '-------- creating map - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)

   cmb_inf = getcmb(nside=nside,cl=cl_inf)
   mrs_convol, cmb_inf, h, cmb
   if k eq 0 then cl = cl_inf*h^2.
   n = randomn(seed,npix)
   ;; -------- bruit non-stationaire -------------
   cmb1 = cmb + n*sigmaNoise
   ;; -------- pas de bruit non-stationaire ------
;;    cmb1 = cmb
   ;; --------------------------------------------
   cmb2 = Mask*cmb1
      
;; ;-------------- Inpainting ------------------------------
   print, '-------- inpainting the mask...'
   rescaledMask = 0
   mrs_write,'in_mask.fits',mask
   mrs_write,'in_data.fits',cmb2
   spawn,'mrs_matmask in_mask.fits out_mat.fits'
   spawn,'mrs_alm_inpainting -m3 -M out_mat.fits  in_data.fits in_mask.fits  out_inp_data.fits'
   cmb3 = rims('out_inp_data.fits')
   print, '-------- inpainting done.'
      
;;    stop
;;    imcl = mrs_powspec(cmb2)
;;    loadct, 15
;;    window,/free
;;    plotcl, cl[0:lmax]
;;    oplotcl, imcl[0:lmax], color=20, line=2   

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

   n = randomn(seed,npix)
   NoiseReal = n*SigmaNoise
   pn = mrs_powspec(NoiseReal) 
   var_noiseBlock = mrs_var_noise_patch(NoiseReal, bkj)
   print, '--- noise power spectrum and noise variance by patch method estimated for one noise realization.'

;; ============================= Compressed Measurements ====================================
   print, '-------------------------- integrating Cl - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)

;; ;-------------- needlet decomposition --------------------
   print, '--- Computing Needlet decomposition of map without mask over', Nscale, ' scales...'
   mrs_wt_coeff, cmb1, needt, bkj
   print, '--- Needlet decomposition done'
   
;; ;-------------- compute the ideal needlet variance  and noise variances ---------------
   if k eq 0 then begin   
      var_noisePn = dblarr(Nscale)
      var_noisePnMC = dblarr(Nscale)
      ellcl = ell*cl[0:lmax]
      ellclPn = ell*pn[0:lmax]
      ellclPnMC = ell*pnMC[0:lmax]
      for i=0,Nscale-1 do begin
         var_needt_ideal(i) = total(bkj[*,i]^2*ellcl)/(4*!dpi)
         var_noisePn(i) = total(bkj[*,i]^2*ellclPn)/(4*!dpi)
         var_noisePnMC(i) = total(bkj[*,i]^2*ellclPnMC)/(4*!dpi)
      endfor
      print, '--- Ideal needlet variance computed'
   endif

;; ;-------------- compute the true needlet variance ---------------
   realcl = mrs_powspec(cmb)
   ellcl = ell*realcl[0:lmax]
   for i=0,Nscale-1 do begin
      var_needt(i,k) = total(bkj[*,i]^2*ellcl)/(4*!dpi)
   endfor
   print, '--- True needlet variance computed'
   
   ;; --------- estimated needlet variance
   for i=0,Nscale-1 do begin
      intcl_pn_ns[i,k] = total(needt.coef[*,i]^2-var_noisePn[i])/npix
      intcl_pn_MC_ns[i,k] = total(needt.coef[*,i]^2-var_noisePnMC[i])/npix
      intcl_block_ns[i,k] = total(needt.coef[*,i]^2-var_noiseBlock[*,i])/npix
      intcl_MC_ns[i,k] = total(needt.coef[*,i]^2-var_noiseMC[*,i])/npix
   endfor

;; ;-------------- needlet decomposition --------------------
   print, '--- Computing Needlet decomposition of inpainted masked map over', Nscale, ' scales...'
   mrs_wt_coeff, cmb3, needt, bkj
   print, '--- Needlet decomposition done'
   
;; --------- estimated needlet variance
  for i=0,Nscale-1 do begin
      intcl_pn_mask[i,k] = total(Mask*(needt.coef[*,i]^2-var_noisePn[i]))/nmask
      intcl_pn_MC_mask[i,k] = total(Mask*(needt.coef[*,i]^2-var_noisePnMC[i]))/nmask
      intcl_block_mask[i,k] = total(Mask*(needt.coef[*,i]^2-var_noiseBlock[*,i]))/nmask
      intcl_MC_mask[i,k] = total(Mask*(needt.coef[*,i]^2-var_noiseMC[*,i]))/nmask
   endfor

;; stop

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

;; ============================= Power Spectrum Reconstruction ====================================
   print, '-------------------------- reconstructing Cl - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)  
   Bmat = mrs_compute_Bmat(bkj)

;; ;; power spectrum from tousi
   print, '------------------- tousi no mask - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)
   mastercl = mrs_powspec(cmb1,lmax=lmax) ;????? lmax = 3*nside ?
   tousicl = mrs_tousi(mastercl, NoisePs=pn)
   tousi_cl_pn_ns[*,k] = tousicl
   tousicl = mrs_tousi(mastercl, NoisePs=pnMC)
   tousi_cl_pn_MC_ns[*,k] = tousicl

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

;; ;; power spectrum from toucan
   print, '------------------- toucan no mask by patches - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)
   toucancl = test_cl_recons(cmb1, sigmaNoise, intcl_block_ns[*,k], filters=bkj, Nscale=Nscale, Firstl=Firstl, lmax=lmax, Niter=Niter, Cl=Cl, Mcl=masterCl, tousil=tousil, part_nMSE = part_nMSE, true_resi=true_resi, uresult=toucanucl)
   toucan_cl_block_ns[*,k] = toucancl
   toucan_ucl_block_ns[*,k] = toucanucl
   toucan_part_nMSE_block_ns[*,k] = part_nMSE
   toucan_true_resi_block_ns[*,k] = true_resi

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

   print, '------------------- toucan no mask by MC - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)
   toucancl = test_cl_recons(cmb1, sigmaNoise, intcl_MC_ns[*,k], filters=bkj, Nscale=Nscale, Firstl=Firstl, lmax=lmax, Niter=Niter, Cl=Cl, Mcl=masterCl, tousil=tousil, part_nMSE = part_nMSE, true_resi=true_resi, uresult=toucanucl)
   toucan_cl_MC_ns[*,k] = toucancl
   toucan_ucl_MC_ns[*,k] = toucanucl
   toucan_part_nMSE_MC_ns[*,k] = part_nMSE
   toucan_true_resi_MC_ns[*,k] = true_resi

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

;; ;; power spectrum from tousi
   print, '------------------- tousi inpainted mask - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)
   mastercl = mrs_master_powspec(cmb3,Mask,lmax=3l*nside) ; cannot put 2l*nside here
   tousicl = mrs_tousi(mastercl, NoisePs=pn)
   tousi_cl_pn_mask[*,k] = tousicl[0:lmax]
   tousicl = mrs_tousi(mastercl, NoisePs=pnMC)
   tousi_cl_pn_MC_mask[*,k] = tousicl[0:lmax]

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

;; ;; power spectrum from toucan
   print, '------------------- toucan inpainted mask by patches - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)
   toucancl = test_cl_recons(cmb3, sigmaNoise, intcl_block_mask[*,k], filters=bkj, Nscale=Nscale, Firstl=Firstl, Mask=Mask, lmax=lmax, Niter=Niter, Cl=Cl, Mcl=masterCl, tousil=tousil, part_nMSE = part_nMSE, true_resi=true_resi, uresult=toucanucl)
   toucan_cl_block_mask[*,k] = toucancl
   toucan_ucl_block_mask[*,k] = toucanucl
   toucan_part_nMSE_block_mask[*,k] = part_nMSE
   toucan_true_resi_block_mask[*,k] = true_resi

   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

   print, '------------------- toucan inpainted mask by MC - iteration number:' + strcompress(k+1) + ' out of' + strcompress(Ntest)
   toucancl = test_cl_recons(cmb3, sigmaNoise, intcl_MC_mask[*,k], filters=bkj, Nscale=Nscale, Firstl=Firstl, Mask=Mask, lmax=lmax, Niter=Niter, Cl=Cl, Mcl=masterCl, tousil=tousil, part_nMSE = part_nMSE, true_resi=true_resi, uresult=toucanucl)
   toucan_cl_MC_mask[*,k] = toucancl
   toucan_ucl_MC_mask[*,k] = toucanucl
   toucan_part_nMSE_MC_mask[*,k] = part_nMSE
   toucan_true_resi_MC_mask[*,k] = true_resi

   if ((k mod 5) eq 4) then save, nside, Nscale, width, Ntest, lin, mex, tousil, bkj, cl, Niter, var_needt, var_needt_ideal, intcl_pn_ns, intcl_pn_MC_ns, intcl_block_ns, intcl_MC_ns, intcl_pn_mask, intcl_pn_MC_mask, intcl_block_mask, intcl_MC_mask, tousi_cl_pn_ns, tousi_cl_pn_MC_ns, toucan_cl_block_ns, toucan_ucl_block_ns, toucan_part_nMSE_block_ns, toucan_true_resi_block_ns, toucan_cl_MC_ns, toucan_ucl_MC_ns, toucan_part_nMSE_MC_ns, toucan_true_resi_MC_ns, tousi_cl_pn_mask, tousi_cl_pn_MC_mask, toucan_cl_block_mask, toucan_ucl_block_mask, toucan_part_nMSE_block_mask, toucan_true_resi_block_mask, toucan_cl_MC_mask, toucan_ucl_MC_mask, toucan_part_nMSE_MC_mask, toucan_true_resi_MC_mask, filename=savename

   Testtime = systime(1) - zerotime
   Testtime = Testtime/((k+1)*60)*(Ntest-k-1)
   currenttime = STRMID(SYSTIME(0), 11, 5)
   print, '--- estimated remaining time at ' + currenttime + ' ='+strcompress(testtime)+' minutes.'

endfor

;; save, nside, Nscale, Ntest, cl, cutoffreq, var_needt, intcl_pn,
;; intcl_block, intcl_MC, filename=savename
;; save, nside, Nscale, Ntest, width, cl, var_needt, intcl_pn, intcl_block, filename=savename


;;============================== Plots ====================================

;;--------------- Compare noise variances ------------------
;; ave_var_noiseBlock = fltarr(Nscale)
;; for i=0, Nscale-1 do ave_var_noiseBlock(i) = mean(var_noiseBlock[*,i])
;; ave_var_noiseMC = fltarr(Nscale)
;; for i=0, Nscale-1 do ave_var_noiseMC(i) = mean(var_noiseMC[*,i])

;; window,/free
;; plot, var_noise
;; oplot, var_noise2, line=2, color=20
;; oplot, ave_var_noiseBlock, color=100
;; oplot, ave_var_noiseMC, line=2, color=200

; cosmic variance
vcos = dblarr(Nscale)
ellcl2 = ell*cl[0:lmax]^2
for i=0,Nscale-1 do vcos[i] = 2/(4*!dpi)*total(bkj[*,i]^4*ellcl2)

MSE:
;;------------------ Compute the average error on one example ---------------------
k = 0
print, '--------'
print, 'est. intcl by powspec - no mask', sigma(var_needt_ideal-intcl_pn_ns[*,k])
print, 'est. intcl by MC powspec - no mask', sigma(var_needt_ideal-intcl_pn_MC_ns[*,k])
print, 'est. intcl by patch - no mask', sigma(var_needt_ideal-intcl_block_ns[*,k])
print, 'est. intcl by MC - no mask', sigma(var_needt_ideal-intcl_MC_ns[*,k])
print, '--------'
print, 'est. intcl by powspec - with mask', sigma(var_needt_ideal-intcl_pn_mask[*,k])
print, 'est. intcl by MC powspec - with mask', sigma(var_needt_ideal-intcl_pn_MC_mask[*,k])
print, 'est. intcl by patch - with mask', sigma(var_needt_ideal-intcl_block_mask[*,k])
print, 'est. intcl by MC - with mask', sigma(var_needt_ideal-intcl_MC_mask[*,k])
print, '--------'
print, 'est. cl by powspec with tousi - no mask', sigma(cl[0:lmax] - tousi_cl_pn_ns[*,k])
print, 'est. cl by MC powspec with tousi - no mask', sigma(cl[0:lmax] - tousi_cl_pn_MC_ns[*,k])
print, '--------'
print, 'est. cl by patch with toucan - no mask', sigma(cl[0:lmax] - toucan_cl_block_ns)
print, 'est. cl by MC with toucan - no mask', sigma(cl[0:lmax] - toucan_cl_MC_ns)
print, '--------'
print, 'est. cl by powspec with tousi - with mask', sigma(cl[0:lmax] - tousi_cl_pn_mask[*,k])
print, 'est. cl by MC powspec with tousi - with mask', sigma(cl[0:lmax] - tousi_cl_pn_MC_mask[*,k])
print, '--------'
print, 'est. cl by patch with toucan - with mask', sigma(cl[0:lmax] - toucan_cl_block_mask)
print, 'est. cl by MC with toucan - with mask', sigma(cl[0:lmax] - toucan_cl_MC_mask)
print, '--------'

stop

;; ;; --------------------------------------- Best iteration -----------------------------------------------------
BESTIT:

print, '--------'
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_part_nMSE_block_ns[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit nMSE by patch with toucan - no mask =', bestit
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_true_resi_block_ns[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit resi by patch with toucan - no mask =', bestit
print, '--------'
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_part_nMSE_MC_ns[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit nMSE by MC with toucan - no mask =', bestit
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_true_resi_MC_ns[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit resi by MC with toucan - no mask =', bestit
print, '--------'
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_part_nMSE_block_mask[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit nMSE by patch with toucan - with mask =', bestit
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_true_resi_block_mask[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit resi by patch with toucan - with mask =', bestit
print, '--------'
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_part_nMSE_MC_mask[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit nMSE by MC with toucan - with mask =', bestit
bestit = 0
for k=0, Ntest-1 do begin
   minvar = min(toucan_true_resi_MC_mask[*,k], ind)
   bestit = bestit+ind
endfor
bestit = bestit/Ntest
print, 'bestit resi by MC with toucan - with mask =', bestit

stop

;; ;; --------------------------------------- Plots -----------------------------------------------------
PLOTS:

;; loadct, 15

;; ;;----------------- Plot the estimated integrated cl for one width ------------------
;; window,/free
;; xtit = 'wavelet filters'
;; ytit = 'wavelet coef variance'
;; plot, var_needt, xtitle=xtit, ytitle=ytit
;; oplot, intcl_pn, color=20
;; oplot, intcl_block, color=60, linestyle=2
;; legend, ['Cl integrated over scale j','needlet est. from noise powspec','needlet est. block'],linestyle=[0,0,2],colors=[0,20,60],/righ


;; ;; ;;--------------------- Log plot of error for each width ------------------

;; error_pn = abs(var_needt-intcl_pn)
;; error_block = abs(var_needt-intcl_block)
;; window,/free
;; tit = 'Error for filter width = ' + strcompress(width)
;; xtit = 'wavelet filters'
;; ytit = 'error on wavelet coef variance'
;; plot, error_pn, title=tit, xtitle=xtit, ytitle=ytit, yrange=[min([error_pn,error_block]), max([error_pn,error_block])], /ylog ,/nodata ;, BACKGROUND = 255, COLOR = 0
;; oplot, error_pn, color=20
;; oplot, error_block, color=60, linestyle=2
;; legend, ['needlet est. w. ns noise','needlet est. block'],linestyle=[0,2],colors=[20,60],/right

;; ;;--------------------- Recons best iteration number ------------------
;; bestit = 0
;; minvar = min(ntousi_true_resi, ind)
;; ;;    print, 'test', i+1, ' : min MSE = ', minvar, ' at iteration', ind
;; ;;    print, 'ntousi_true_resi = ', ntousi_true_resi[ind,i], ' at iteration', ind
;; bestit = ind
;; print, '--- best iteration for width'+ strcompress(width) + ' =' + strcompress(bestit)

;; bestit = 0
;; minvar = min(ntousi_part_nMSE, ind)
;; ;; print, 'test', i+1, ' : min MSE = ', minvar, ' at iteration', ind
;; ;; print, 'ntousi_part_nMSE = ', ntousi_part_nMSE[ind,i], ' at iteration', ind
;; bestit = ind
;; print, '--- best iteration for width'+ strcompress(width) + ' =' + strcompress(bestit)

;; ;----------------------------- error --------------------------------------

;; toucan_error = fltarr(lmax+1,Ntest)

;; for j=2,lmax do begin
;;    toucan_error[j] = abs(cl[j] - ntousi_cl_test[j])
;; endfor

;; ll = findgen(lmax+1)
;; ll = ll*(ll+1.)/(2*!pi)

;; ;; for i=0, Nbscale-1 do begin
;; ;;    toucan_error[*,i] = toucan_error[*,i]*ll
;; ;; endfor

;; coloring = floor(200*findgen(Ntest)/(Ntest-1))
;; window,/free
;; tit = 'Abs Error'
;; xtitle = '!6 Multipole !12 l !6'
;; ytitle = '!12  l (l+1) !6RMS!n /!6 2!7p!n!6 '
;; plot, toucan_error*ll, title=tit, xtitle=xtitle, ytitle=ytitle, charsize=1.3, xrange=[0, lmax], xstyle=1 ;, /ylog ;, xrange=[0,500];, yrange=[0,300]
;; ;; for i=1, Ntest-1 do begin
;; ;;    oplot, toucan_error[*,i]*ll, line=2, color=coloring[i]
;; ;; endfor
;; ;; legend, 'width ='+strcompress(width[0:Ntest-1]),linestyle=[0,2+intarr(Ntest-1)],colors=coloring,/righ, charsize=1.3

;; ;----------------------------- Cl --------------------------------------

;; ;; coloring = 20 + floor(180*findgen(Ntest)/(Ntest-1))
;; window,/free
;; plotcl, cl[0:lmax], xrange=[0, lmax] ;, yrange = [0, 10000]
;; ;; for i=0, Ntest-1 do begin
;;    oplotcl, ntousi_cl_test, line=2, color=20 ;, color=coloring[i]
;; ;; endfor
;; legend, ['True Cl','width ='+strcompress(width)],linestyle=[0,2],colors=[0,20],/right, charsize=1.3 ;,/bottom
;; ;; legend, ['True Cl','width ='+strcompress(width[0:Ntest-1])],linestyle=[0,2+intarr(Ntest)],colors=[0,coloring],/right, charsize=1.3 ;,/bottom


end
