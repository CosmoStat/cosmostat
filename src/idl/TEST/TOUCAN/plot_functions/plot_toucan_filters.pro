;; Make some of the plots for the TOUCAN paper
; .compile mrs_get_cl_theo_powspec

;; ============ Parameters of Toucnan =============
width = 140
nside = 512 ; 512; 256
Nscale = 60
Firstl = 2 ;800
lmax = 2l*nside
Niter = 240

;; ============ Parameters of cmb =============
npix = nside2npix(nside)
ell = findgen(lmax+1)*2.+1.
fsky = 1

;; ============ Generate the wavelet filters =============
;; bkj = toucan_filters(Nscale, Firstl=Firstl, lmax=lmax, filtwidth=width, Bmat=Bmat) ; could add lin or log, mod, width etc...

;; =================== Plots =======================
loadct, 39
Ncolors = 5
;; colors = 50+(254-50)*(findgen(Ncolors)/Ncolors)
colors = [50,100,150,200,250]

N1 = 0 & N2 = Nscale-1

imagename = '~/images/toucan_filters_'+strcompress(nside,/remove_all)+'_'+strcompress(nscale,/remove_all)+'.eps'
;; ops,file=imagename, /COLOR,/landscape
window, /free
plot, bkj[*,0], yrange=[0,1], xrange=[0,1200], xtitle='!6 Multipole !12 l !6',/nodata, charsize=1.3
for j=N1, N2 do oplot, bkj[*,j], color=colors((j mod Ncolors))
;; cps

end

