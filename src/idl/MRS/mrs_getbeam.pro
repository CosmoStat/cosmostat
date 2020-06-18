function mrs_getbeam, Fwhm=Fwhm, lmax=lmax, Sigma=Sigma
if not keyword_set(lmax) then lmax=4000
if not keyword_set(Fwhm) then begin
   Fwhm=10.  ; 10 arc minutes
   if keyword_set(sigma) then Fwhm= 2. * Sigma * sqrt(2.*alog(2.))
end

F = Fwhm / 60. * !dtor
l = findgen(lmax+1)
ell = l*(l+1)
 bl = exp(-ell*F^2. /16./alog(2.))
  
;plot, bl
print, bl[0:10]
return, bl
end


;=============================================

function mrs_bandpassfilter, Fwhm1=Fwhm1, lmax=lmax, Sigma1=Sigma1, Fwhm2=Fwhm2,  Sigma2=Sigma2
if not keyword_set(lmax) then lmax=4000
if not keyword_set(Fwhm1) then begin
   Fwhm1=10.  ; 10 arc minutes
   if keyword_set(sigma1) then Fwhm1= 2. * Sigma1 * sqrt(2.*alog(2.))
end

if not keyword_set(Fwhm2) then begin
   Fwhm2=Fwhm1*2.
   if keyword_set(sigma2) then Fwhm2= 2. * Sigma2 * sqrt(2.*alog(2.))
end

F = Fwhm1 / 60. * !dtor
l = findgen(lmax+1)
ell = l*(l+1)
bl = exp(-ell*F^2. /16./alog(2.))

F2 = Fwhm2 / 60. * !dtor
bl2 = exp(-ell*F2^2. /16./alog(2.))
 
bl = bl - bl2
bl = bl / max(bl) 
;plot, bl
print, bl[0:10]
return, bl
end

;=============================================

function mrs_beamconvol, data, beam
vs = size(beam)
lmax = vs[1] - 1
mrs_almtrans, data, T_ALM, lmax=lmax
alm = T_ALM.alm
alm_product2, alm, beam, res
T_ALM.alm = res
mrs_almrec, T_ALM, map
return, map
end

;=============================================

function testconv, map, almwmap=almwmap, sig1, sig2
if not keyword_set(almwmap) then mrs_almtrans, map, almwmap
if not keyword_set(sig1) then sig1=250.
if not keyword_set(sig2) then sig2=500.

T_ALM = almwmap
beam =  mrs_bandpassfilter( sigma1=sig1, sigma2=sig2, lmax=almwmap.lmax)
alm = T_ALM.alm
alm_product2, alm, beam, res
T_ALM.alm = res
mrs_almrec, T_ALM, map
return, map
end

function wmapthreshold, map, nsig
signoise = mrs_mad(map)
ThresholdLevel = signoise * nsig
fil =  mrs_absthreshold(map, ThresholdLevel)
tvs, fil
return, fil
end

pro test_wmap
map = mrs_read('/dsm/sappcwf/planck/wmap_3years/wmap_ilc_3yr_v2.fits')
map = map[*,0]
mrs_almtrans, map, ALmwmap
sig1=250.
sig2=500.
mapconv = testconv(map, almwmap=almwmap, sig1, sig2)
nsig=5
t = wmapthreshold(mapconv, nsig)

g = mrs_healpix2glesp(map, /alm)
g.sky = g.t_sky / 500.
mrs_mca, residual=r, g, tabCw, SelectTrans=[4,2], niter=100, /mom,  LastThreshold=LastThreshold, SigmaNoise=1, NbrScale=8

end

;=============================================
