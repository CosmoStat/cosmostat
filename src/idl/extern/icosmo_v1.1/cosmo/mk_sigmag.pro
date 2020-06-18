function mk_sigmag,pk,r_ran=r_ran,n_r=n_m,plotit=plotit

; July 2008 - Modified by AR to turn it into a function and use pk
;             structure as input
; May 2003 - written by A. Refregier
;
; PURPOSE: compute the linear rms density contrast sigma_g(R) with a gaussian
; window function of size R. This is to be input in the computation
; of the non-linear power spectrum using the Smith et al. 
; (astro-ph/0207664) fitting function.
; INPUT: pk: (linear) power spectrum from mk_pk_lin.pro (or mk_pk_nl.pro)
; OPTIONAL INPUT: r_ran: radius R range [h^-1 Mpc]
;                 n_r: number of r values (default=50)
;                 plotit: plot results if set
; OUTPUT: sigmag structure:
;           r: R(M) [h^-1 Mpc]
;           sigma: sigma(M)    [1]

; declaration
if not keyword_set(r_ran) then r_ran=[.005,50.]               
                                ; radius R range [h^-1 Mpc]
if not keyword_set(n_r) then n_r=50                          
                                ; number of radii to consider
r=xgen(r_ran(0),r_ran(1),n=n_r,/log)
n_k=n_elements(pk.k)        ; number of k values to consider
k=pk.k

; check if range of k values is sufficient
if k(0) gt .001 or k(n_k-1) lt 100. then $
  stop,'ERROR: mk_sigmag.pro: k range insufficient' 

;; compute linear power spectrum 
;if not keyword_set(z) then z=0.
;del2_l=del2_lin(cosmo,k,z)
underflow=check_math(mask=32,/noclear) ; check floating underflow flag
                                ; (to supress floating undeflow error message)

; compute sigma_g(R) and related quantities for each R (see eq 60-61
; and appendix C in Smith et al.)
del2_l=pk.del2_l
sigma=fltarr(n_r)    
neff=fltarr(n_r)   
curv=fltarr(n_r)
for i=0,n_r-1 do begin
  yi=k*r(i)
  sigma[i]=sqrt(int_tabulated(alog(k),del2_l*exp(-yi^2)))
  neff[i]=-3.+2./sigma[i]^2*$
          int_tabulated(alog(k),del2_l*yi^2*exp(-yi^2))
  curv[i]=(3.+neff(i))^2+4./sigma[i]^2*$
          int_tabulated(alog(k),del2_l*(yi^2-yi^4)*exp(-yi^2))
;  print,'r,sigma,neff,curv:',r(i),sigma(i),neff(i),curv(i)
;  read,test
;  if test eq 1 then stop
endfor

; compute k_sigma defined as sigma(k_sigma^-1)=1, and related
; quantities
r_sigma=exp(interpol(alog(r),alog(sigma),0.))
if r_sigma lt r_ran(0) or r_sigma gt r_ran(1) then begin
  print,'mk_sigmag: r_sigma is out of range'
  stop
endif
if underflow eq 0 then underflow=check_math(mask=32) 
                                ; clear floating underflow flag if it
                                ; was originally clear
k_sigma=1./r_sigma
neff_sigma=interpol(neff,alog(sigma),0.)
curv_sigma=interpol(curv,alog(sigma),0.)

; plot results if requested
if keyword_set(plotit) then begin
  !p.multi=[0,1,3]
  old=!p.charsize & !p.charsize=2
  plot,r,sigma,/xtype,/ytype,/xstyle,/ystyle,$
    xtitle='R (h^-1 Mpc)',ytitle='sigma_g(R)'
  oplot,r_sigma*[1.,1.],[.001,1000.],lines=1
  plot,r,neff,/xtype,/xstyle,/ystyle,$
    xtitle='R (h^-1 Mpc)',ytitle='neff(R)'
  oplot,r_sigma*[1.,1.],[-10.,10.],lines=1
  plot,r,curv,/xtype,/xstyle,/ystyle,$
    xtitle='R (h^-1 Mpc)',ytitle='curv(R)'
  oplot,r_sigma*[1.,1.],[-10.,10.],lines=1
;  plot,sigma,neff,/xtype,xran=[.04,10.],yran=[-2.6,0.],/xsty,/ysty,$
;    xtitle='sigma_g(R)',ytitle='n_eff(R)'
;  plot,sigma,curv,/xtype,xran=[.04,10.],yran=[0.1,0.8],/xsty,/ysty,$
;    xtitle='sigma_g(R)',ytitle='Curv(R)'
  !p.charsize=old & !p.multi=0
  print,'r_sigma, k_sigma [h^-1 Mpc,h Mpc^-1]:',r_sigma,k_sigma
  print,'n_sigma,c_sigma:',neff_sigma,curv_sigma
endif

;; compute d(ln sigma)/d(ln M) numerically
;dlnsdlnm=deriv(alog(m),alog(sigma))

; store in structure
sigmag={r:r,sigma:sigma,n:neff,c:curv,$
  k_sigma:k_sigma,n_sigma:neff_sigma,c_sigma:curv_sigma}

return,sigmag
end


