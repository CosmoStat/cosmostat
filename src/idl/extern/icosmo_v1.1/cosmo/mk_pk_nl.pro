function mk_pk_nl,fid,cosmo_new,pk_lin,z=z

; July 08 - written by AR from mk_del2.pro
;
; PURPOSE: Compute the nonlinear matter power spectrum given the
; linear power spectrum.This is done using the from the prescrition of
; Peacock and Dodds (1996, mnras, 280, L19), as summarized by Peacock
; (1997, mnras, 284, 885). Optionally, the Ma et al. (1999, ApJ, 521,
; L1) or Smith et al. (2002, astro-ph/0207664) formula is used.
; INPUT: fid: fiducial parameters from set_fiducial.pro
;        cosmo: cosmological parameter structure from mk_cosmo.pro
; OPTIONAL INPUT:
;        pk_l: linear power spectrum from mk_pk_lin.pro (if not
;              provided it is computed)
;        z: redshift (default: z=0)
; FIUCIAL PARAMETERS USED:  
;                     kl_ran: range of linear k values to consider 
;                             [h Mpc^-1] (h-corrected physical units;
;                             default= [.01,100.] h Mpc^-1). Note that
;                             the non-linear k range will be shifted
;                             from that. 
;                     n_k: number of k values to consider (default=100)
;                     fit_nl: fitting function for non-linear corrections
;                             to use: 0: Peacock and Dodds, 1: Ma et
;                             al., 2: Smith et al. 
;                     fit_tk: fitting function for transfer function
;                     to use (Note: This is not used by mk_del2
;                     but is passed to del2_lin) : 
;			0: E&H without wiggles (default)	
;			1: E&H with wiggles (UNSTABLE w/ NON-LINEAR CORRECTION)
;			2: BBKS as summarised by Peacock & Dodds (1997)
; OUTPUT: pk_nl structure containing:
;            k: (nonlinear) comoving wave number in h-corrected physical units
;                i.e. k=k_comoving/R_0/h   [h Mpc^-1]
;            del2: Delta^2(k) (nonlinear) [1]
;            pk: Power spectrum P(k) (nonlinear) [h^-3 Mpc^3]
;            del2_l: linear Delta^2(k_nl)  [1]
;            pk_l: linear power spectrum [h^-3 Mpc^3]
;
; NOTE: Peacock and Dodds and Ma et al. not supported yet as the
; linear wavenumbers k_l need to be mapped to final k values

; declarations: 
if not keyword_set(z) then z=0.    ; redshift
aa=1./(1.d +z)
cosmo=cosmo_new.const  ; map to old structures
evol=cosmo_new.evol
;kl_ran=fid.calc.k_ran  ; Range of k (linear k in h-corrected physical
;                       ; coordinates, i.e. k=k_comoving/R_0/h   [h Mpc^-1]) 
;n_k=fid.calc.n_k       ; Number of k values for integration
fit_nl=fid.calc.fit_nl     ; non-linear correction (0: P&D, 1: Ma et al,
                       ; 2:Smith et al)

; compute linear power spectrum if not provided
if not keyword_set(pk_lin) then pk_lin=mk_pk_lin(fid,cosmo_new,z=z)
kl=pk_lin.k          ; (linear) wave number [h Mpc^-1]
n_k=n_elements(kl)
del2_l=pk_lin.del2_l   ; linear variance and power spectrum
pk_l=pk_lin.pk_l

; compute non-linear correction for the choice of fitting functions
case fit_nl of

  0: begin ; Peacock & Dodds

;     stop,'ERROR: mk_pk_nl: Peacock & Dodds not implemented'
  
    ; compute the effective spectral index by computing the
    ; logarithmic derivative of the power spectrum (see eq. 25 in 
    ; Peacock 97)
    neff=dblarr(n_k-2)
    for i=1,n_k-2 do $
      neff(i-1)=-3.+alog(del2_l(i+1)/del2_l(i-1)) / alog(kl(i+1)/kl(i-1))
    kl_neff=2.d*kl(1:n_k-2)

    ; parameters from Peacock (1997) eqs 20-24
    n=interpol(neff,kl_neff,kl)
    a=0.482d *(1.d +n/3.d)^(-.947d)
    b=0.226d*(1.d +n/3.d)^(-1.778d)
    alpha=3.310d*(1.d +n/3.d)^(-.244d)
    beta=0.862d*(1.d +n/3.d)^(-0.287d)
    v=11.55d*(1.d +n/3.d)^(-0.423d)

    ; compute the nonlinear variance spectrum
    x=del2_l
;;    mk_evol,cosmo,evol,z=z
;;    g=evol.g
    d=growth_a_ode(cosmo,aa,/unnorm)
;    d=interpol(cosmo_new.evol.d,cosmo_new.evol.a,aa)
    g=d/aa
    del2_nl=x*( (1.d +b*beta*x+(a*x)^(alpha*beta)) /$
                (1.d +( (a*x)^alpha*g^3.d /(v*sqrt(x)) )^beta ) )^(1.d /beta)
    
                                ; compute the non-linear wave number
    k_nl=(1.d +del2_nl)^(1.d /3.d )*kl


 end                            ; (Peacock and Dodds)
  
  1: begin                      ; Ma et al. (199)

     stop,'ERROR: mk_pk_nl: Ma et al. not implemented'
  
    ; check if model is supported
    if cosmo.wa ne 0. then $
      message,'evolving w models not supported with Ma et al. fitting function'

    ; parameters for non-quintessence models
    c1=1.08d-4
    c2=2.10d-5  
    beta=0.7d
    mk_evol,cosmo,evol,z=z
    g_0=evol.g0
    g=evol.g
    if cosmo.norm eq 1 then sig8=cosmo.sigma8 else sig8=sigma8(cosmo)
    sig8=sig8*evol.d    ; sigma8 at z of interest
  
    ; corrections for quintessence models
    if cosmo.omega_l ne 0. and cosmo.w0 ne -1. then begin
      if cosmo.omega ne 1. then begin
        print,'Non flat quintessence models not supported'
        stop
      endif
      beta=0.83d
      g_0=abs(cosmo.w_l)^(1.3d *abs(cosmo.w0)-.76d )*evol.g0
    endif

    ; compute non-linear variance and power spectrum
    x=del2_l/g_0^1.5d /sig8^beta
    del2_nl=del2_l*(1.d +alog(1.+0.5d *x))*(1.d +0.02d *x^4.d +c1*x^8.d /g^3.d )/$
                (1.d +c2*x^7.5d)   
    ; compute the non-linear wave number
    k_nl=(1.d +del2_nl)^(1.d /3.d )*kl

     end  ; (Ma et al.)

  2: begin  ; Smith et al. (2002)

      ; compute non-linear scale and related quantities
;      mk_sigmag,cosmo,sigmag,z=z         
      sigmag=mk_sigmag(pk_lin)    
      n=sigmag.n_sigma    ; effective slope at k_sigma
      c=sigmag.c_sigma    ; curvature at k_sigma
      k_nl=kl
      y=k_nl/sigmag.k_sigma  ; defined as sigma_gaussian(k_sigma^-1)=1      

      ; parameters from the appendix of smith et al.
      a_n=10.d ^(1.4861d +1.8369d *n+1.6762d *n^2.d +0.7940d *n^3.d +0.1670d *n^4.d -0.6206d*c)
      b_n=10.d ^(0.9463d +0.9466d *n+0.3084d *n^2.d -0.9400d *c)
      c_n=10.d ^(-0.2807d +0.6669d *n+0.3214d *n^2.d -0.0793d *c)
      gamma_n=0.8649d +0.2989d *n+0.1631d *c
      alpha_n=1.3884d +0.3700d *n-0.1452d *n^2.d
      beta_n=0.8291d +0.9854d *n+0.3401d *n^2.d 
      mu_n=10.d ^(-3.5442d +0.1908d *n)
      nu_n=10.d ^(0.9589d +1.2857d *n)

      ; cosmology dependent functions
;      if cosmo.w_l ne -1. or cosmo.omega gt 1. or cosmo.w1_l ne 0. then begin
      if cosmo.omega gt 1. then begin
;        print,'mk_del2: Smith et al. formula not defined for w<>-1 or Omega>1'
;         print,'NEED to be careful about closed models!!!'
        ;stop
      endif      
      om_m=interpol(evol.omega_m_a,evol.z,z)   ; Omega_m(z)
      om_v=interpol(evol.omega_l_a,evol.z,z)   ; Omega_l(z)
      om_k=interpol(evol.omega_k_a,evol.z,z)   ; Omega_k(z)
      w_de=interpol(evol.w_a,evol.z,z)         ; w_de(z)
;      stop
      if om_m eq 1. then begin
         f1=1.d & f2=1.d & f3=1.d
      endif else begin
        f1a=om_m^(-0.0732d)
        f2a=om_m^(-0.1423d)
        f3a=om_m^(0.0725d)
        f1b=om_m^(-0.0307d)
        f2b=om_m^(-0.0585d)
        f3b=om_m^(0.0743d)       
 
      ; begin with Omega_m+Omega_de+Omega_k=1
      ;1 - define frac1 as frac1 = Omega_de/(Omega_de+Omega_k)
      ;    (i.e. frac1=1 means Omega_k=0 and f1=0 is Omega_de=0)
      ;2 - define effective w as we= frac1*w1 +(1-frac1)*w2
      ;    (where w1 is DE w and w2 is curvature w=-1./3.)
      ;3 - define frac2 as frac2 = -(3we +1)/2
      ;    this is used to interpolate between the two 
      ;    fitting funtions given by Smith et al.
        frac1=om_v/(om_v+om_k) 
        we=(frac1*w_de)+((1.d -frac1)*(-1.d/3.d))
        frac2=(-1.d)*((3.d *we)+1.d )/2.d


       f1=frac2*f1b+(1.d -frac2)*f1a
       f2=frac2*f2b+(1.d -frac2)*f2a
       f3=frac2*f3b+(1.d -frac2)*f3a
      endelse

      ; compute non-linear power spectrum
      underflow=check_math(mask=32,/noclear)  ; check underflow flag
      del2_q=del2_l*(1.d +del2_l)^beta_n/(1.d +alpha_n*del2_l)*exp(-y/4.d -y^2.d /8.d)
      if underflow eq 0 then underflow=check_math(mask=32) 
                    ; supress underflow error message if necessary
      del2_h=a_n*y^(3.d *f1)/(1.d +b_n*y^f2+(f3*c_n*y)^(3.d -gamma_n))
      del2_h=del2_h/(1.d +mu_n/y+nu_n/y^2.d)
      del2_nl=del2_q+del2_h

    end  ; (smith et al.)

endcase

; if needed map linear k's onto nonlinear k's
if fit_nl ne 2 then del2_nl=interpol(del2_nl,alog(k_nl),alog(kl))

; compute nonlinear power spetrum
pk_nl=2.*!pi^2/kl^3*del2_nl

; store result in a structure
pk_nl={z:z,k:kl,pk:pk_nl,del2:del2_nl,pk_l:pk_l,del2_l:del2_l}

return,pk_nl
end



