function sigma8,pk

; July 2008 - Modified by AR to allow computation from provided pk
;             structure and to reflect new conventions for the cosmo structure
; April 2000 - Written by A. Refregier
; May 2001 - modified by AR to act as a function and only consider the
;            linear power spectrum
; PURPOSE: compute sigma8, the rms density contrast in 8 h^-1 Mpc
; from an input (linear) power spectrum
; INPUT: cosmo: cosmological parameter structure read by rd_cosmo.pro
; OUTPUT: sigma8 is returned

; declarations
kl=pk.k     ; wavenumber k for linear power spectrum [h Mpc^-1]
n_kl=n_elements(kl)
k_ran_sigma8=[.01,10.] ; minimum range of k-values to compute sigma8 [h Mpc^-1]

; check that k-range is sufficient
if kl(0) gt k_ran_sigma8(0) or kl(n_kl-1) lt k_ran_sigma8(1) then $
       message,'Error: mk_pk.pro:k range insufficient to compute sigma8'

; compute power spectrum
;del2_l=del2_lin(cosmo,kl,0.)
;pk_l=2.*!pi^2/kl^3*del2_l
;del2={del2_l:del2_l,pk_l:pk_l,k_l:kl}

; compute sigma8
r8=8.                      ; [h^-1 Mpc]
w8=3./(kl*r8)^2*( sin(kl*r8)/(kl*r8) - cos(kl*r8) )
sig8=sqrt( 1./(2.*!pi^2) *$
   int_tabulated(alog(kl),kl^3*pk.pk_l*w8^2) )

return,sig8
end
