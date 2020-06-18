function hubble_a,cosmo,a
   
; May 06 - Modified by AA to be double precision
; May 05 - output turned into a structure by AR
; Jan 05 - Written by A. Refregier
;
; PURPOSE: compute the hubble constant H(a) as a function of the
; expansion parameter a and related basic evolution parameters. This routine
; is meant to centralise all the basic dark energy evolution dependence.
; INPUT: a: expansion parameter
;        cosmo: cosmological parameter structure (This is the
;        cosmo.cont strucutre created by mk_cosmo)
; OUTPUT: hubble_a: structure containing:
; 		H(a): Hubble constant [km/s/Mpc]
;		w_a: dark energy equation of state parameter w(a) [1]
;             	qevol: evolution of dark energy density 
;                                  rho_l(a)/rho_l(0)   [1]
;             	h2evol: (H(a)/H_0)^2  [1]
;            	dlnh_dlna:  d(ln H)/d(ln a)  [1]

w_a=cosmo.w0+(1.d -a)*cosmo.wa    ; w(a) in Linder-Polarsky parametrisation

qevol=a^(-3.d*(1.d +cosmo.w0+cosmo.wa))*exp(-3.d*cosmo.wa*(1.d - a))  
                                                    ; rho_l(a)/rho_l(0)

h2evol=cosmo.omega_m*a^(-3.d)+cosmo.omega_l*qevol+cosmo.omega_k*a^(-2.d)
                           ; (H(a)/H_0)^2  [1]

hc=cosmo.h0*sqrt(h2evol)   ; H(a) [km/s/Mpc]   

dlnh_dlna=-.5d/h2evol*(3.d*cosmo.omega_m*a^(-3.d)+$
               3.d*cosmo.omega_l*qevol*(1.+w_a)-2.d*cosmo.omega_k*a^(-2.d))
                                 ; d(ln H)/d(ln a)  [1]

hubble={w_a:w_a,hc:hc,qevol:qevol,h2evol:h2evol,dlnh_dlna:dlnh_dlna}

return,hubble
end

