;function chi_a,cosmo,a
;
; Jul 08 - modified by AA to be consistent with v0.11
; May 06 - modified by AA to be double precision
; Jan 05 - modified by AR to call hubble_a.pro
; Dec 04 - modified by AR to include varying dark energy eq. of state
; Mar 99 - Written by A. Refregier
;
; PURPOSE: compute the comoving distance Chi(a) as a function of the
; expansion parameter a, for a given cosmology. This is done by numerical
; integration
; INPUT: a: expansion parameter
;        cosmo: cosmological parameter structure
; OUTPUT: chi(a): comoving distance in unit of the hubble radius
;                 i.e. chi/R_0 is given rather than chi

function chi_a_int,ap    ; integrand for qromb.pro

common cosmology,cosmo2

;qevol=ap^(-3.*(1.+cosmo2.w_l+cosmo2.w1_l))*exp(-3*cosmo2.w1_l*(1-ap))  
;                                ; rho_l(a)/rho_l(0)
;return,-cosmo2.sqrtk/$
;  sqrt(cosmo2.omega_m*ap+cosmo2.omega_l*ap^4*qevol+$
;  cosmo2.omega_k*ap^2)

hubble=hubble_a(cosmo2,ap)
return,-cosmo2.sqrtk*cosmo2.h0/(ap^2.d*hubble.hc)

end


function chi_a,cosmo,a

common cosmology,cosmo2        ; to pass to qromb
cosmo2=cosmo
;stop
; complete unspecified parameters if needed
if not tag_exist(cosmo2,'w0') then cosmo2=create_struct(cosmo2,'w0',-1.d)
if not tag_exist(cosmo2,'wa') then cosmo2=create_struct(cosmo2,'wa',0.d)

; compute Chi(a) by numerical integration
chi=qromb('chi_a_int',1.d,a,/double)
;stop
return,chi
end
