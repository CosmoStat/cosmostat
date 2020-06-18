;function t_a,cosmo,a
;
; Jan 05 - modified by AR to call hubble_a.pro
; Dec 04   - modified by AR to include varying dark energy eq. of state
; March 99 - Written by A. Refregier
;
; PURPOSE: compute the look back time as a function of the
; expansion parameter a, for a given cosmology. This is done by numerical
; integration of t=int_1^a da'/(a'H(a')) 
; INPUT: a: expansion parameter
;        cosmo: cosmological parameter structure
; OUTPUT: t(a): look back time as a function of a in units of [Ho^-1]


function t_a_int,ap    ; integrand for qromb.pro

common cosmology,cosmo2

;qevol=ap^(-3.*(1.+cosmo2.w_l+cosmo2.w1_l))*exp(-3*cosmo2.w1_l*(1-ap))  
;                                                    ; rho_l(a)/rho_l(0)

;return,1./( ap*sqrt(cosmo2.omega_m*ap^(-3)+cosmo2.omega_l*qevol+$
;  cosmo2.omega_k*ap^(-2)))                     ; in units of [Ho^-1]
h=hubble_a(cosmo2,ap)
return,cosmo2.h0/(ap*h.hc)   ; ; in units of [Ho^-1]

end


function t_a,cosmo,a

common cosmology,cosmo2        ; to pass to qromb
cosmo2=cosmo

; complete unspecified parameters if needed
if not tag_exist(cosmo2,'w0') then cosmo2=create_struct(cosmo2,'w0',-1.)
if not tag_exist(cosmo2,'wa') then cosmo2=create_struct(cosmo2,'wa',0.)

; compute Chi(a) by numerical integration
;stop
t=-qromb('t_a_int',1.,a)

return,t
end
