function mk_sh, cosmo_param

;Jul 08 - modified by AA - made compatible with call in mk_cosmo
;Jul 08 - modified by AA - name change 
;Written by Anais Rassat, July 2008
; ***************// HELP BEGIN //**************
; PURPOSE: Computes the sound horizon for a given cosmology
;          Uses equation (A21) of Parkinson et al. 2007 (MNRAS, Volume
;          377, Issue 1,;pp. 185-197) or equation (A22) in astro-ph version.
; INPUT: cosmo structure output from cosmo = mk_cosmo(fid) (icosmo_v0.11 and above)
; OPTIONAL INPUT: None
; OUTPUT: sound horizon in Mpc.
; OPTIONAL OUTPUT: None
; ***************// HELP END //**************

w_b = cosmo_param.omega_b*cosmo_param.h^2
w_m = cosmo_param.omega_m*cosmo_param.h^2

eta_nu = 0.68130d0 ;ratio of the energy density in neutrinos 
                   ;to the energy in photons

a_eq=1.d0/(24185.d0*(1.6813d/(1.d0+eta_nu))*w_m)
g1=0.0780d*w_b^(-0.2380d)*(1.d0+39.5d0*w_b^0.7630d0)^(-1)
g2=0.560d*(1.d0+21.1*w_b^1.81d0)^(-1)

r_a_const=30496.d*cosmo_param.omega_b*cosmo_param.h^2

;redshift of recomination
z_r= 1048.d*(1.d0+0.00124d*w_b^(-0.7380d))*(1.d0+g1*w_m^g2)
R_zr = r_a_const/(1.d0+z_r) ;Equation (A22) of Parkinson et al. 
                          ;MNRAS, Volume 377, Issue 1, pp. 185-197

;redshift of equality
z_eq = 1.d0/a_eq - 1.d0
R_eq= r_a_const/(1.d0+z_eq)

Frac = (sqrt(1.d0+R_zr)+sqrt(R_zr+R_eq))/(1.d0+sqrt(R_eq))

sh = 4000.d0/sqrt(w_b)*sqrt(a_eq)/sqrt(1.d0+eta_nu)$
     *alog(frac)                ; in units of Mpc

return, sh
end
