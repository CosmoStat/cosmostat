;function tk,k,cosmo,fit_tk=fit_tk
;
;
; Nov 07 - Updated by AR to compute the Eisenstein & Hu transfer function
;          for a treatement of the baryons
; March 99 - Written by A. Refregier
;
; PURPOSE: compute the (unnormalised) linear transfer function T(k) for a
; given cosmological model and including baryonic corections. This is computed ; using the analytical fits of Eisenstein and Hu (1998, ApJ, 496, 605)
; INPUT: k: comoving wave number [h Mpc^-1]
;        cosmo: cosmological parameter structure
; OPTIONAL INPUT: fit_tk: fitting function to use: 0: E&H without wiggles (default), 
;                         1:E&H with wiggles, 2: BBKS as summarised 
;                           by Peacock & Dodds (1997) 
; OUTPUT: T(k) from Eisenstein & Hu (1998)


;-----------------

function t0_tilde,q,alpha,beta  
; intermediate function To_tilde(q,alpha,beta) in Eq. 19 in E&H (1998)

c=14.2d/alpha+386.d/(1.d0+69.9d*q^1.08)
ee=exp(1.)

return,alog(ee+1.8d*beta*q)/(alog(ee+1.8d*beta*q)+c*q^2)
end

;------------------

function tk,k,cosmo,fit_tk=fit_tk

; declarations
if not keyword_set(fit_tk) then begin
	fit_tk=0 ; choice of fitting function
;		print, 'Default option has changed to Eisenstein and Hu without wiggles.'
;stop
endif

; compute transfer functions for requested fitting function
if fit_tk eq 2 then begin   ; use P&D
;print, 'BBKS/P&D'
  ; transfer function from Peacock and Dodds (1996, mnras, 280, L19)
  ; as summarized by Peacock (1997, mnras, 284, 885)
  q_pd=k/cosmo.gamma
  tk_out=alog(1.+2.34*q_pd)/(2.34*q_pd)*$
     ( 1.+3.89*q_pd+(16.1*q_pd)^2+(5.46*q_pd)^3+(6.71*q_pd)^4 )^(-.25)

endif else begin       ; use E&H

  ; declarations
  t_cmb=2.728   ; [K] from Fixsen et al. (1996)
  theta_27=t_cmb/2.7d
  omm=cosmo.omega_m
  omb=cosmo.omega_b
  h=cosmo.h

  ; compute various scales
  kk=k*h                 ; k in [Mpc^-1]
  z_eq=2.50d4*omm*h^2*theta_27^(-4)  ; matter/radiation equ. [1]
  k_eq=7.46d-2*omm*h^2/theta_27^2    ; [Mpc^-1]
  b1=0.313d*(omm*h^2)^(-0.419)*(1.d0+0.607d*(omm*h^2)^0.674)
  b2=0.238d*(omm*h^2)^0.223
  z_d=1291.d*(omm*h^2)^0.251d$       ; decoupling
      /(1.d0+0.659d*(omm*h^2)^0.828)*(1.d0+b1*(omb*h^2)^b2)
  r_eq=31.5d*omb*h^2/theta_27^4*1e3/z_eq
  r_d=31.5d*omb*h^2/theta_27^4*1e3/z_d
  s=2.d/(3.d*k_eq)*sqrt(6.d/r_eq)*$     ; sound horizon
           alog((sqrt(1.d0+r_d)+sqrt(r_d+r_eq))/(1.d0+sqrt(r_eq)))
  k_silk=1.6d*(omb*h^2)^0.52*(omm*h^2)^0.73*$
           (1.d0+(10.4d*omm*h^2)^(-0.95)) ; [Mpc^-1]

  if fit_tk eq 1 then begin  ; E&H with wiggles

    ; transfer functions parameters
    a1=(46.9d*omm*h^2)^0.670*(1.d0+(32.1d*omm*h^2)^(-0.532))
    a2=(12.0d*omm*h^2)^0.424*(1.d0+(45.0d*omm*h^2)^(-0.582))
    alpha_c=a1^(-omb/omm)*a2^(-(omb/omm)^3)
    bb1=0.944d/(1.d0+(458.d*omm*h^2)^(-0.708))
    bb2=(0.395*omm*h^2)^(-0.0266)
    beta_c=1.d/(1.d0+bb1*((cosmo.omega_dm/omm)^bb2-1.d0))

    ; compute transfer function
    q=kk/(13.41*k_eq)   ; dimensionless wave number
    ; CDM part of T(k)
    f=1.d/(1.d0+(kk*s/5.4d)^4)  ; interpolation factor
    tk_c=f*t0_tilde(q,1.d,beta_c)+(1.-f)*t0_tilde(q,alpha_c,beta_c)
    ; Baryonic part of T(k)
    beta_mode=8.41d*(omm*h^2)^0.435
    s_tilde=s/(1.d0+(beta_mode/(kk*s))^3)^(1./3.)
    y=(1.d0+z_eq)/(1.d0+z_d)
    gy=y*(-6.d*sqrt(1.d0+y)+(2.d0+3.d*y)*$
          alog((sqrt(1.d0+y)+1.d)/(sqrt(1.d0+y)-1.d)))
    alpha_b=2.07d*k_eq*s*(1.d0+r_d)^(-3./4.)*gy
    beta_b=0.5d0+omb/omm+(3.d0-2.d*omb/omm)*sqrt((17.2*omm*h^2)^2+1.d) 
    tk_b=(t0_tilde(q,1.d,1.d)/(1.d0+(kk*s/5.2d)^2)+$
          alpha_b/(1.d0+(beta_b/kk/s)^3)*exp(-(kk/k_silk)^1.4))*$
          sin(kk*s_tilde)/(kk*s_tilde)
    ; total transfer function
    tk_out=omb/omm*tk_b+cosmo.omega_dm/omm*tk_c

  endif else begin  ; E&H without wiggles

    ee=exp(1.)
    alpha_gamma=1.d0-0.328d*alog(431.d*omm*h^2)*omb/omm+$
                0.38d*alog(22.3d*omm*h^2)*(omb/omm)^2
    gamma_eff=omm*h*(alpha_gamma+(1.d0-alpha_gamma)/(1d0+(0.43d*kk*s)^4))
    q=k*theta_27^2/gamma_eff    ; rescaled wave number
    c0=14.2d0+731.d/(1.d0+62.5d*q)
    l0=alog(2.d*ee+1.8d*q)
    tk_out=l0/(l0+c0*q^2)
	
  endelse   ; fit_tk=0 or 1?
  
endelse  ; fit tk=2?

return,tk_out
end
