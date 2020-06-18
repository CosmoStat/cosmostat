;function growth_a_ode,cosmo,a,dlnddlna,unnorm=unnorm
;
; Jul 08 - modified by AA - cleaned and rearranged
; Jul 08 - modified by An.R. (taken from previous modification in
;           icosmo_v8 by A.R.)
; May 06 - modified by AA to be double precision
; Jan 05 - Written by A. Refregier
;
; PURPOSE: compute the linear growth factor D as a function of the
;  expansion parameter a, for a given cosmology. This is done by numerical
; integration of the ODE for the growth factor. The default
; normalisation is such that D(z=0)=1, but it can be set such that
; D(a)=a in a matter dominated universe (omega_m=omega=1)
; INPUT: a: expansion parameter
;        cosmo: cosmological parameter structure
; KEYWORD: unnorm: normalise such that D(a)=a for omega_m=omega=1
;                  universe instead default of D(z=0)=1
; OUTPUT: D(a): growth factor as a function of a [1]
;         dlnddlna: derivative of growth factor dlnD/dlnA

function growth_a_derivs,x,y    ; integrand for qromb.pro

common cosmology,cosmo2

; compute evolution quantities
a=exp(x)                    ; x=ln(a) 
hubble=hubble_a(cosmo2,a)  
; structure contains: hc: H; h2evol: (H(a)/H0)^2; dlnh_dlna: d(ln H)/d(ln a)

omega_m_a=cosmo2.omega_m*a^(-3.d)/hubble.h2evol  ; Omega_m(a)

; compute derivatives of y1=G and y2=dG/dlna
f1=y(1)                   ; f1=dy1/dx=y2=dG/dlna
f2=-(4.d +hubble.dlnh_dlna)*y(1)-(3.d +hubble.dlnh_dlna-1.5d *omega_m_a)*y(0)

return,[f1,f2]
end


function growth_a_ode,cosmo,a,dlnddlna,unnorm=unnorm

common cosmology,cosmo2       ; to pass to odeint
cosmo2=cosmo

; declarations
a_start=double(min([a,1.1d-3])) ; initial condition for integration
x_start=alog(a_start)
y_start=[1.d,0.d]               ; initial conditions: y1=G=1, y2=dG/dlna=0.
eps=.05d                        ; integration accuracy
hwant=.2d                       ; desired (i.e.) maximum step in x
h1=hwant                        ; guess for first step in x


; complete unspecified parameters if needed
if not tag_exist(cosmo2,'w_l') then cosmo2=create_struct(cosmo2,'w_l',-1.d)
if not tag_exist(cosmo2,'w1_l') then cosmo2=create_struct(cosmo2,'w1_l',0.d)

; compute G(a)=D(a)/a from a=a_start to a=1 by numerical integration for 
; the ODE system:
; x=ln(a), y1=G, y2=dG/dlna
; dy1/dx=dG/dlna=y2
; dy2/dx=d^2G/d(lna)^2=-(4+dlnH/dlna)y2-(3+dlnH/dlna-1.5*Omega_m(a))y1
; with initial conditions: y1=G=1 and y2=dG/dlna=0 for a->0
odeint,y_start,x_start,0.d,eps,h1,0.d,nok,nbad,'growth_a_derivs',$
       xp,yp,hwant=hwant

; interpolate to desired values and set normalisation
d=interpol(yp(0,*),xp,alog(a))*a             ; D(a) normalised so that D=a
                                ; in a matter dominated universe
                                ; (omega_m=omega=1)
dlnddlna=interpol(yp(1,*),xp,alog(a))*a/d+1. ; d(ln D)/d(ln a) 

if not keyword_set(unnorm) then d=d/y_start(0)

return,d
end
