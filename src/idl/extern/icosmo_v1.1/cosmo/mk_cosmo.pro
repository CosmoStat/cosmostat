function mk_cosmo,fiducial,step_param=step_param,step_size=step_size,nopk=nopk

; Jul 08 modified by AA to include volume and volume element (from ARe)
; Jul 08 modified by AA add calculation of sound horizon  
; Jul 08 modified by AA switch on look back time calculations (should
; not do this in speed mode) !!!!!
; July 08 modified by AA to include comoving radius in [Mpc]
; July 08 modified by AA to merge with mk_evol
; Written by Adam Amara May 2008.  
; Replaces an earlier routine (rd_cosmo) which was originally written
; by Alex Refregier 03/1999
; ***************// HELP BEGIN //**************
; PURPOSE: This routine creates a structure with all the relevant
; static (e.g. without z evolution) comological quantities. This 
;
; INPUTS: fiducial: structure containing the core cosmology parameters
; OPTIONAL INPUT: step_param: name of parameter which will be used to
;                   step away from the fiducial values. This keyword
;                   allows the calculation of an cosmo structure for a
;                   cosmology model where one of the parameters is a
;                   small step  away from the one set in the fiducial
;                   structure.  This is useful when doing a Fisher
;                   matrix calculation. See the routine
;                   set_cosmo_fiducial for a list of possible tags.
;                 step_size: this sets the step size of the step taken
;                   by the step_param parameter, away from its
;                   fiducial value.          
; OUTPUT: returns a structure with the static cosmology parameters.
; ----
; Example 1:cosmo = mk_cosmo(fiducial)
; This will take in the structure fiducial (which is produced using
; the function set_fiducial) and returns a cosmo structure.
;
; Example 2: cosmo = mk_cosmo(fiducial,step_param='omega_m')
; This will produce a cosmo structure where the value of Omega matter
; is different from that of the fiducial model (all other parameter
; use the fiducial values). step_size will use the default value of
; step_size=0.006d.
;
; Example 3: cosmo = mk_cosmo(fiducial,step_param='w0', step_size=0.01d)
; This will produce a cosmo structure where the value of w0 is
; different from that of the fiducial model (all other parameter 
; use the fiducial values). step_size will is set to 0.01d
; ----
; ***************// HELP END //**************


;Option to use fiducial value or take a step away from fiducial value 
if keyword_set(step_param) then cosmo_param=mk_cosmo_step(fiducial.cosmo,step_param,step_size) else cosmo_param=fiducial.cosmo
;NEED TO MAKE SURE THAT THE KEYWORDS STRUCTURE IS CARRIED THROUGH TO
;MK_COSMO_STEP !!!

;***************** Begin: set up constant Structure *********************  
c = 299792458.d                 ; Speed of light [m s^-1]  
tcmb = 2.726d *1d6              ; CMB temperature [micro Kelvin]
h0 = 100.0d *cosmo_param.h      ; Hubble constant [Km/s/Mpc]
rh = c/1.d3/h0                  ; Hubble radius
omega_dm = cosmo_param.omega_m-cosmo_param.omega_b ; Dark matter density
omega=cosmo_param.omega_m+cosmo_param.omega_l      ; Total density
omega_k = 1.0d - omega                             ; Curvature
gamma=cosmo_param.omega_m*cosmo_param.h*exp(-cosmo_param.omega_b*(1.d + sqrt(2.d*cosmo_param.h)/cosmo_param.omega_m)) ;Sugiyama (1995, APJS, 100, 281)

norm=1. ; This flag decides whether to use sigma8 or CMB normalisation. for the moment sigma8 is used. !! need to change this!!

norma=-1.d   ; normalisation of the power spectrum; set to -1. if not normalised yet

; useful factors
  case 1 of
     omega gt 1.: begin         ; closed
        k=1.d                         
        sqrtk=sqrt(omega-1.d)   ; useful factor
     end                        ; with  a=r/r0
     omega eq 1.: begin
        k=0.d                   ; flat
        sqrtk=1.d
     end
     omega lt 1.: begin         ; open
        k=-1.d   
        sqrtk=sqrt(1.d -omega)       
     end
  endcase
  r0=rh/sqrtk                   ; scale radius at z=0 [Mpc] 
  
;compute the sound horizon

sh = mk_sh(cosmo_param)                ; in units of Mpc

; construct cosmo structure
const={h:cosmo_param.h,h0:h0,rh:rh,omega_dm:omega_dm,omega_b:cosmo_param.omega_b,omega_m:cosmo_param.omega_m,omega_l:cosmo_param.omega_l,omega:omega,omega_k:omega_k,w0:cosmo_param.w0,wa:cosmo_param.wa,k:k,sqrtk:sqrtk,r0:r0,n:cosmo_param.n,sigma8:cosmo_param.sigma8,tau:cosmo_param.tau,gamma:gamma,c:c,tcmb:tcmb,norm:norm,norma:norma,sh:sh}
;***************** End: set up constant Structure *********************  

;***************** Begin: set up evolve Structure *********************  
; ---- START: Keyword declarations ----
void=tag_check(fiducial.calc,'nz_fn',val=nz_fn) 
void=tag_check(fiducial.calc,'ran_z',val=ran_z)
; ---- END: Keyword declarations ----

; redshift and expansion parameter a=R/Ro=1/(1+z)   [1]
;stop
z=xgen(ran_z(0),ran_z(1),np=nz_fn,/doub)
a=1.d/(1.d +z)

; hubble constant and dark energy quantities
hubble=hubble_a(const,a)  ; structure contains:
	      ; hc: hubble constant H [km/s/Mpc]
              ; qevol: dark energy density evolution rho_l(a)/rho_l(0) [1]
              ; w_a: dark energy equation of state parameter w(a)   [1]
w_a=hubble.w_a

; evolution of the density parameters
h_h0=hubble.hc/const.h0    ; H/H_0
omega_m_a=const.omega_m*a^(-3.d)*h_h0^(-2.d)
omega_l_a=const.omega_l*hubble.qevol*h_h0^(-2.d)
omega_a=omega_m_a+omega_l_a
omega_k_a=1.d -omega_a

; comoving distance in units of R_0  [1]
chi=chi_a(const,a)

; look back time (commented out to save computing time)
if keyword_set(fiducial.calc.speed) then begin
th0=-1 & t=-1
endif else begin
   th0=t_a(const,a)                         ; t*Ho [1]
   t=th0/const.h0*3.085678e19/3.15581498e7/1e9 ; t [Gyr]
endelse
; sk(chi): radial comoving distance in units of R_0  [1],
case 1 of
  const.omega eq 1.: sk=chi
  const.omega gt 1.: sk=sin(chi)
  const.omega lt 1.: sk=sinh(chi)
endcase

; angular-diameter and luninosity distance  [Mpc]
da=const.r0*sk/(1.d +z)
dl=da*(1.d +z)^2.d

; compute the linear growth factor by numerical ODE solving
d=growth_a_ode(const,a)           ; normalised to D(z=0)=1 
d_u=growth_a_ode(const,a,/unnorm) ; normalised to D(a)=a in om=om_m=1 model

;dzdr = 1/const.rh*hubble.hc/const.h0 ;test by ar
dzdr = hubble.hc/(const.c/1000.d)

;Volume universe [Mpc^3] to a redshift z. For equation of Volumes, see
;Hogg 1999 (astroph only), Equation (29). 
CASE 1 OF
   const.omega gt 1.: begin     ; closed
      vol = (4.d0*!dpi*const.rh^3/(2.d0*const.omega_k))*$
                (sk/const.rh*sqrt(1.d0+const.omega_k*sk^2/const.rh^2)-$
                 1.d/sqrt(abs(const.omega_k))*$
                 asinh(sqrt(abs(const.omega_k))*sk/const.rh))
   end                        
   const.omega lt 1.: begin     ;open
      vol = (4.d*!dpi*const.rh^3/(2.d0*const.omega_k))*$
                (sk/const.rh*sqrt(1.d0+const.omega_k*sk^2/const.rh^2)-$
                 1.d/sqrt(abs(const.omega_k))*$
                 asin(sqrt(abs(const.omega_k))*sk/const.rh))
   end
   const.omega eq 1.: begin     ;flat
      vol = 4.d0/3.d0*!dpi*(sk*const.r0)^3             
   end
ENDCASE   

; store in structure
evol={z:z,a:a,hc:hubble.hc,chi:chi,sk:sk,da:da,dl:dl,d:d,$
      w_a:w_a,omega_a:omega_a,omega_m_a:omega_m_a,$
      omega_l_a:omega_l_a,omega_k_a:omega_k_a, $
      dzdr: dzdr,th0:th0,t:t,vol:vol,ran_z:ran_z}
;***************** End: set up evolve Structure *********************  

;***************** Begin: set up Pk Structure *********************  
cosmo={const:const,evol:evol}   ; temporary structure
if keyword_set(nopk) then pk=-1 else begin
; compute power spectrum grid P(k,z)
   pk=mk_pk_grid(fiducial,cosmo)
endelse
;***************** End: set up Pk Structure *********************  

cosmo={const:cosmo.const,evol:cosmo.evol,pk:pk}
;stop
return,cosmo
end





