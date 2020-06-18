;pro mk_cov_blake, cosmo, x,sv,bao_type
function mk_cov_blake,cosmo,sv,bao_type

;Aug 08 - modified by AnR, now check that x_t>x and stop the code if
;         this is not the case
;Aug 08 - modified by AnR, small bug reported by Tom Kitching fixed
;Jul 08 - modified by AA use volume now calculated in mk_cosmo
;Jul 08 - modified by AA make compatible with v0.11  
;Written by Anais Rassat, July 2008
;Calculates the transverse and radial errors for the BAO scale.
;INPUT:  cosmo (structure, version 0.10.1 and above)
;        x: either radial or transverse BAO scale (i.e. ouput of either
;        mk_bao_rad or mk_bao_trans
;OPTIONAL INPUT:
; sv: survey (structure, version 0.10.1 and above)
; bao: type of BAO Calculation
;     0: Spectroscopic Tangential
;     1: Spectroscopic Radial
;     2: Photometric Tangential
;OUTPUT: x = Delta y/y or  xp = Delta yp/yp
  
;TO DO:
;4) Photometric redshift errors

if not keyword_set(cosmo) then stop,'Error: must input cosmo structure - stopped in mk_cov_blake'
if not keyword_set(sv) then stop,'Error: must input survey structure - stopped in mk_cov_blake'
;if not keyword_set(bao_type) then stop,'Error: must input  bao_type - stopped in mk_cov_blake'


;Inputs from sv
  zmin = sv.z_min
  zmax = sv.z_max
  n = get_neff(cosmo, sv) ;translate galaxies /arcmin2 to galaxies / Mpc^3
  b_0 = sv.bias
;  print, 'mk_cov_blake.pro: Check output from survey.'
  
;General Coefficients taken from Blake et al. 
;MNRAS Volume 365, Issue 1, pp. 255-264
;----------------------------------------------------------
  
;----------------------------------------------------------
  If not keyword_set(bao_type) then bao = 0
  CASE BAO_type OF
;For spectro-z tangential
     0: begin
        blake_p=[0.85d0, 0.82d*1d-3, 1.4d0, 0.5d0, 0.52d0, 2.d0,7.3d0,0.26d0,0.27d0]
     end
;For sepctro-z radial
     1: begin
        blake_p=[1.48d0, 0.82*1d-3, 1.4d0, 0.5d0, 0.52d0, 2.d0, 10.6d0, 0.49d0, 1.d0]
     end
;For Photo-z tangential
     2:begin
        blake_p=[0.85d0, 0.82*1d-3, 1.4d0, 0.5d0, 0.52d0, 4.d0, 4.2d0, 0.11d0, 0.42d0]
     end
  ENDCASE
  
;Rewrite params for Blake et al. in format of paper
  x_0 = blake_p[0]                 ;per cent
  n_0 = blake_p[1]*cosmo.const.h^3 ;*1e-3 Mpc^{-3}
  z_m = blake_p[2]                 ;
  gam = blake_p[3]
  b = blake_p[4]
  
  
  
;Other terms needed for General Equation:
;Assume Flat Universe for the moment
;  V = mk_volume(cosmo, sv)
  Vmax = interpol(cosmo.evol.vol, cosmo.evol.z,sv.z_max)
  Vmin = interpol(cosmo.evol.vol, cosmo.evol.z,sv.z_min)
  V =  (vmax-vmin)*sv.f_sky
;  V=cosmo.evol.vol
  zmean = (zmin+zmax)/2.d0 ;Only true for histogram!!!

  growth = interpol(cosmo.evol.d, cosmo.evol.z, zmean)

  z_0 = 1.d0                       ;random fiducial redshift (how to chose this one?)
;  V_0 = mk_volume(cosmo, {z_min:0.6d0, z_max:1.4d0, f_sky:1.d0/40.d0}) ;last term for 1000deg^2 sky survey)
  Vmax_0 = interpol(cosmo.evol.vol, cosmo.evol.z, 1.4d)
  Vmin_0 = interpol(cosmo.evol.vol, cosmo.evol.z, 0.6d)
  fsky_0=1.d0/40.d0

  v_0 =  (vmax_0-vmin_0)*fsky_0

  D_0 = interpol(cosmo.evol.d, cosmo.evol.z, z_0)
;stop 

;Define n_eff. "Given that the amplitude of the power spectrum
;decreases with increaing k, the variation in the extent of the linear
;regime with z changes the average amplitude of P(k) included in the
;analysis, and hence the value of n_eff. (cf Equation (7) of Blake et
;al.)

  n_eff = make_array(n_elements(zmean),value=n_0)
  ind = where(zmean le z_m, count)
  if count ne 0 then n_eff[ind] = n_0*[1.d0-b*(1.d0-zmean[ind]/z_m)] 

;Equation (6a) General Equation for x for z>z_m for tangential part:
  x = x_0*sqrt(V_0/V)*(1.d0+n_eff/n*D_0^2/(b_0^2*growth^2))
;Include redshift errors if there are any;
;TODO
; if bao eq 2 then begin
;photoz=
;sig_r = 
;sig_0 =    
;    x = x_0*sqrt(V_0/V)*sqrt(sig_r/sig_0)*(1.d0+n_eff/n*D_0^2/(b_0^2*growth^2))
; endif



;---------------------------------------------------
;Equation (6b)
;---------------------------------------------------
  if count ne 0 then begin 
     x[ind] = x[ind]*(z_m/zmean[ind])^gam
  endif
;---------------------------------------------------
;Equation (9)
;---------------------------------------------------
  p = blake_p[5]
  a = blake_p[6]
  alpha=blake_p[7]
  beta =blake_p[8]
  
  x_t = a*(v/v_0)^alpha

  if count ne 0 then x_t[ind]=x_t[ind]*(z_m/zmean[ind])^beta
;---------------------------------------------------
;Equation (8)
;---------------------------------------------------

;New Test, for Equation (8) to be applied, x_t must be > than x, or
;negative errors will be returned
;Here we check that x_t > x and stop the code if this is not the case.
test = where(x_t le x, count_test)
if count_test eq 0 then begin
   x = x/(1.d0-(x/x_t)^p)
endif else begin
   print, 'WARNING: Acoustic features are not well resolved with this survey.  BAO calculation cannot be performed.'
   print, '(See section 3 of Blake et al. 2005, astro-ph/0510239)'
   print, 'WARNING: Try reducing the number of redshift bins.'
   stop
endelse

return,x  
end
