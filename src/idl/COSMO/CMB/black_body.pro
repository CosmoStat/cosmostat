
;=============================================


pro plot_dust, nu=nu
T_Dust=15. 
units='GhZ'
nuGHz = findgen(10000) / 10000. * 3500. + 10.
c=2.99792458e10                 ; cm / s

units='Microns'
numicron = c*10^(-2.) / (nuGHz*10.^9.) * 10.^6 
info, numicron
beta=2.
if keyword_set(nu) then numicron = nu
P2 = planck_function(T_Dust, numicron, 0, /MJy, units=units) * (numicron *1.e-6)^(-beta)
if not keyword_set(nu) then plot, nuGHz,  (p2+1)  $
else print, p2
end

;=============================================

pro plot_cmb
T_cmb=2.726
T_Dust=15. 
units='GhZ'
nuGHz = findgen(10000) / 10000. * 1000. + 10.
c=2.99792458e10                 ; cm / s

units='Microns'
numicron = c*10^(-2.) / (nuGHz*10.^9.) * 10.^6 
beta=2.
; P2 = planck_function(T_cmb, numicron, 0, /MJy, units=units)  ; * (numicron *1.e-6)^(-beta)
P2 = planck_function(T_cmb, nuGHz, 0, /MJy, units='GhZ')  
;info, p2-p3
plot, nuGHz, p2
end

;=============================================

function bb_flux_at_lambda, Temperature, lambda, beta=beta, dustmodel=dustmodel
; Return the flux at a given wavelenth for a black body at
; a given temperature 
; Flux in W / (m^3.sr)
; Ex: T_cmb=2.726
;     F_cmb = bb_flux_at_lambda(T_cmb)
; if beta is different from zero, then it 
; corresponds to the model dust emission:
;      F_lambda =  B_lambda(T) lambda^-beta
T = Temperature
if not keyword_set(beta) then begin
   beta=0
   if keyword_set(dustmodel) then begin
     p1 = 0.4
     p2 = 0.008
    beta = 1. / (P1 + P2*Temperature)
  end
end
;units='GhZ'
units='Microns'
P2 = planck_function(T, lambda, 0, /MJy, units=units) * (lambda*1.e-6)^(-beta)
return, p2 
end

;=============================================

function bb_wien_flux_at_lambda, T, lambda, beta=beta
; Return the flux  in W / (m^3.sr) at a given wavelenth for a black body at
; a given temperature using the Wien approximation
if not keyword_set(beta) then beta=0
h = 6.62617e-34  ; J.s (constante de Planck
c = 299792458.     ;  m/s (vitesse de la lumière)
c1 =  3.7418310e-16 ; = 2 PI h c^2  ( W.m2  )
c2 = 1.43879e-2 ; hc /k  (m.K)
L = 1/ !PI * c1 * lambda^(-5) / (exp(c2 / (lambda*T))) * lambda^(-beta)
return, L  
end

;=============================================

function bb_get_temp, lambda, Flux, beta=beta
; Return the temperature (in K) of a blackbody observed at the
; wavelength lambda, with a given flux (in W / (m^3.sr) )
L = flux
if not keyword_set(beta) then beta=0
L_angstrom = lambda * 10.^10
c1 =  3.7418310e-16 ; = 2 PI h c^2  ( W.m2  )
c2 = 1.43879e-2 ; hc /k  (m.K)

Val =  c1 * lambda^(-(5+beta)) / (L*!PI) + 1.
LogVal = alog(Val)
T = c2 / (LogVal * lambda)
return, T
end
;=============================================

function bb_get_temp_from_two_obs, lambda1, Lambda2,  Flux1, Flux2, beta=beta
; Return the temperature of a blackbody observed at two 
; wavelengths lambda1 and Lambda2
; As only ratio lambda1/Lambda2 and Flux1/Flux2, unit are not important 
L1 = Flux1
L2 = Flux2
if not keyword_set(beta) then beta=0
c1 =  3.7418310e-16 ; = 2 PI h c^2  ( W.m2  )
c2 = 1.43879e-2 ; hc /k  (m.K)

R = L1 / L2
LL =  (lambda1/ Lambda2)^(-(5+beta))
R = alog(R / LL)
T = c2 * (1./lambda2 - 1./ lambda1) / R
return, T
end

;=============================================

function bb_get_beta, lambda, Temperature, Flux
; Return the spectral index beta in the emission model (dust)
;      F_lambda =  B_lambda(T) lambda^-beta
T = Temperature
L = Flux
if not keyword_set(beta) then beta=0
L_angstrom = lambda * 10.^10
c1 =  3.7418310e-16 ; = 2 PI h c^2  ( W.m2  )
c2 = 1.43879e-2 ; hc /k  (m.K)

den = !PI * ( exp(c2/(lambda*T))-1.)
Val =  den / c1 * L  
LogVal = alog(Val)
beta = -5. - LogVal / alog(lambda)
return, beta
end

;=============================================

function bb_get_beta_from_two_obs, lambda1, lambda2, Temperature, Flux1, Flux2
; Return the spectral index beta in the emission model (dust) observed at two 
; wavelengths lambda1 and Lambda2
; As only ratio lambda1/Lambda2 and Flux1/Flux2, unit are not important   

T = Temperature
L1 = Flux1
L2 = Flux2
c1 =  3.7418310e-16 ; = 2 PI h c^2  ( W.m2  )
c2 = 1.43879e-2 ; hc /k  (m.K)

R = L1 / L2
e1 = exp(c2/(lambda1*T))-1.
e2 = exp(c2/(lambda2*T))-1.
Val = R * e1 / e2
LogVal = alog(Val)
beta = -5. - LogVal / alog(lambda1/lambda2)
return, beta
end

;=============================================

function bb_get_flux_from_other_obs, lambda1, lambda2, Temperature, Flux2, beta=beta
; Return the flux at wavelength lambda1 from an observation at lambda2 with flux Flux2
; for a component at Temperature
if not keyword_set(beta) then beta=0
T = Temperature
ind = where(t lt 5, count)
if count GT 0 then t[ind]= 5
c1 =  3.7418310e-16 ; = 2 PI h c^2  ( W.m2  )
c2 = 1.43879e-2 ; hc /k  (m.K)

LL =  (lambda1/ Lambda2)^(-(5+beta))
e1 = exp(c2/(lambda1*T))-1.
e2 = exp(c2/(lambda2*T))-1.

F = LL * e2 / e1 *  Flux2
return, f
end

;=============================================

function planck_function, Ttemp, nu_or_lambda, dBdT, units=units, mjy=mjy, wcm2=wcm2

;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     PLANCK returns the spectral radiance of a blackbody.
;
;DESCRIPTION:  
;    IDL function to return the spectral radiance of a blackbody,
;    i.e. the Planck curve, in units of either MJy/steradian (I_nu)
;    or watts/cm^2/steradian (nu_I_nu).
;    The blackbody temperature and either frequency (in icm or GHz)
;    or wavelength (in microns) are inputs to the function.  The
;    routine also optionally returns the derivative with respect to 
;    temperature, in units of MJy/sr/K or W/cm^2/sr/K.
;
;CALLING SEQUENCE:  
;     RESULT = PLANCK (temperature, nu_or_lambda [,dBdT] $
;              [,UNITS=units], [/MJY], [/WCM2])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     RESULT        O   flt [arr]  Spectral radiance at each wavelength. 
;                                  Units: W/cm^2/sr/K if /WCM2 specified
;                                         MJy/sr      if /MJY specfied
;     TEMPERATURE   I   flt        Temperature of blackbody, in K.
;     NU_OR_LAMBDA  I   flt        Frequency or wavelength at which to 
;                                  calculate spectrum. Units are as 
;                                  specified with UNITS keyword.
;     dBdT         [O]  flt [arr]  Derivative of Planck with respect to 
;                                  temperature. 
;     UNITS        [I]  str        'Microns', 'icm', or 'GHz' to 
;                                  identify units of NU_OR_LAMBDA. Only 
;                                  first character is required.  If 
;                                  left out, default is 'microns'.
;     /MJY          I   key        Sets output units to MJy/sr
;     /WCM2         I   key        Sets output units to W/cm^2/sr
;
;WARNINGS:
;     1.  One of /MJY or /WCM2 MUST be specified.  
;     2.  Routine gives incorrect results for T < 1 microKelvin and
;            wavelengths shortward of 1.e-10 microns.  (So sue me).
;
;EXAMPLE:
;     To produce a 35 K spectrum in MJy/sr at 2, 4, 6, 8, 10 microns:
;
;       wavelength = 2. + 2.*findgen(5)
;       temp = 35.
;       blackbody = planck(temp, wavelength, units='micron', /mjy)
;
;     One could also get back the derivative by including it in the
;     call:
;       blackbody = planck(temp, wavelength, deriv, units='m', /mjy)
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Identifies units using the UNITS keyword, then converts the 
;     supplied independent variable into microns to evaluate the 
;     Planck function.  Uses Rayleigh-Jeans and Wien approximations 
;     for the low- and high-frequency end, respectively.  Reference: 
;     Allen, Astrophysical Quantities, for the Planck formula.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  17 April 1991;     UNITS        [I]  str        'Microns', 'icm', or 'GHz' to 
;                                  identify units of NU_OR_LAMBDA. Only 
;                                  first character is required.  If 
;    Revised by RBI 27 Jan 1992 to use updated fundamental constants 
;         (SPR 9449)
;    Revised by RBI 29 Jan 1992 to calculate derivatives only when 
;         necessary
;    Revised by RBI 3 Feb 1992 to redimension output to a scalar if only 
;       a single wavelength is supplied  (SPR 9459)
;    Revised by RBI 6 Mar 92 to return either MJy/sr or (!) W/cm^2/sr
;    Revised by RBI 1 Jun 92 to fix single-wavelength failure when no
;       derivative is requested (SPR 9738), and to use MESSAGE.
;    RBI corrected error in derivative calculation SPR 9817 (17 Jul 92)
;    RBI corrected error in Wien and RJ tails SPR 10392 (24 Dec 92)
;	 but didn't get it quite right (Piper/Kryszak, 28-Dec-92)
;    J.P. Bernard Changed the name from planck.pro into planck_function.pro, to
;        avoid confusion with Astron library planck.pro
;
; SPR 9616
;.TITLE
; Routine PLANCK
;-
;
; Check on input parameters
;

IF n_elements(nu_or_lambda) lt 1 or n_elements(Ttemp) lt 1 then BEGIN
    message,'CALLING SEQUENCE: bbflux = planck_function(temp,wavelength[,dbdt=][,units=][,/mjy][,/wcm2])',/info
    val=0.
    goto,sortie
ENDIF

if not keyword_set(mjy) and not keyword_set(wcm2) then $
     message, 'Either /MJy or /Wcm2 must be specified!'
if keyword_set(mjy) and keyword_set(wcm2) then $
     message, 'Only one of /MJy or /Wcm2 may be specified!'
;
makederiv = n_params() eq 3           ; see whether dBdT is requested
if n_elements(units) lt 1 then begin
   units = 'm'
   message, /continue, "Wavelength units are assumed to be microns."
endif
T = float(ttemp) > 1.e-06                ;force temperature to be real*4
;
; Define some necessary constants
;
c = double(299792.458d0)                        ;speed of light, Physics Today Aug 90
hck = double(14387.69d0 )                        ;h*c/k             "      "
thcc = double(1.1910439d0)                       ;2*h/c^2           "     "  
coeff = thcc/hck^4 * 1.e04               ;Stephan-Boltzmann * 15/pi^5
;
; Convert nu_or_lambda into lambda (in microns) depending on specified units
;
units0 = strupcase (strmid (units,0,1))
case units0 of
   'M':   lambda = nu_or_lambda > 1.e-10                ; microns
   'I':   lambda = 1.e04 / (nu_or_lambda > 1.e-10)      ; icm, avoiding nu=0.
   'G':   lambda = c / (nu_or_lambda > 1.e-10)          ; GHz, avoiding nu=0.
   else:  message, "Units must be specified as 'microns', 'icm', or 'GHz'"
endcase
; 
;  Variable fhz is a scale factor used to go from units of nu_I_nu units to 
;  MJy/sr if the keyword is set.  In that case, its value is
;  frequency in Hz * w/cm^2/Hz ==> MJy
;
if keyword_set(wcm2) then fhz = 1. + fltarr(n_elements(lambda))
if keyword_set(mjy) then fhz = c/lambda * 1.e-15         
;
;  Introduce dimensionless variable chi, used to check whether we are on 
;  Wien or Rayleigh Jeans tails
;
chi = hck / lambda / T
val = fltarr(n_elements(chi))
if makederiv then dBdT = fltarr(n_elements(chi))
;
;  Start on Rayleigh Jeans side
;
rj = where (chi lt 0.001)
if rj(0) ne -1 then begin
    val(rj) = coeff * T^4 * chi(rj)^3 / fhz(rj)
    if makederiv then dBdT(rj) = val(rj) / T
endif
;
;  Now do nonapproximate part
;
exact = where (chi ge 0.001 and chi le 50)
if exact(0) ne -1 then begin
    chi_ex = chi(exact)
    val(exact) = double(coeff) * T^4 * chi_ex^4 / (exp(chi_ex) - 1.) / fhz(exact)
    if makederiv then dBdT(exact) = $
	val(exact) * chi_ex / T / (1. - exp(-chi_ex))
endif
;
;  ...and finally the Wien tail
;
wien = where (chi gt 50.)
if wien(0) ne -1 then begin
    chi_wn = chi(wien)
    val(wien) = coeff * T^4 * chi_wn^4 * exp(-chi_wn) / fhz(wien)
    if makederiv then dBdT(wien) = $
	val(wien) * chi_wn / T / (1. - exp(-chi_wn))
endif
;
;  Redimension to a scalar if only 1 wavelength supplied, then return
;
if n_elements(nu_or_lambda) eq 1 then begin
    if makederiv then dBdT = dBdT(0)
    val = val(0)
endif

sortie:

return, val
;
end


;=============================================

function planck_ratio, lambda1, lambda2, Temperature, beta=beta, units=units
; Return the flux at wavelength lambda1 from an observation at lambda2 with flux Flux2
; for a component at Temperature
;     UNITS        [I]  str        'Microns', 'icm', or 'GHz' to 
;                                  identify units of NU_OR_LAMBDA. Only 
;                                  first character is required.  If 

if not keyword_set(beta) then beta=0
if not keyword_set(beta) then units='GhZ'

T = Temperature
LL =  (lambda1/ Lambda2)^(-beta)
P1 = planck_function(T, lambda1, 0, /mjy, units=units)
P2 = planck_function(T, lambda2, 0, /mjy, units=units)

F = P1/ P2 * LL
return, f
end

;=============================================
