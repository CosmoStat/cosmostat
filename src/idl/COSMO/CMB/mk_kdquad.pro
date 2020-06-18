function mk_kdquad, nside=nside, muk=muk

; Written Jan 2011 by Anais Rassat
; PURPOSE: Make a map of the kinetic Doppler quadrupole
; See equation in section 2.4 of Copi et al 2006, and email from John
; Peacock
; INPUT: nside (optional, default = 64)
; OUTPUT: map of quadrupole [unitless]. I.e. multiply output by
; 2.725*1d6 to have the map in units of microK (muK)

if not keyword_set(nside) then nside=64

; ========================================
; Usefule numbers
; ========================================
  v = 370.d                     ; km/s - velocity of dipole
  c = 299792.46                 ; km/s - speed of light


; ========================================
; Galactic coordinates of Dipole direction
; ========================================
ldip = 263.85d ; degrees
bdip = 48.25d  ; degrees
thetadip = 90.d - bdip ; translate b into an angle theta in range [0,180], where 0 corresponds to the north pole

; ========================================
; Create Map 
; ========================================
  kdquad = dblarr(nside2npix(nside))

; ========================================
; Calculate dT/T value for each angle on the sky
; ========================================
  for i = 0L, n_elements(kdquad)-1L do begin
     pix2ang_nest_oldv2p20, nside, i, theta, phi
     theta = theta * 180.d / !dpi ; translate theta from rad to deg
     phi = phi * 180.d / !dpi ; translate phi from rad to deg

     cosangle = dotproduct_sphere(phi, theta, ldip, thetadip, degrees=1)
 ; calculate cosine of angle between position on sky and direction of dipole

     dt_t  = (v/c)^2 * (cosangle^2 - 1.d/3.d)

     kdquad(i) = dt_t
  endfor

if keyword_set(muk) then kdquad = kdquad * 2.725*1d6
return, kdquad

end
