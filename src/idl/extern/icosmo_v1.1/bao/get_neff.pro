function get_neff, cosmo, sv
;Sept 08 modified by AnR.: fixed bug found by AnR. 
;Jul 08 modified by AA use the volume in cosmo structure
;Written by Anais Rassat, July 2008
;PURPOSE: Computes the effective number of galaxies per Mpc^3, for a
;given survey configuration (i.e. number of galaxies per arcmin^2 and
;redshift slices). 
;INPUT: cosmo
;     : sv

;if (sv.dndztype ne 'hist') then stop,'Only valid for hist setting'

neff = sv.n_g

dradeg = 180.d0/!dpi ;to convert from radians to degrees in double precision
sterad_to_arcmin2 = (180.d0*60.d0/!dpi)^2

;convert n_g from galaxies per arcmin^2/deltaz to galaxies per steradian/deltaz
neff = neff*sterad_to_arcmin2

;convert to galaxies per bin
neff = neff*4.*!dpi*sv.f_sky

;calculate volume of each bin
Vmax = interpol(cosmo.evol.vol, cosmo.evol.z, sv.z_max)
Vmin = interpol(cosmo.evol.Vol, cosmo.evol.z, sv.z_min)

v_bin = Vmax-Vmin

;calculate number of galaxies per Mpc^3
neff = neff/v_bin
return, neff
End


;pro get_narcmin,cosmo,  neff, sv
;;Get number of galaxies per arcmin, (inverse of get_neff)
;
;;calculate volume of each bin
;Vmax = interpol(cosmo.evol.vol, cosmo.evol.z, sv.z_max)
;Vmin = interpol(cosmo.evol.Vol, cosmo.evol.z, sv.z_min)
;
;v_bin = Vmax-Vmin
;
;;calculate number of galaxies per redshift bin
;neff = neff*v_bin
;
;return, narcmin
;end
