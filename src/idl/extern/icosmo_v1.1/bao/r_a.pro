function r_a, cosmo, z
;Written by Anais Rassat, July 2008.
;Equation (A22) of Parkinson et al. MNRAS, Volume 377, Issue 1, pp. 185-197
a = 1.d0/(1.d0+z)
r_a = 30496.d*cosmo.const.omega_b*cosmo.const.h^2*a
return, r_a
end

