pro wt_params_camb, cosmo, file=file, verbose=verbose, tk=tk

; Nov 08 - Modified by Anais Rassat - fixed bug where wa value was omitted
; Oct 08 - Modified by Anais Rassat - added sigma8 as output
; Oct  08 - Modified by Anais Rassat to be compatible with
;           Jochen's Planck Fisher matrix code
; Oct  08 - Modified by Anais Rassat to be compatible with icosmo_v1.0
;           and onwards
; June 08 - Written by Anais Rassat
;
; PURPOSE: Write a 'params.ini' file that CAMB (http://camb.info/)
; will read. The file has a unique filename which depends on the date (YYYYMMDDhhmmss_cambtemp)
; it was created unless the user specifies a filename.

;Now permits you to specify what file it will be written to
;make a filename that is unique
if not keyword_set(file) then get_filedate, file

;sort out the output file names
;output name for scalar Cl's
scl_temp   = strcompress('scalCls_'+string(file)+'.dat',/remove_all)
scl_temp   = strcompress('scalar_output_file = '+scl_temp)
;output name for tk
tk_temp = strcompress('tk_'+string(file)+'.dat',/remove_all)
tk = tk_temp
tk_temp = strcompress('transfer_filename(1)    = '+tk_temp)

if keyword_set(verbose) then print, scl_temp

;prepare what to ask CAMB to calcualte:

if keyword_set(tk) then begin
   do_tk = strcompress('get_transfer    = T')
endif else begin 
   do_tk = strcompress('get_transfer    = F') 
endelse

   
if keyword_set(verbose) then begin
   print, tk 
   print, do_tk
 endif
;prepare variables
;for 'physical parameters'
ombh2  = strcompress('ombh2          = '+string(cosmo.const.omega_b*cosmo.const.h^2))
omch2  = strcompress('omch2          = '+string(cosmo.const.omega_dm*cosmo.const.h^2))
omk    = strcompress('omk            = '+string(cosmo.const.omega_k))
;for physical=F 
om_cdm = strcompress('omega_cdm      = '+string(cosmo.const.omega_dm))
omb    = strcompress('omega_baryon   = '+string(cosmo.const.omega_b))
om_de  = strcompress('omega_lambda   = '+string(cosmo.const.omega_l))
;Others
bigh   = strcompress('hubble         = '+string(cosmo.const.h*100.))
w0     = strcompress('w0             = '+string(cosmo.const.w0))
wa     = strcompress('wa             = '+string(cosmo.const.wa))
n_s     = strcompress('scalar_spectral_index(1)  = '+string(cosmo.const.n))
transfer_kmax = strcompress('transfer_kmax           = '+string(3./cosmo.const.h))
if keyword_set(verbose) then begin
   print, om_cdm
   print, omb
   print, om_de
   print, omk
   print, bigh
   print, w0
   print, wa
endif

;write the params.ini file
openw, 3, file
printf, 3, '#Parameters for CAMB'
printf, 3,  ''
printf, 3,  ''

printf, 3, '#output_root is prefixed to output file names'
printf, 3,  'output_root = test2'

printf, 3,  '#What to do'
printf, 3,  'get_scalar_cls = T'
printf, 3,  'get_vector_cls = F'
printf, 3,  'get_tensor_cls = F'
printf, 3,   do_tk
printf, 3, ''

printf, 3, '#if do-lensing then scalar_outputfile contains additional columns of l^4C_l^{pp} and l^3C_l^{pT}'
printf, 3, '#where p is the projected potential.  Output lensed CMB Cls (without tensors) are in lensed_output_file below.'
printf, 3, 'do_lensing = F'
printf, 3, '# 0: linear, 1:non-linear matter power (HALOFIT), 2: non-linear CMB lensing (HALOFIT)'
printf, 3, 'do_nonlinear   = 0'
printf, 3, ''
printf, 3, ''

printf, 3,  '#Maximum multipole and k*eta. '
printf, 3,  '#  Note that C_ls near l_max are inaccurate (about 5%), go to 50 more than you need'
printf, 3,  '#  Lensed power spectra are computed to l_max_scalar-250 where accurate at %-level'
printf, 3,  '#  For high accuracy lensed spectra set l_max_scalar = (l you need) + 500'
printf, 3,  '#  To get accurate lensed BB need to have l_max_scalar>2000, k_eta_max_scalar > 10000'
printf, 3,  '#  Otherwise k_eta_max_scalar=2*l_max_scalar usually suffices'
printf, 3,  'l_max_scalar      = 2000'
printf, 3,  'k_eta_max_scalar  = 4000'
printf, 3, ''


printf, 3,  '#  Tensor settings should be less than or equal to the above'
printf, 3,  'l_max_tensor      = 1500'
printf, 3,   'k_eta_max_tensor  = 3000'
printf, 3, ''

printf, 3,   '#Main cosmological parameters, neutrino masses are assumed degenerate'
printf, 3,   '# If use_phyical set phyiscal densities in baryone, CDM and neutrinos + Omega_k'
printf, 3,   'use_physical   = F'
;printf, 3,   ombh2
;printf, 3,   omch2
;printf, 3,   strcompress('omnuh2         = '+string(0.))
;printf, 3,   omk
printf, 3,   bigh
printf, 3,   '#effective equation of state parameter for dark energy, assumed constant'
printf, 3,   w0
printf, 3,   wa
printf, 3,   '#constant comoving sound speed of the dark energy (1=quintessence)'
printf, 3,   'cs2_lam        = 1'
printf, 3, ''

printf, 3,   '#if use_physical = F set parameters as here'
printf, 3,   omb
printf, 3,   om_cdm
printf, 3,   om_de
printf, 3,   strcompress('omega_neutrino         = '+string(0.))

;printf, 3,   '#omega_baryon   = 0.0462'
;printf, 3,   '#omega_cdm      = 0.2538'
;printf, 3,   '#omega_lambda   = 0.7'
;printf, 3,   '#omega_neutrino = 0'
printf, 3, ''

printf, 3,   '#massless_neutrinos is the effective number (for QED + non-instantaneous decoupling)'
printf, 3,   'temp_cmb           = 2.726'
printf, 3,   'helium_fraction    = 0.24'
printf, 3,   'massless_neutrinos = 3.04'
printf, 3,   'massive_neutrinos  = 0'
printf, 3, ''

printf, 3,   '#Neutrino mass splittings'
printf, 3,   'nu_mass_eigenstates = 1'
printf, 3,   '#nu_mass_degeneracies = 0 sets nu_mass_degeneracies = massive_neutrinos'
printf, 3,   '#otherwise should be an array'
printf, 3,   '#e.g. for 3 neutrinos with 2 non-degenerate eigenstates, nu_mass_degeneracies = 2 1'
printf, 3,   'nu_mass_degeneracies = 0  '
printf, 3,   '#Fraction of total omega_nu h^2 accounted for by each eigenstate, eg. 0.5 0.5'
printf, 3,   'nu_mass_fractions = 1'
printf, 3, ''

printf, 3,   '#Initial power spectrum, amplitude, spectral index and running. Pivot k in Mpc^{-1}.'
printf, 3,   'initial_power_num         = 1'
printf, 3,   'pivot_scalar              = 0.05'
printf, 3,   'pivot_tensor              = 0.05'
printf, 3,   'scalar_amp(1)             = 2.5e-9'
printf, 3,   n_s 
printf, 3,   'scalar_nrun(1)            = 0'
printf, 3,   'tensor_spectral_index(1)  = 0'
printf, 3,   '#ratio is that of the initial tens/scal power spectrum amplitudes'
printf, 3,   'initial_ratio(1)          = 1'
printf, 3,   '#note vector modes use the scalar settings above'
printf, 3,   ''
printf, 3, ''

printf, 3,   '#Reionization, ignored unless reionization = T, re_redshift measures where x_e=0.5'
printf, 3,   'reionization         = T'
printf, 3,   ''

printf, 3,   're_use_optical_depth = T'
printf, 3,   're_optical_depth     = 0.09'
printf, 3,   '#If re_use_optical_depth = F then use following, otherwise ignored'
printf, 3,   're_redshift          = 12'
printf, 3,   '#width of reionization transition. CMBFAST model was similar to re_delta_redshift~0.5.'
printf, 3,   're_delta_redshift    = 1.5'
printf, 3,   '#re_ionization_frac=-1 sets to become fully ionized using YHe to get helium contribution'
printf, 3,   '#Otherwise x_e varies from 0 to re_ionization_frac'
printf, 3,   're_ionization_frac   = 1'
printf, 3, ''

printf, 3,   ''
printf, 3,   '#RECFAST 1.4 recombination parameters'
printf, 3,   'RECFAST_fudge = 1.14'
printf, 3,   'RECFAST_fudge_He = 0.86'
printf, 3,   'RECFAST_Heswitch = 6'
printf, 3,   ''

printf, 3,  '#Initial scalar perturbation mode (adiabatic=1, CDM iso=2, Baryon iso=3, '
printf, 3,   '# neutrino density iso =4, neutrino velocity iso = 5) '
printf, 3,   'initial_condition   = 1'
printf, 3,   '#If above is zero, use modes in the following (totally correlated) proportions'
printf, 3,   '#Note: we assume all modes have the same initial power spectrum'
printf, 3,   'initial_vector = -1 0 0 0 0'
printf, 3, ''

printf, 3,   '#For vector modes: 0 for regular (neutrino vorticity mode), 1 for magnetic'
printf, 3,   'vector_mode = 0'
printf, 3, ''

printf, 3,   '#Normalization'
printf, 3,   'COBE_normalize = F'
printf, 3,   '##CMB_outputscale scales the output Cls'
printf, 3,   '#To get MuK^2 set realistic initial amplitude (e.g. scalar_amp(1) = 2.3e-9 above) and'
printf, 3,   '#otherwise for dimensionless transfer functions set scalar_amp(1)=1 and use'
printf, 3,   '#CMB_outputscale = 1'
printf, 3,   'CMB_outputscale = 7.4311e12'
printf, 3, ''
printf, 3, 'sigma8 = '+string(cosmo.const.sigma8)

printf, 3,   '#Transfer function settings, transfer_kmax=0.5 is enough for sigma_8'
printf, 3,   '#transfer_k_per_logint=0 sets sensible non-even sampling; '
printf, 3,   '#transfer_k_per_logint=5 samples fixed spacing in log-k'
printf, 3,   '#transfer_interp_matterpower =T produces matter power in regular interpolated grid in log k; '
printf, 3,   '# use transfer_interp_matterpower =F to output calculated values (e.g. for later interpolation)'
printf, 3,   'transfer_high_precision = F'
printf, 3,   'transfer_kmax           = 4.28571'
printf, 3,   'transfer_k_per_logint   = 0'
zbin = 1
printf, 3,   strcompress('transfer_num_redshifts  = '+string(zbin))
printf, 3,   'transfer_interp_matterpower = T'


if zbin ne 1 then begin
   z = reverse(findgen(zbin)/zbin*2.) ; just a random distribution of redshift bins
   for i = 1, zbin do begin
      printf, 3, strcompress('transfer_redshift('+string(i)+')='+string(z[i-1]),/remove_all)
   endfor
endif else begin
   printf, 3,   'transfer_redshift(1)    = 0'
endelse

printf, 3,   tk_temp
printf, 3,   '#Matter power spectrum output against k/h in units of h^{-3} Mpc^3'
;output name for tk
pk_temp = strcompress('pk_'+string(file)+'.dat',/remove_all)
pk_temp = strcompress('transfer_matterpower(1)    = '+pk_temp)
printf, 3, pk_temp
;printf, 3,   'transfer_matterpower(1) = matterpower.dat'
printf, 3, ''

printf, 3,   '#Output files not produced if blank. make camb_fits to use use the FITS setting.' 
printf, 3, scl_temp
printf, 3,   'vector_output_file = vecCls.dat'
printf, 3,   'tensor_output_file = tensCls.dat'
printf, 3,   'total_output_file  = totCls.dat'
printf, 3,   'lensed_output_file = lensedCls.dat'
printf, 3,   'lensed_total_output_file  =lensedtotCls.dat'
printf, 3,   'FITS_filename      = scalCls.fits'
printf, 3, ''

printf, 3,   '##Optional parameters to control the computation speed,accuracy and feedback'
printf, 3, ''

printf, 3,   '#If feedback_level > 0 print out useful information computed about the model'
printf, 3,   'feedback_level = 1'
printf, 3, ''

printf, 3,   '# 1: curved correlation function, 2: flat correlation function, 3: inaccurate harmonic method'
printf, 3,   'lensing_method = 1'
printf, 3,   'accurate_BB = F'
printf, 3, ''


printf, 3,   '#massive_nu_approx: 0 - integrate distribution function'
printf, 3,   '#                   1 - switch to series in velocity weight once non-relativistic'
printf, 3,   '#                   2 - use fast approximate scheme (CMB only- accurate for light neutrinos)'
printf, 3,   '#                   3 - intelligently use the best accurate method'
printf, 3,   'massive_nu_approx = 3'
printf, 3, ''

printf, 3,   '#Whether you are bothered about polarization. '
printf, 3,   'accurate_polarization   = T'
printf, 3, ''

printf, 3,   '#Whether you are bothered about percent accuracy on EE from reionization'
printf, 3,   'accurate_reionization   = T'
printf, 3, ''

printf, 3,   '#whether or not to include neutrinos in the tensor evolution equations'
printf, 3,   'do_tensor_neutrinos     = F'
printf, 3, ''

printf, 3,   '#Whether to turn off small-scale late time radiation hierarchies (save time,v. accurate)'
printf, 3,   'do_late_rad_truncation   = T'
printf, 3, ''

printf, 3,   '#Computation parameters'
printf, 3,   '#if number_of_threads=0 assigned automatically'
printf, 3,   'number_of_threads       = 0'
printf, 3, ''

printf, 3,   '#Default scalar accuracy is about 0.3% (except lensed BB). '
printf, 3,   '#For 0.1%-level try accuracy_boost=2, l_accuracy_boost=2.'
printf, 3, ''


printf, 3,   '#Increase accuracy_boost to decrease time steps, use more k values,  etc.'
printf, 3,   '#Decrease to speed up at cost of worse accuracy. Suggest 0.8 to 3.'
printf, 3,   'accuracy_boost          = 1'
printf, 3, ''

printf, 3,   '#Larger to keep more terms in the hierarchy evolution. '
printf, 3,   'l_accuracy_boost        = 1'
printf, 3, ''

printf, 3,   '#Increase to use more C_l values for interpolation.'
printf, 3,   '#Increasing a bit will improve the polarization accuracy at l up to 200 -'
printf, 3,   '#interpolation errors may be up to 3%'
printf, 3,   '#Decrease to speed up non-flat models a bit'
printf, 3,   'l_sample_boost          = 1'
printf, 3, ''

close, 3

end
