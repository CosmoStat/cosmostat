function mk_fisher_planck,fid, cosmo, verbose=verbose, cambpath=cambpath, params=params, file=file, time=time
;Nov 08 - modified by An.R. - added time keyword
;Nov 08 - modified by An.R. - fid now an entry. cambpath now read
;         directly from fid structure
;Oct 08 - Modified by Anais Rassat - now does call in CAMB directory
;         and comes back to the current directory
;Oct 08 - Written by Anais Rassat
;Based on call_camb written in June 08 by An.R.
;PURPOSE: Call Planck fisher matrix code written by Jochen Weller and
;read in the output file
;NB: You can run this routine in any directory and it will
;automatically call the CAMB routines. The default CAMB directory is
;defined in fid structure (created by set_fiducial.pro)


;if not keyword_set(path) then cambpath = '/Users/arassat/Work/icosmo_v1.1/camb/CAMB_jw_Planck/CAMB/' else path=path
if not keyword_set(params) then params = "cambtemp.ini"
if not keyword_set(file) then file = '../DATA/fisher_planck_EUCLID.dat'

;Write out the params.ini file need by CAMB
wt_params_camb, cosmo, file=params, verbose=verbose

;Find out current directory & change to directory where CAMB stuff is
cambpath = fid.calc.cambpath; '~/Work/CAMB_for_iCosmo/CAMB/'
spawn, 'cp -f '+string(params)+' '+cambpath
cd, cambpath, current=old_dir   ;old_dir is the directory in which we are before changing to the CAMB directory

;CAMB command directly in shell
if not keyword_set(time) then camb = strcompress(cambpath+'./fisher_EUCLID',/remove_all) else camb = 'time '+strcompress(cambpath+'./fisher_EUCLID', /remove_all)

command =strcompress(camb+' '+string(params)) ; space between path and string(f) is important
spawn, command

;Read in the output planck fisher matrix
rd_planck_fisher, cosmo,planck, file=file
;Return to current directory
cd, old_dir
return, planck

end
