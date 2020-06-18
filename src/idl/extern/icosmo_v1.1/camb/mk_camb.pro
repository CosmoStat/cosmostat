function mk_camb,fid,cosmo, verbose=verbose, params=params, filetk=filetk, filepk=filepk, filecl=filecl, time=time

; Dec 08 - Written by Anais Rassat
;
; Based on mk_fisher_planck  written in Oct. 08 by An.R.
; PURPOSE: Call CAMB 

; NB: You can run this routine in any directory and it will
; automatically call the CAMB routines. The default CAMB directory is
; defined in fid structure (created by set_fiducial.pro)


;Define here the names of various files to use
if not keyword_set(params) then params = "cambtemp.ini"
if not keyword_set(filetk) then filetk = 'transfer_out.dat'
if not keyword_set(filepk) then filepk = 'matterpower.dat'
if not keyword_set(filecl) then filecl = 'scalCls.dat'

;Write out the params.ini file need by CAMB

wt_params_camb_general, fid, cosmo, fileparams=params, filetk=filetk, filepk=filepk, filecl=filecl

;Find out current directory & change to directory where CAMB stuff is
cambpath = fid.calc.cambpath
;Copy the params.ini file to the directory where CAMB is located
spawn, 'cp -f '+string(params)+' '+cambpath
cd, cambpath, current=old_dir   ;old_dir is the directory in which we are before changing to the CAMB directory

;CAMB command directly in shell
if not keyword_set(time) then camb = strcompress(cambpath+'./camb',/remove_all) else camb = 'time '+strcompress(cambpath+'./camb', /remove_all)

command =strcompress(camb+' '+string(params)) ; space between path and string(f) is important
spawn, command

;Eventually add here a command that reads in the T(k) in required format
; rd_camb_tk, cosmo, tk, file=filetk

;Eventually add here a command that reads in the P(k) in required format
; rd_camb_pk, cosmo, pk, file = filepk

;Eventually add here a command that reads in the C(l) in required format
rd_camb_cl, fid, cosmo, cl, file = filecl

;Return to current directory
cd, old_dir
return, cl

end
