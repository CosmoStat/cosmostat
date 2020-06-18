pro rd_planck_fisher, cosmo, planck, file=file
; March 09 - Modified by Anais Rassat: Fixed bug which had sigma8 and
;            Omega_b parameters inverted (from icosmo_v1.0 onwards).
; Oct 08 - Written by Anais Rassat
; PURPOSE: Reads output from Jochen Weller's CAMB Planck Fisher
; Matrix code.

;Define the Planck matrix, all in double precision
planck = dblarr(8,8)
fish = double(0)

;Path to the Planck fisher matrix file
if not keyword_set(file) then file='/Users/arassat/Work/icosmo_v1.1/camb/CAMB_jw_Planck/Data/fisher_planck_EUCLID.dat'

;Open file to read
openr, inunit, file, /get_lun

;index changes parameter order to the one in mk_fisher_lens
n = 35 & index = [0, 7, 3, 5, 4, 1, 2, 6]
for i = 0, 36-1 do begin
   readf, inunit, fish, a, b;, format='(F25.14,I5,I5)'
   planck(index[a-1], index[b-1]) = fish
   planck(index[b-1], index[a-1]) = fish ;symmetric
endfor

;Close file and Deallocate unit
free_lun, inunit

;Create rest of 'fisher' structure for Planck fisher matrix: 
;...................................
;Names of parameters varied
name =  ["!7X!6!im!n",  "w0", "wa", "h",  "!7X!6!ib!n", "!4o!6!d8!n", "n", "!7X!6!il!n"]

;Central values of parameters varied
p=[cosmo.const.omega_m,cosmo.const.w0,cosmo.const.wa,$
   cosmo.const.h,cosmo.const.omega_b,cosmo.const.sigma8,$
   cosmo.const.n,cosmo.const.omega_l]

;Covariance matrix: 
cov = la_invert(planck, /double)

;Marginalised 1 sigma errors: 
sigma_marg = sqrt(diag_matrix(cov))

;Fixed 1 sigma errors: 
sigma_fix = 1.d0/(sqrt(diag_matrix(planck)))

;Check status of planck and covariance matrices
f_check=check_matrix(planck)
cov_check=check_matrix(cov)
status={f_d:f_check.CHECK_D,f_ev:f_check.CHECK_ev,$
        cov_d:cov_check.CHECK_D,cov_ev:cov_check.CHECK_ev}

;Create planck fisher matrix
planck={pname:name,p:p,f:planck,cov:cov,sigma_marg:sigma_marg,sigma_fix:sigma_fix,status:status}
end
