function check_matrix, m, d=d, eps=eps, _extra=extra

;Aug 08 modified by AA, machine precision calculated for double 
;Aug 08 modified by AnR, now also checks if evalues and diagonal
;elements are positive and below precision
;Aug 08 modified by AnR, re-added _extra keyword, required to find precision
;Aug 08 modified by AnR, changed 'print' outputs
;Jul 08 modified by AA merged with check_diag and made generic  
;Written by Anais Rassat, July 2008
;Modified by An.Ra. 24th July 2008: Added double precision to
;calculation of eigenvalues. Also speeded up e-values calculation by
;not calculating e-vectors. Added eps keyword (outputs precision).
;PURPOSE: Checks that the input matrix has positive diagonals and
;positive eigen values.
;definite
;
;INPUT: f: Matrix (not structure)
;OPTIONAL INPUT: _extra: any keyword accepted by machar() (e.g.: double)
;OUTPUT: check: 0 no problem
;        check: 1 negative but very small eigenvalues - or positive
;                                                       and very small
;evalues. Problem is
;        ignored.
;        check: 2 negative and not that small. 
; If check = 2, you should not trust the Fisher matrix!!
;OPTIONAL OUTPUT: d: eigenvalues
;                 eps: precision


  mtemp = m                      ;so as to keep input matrix intact

;Find out what the accuracy is
  Eps = (machar(_extra=extra,/double)).eps


general_warning = ' Fisher and Covariance matrices should be positive definite.'

;*** begin: check diagonals ***
check_d=0
ind = where(diag_matrix(mtemp) lt 0.d0, count)
if count ne 0 then if max((abs(diag_matrix(mtemp)))(ind)) le eps then begin
   check_d = 1
   print,'Warning: Matrix has (negative) very small diagonal elements. These are smaller than the precision.'+ general_warning
endif else if count ne 0 then if  max((abs(diag_matrix(mtemp)))(ind))  gt eps then begin
   check_d=2
   print,'Warning: Matrix has negative diagonal elements. These may be small, but not smaller than the precision.'+ general_warning
endif

;Find any positive diaogonal elements that are very close to zero 
ind = where(diag_matrix(mtemp) lt eps and diag_matrix(mtemp) ge 0.d0, count)
if count ne 0 then begin 
   if check_d eq 0 then check_d = 1 ; only do this if check_d = 0 i.e. if there are no other problems
   print, 'Warning: Matrix has positive diagonal elements that are smaller than the precision.'+ general_warning
endif
                           ;*** end: check diagonals ***

;*** begin: check eigen values ***
  check_ev=0
;Compute the triddiagonal form of f:
  trired, mtemp, d, e,/double
;Compute the eigenvalues (returned in vector D) and the eigenvectors
;(returned in the rows of the array mtemp):
  triql_novec, d, e, /double

;Find negative evalues if any
  ind = where (d lt 0.d0, count)
;For negative eigenvalues that are very close to zero
  if count ne 0 then if max(abs(d[ind])) le eps then begin
     check_ev = 1
     print, 'Warning: Matrix has (negaitve) very small eigenvalues. Problem is probably due to precision. Problem is ignored.'+ general_warning
  endif else if count ne 0 then if max(abs(d[ind])) gt eps then begin
     check_ev=2
     print,'Warning: Matrix has negative eigenvalues. These may be small, but not smaller than the precision.'+ general_warning
  endif

;Find any positive evalues that are very close to zero 
ind = where(d lt eps and d ge 0.d0, count)
if count ne 0 then begin 
   if check_ev eq 0 then check_ev = 1 ; only do this if check_ev = 0 i.e. if there are no other problems
   print, 'Warning: Matrix has positive eigenvalues that are smaller than the precision.'+ general_warning
endif

;*** end: check eigen values ***
  
; construct structure:  
check={check_d:check_d,check_ev:check_ev}
return, check
end
