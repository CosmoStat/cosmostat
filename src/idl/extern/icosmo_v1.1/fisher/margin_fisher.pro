pro margin_fisher,fishin,fishout,pmarg

; Jan 03 - Written by A. Refregier
;
; PURPOSE: Marginalise a fisher matrix over unwanted parameters
; INPUT: fishin: input fisher matrix structure given by mk_fisher.pro
;        pmarg: vector of flags to indicate whether a parameter is kept (1)
;               or whether it is marginalised over (0). Its size must
;               equal the number of parameters in fishin and it must
;               contain only 0's and 1's.
; OUTPUT: fishout: new (smaller) fisher matrix with only interesting
;            parameters. the other parameters have been marginalised over

; declarations
out=where(pmarg eq 1,n_out)        ; output parameters
marg=where(pmarg eq 0,n_marg)      ; marginalised parameters        
n_in=n_elements(fishin.p)
if not(n_in eq (n_out+n_marg)) then begin
  print,'margin_fisher: incompatible parameter lists'
  stop
endif

; construct output vectors
pname_out=fishin.pname(out)
p_out=fishin.p(out)

; construct output fisher and covariance matrices by marginalising
; over unwanted parameters
cov_out=dblarr(n_out,n_out)
for i=0,n_out-1 do begin
  for j=0,n_out-1 do begin
    cov_out(i,j)=fishin.cov(out(i),out(j))
  endfor
endfor
fish_out=invert(cov_out,/double)

; construct error vectors
sigma_marg_out=dblarr(n_out)
sigma_fix_out=dblarr(n_out)
for i=0,n_out-1 do begin
  sigma_fix_out(i)=sqrt(1.d/fish_out(i,i))
  sigma_marg_out(i)=sqrt(cov_out(i,i))
endfor

; store in structure
fishout={pname:pname_out,p:p_out,f:fish_out,cov:cov_out,$
  sigma_marg:sigma_marg_out,sigma_fix:sigma_fix_out}

end
