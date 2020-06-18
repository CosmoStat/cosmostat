function comb_fisher,fishera,fisherb,printit=printit

; Aug 08 - modified by AA - able to combine Fisher with different dimensions
; Jun 2002 - Written by A. Refregier
;
; PURPOSE: compute the fisher matrix for a combination of 2 measurements
; INPUT: fisher1,2: fisher matrices for each measurements from
;                   mk_fisher.pro. The parameters and their central values
;                   are assumed to be the same
; OPTIONAL INPUT: printit: print results if requested
; OUTPUT: the fisher matrix structure for the combined measurement is
;         returned
;
; !!! NEED TO CHECK THAT THE ADDITION CAN BE DONE !!!

;make fisher1 the one with the most elements
if (n_elements(fishera.p) ge n_elements(fisherb.p)) then begin
   fisher1=fishera
   fisher2=fisherb
endif else begin
   fisher2=fishera
   fisher1=fisherb
endelse

; declarations
n_p1=n_elements(fisher1.p)
n_p2=n_elements(fisher2.p)
ind1=make_array(n_p2,/int)
ind2=make_array(n_p2,/int)

;Find where elements of fisher matricies are the same
num1=0L
num2=0L
for i=0,n_p2-1 do begin
   temp=where(fisher1.pname eq fisher2.pname(i))
   if (temp(0) ne -1) then begin
      ind1(num1)=i
      num1=num1+1L
 ;     stop
   endif else begin
      ind2(num2)=i
      num2=num2+1L
   endelse
endfor
if (num1 eq 0) then ind1=-1 else ind1=ind1(0:num1-1)
if (num2 eq 0) then ind2=-1 else ind2=ind2(0:num2-1)

;total number of cosmology elements
n_pt=n_p1+n_p2-num1

; NEED TO ADD ERROR CHECKS TO MAKE SURE THAT COSMOLOGIES FOR THE
; FISHER MATRICIES ARE THE SAME.

;construct new pname array
pname=make_array(n_pt,/string)
pname(0:n_p1-1)=fisher1.pname
for i=0,num2-1 do pname(n_p1+i)=fisher2.pname(ind2(i))

;construct p array
p=make_array(n_pt)
p(0:n_p1-1)=fisher1.p
for i=0,num2-1 do p(n_p1+i)=fisher2.p(ind2(i))

;construct fisher matrix
ft1=make_array(n_pt,n_pt,/double)
ft2=make_array(n_pt,n_pt,/double)
ft1(0:n_p1-1,0:n_p1-1)=fisher1.f

ind3=make_array(n_p2,/int)
for i=0,n_p2-1 do begin
   ind3(i)=where(pname eq fisher2.pname(i))
endfor

for i=0,n_p2-1 do begin
   for j=0,n_p2-1 do begin
      ft2(ind3(i),ind3(j))=fisher2.f(i,j)
   endfor
endfor

f=ft1+ft2
;stop

; construct the rest of the fisher strucutre
cov=invert(f)
sigma_fix=dblarr(n_pt)
sigma_marg=dblarr(n_pt)
for i=0,n_pt-1 do begin
  sigma_fix(i)=sqrt(1./f(i,i))
  sigma_marg(i)=sqrt(cov(i,i))
endfor

; print results if requested
if keyword_set(printit) then begin
  print,'combined fisher matrix: F_ij:',f
  print,'sigma (fixed):',sigma_fix
  print,'cov(pi,pj):',cov
  print,'sigma (marginalised):',sigma_marg
endif

; store in structure
fish={pname:pname,p:p,f:f,cov:cov,sigma_marg:sigma_marg,$
  sigma_fix:sigma_fix}

return,fish
end
