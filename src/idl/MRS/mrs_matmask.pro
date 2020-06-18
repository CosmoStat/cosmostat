;+
; NAME:
;        mrs_matmask
;
; PURPOSE:
;   Return the matrix Mat related to a mask.
;   For a given mask M, the power spectrum of an image I is equal to:
;         Mat # PowSpec(I) = PowSpec ( M . I ) 
;
; CALLING:
;     Mat = mrs_matmask(Mask)
;
; INPUTS:
;     Mask -- IDL array of healpix map: Input mask 
;    
; OUTPUTS:
;     Mat -- 2D IDL fltarr: Matrix related to the mask.
;
; EXAMPLE:
;       Compute the power spectrum of an image. 
;               Mat = matmask(mask) 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2009
;	March, 2009 File creation
;--------------------------------------------------------------------------------------------------- 


;================================================================

function mrs_mwigner3j2,il1,il2,il3,lnwpro
; This routine is extract from the MASTER package
l1=long(il1)
l2=long(il2)
l3=long(il3)
L=l1+l2+l3
L_2=L/2
min=abs(l1-l2)
max=l1+l2
c=l3*0d
w=where((L_2*2-L) eq 0 and l3 ge min and l3 le max)
if w(0) ne -1 then begin
    lnw1=lnwpro(L_2(w)-l1)
    lnw2=lnwpro(L_2(w)-l2)
    lnw3=lnwpro(L_2(w)-l3(w))
    lnwl=lnwpro(L_2(w))
    lnc=-alog(L(w)+1d)-lnwl+lnw1+lnw2+lnw3
    c(w)=exp(lnc)
endif
return,c
end

;================================================================

function mrs_mmake_mll,ell,well,lmax
; This routine is extract from the MASTER package
ellm=ell(0:lmax)
; master_wigner_init,n_elements(ell),lnwpro
k=indgen(2+n_elements(ell)*3)
lnwpro=lngamma(2*k+1d)-2*lngamma(k+1d)

m=dblarr(n_elements(ellm),n_elements(ellm))

 ;  title='Computing the Mll matrix ...',/xsty,/ysty
for l1=0,n_elements(ellm)-1 do begin 
   c=(2*ell+1)*well  
  for l2=0,n_elements(ellm)-1 do begin 
     b=c*(mrs_mwigner3j2(ellm(l1),ellm(l2),ell,lnwpro))  
     a=total(b)   
     m(l1,l2)=(2*ellm(l2)+1)*a/(4.*!pi)  
  end 
end
return,m
end

;=====================================================================

function mrs_matmask, Mask, lmax=lmax, PSMask=PSMask 
COMMON MR1ENV

if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: Mat = matmask(Mask)'
	    MatMask =-1
        goto, DONE
        end

if not keyword_set(PSMask ) then pm = double ( mrs_powspec(mask, /nonorm, lmax=lmax) ) $
else pm = double ( PSMask )
ell = dindgen(n_elements(pm))
lmax = n_elements(pm) -1
MatMask= mrs_mmake_mll(ell, pm, lmax) 


DONE: 
return, MatMask
end

;====================================================================
