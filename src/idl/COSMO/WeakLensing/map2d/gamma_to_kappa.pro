;+
; NAME: 
;				GAMMA_TO_KAPPA
;
; PURPOSE: 
;				to compute the convergence kappa from the 
;				shear gamma
; 				simple output structure
;
; CALLING:
;				gamma_to_kappa, gamma, kappa, /cat, /bmode
;
; INPUTS: 
;				gamma ---  shear structure 
;
; KEYWORD:			bmode --- to compute the B-mode from the 
;				shear gamma
;				cat --- to use a structure with more 
;				parameters about the catalogue (real data) 
; OUTPUT: 
;				kappa --- kappa structure 
;
; HISTORY:
;				Written by A. Refregier & S. Pires Nov 2005
;
;-
;-------------------------------------------------------------------------------
pro gamma_to_kappa, gamma, kappa, cat = cat, bmode = bmode

; check dimensions
n1=gamma.n1
n2=gamma.n2
theta1=gamma.theta1
theta2=gamma.theta2
delta1=gamma.delta1
delta2=gamma.delta2


; compute B-mode by rotating shears by 45deg if requested
if not keyword_set(bmode) then begin   
  gamma1=gamma.gamma1
  gamma2=gamma.gamma2
endif else begin
  gamma1=-gamma.gamma2
  gamma2=gamma.gamma1
endelse


; take FFT of kappa map
;print,'Taking FFT of gamma'
gamma1_tilde=fft(gamma1)
gamma2_tilde=fft(gamma2)

; compute corresponding l values
;Computing l1 values
l1=indgen(n1)#replicate(1,n2)-(n1-1)/2
l1=shift(l1,-(n1-1)/2,-(n1-1)/2)


;Computing l2 values
l2=replicate(1,n1)#indgen(n2)-(n2-1)/2
l2=shift(l2,-(n2-1)/2,-(n2-1)/2)

l1=2.*!pi/theta1*float(l1)
l2=2.*!pi/theta2*float(l2)


; compute FT Psi_tilde of lensing potentional Psi, as well
; as other quantities
kappa_tilde=((l1^2-l2^2)*gamma1_tilde+(2.*l1*l2)*gamma2_tilde)/(l1^2+l2^2)
kappa_tilde(0,0)=0.       ; set k=0 mode to 0.
sz=size(kappa_tilde)
;print, sz[1]
;print, sz[2]
; take inverse FFT
k=float(fft(kappa_tilde,/inverse))

; store in structure

kappa={n1:gamma.n1,theta1:gamma.theta1,delta1:gamma.delta1,$
       n2:gamma.n2,theta2:gamma.theta2,delta2:gamma.delta2,$
       kappa:k,$
       gamma1:gamma.gamma1,gamma2:gamma.gamma2}

if keyword_set(cat) then begin
kappa={n1:gamma.n1,theta1:gamma.theta1,delta1:gamma.delta1,$
       n2:gamma.n2,theta2:gamma.theta2,delta2:gamma.delta2,$
       kappa:k,$
       gamma1:gamma.gamma1,gamma2:gamma.gamma2,$
       ng:gamma.ng, wtot:gamma.wtot, mask:gamma.mask,$
       gamma_err:gamma.gamma_err,kappa_err:gamma.kappa_err,ng_eff:gamma.ng_eff,$
       sigma_gamma:gamma.sigma_gamma, x1_ran:gamma.x1_ran,$
       x2_ran:gamma.x2_ran, x1_m:gamma.x1_m, x2_m:gamma.x2_m}
end

end
