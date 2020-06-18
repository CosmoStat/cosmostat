;+
; NAME: 
;				KAPPA_TO_GAMMA
;
; PURPOSE: 
;				to compute the shear gamma from the 
;				projected mass map kappa
;
; CALLING:
;				kappa_to_gamma, kappa, gamma, /cat
;
; INPUTS: 
;				kappa --- kappa structure 
;
; KEYWORD:			cat --- to use a structure with more 
; 				parameters (about the catalogue) 
;
; OUTPUT: 
;				gamma --- shear structure 
;
; HISTORY:
;				Written by A. Refregier & S. Pires Nov 2005
;-
;-------------------------------------------------------------------------------

pro kappa_to_gamma,kappa,gamma, cat = cat

; check dimensions
n1=kappa.n1
n2=kappa.n2
theta1=kappa.theta1
theta2=kappa.theta2
delta1=kappa.delta1
delta2=kappa.delta2

; take FFT of kappa map
kappa_tilde=fft(kappa.kappa)

; compute corresponding l values
l1=indgen(n1)#replicate(1,n2)-(n1-1)/2
l1=shift(l1,-(n1-1)/2,-(n1-1)/2)

l2=replicate(1,n1)#indgen(n2)-(n2-1)/2
l2=shift(l2,-(n2-1)/2,-(n2-1)/2)

l1=2.*!pi/kappa.theta1*float(l1)
l2=2.*!pi/kappa.theta2*float(l2)

; compute FT Psi_tilde of lensing potentional Psi, as well
; as other quantities
psi_tilde=-2./(l1^2+l2^2)*kappa_tilde
psi_tilde(0,0)=0.        ; set k=0 mode to 0.
gamma1_tilde=-.5*(l1^2-l2^2)*psi_tilde
gamma2_tilde=-l1*l2*psi_tilde

; take inverse FFT
gamma1=float(fft(gamma1_tilde,/inverse))
gamma2=float(fft(gamma2_tilde,/inverse))

; store in structure
gamma={n1:kappa.n1,n2:kappa.n2,theta1:kappa.theta1,theta2:kappa.theta2,delta1:kappa.delta1,delta2:kappa.delta2,kappa:kappa.kappa,$
       gamma1:gamma1,gamma2:gamma2}

if keyword_set(cat) then begin
kappa={n1:kappa.n1,theta1:kappa.theta1,delta1:kappa.delta1,$
       n2:kappa.n2,theta2:kappa.theta2,delta2:kappa.delta2,$
       kappa:kappa.kappa,$
       gamma1:gamma1,gamma2:gamma2,$
       ng:kappa.ng, wtot:kappa.wtot, mask:kappa.mask,$
       gamma_err:kappa.gamma_err,kappa_err:kappa.kappa_err,ng_eff:kappa.ng_eff,$
       sigma_gamma:kappa.sigma_gamma, x1_ran:kappa.x1_ran,$
       x2_ran:kappa.x2_ran, x1_m:kappa.x1_m, x2_m:kappa.x2_m}
end


end
