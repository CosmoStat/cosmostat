;+
; NAME: 
;				G2K
;
; PURPOSE: 
;				Calculate weak lensing E mode from a shear data set.
;
; CALLING:
;				Emode =  g2k(Gamma, psize=psize, sparse=sparse)
;
; INPUT: 
;				Gamma ---- IDL array [*,*,0:1] : shear data set
;
; OUTPUT:
;                 Emode    ---- IDL array [*,*] :  E mode  
;
; INPUT KEYWORD:			
;              psize --- float : pixel size 
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

function g2k, G, psize=psize
; psize in degree

if not keyword_set(psize) then psize=0.98
vs = size(G)
n1 = vs[1]
n2 = vs[2]
k = dblarr(n1, n2)
mk_kappa, k, n1, n2, psize, psize, amin=amin, kappa

gamma={n1:kappa.n1,n2:kappa.n2,theta1:kappa.theta1,theta2:kappa.theta2,delta1:kappa.delta1,delta2:kappa.delta2,kappa:kappa.kappa,$
       gamma1:G[*,*,0],gamma2:G[*,*,1]}

gamma_to_kappa, gamma, kappa,  bmode=bmode
return, kappa.Kappa
end

