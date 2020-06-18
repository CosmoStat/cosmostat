;+
; NAME: 
;				G2EB
;
; PURPOSE: 
;				Calculate weak lensing E and B mode from a shear data set.
;
; CALLING:
;				EB =  g2eb(Gamma, psize=psize)
;
; INPUT: 
;				Gamma ---- IDL array [*,*,0:1] : shear data set
;
; OUTPUT:
;                 EB    ---- IDL array [*,*,0:1] :  E and B modes. E mode = EB[*,*,0] and B mode = EB[*,*,1] 
;
; INPUT KEYWORD:			
;              psize --- float : pixel size 
;              sparse --- scalar : if set, an iterative method is used, with sparse constraints on the border 
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

function g2eb, G, psize=psize, sparse=sparse
; psize in degree
if not keyword_set(psize) then psize=0.98

if not keyword_set(sparse) then begin
	vs = size(G)
	n1 = vs[1]
	n2 = vs[2]
	k = dblarr(n1, n2)
	mk_kappa, k, n1, n2, psize, psize, amin=amin, kappa

	gamma={n1:kappa.n1,n2:kappa.n2,theta1:kappa.theta1,theta2:kappa.theta2,delta1:kappa.delta1,delta2:kappa.delta2,kappa:kappa.kappa,$
       gamma1:G[*,*,0],gamma2:G[*,*,1]}

	gamma_to_kappa, gamma, kappa,  bmode=bmode
	EB = dblarr(n1, n2, 2)
	EB[*,*,0] = kappa.Kappa
	gamma_to_kappa, gamma, kappa, /bmode
	EB[*,*,1] = kappa.Kappa
end else begin
   EB = cst_g2eb(G, BZero=BZero, BBorderZero=BBorderZero, /sparse, soft=soft, Niter=Niter, TrueSolEB=TrueSolEB, psize = psize, nscale=nscale)
end
   

return, EB
end
