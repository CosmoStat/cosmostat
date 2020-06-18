;+
; NAME: 
;				EB2G
;
; PURPOSE: 
;				Calculate weak lensing shear data  from E and B mode.
;
; CALLING:
;				G =  g2eb(EB, psize=psize, sparse=sparse)
;
; INPUT: 
;				EB ---- IDL array [*,*,0:1] :  E and B modes. E mode = EB[*,*,0] and B mode = EB[*,*,1] 
;
; OUTPUT:
;                 G    ---- IDL array [*,*,0:1] :   shear. gamma 1  = G[*,*,0] and gamma 2 = G[*,*,1] 
;
; INPUT KEYWORD:			
;              psize --- float : pixel size 
;              sparse --- scalar : if set, an iterative method is used, with sparse constraints on the border 
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

function EB2G, EB, psize=psize, sparse=sparse

if not keyword_set(psize) then psize=0.98

if not keyword_set(sparse) then begin
	vs = size(EB)
	n1 = vs[1]
	n2 = vs[2]
	k1 = EB[*,*,0]
	k2 = EB[*,*,1]
; mk_kappa, k, n1, n2, psize, psize, amin=amin, m
	mk_kappa, k1, n1, n2, psize, psize, m1
	mk_kappa, k2, n1, n2, psize, psize, m2
	kappa_to_gamma, m1, g1
	kappa_to_gamma, m2, g2
   ; kappa_to_gamma, M,  gamma, cat = cat
	G = dblarr(n1,n2,2)
	G[*,*,0] = g1.gamma1 + g2.gamma2
	G[*,*,1] = g1.gamma2 - g2.gamma1
end else begin
   G = cst_eb2g(EB, BZero=BZero, BBorderZero=BBorderZero, /sparse, soft=soft, Niter=Niter, TrueSolGamma=TrueSolGamma, psize=psize, nscale=nscale)
 end
   
return, G
end


