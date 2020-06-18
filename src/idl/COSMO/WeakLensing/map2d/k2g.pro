;+
; NAME: 
;				K2G
;
; PURPOSE: 
;				Calculate weak lensing shear data  from the E mode.
;
; CALLING:
;				G =  k2g(Emode, psize=psize)
;
; INPUT: 
;				Emode ---- IDL array [*,*] :  E  mode 
;
; OUTPUT:
;               G  ---- IDL array [*,*,0:1] :   shear. gamma 1  = G[*,*,0] and gamma 2 = G[*,*,1] 
;
; INPUT KEYWORD:			
;              psize --- float : pixel size 
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

function k2g, K, psize=psize
; psize in degree

if not keyword_set(psize) then psize=0.98
vs = size(K)
n1 = vs[1]
n2 = vs[2]
mk_kappa, k, n1, n2, psize, psize, amin=amin, m
kappa_to_gamma, M,  gamma, cat = cat
G = dblarr(n1,n2,2)
G[*,*,0] = Gamma.Gamma1
G[*,*,1] = Gamma.Gamma2
; tvxy, g, /pola, skip=15, ima=k, leng=80
return, G
end


