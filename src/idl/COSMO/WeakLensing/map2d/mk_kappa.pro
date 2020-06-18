;+
; NAME: 
;				MK_KAPPA
;
; PURPOSE: 
;				to create a structure from mass map
;
; CALLING:
;				mk_kappa, k, n1, n2, theta1, theta2, /amin, m
;
; INPUTS: 
;				k --- mass map
;				n1, n2 --- map size in pixel
;				theta1, theta2 --- map size in amin or degree 
; KEYWORD:
;				amin --- if map size in amin
;
; OUTPUTS: 
;				m --- mass map structure 
; CALLING:
;				mk_kappa, map, 1024, 1024,200, 200, amin='amin', m
;
; HISTORY:
;				Written: S. Pires Nov. 05
;-
;-------------------------------------------------------------------------------
pro mk_kappa, k, n1, n2, theta1, theta2, amin=amin, m

;if map size in degrees
t1=theta1/!radeg   ; map size [rad]
t2=theta2/!radeg   ; map size [rad]
delta1=t1/float(n1)    ; pixel size [rad]
delta2=t2/float(n2)    ; pixel size [rad]

;if map size in arcmin
if keyword_set(amin) then begin
  t1=theta1/60./!radeg   ; map size [rad]
  t2=theta2/60./!radeg   ; map size [rad]
  delta1=t1/float(n1)    ; pixel size [rad]
  delta2=t2/float(n2)    ; pixel size [rad]
endif

; store in structure
m={n1:n1,n2:n2,theta1:t1,theta2:t2,delta1:delta1,delta2:delta1,kappa:k}
end
