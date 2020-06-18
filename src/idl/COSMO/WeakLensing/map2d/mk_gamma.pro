;+
; NAME: 
;				MK_GAMMA
;
; PURPOSE: 
;				to create a structure from shear maps
;
; CALLING:
;				mk_gamma, gamma1, gamma2, n1, n2, theta1, $
;				theta2, g, /amin
;
; INPUTS: 
;				gamma1, gamma2 --- shear maps
;				n1, n2 --- map size in pixel
;				theta1, theta2 --- map size in amin or degree 
; KEYWORD:
;				amin --- if map size in amin
;
; OUTPUTS: 
;				g --- gamma structure 
; CALLING:
;				mk_gamma, g1, g2, 1024, 1024,200, 200, $
;				amin='amin', g
;
; HISTORY:
;				Written by S. Pires Nov. 05
;-
;-------------------------------------------------------------------------------
pro mk_gamma, gamma1, gamma2, n1, n2, theta1, theta2, amin = amin, g

;if map size in degrees
t1=theta1/!radeg   ; map size [rad]
t2=theta2/!radeg   ; map size [rad]
delta1=t1/float(n1)    ; pixel size [rad]
delta2=t2/float(n2)    ; pixel size [rad]

;if map size in seconds
if keyword_set(amin) then begin
  t1=theta1/60./!radeg   ; map size [rad]
  t2=theta2/60./!radeg   ; map size [rad]
  delta1=t1/float(n1)    ; pixel size [rad]
  delta2=t2/float(n2)    ; pixel size [rad]
endif

; store in structure
g={n1:n1,n2:n2,theta1:t1,theta2:t2,delta1:delta1,delta2:delta1,gamma1:gamma1, gamma2:gamma2}
end
