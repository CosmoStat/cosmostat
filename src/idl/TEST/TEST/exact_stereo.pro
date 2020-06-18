;+
; NAME: 
;      EXACT_GNOME
;
; PURPOSE: 
;      Do the gnomonic projection
;
; CALLING:
;      exact_gnome

;
; INPUT:
;     theta, phi --- coordinates of the pixel to be projected
;     theta0, phi0 --- center coordinates
; 
; OUTPUT:
;     x, y --- pixel coordinates in the patch 
;
; HISTORY:
;	Written: Sandrine Pires Sept 2009
;-
pro exact_stereo, theta, phi, theta0, phi0, x, y

temp = 1+(cos(theta)*cos(theta0)+sin(theta)*sin(theta0)*cos(phi-phi0))
x =2*((cos(theta)*sin(theta0)-sin(theta)*cos(theta0)*cos(phi-phi0)))/temp
y =2*(sin(theta)*sin(phi-phi0))/temp

end

