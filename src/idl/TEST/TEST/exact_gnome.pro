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
pro exact_gnome, theta, phi, theta0, phi0, x, y

temp = (cos(theta)*cos(theta0)+sin(theta)*sin(theta0)*cos(phi-phi0))
x =(cos(theta)*sin(theta0)-sin(theta)*cos(theta0)*cos(phi-phi0))/temp
y =(sin(theta)*sin(phi-phi0))/temp

end

