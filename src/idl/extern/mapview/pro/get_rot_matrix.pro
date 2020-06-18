;+
; NAME:
;     GET_ROT_MATRIX
;
; PURPOSE:
;     Function to return the appropriate rotation matrix that will convert
;     between ecliptic, celestial and Galactic coordinates in J2000.
;     
; CALLING SEQUENCE:
;     matrix = get_rot_matrix(rot_code)
;
; INPUTS:
;     code - char string - Character string specifying the desired
;                         coordinate conversion.  Only ONE code
;                         may be requested at a time.  Valid codes are:
;                         'c2e' - from = celestial   to = ecliptic
;                         'e2c' - from = ecliptic    to = celestial
;                         'g2e' - from = galactic    to = ecliptic
;                         'e2g' - from = ecliptic    to = galactic
;                         'c2g' - from = celestial   to = galactic
;                         'g2c' - from = galactic    to = celestial
;
;
; OUTPUTS:  
;     matrix - dblarr -   The 3x3 rotation matrix appropriate to
;                         the specified code.  The rotation matrix operates
;                         on a unit vector with 3 ROWS: x,y,z.
;
; KEYWORDS:
;     None.
;
; ROUTINES CALLED:
;     rotx(), rotz()
;
; EXAMPLE:
;     Convert a single unit vector in celestial to galactic.  IDL's
;     (column x row) convention makes this a little tricky:
;       cel_uvec   = [1,0,0]
;       rot_matrix = get_rot_matrix('c2g')
;       gal_uvec   = transpose(rot_matrix) # cel_uvec
;     
; COMMENTS:
;     Intended as lower-level routine used by COORTRANS.PRO.
;     
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 04 May 1999
;
;-
;======================================================================
function get_rot_matrix,rot_code
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
if (N_params() LT 1 ) then begin
    print,'Syntax: result = get_rot_matrix(rot_code)'
    return,-1
  endif
;
d2r = !DPI/180.d0    ; degrees-to-radians conversion factor
;
; angle is mean obliquity of the ecliptic at J2000 which is exactly
; 23d 26' 21.448" (Table 1.2.2 Hipparcos Explanatory Guide) 
;
  angle = 23.4392911111d0 * d2r
;
; rotation angles for galactic/celestial conversion
; see j2000_coor_trans.tex
;
  psi   = 282.85948d0 * d2r
  theta =  62.871750d0 * d2r
  phi   = 327.06808d0 * d2r

;
; The rotation matrix for equatorial to Galactic in J2000
; and vice-versa. are from the Hipparcos Explanatory Guide
; The rotation matrix for ecliptic to Galactic and vice-versa is computed
; by combining rotation matrices above
;  i.e., 'g2e' = 'g2c'#'c2e'  & 'e2g' = 'e2c'#'c2g'
;

 case strlowcase(strtrim(rot_code,2)) of 
 'c2e':  return, rotx(angle)
 'e2c':  return, rotx(-angle)
 'c2g':  return, rotz(psi)#rotx(theta)#rotz(phi)
 'g2c':  return, rotz(-phi)#rotx(-theta)#rotz(-psi)
 'e2g':  return, rotx(-angle)#rotz(psi)#rotx(theta)#rotz(phi)
 'g2e':  return, rotz(-phi)#rotx(-theta)#rotz(-psi)#rotx(angle) 
  else:  begin
         message,'Invalid rotation code string: ', rot_code,/CON
	 return,-1
	 end
 endcase	
;
end







