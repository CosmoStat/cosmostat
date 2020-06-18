;+
; NAME:
;     rotz
;
; PURPOSE:
;     Computes the rotation matrix for rotations about the z axis.
;     
; CALLING SEQUENCE:
;     rot_matrix = rotz(angle)
;
; INPUTS:
;     angle - scalar, float - The rotation angle, in RADIANS.
;
; OUTPUTS:  
;     rot_matrix - dblarr - The 3x3 rotation matrix.
;
; KEYWORDS:
;     None.
;
; COMMENTS:
;     Specify a single angle only.
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 04 May 1999
;
;-
;======================================================================
;
function rotz,angle
;
; rotation matrix for rotation of angle 'angle' about the z axis
;    input: angle in radians
;    output: 3x3 rotation matrix


matrix = [[ cos(angle) , sin(angle),  0.0d0], $
          [-sin(angle) , cos(angle),  0.0d0], $
          [0.0d0       , 0.0d0     ,  1.0d0]    ]

return, matrix
end
