;+
; NAME:
;     ROTX
;
; PURPOSE:
;     Computes the rotation matrix for rotations about the x axis.
;     
; CALLING SEQUENCE:
;     rot_matrix = rotx(angle)
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
function rotx,angle
;
; rotation matrix for rotation of angle 'angle' about the x axis
;    input: angle in radians
;    output: 3x3 rotation matrix


matrix = [[1.0d0,  0.0d0      , 0.0d0     ], $
          [0.0d0,  cos(angle) , sin(angle)], $
          [0.0d0, -sin(angle) , cos(angle)]    ]

return, matrix
end
