;+
; NAME:
;     q2m
;
; PURPOSE:
;     Converts quaternions to rotation matrices.
;     
; CALLING SEQUENCE:
;     q2m,q,a
;
; INPUTS:
;     q - Quaternions.  May be a 2-D array dimensioned 4xN or
;         simply a vector dimensioned 4.  If q is double, then
;         the computation is done in double precision and
;         a double-precision result is returned.
;
; OUTPUTS:  
;     a - Cube of attitude rotation matrices, 3x3xN (or 3x3
;         if only one input quaternion).
;
; KEYWORDS:
;     SINGLE - Set this flag to force a single-precision result,
;              e.g., if storage is an issue.
;
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     None.
;
; COMMENTS:
;     Operates on a list of quaternions or a single quaternion.
;     
; MODIFICATION HISTORY:
;     Initial version, G. Hinshaw, 30 Sep 1997
;     Vectorized; double precision option.  
;         R. S. Hill, RITSS, 16 Apr 2001
;-
;======================================================================
;
PRO Q2M,Q,A, SINGLE=single

; Converts quaternions q into rotation matrices a.

; Convert to CSC style quaternion:  conversion uses V(cel)=A*V(sc) 
; but CSC uses V(sc)=A*V(cel), so use -(spatial elements).

noforce = 1b - keyword_set(single)
dub = (size(q, /type) EQ 5) AND noforce

message_not_displayed = 1

q1=-reform(q[0,*])
q2=-reform(q[1,*])
q3=-reform(q[2,*])
q4= reform(q[3,*])

NORMALIZATION_LOOP:
q11=q1*q1
q22=q2*q2
q33=q3*q3
q44=q4*q4
s=q11+q22+q33+q44
w = where(abs(s-1.0) GT 1.0e-5, nw)
IF nw GT 0 THEN BEGIN
    IF message_not_displayed THEN $
        print,'Normalization error in q2m'
    message_not_displayed = 0
    s=sqrt(s)
    q1=q1/s
    q2=q2/s
    q3=q3/s
    q4=q4/s
    GOTO, NORMALIZATION_LOOP
ENDIF

q12=q1*q2
q13=q1*q3
q14=q1*q4
q23=q2*q3
q24=q2*q4
q34=q3*q4

IF dub THEN a = dblarr(n_elements(q1),3,3) $
       ELSE a = fltarr(n_elements(q1),3,3)

a[0,0,0] = q11 - q22 - q33 + q44
a[0,0,1] = 2. * ( q12 + q34 )
a[0,0,2] = 2. * ( q13 - q24 )
a[0,1,0] = 2. * ( q12 - q34 )
a[0,1,1] = -q11 + q22 - q33 + q44
a[0,1,2] = 2. * ( q23 + q14 )
a[0,2,0] = 2. * ( q13 + q24 )
a[0,2,1] = 2. * ( q23 - q14 )
a[0,2,2] = -q11 - q22 + q33 + q44

a = transpose(a, [1,2,0])

RETURN
END
